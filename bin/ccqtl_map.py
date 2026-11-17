#!/usr/bin/env python
"""
ccQTL association scan for one chromosome with tensorQTL's map_trans: every variant
(MAF >= --maf in the analysis' samples) against every cell-type composition phenotype
of every analysis (<name>.phenotype.bed + <name>.covariates.tsv from ccqtl_composition.py).
Genotype dosages are read from the plink2 pgen in chunks (pgenlib), so the whole
chromosome is never held as a DataFrame.

Scan mode (--n-perm 0), outputs in --outdir:
  assoc/<analysis>/<celltype>/<celltype>__chr<chrom>.tsv.gz
        every test, unfiltered: variant_id chrom pos ref alt af n_samples b b_se t pval
  stats/<analysis>__chr<chrom>.hist.tsv     -log10(p) histogram per cell type (QQ / lambda)
  stats/<analysis>__chr<chrom>.thin.tsv.gz  p < 1e-3 plus a 1% random sample of the rest (Manhattan)
  stats/<analysis>__chr<chrom>.lookup.tsv   rows for --lookup-variants

Permutation mode (--n-perm N > 0): each cell type is residualised on the covariates,
samples are shuffled (one shuffle per permutation, shared by all cell types of an
analysis so their correlation is kept), and the chromosome is rescanned. Outputs:
  perm/<analysis>__chr<chrom>.perm_stats.tsv  per permutation x cell type: n_tests, lambda,
                                              min_p, n_p_lt_1e-5, n_p_lt_5e-8
  perm/<analysis>__chr<chrom>.perm_hist.tsv   -log10(p) histogram pooled over permutations
"""
import argparse
import glob
import gzip
import inspect
import os
import time

import numpy as np
import pandas as pd
from scipy import stats

# -log10(p) bins: fine around the median (lambda), coarser in the tail
HIST_EDGES = np.concatenate([np.arange(0, 2, 0.001), np.arange(2, 10, 0.01), np.arange(10, 50, 0.5), [50, np.inf]])
CHI2_MEDIAN = stats.chi2.ppf(0.5, 1)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pgen-prefix", required=True, help="plink2 prefix (.pgen/.pvar/.psam)")
    p.add_argument("--chrom", required=True, help="Chromosome as in the pvar (with or without 'chr')")
    p.add_argument("--analysis-dir", default=".", help="Directory with <name>.phenotype.bed / <name>.covariates.tsv")
    p.add_argument("--maf", type=float, default=0.05)
    p.add_argument("--chunk-size", type=int, default=100000, help="Variants read per chunk")
    p.add_argument("--batch-size", type=int, default=20000, help="map_trans batch size")
    p.add_argument("--lookup-variants", default="", help="Comma-separated variant IDs to report separately")
    p.add_argument("--n-perm", type=int, default=0)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--outdir", default=".")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Genotypes
# ---------------------------------------------------------------------------
def read_pvar_chrom(pvar, chrom):
    """Variants of one chromosome with their row index in the pvar."""
    n_meta = 0
    with open(pvar) as f:
        for line in f:
            if line.startswith("##"):
                n_meta += 1
            else:
                break
    want = {str(chrom), str(chrom).replace("chr", ""), "chr" + str(chrom).replace("chr", "")}
    parts, offset = [], 0
    for chunk in pd.read_csv(pvar, sep="\t", skiprows=n_meta, usecols=[0, 1, 2, 3, 4],
                             dtype={0: str}, chunksize=2_000_000):
        chunk.columns = ["chrom", "pos", "variant_id", "ref", "alt"]
        chunk.index = np.arange(offset, offset + len(chunk))
        offset += len(chunk)
        parts.append(chunk.loc[chunk["chrom"].isin(want)])
    v = pd.concat(parts)
    v["idx"] = v.index
    if len(v) == 0:
        raise ValueError(f"No variants on chromosome {chrom} in {pvar}")
    if not np.all(np.diff(v["idx"].values) == 1):
        raise ValueError(f"Chromosome {chrom} variants are not contiguous in {pvar}")
    return v.reset_index(drop=True), offset


def read_psam(psam):
    df = pd.read_csv(psam, sep="\t", dtype=str)
    col = "IID" if "IID" in df.columns else ("#IID" if "#IID" in df.columns else df.columns[0])
    return df[col].astype(str).tolist()


class Dosages:
    """Chunked dosage reader for a subset of samples (pgenlib)."""

    def __init__(self, prefix, samples, n_variants_total):
        import pgenlib
        psam_ids = read_psam(prefix + ".psam")
        pos = {s: i for i, s in enumerate(psam_ids)}
        missing = [s for s in samples if s not in pos]
        if missing:
            raise ValueError(f"{len(missing)} samples not in the psam, e.g. {missing[:3]}")
        idx = np.array(sorted(pos[s] for s in samples), dtype=np.uint32)
        self.samples = [psam_ids[i] for i in idx]          # pgenlib returns samples in file order
        self.reader = pgenlib.PgenReader(bytes(prefix + ".pgen", "utf8"), raw_sample_ct=len(psam_ids),
                                         variant_ct=n_variants_total, sample_subset=idx)

    def read(self, start, end):
        """Dosages of variants [start, end) as float32 (variants x samples); missing = NaN."""
        arr = np.empty((end - start, len(self.samples)), dtype=np.float32)
        self.reader.read_dosages_range(start, end, arr)
        arr[arr == -9] = np.nan
        return arr


def impute_and_filter(g, maf):
    """Mean-impute missing dosages per variant; return (imputed, keep mask by MAF)."""
    mean = np.nanmean(g, axis=1)
    nan = np.isnan(g)
    if nan.any():
        g = np.where(nan, mean[:, None], g)
    af = mean / 2
    keep = (np.minimum(af, 1 - af) >= maf) & np.isfinite(af)
    return g, keep


# ---------------------------------------------------------------------------
# Phenotypes / tensorQTL
# ---------------------------------------------------------------------------
def load_analyses(analysis_dir):
    out = {}
    for bed in sorted(glob.glob(os.path.join(analysis_dir, "*.phenotype.bed"))):
        name = os.path.basename(bed)[:-len(".phenotype.bed")]
        ph = pd.read_csv(bed, sep="\t", index_col=3).drop(columns=["#chr", "start", "end"])
        cov = pd.read_csv(os.path.join(analysis_dir, f"{name}.covariates.tsv"), sep="\t", index_col=0).T
        cov = cov.loc[ph.columns]
        ph.columns = ph.columns.astype(str)
        cov.index = cov.index.astype(str)
        out[name] = (ph.astype(np.float32), cov.astype(np.float32))
    if not out:
        raise ValueError(f"No *.phenotype.bed in {analysis_dir}")
    return out


def make_map_trans():
    from tensorqtl import trans
    sig = inspect.signature(trans.map_trans).parameters

    def run(genotype_df, phenotype_df, covariates_df, batch_size):
        kw = dict(covariates_df=covariates_df, return_sparse=True, pval_threshold=1.0,
                  maf_threshold=0, batch_size=batch_size)
        if "verbose" in sig:
            kw["verbose"] = False
        res = trans.map_trans(genotype_df, phenotype_df, **{k: v for k, v in kw.items() if k in sig})
        if "af" not in res.columns and "maf" in res.columns:
            res = res.rename(columns={"maf": "af"})
        for c in ["variant_id", "phenotype_id", "pval", "b", "b_se"]:
            if c not in res.columns:
                raise RuntimeError(f"map_trans output has no '{c}' column (columns: {list(res.columns)}); "
                                   "tensorQTL version too old?")
        return res
    return run


def residualise(ph, cov):
    """Residuals of each phenotype (rows) on covariates + intercept (samples are columns)."""
    C = np.column_stack([np.ones(cov.shape[0]), cov.values.astype(np.float64)])
    Y = ph.values.T.astype(np.float64)
    beta, *_ = np.linalg.lstsq(C, Y, rcond=None)
    return pd.DataFrame((Y - C @ beta).T.astype(np.float32), index=ph.index, columns=ph.columns)


def hist_counts(p):
    with np.errstate(divide="ignore"):
        x = -np.log10(np.asarray(p, dtype=np.float64))
    return np.histogram(x, bins=HIST_EDGES)[0]


def hist_table(hists, extra):
    rows = []
    for ct, counts in hists.items():
        nz = np.nonzero(counts)[0]
        rows.append(pd.DataFrame({"celltype": ct, "bin_lo": HIST_EDGES[nz], "bin_hi": HIST_EDGES[nz + 1],
                                  "count": counts[nz]}))
    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["celltype", "bin_lo", "bin_hi", "count"])
    for k, v in extra.items():
        out.insert(0, k, v)
    return out


# ---------------------------------------------------------------------------
# Modes
# ---------------------------------------------------------------------------
def scan(args, variants, dosages, analyses, map_trans, chrom_label):
    lookup = set(v for v in args.lookup_variants.split(",") if v)
    rng = np.random.default_rng(args.seed)
    for d in ["assoc", "stats"]:
        os.makedirs(os.path.join(args.outdir, d), exist_ok=True)
    handles, hists, thin, looked = {}, {}, {n: [] for n in analyses}, {n: [] for n in analyses}
    for name, (ph, _) in analyses.items():
        for ct in ph.index:
            d = os.path.join(args.outdir, "assoc", name, ct)
            os.makedirs(d, exist_ok=True)
            handles[(name, ct)] = gzip.open(os.path.join(d, f"{ct}__chr{chrom_label}.tsv.gz"), "wt", compresslevel=4)
            handles[(name, ct)].write("variant_id\tchrom\tpos\tref\talt\taf\tn_samples\tb\tb_se\tt\tpval\n")
            hists[(name, ct)] = np.zeros(len(HIST_EDGES) - 1, dtype=np.int64)
    col_idx = {name: [dosages.samples.index(s) for s in ph.columns] for name, (ph, _) in analyses.items()}

    start0, n = variants["idx"].iloc[0], len(variants)
    n_tested = {name: 0 for name in analyses}
    for off in range(0, n, args.chunk_size):
        vchunk = variants.iloc[off:off + args.chunk_size]
        g_all = dosages.read(start0 + off, start0 + off + len(vchunk))
        for name, (ph, cov) in analyses.items():
            g, keep = impute_and_filter(g_all[:, col_idx[name]], args.maf)
            if not keep.any():
                continue
            vk = vchunk.loc[keep]
            gdf = pd.DataFrame(g[keep], index=vk["variant_id"].values, columns=ph.columns)
            res = map_trans(gdf, ph, cov, args.batch_size)
            res = res.merge(vk[["variant_id", "chrom", "pos", "ref", "alt"]], on="variant_id", how="left")
            res["n_samples"] = ph.shape[1]
            res["t"] = res["b"] / res["b_se"]
            n_tested[name] += int(keep.sum())
            for ct, r in res.groupby("phenotype_id", sort=False):
                r[["variant_id", "chrom", "pos", "ref", "alt", "af", "n_samples", "b", "b_se", "t", "pval"]].to_csv(
                    handles[(name, ct)], sep="\t", header=False, index=False, float_format="%.6g")
                hists[(name, ct)] += hist_counts(r["pval"])
            sel = (res["pval"] < 1e-3) | (rng.random(len(res)) < 0.01)
            thin[name].append(res.loc[sel, ["phenotype_id", "variant_id", "chrom", "pos", "pval", "b"]])
            if lookup:
                looked[name].append(res.loc[res["variant_id"].isin(lookup)])
        log(f"chr{chrom_label}: {min(off + args.chunk_size, n):,} / {n:,} variants")

    for h in handles.values():
        h.close()
    for name in analyses:
        stem = os.path.join(args.outdir, "stats", f"{name}__chr{chrom_label}")
        hist_table({ct: hists[(name, ct)] for (nm, ct) in hists if nm == name},
                   {"chrom": chrom_label, "analysis": name}).to_csv(stem + ".hist.tsv", sep="\t", index=False)
        t = pd.concat(thin[name]) if thin[name] else pd.DataFrame(columns=["phenotype_id", "variant_id", "chrom", "pos", "pval", "b"])
        t.to_csv(stem + ".thin.tsv.gz", sep="\t", index=False, float_format="%.6g", compression="gzip")
        lk = pd.concat(looked[name]) if looked[name] else pd.DataFrame()
        if len(lk):
            lk.insert(0, "analysis", name)
        lk.to_csv(stem + ".lookup.tsv", sep="\t", index=False, float_format="%.6g")
        log(f"{name}: {n_tested[name]:,} variants tested (MAF >= {args.maf}) x {analyses[name][0].shape[0]} cell types")


def permute(args, variants, dosages, analyses, map_trans, chrom_label):
    os.makedirs(os.path.join(args.outdir, "perm"), exist_ok=True)
    start0 = variants["idx"].iloc[0]
    g_all = dosages.read(start0, start0 + len(variants))
    log(f"chr{chrom_label}: loaded {g_all.shape[0]:,} variants for {args.n_perm} permutations")
    for name, (ph, cov) in analyses.items():
        rng = np.random.default_rng(args.seed)
        g, keep = impute_and_filter(g_all[:, [dosages.samples.index(s) for s in ph.columns]], args.maf)
        gdf = pd.DataFrame(g[keep], index=variants.loc[keep, "variant_id"].values, columns=ph.columns)
        resid = residualise(ph, cov)
        stats_rows, hists = [], {ct: np.zeros(len(HIST_EDGES) - 1, dtype=np.int64) for ct in ph.index}
        for b in range(1, args.n_perm + 1):
            perm = rng.permutation(ph.shape[1])
            ph_perm = pd.DataFrame(resid.values[:, perm], index=ph.index, columns=ph.columns)
            res = map_trans(gdf, ph_perm, cov, args.batch_size)
            for ct, r in res.groupby("phenotype_id", sort=False):
                p = r["pval"].values
                stats_rows.append({"analysis": name, "chrom": chrom_label, "celltype": ct, "perm": b,
                                   "n_tests": len(p),
                                   "lambda": np.median(stats.chi2.isf(p, 1)) / CHI2_MEDIAN,
                                   "min_p": p.min(), "n_p_lt_1e-5": int((p < 1e-5).sum()),
                                   "n_p_lt_5e-8": int((p < 5e-8).sum())})
                hists[ct] += hist_counts(p)
            if b % 10 == 0 or b == args.n_perm:
                log(f"{name}: {b} / {args.n_perm} permutations")
        stem = os.path.join(args.outdir, "perm", f"{name}__chr{chrom_label}")
        pd.DataFrame(stats_rows).to_csv(stem + ".perm_stats.tsv", sep="\t", index=False, float_format="%.6g")
        hist_table(hists, {"n_perm": args.n_perm, "chrom": chrom_label, "analysis": name}).to_csv(
            stem + ".perm_hist.tsv", sep="\t", index=False)


def main():
    args = parse_args()
    import torch
    torch.set_num_threads(args.threads)
    chrom_label = str(args.chrom).replace("chr", "")
    analyses = load_analyses(args.analysis_dir)
    log(f"Analyses: {', '.join(f'{k} ({v[0].shape[0]} cell types x {v[0].shape[1]} samples)' for k, v in analyses.items())}")

    variants, n_total = read_pvar_chrom(args.pgen_prefix + ".pvar", args.chrom)
    log(f"chr{chrom_label}: {len(variants):,} variants in the pvar")
    samples = sorted(set().union(*[set(ph.columns) for ph, _ in analyses.values()]))
    dosages = Dosages(args.pgen_prefix, samples, n_total)
    map_trans = make_map_trans()

    if args.n_perm > 0:
        permute(args, variants, dosages, analyses, map_trans, chrom_label)
    else:
        scan(args, variants, dosages, analyses, map_trans, chrom_label)
    log("Done")


if __name__ == "__main__":
    main()
