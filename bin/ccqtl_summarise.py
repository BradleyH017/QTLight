#!/usr/bin/env python
"""
Summarise a ccQTL run from the per-chromosome outputs of ccqtl_map.py.

Inputs (--stats-dir / --perm-dir): <analysis>__chr<c>.{hist.tsv,thin.tsv.gz,lookup.tsv}
and <analysis>__chr<c>[__perm<a>-<b>].{perm_stats.tsv,perm_hist.tsv} (split permutation jobs are pooled).

Outputs (--outdir):
  lambda.tsv                     genome-wide lambda + n tests per analysis x cell type
  hits.tsv                       lead variants (p < 5e-8, distance-clumped per cell type,
                                 --clump-kb), flagged if also < 5e-8 / n cell types in the analysis
  lookup.tsv                     results for the --lookup-variants
  permutation_calibration.tsv    per analysis x cell type on the permuted chromosome(s): observed
                                 vs permutation lambda, null rate of p < 1e-5, family-wise rates
  plots/<analysis>/qq.png        genome-wide QQ, one line per cell type
  plots/<analysis>/qq_permuted_chromosomes.png   observed vs permutation null, permuted chromosomes
  plots/<analysis>/manhattan/<celltype>.png
"""
import argparse
import glob
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

GWS = 5e-8
CHI2_MEDIAN = stats.chi2.ppf(0.5, 1)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stats-dir", required=True)
    p.add_argument("--perm-dir", default="")
    p.add_argument("--clump-kb", type=int, default=1000)
    p.add_argument("--outdir", default="summary")
    return p.parse_args()


def read_many(pattern, **kw):
    dfs = []
    for f in sorted(glob.glob(pattern)):
        try:
            d = pd.read_csv(f, sep="\t", **kw)
        except pd.errors.EmptyDataError:     # e.g. lookup files on chromosomes without lookup variants
            continue
        if len(d):
            dfs.append(d)
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def hist_lambda(h):
    """Lambda from a -log10(p) histogram (median interpolated linearly within its bin)."""
    h = h.groupby(["bin_lo", "bin_hi"], as_index=False)["count"].sum().sort_values("bin_lo")
    n = h["count"].sum()
    cum = h["count"].cumsum().values
    i = int(np.searchsorted(cum, n / 2))
    lo, hi, c = h["bin_lo"].iloc[i], h["bin_hi"].iloc[i], h["count"].iloc[i]
    x = lo + (n / 2 - (cum[i] - c)) / c * (hi - lo)
    return stats.chi2.isf(10 ** -x, 1) / CHI2_MEDIAN, n


def hist_qq(h):
    """QQ points from a histogram: for each bin, observed = bin_lo, expected = -log10(tail fraction)."""
    h = h.groupby(["bin_lo"], as_index=False)["count"].sum().sort_values("bin_lo", ascending=False)
    n = h["count"].sum()
    tail = h["count"].cumsum().values
    return -np.log10(tail / n), h["bin_lo"].values


def chrom_key(c):
    c = str(c).replace("chr", "")
    return int(c) if c.isdigit() else {"X": 23, "Y": 24, "XY": 25, "MT": 26}.get(c, 99)


def clump(d, kb):
    """Greedy distance clumping: lead = smallest p; drop others within kb on the same chromosome."""
    d = d.sort_values("pval")
    leads = []
    taken = {}
    for _, r in d.iterrows():
        pos = taken.setdefault(r["chrom"], [])
        if all(abs(r["pos"] - p) > kb * 1000 for p in pos):
            leads.append(r)
            pos.append(r["pos"])
    return pd.DataFrame(leads)


def manhattan(t, title, path):
    t = t.copy()
    t["ck"] = t["chrom"].map(chrom_key)
    t = t.sort_values(["ck", "pos"])
    offsets, ticks, x0 = {}, [], 0
    for c, d in t.groupby("ck", sort=True):
        offsets[c] = x0
        ticks.append((x0 + d["pos"].max() / 2, str(d["chrom"].iloc[0])))
        x0 += d["pos"].max() + 5e6
    x = t["pos"] + t["ck"].map(offsets)
    fig, ax = plt.subplots(figsize=(11, 3.4))
    ax.scatter(x, -np.log10(t["pval"].clip(lower=1e-300)), s=2,
               c=np.where(t["ck"] % 2 == 0, "#3b6ea5", "#9aa9b8"), rasterized=True)
    ax.axhline(-np.log10(GWS), color="firebrick", lw=0.8, ls="--")
    ax.set_xticks([a for a, _ in ticks])
    ax.set_xticklabels([b for _, b in ticks], fontsize=6)
    ax.set_ylabel("-log10(p)")
    ax.set_title(title, fontsize=9)
    ax.set_xlim(0, x0)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    hist = read_many(os.path.join(args.stats_dir, "*.hist.tsv"))
    lookup = read_many(os.path.join(args.stats_dir, "*.lookup.tsv"))
    if hist.empty:
        raise ValueError(f"No *.hist.tsv in {args.stats_dir}")
    # thin files carry no analysis column: take it from the file name
    thin_parts = []
    for f in sorted(glob.glob(os.path.join(args.stats_dir, "*.thin.tsv.gz"))):
        d = pd.read_csv(f, sep="\t")
        if len(d):
            d["analysis"] = re.sub(r"__chr[^_]+\.thin\.tsv\.gz$", "", os.path.basename(f))
            thin_parts.append(d)
    thin = pd.concat(thin_parts, ignore_index=True) if thin_parts else pd.DataFrame()

    # ---------------- lambda, QQ
    lam_rows = []
    for (an, ct), h in hist.groupby(["analysis", "celltype"]):
        lam, n = hist_lambda(h)
        lam_rows.append({"analysis": an, "celltype": ct, "n_tests": n, "lambda": round(lam, 4)})
    lam = pd.DataFrame(lam_rows)
    n_ct = lam.groupby("analysis")["celltype"].nunique()
    lam["n_celltypes_in_analysis"] = lam["analysis"].map(n_ct)
    lam.to_csv(os.path.join(args.outdir, "lambda.tsv"), sep="\t", index=False)
    print(lam.groupby("analysis")["lambda"].describe()[["count", "min", "50%", "max"]])

    for an, ha in hist.groupby("analysis"):
        pdir = os.path.join(args.outdir, "plots", an)
        os.makedirs(os.path.join(pdir, "manhattan"), exist_ok=True)
        fig, ax = plt.subplots(figsize=(5, 5))
        for ct, h in ha.groupby("celltype"):
            e, o = hist_qq(h)
            ax.plot(e, o, lw=0.8, alpha=0.7)
        lim = max(ax.get_xlim()[1], 8)
        ax.plot([0, lim], [0, lim], color="grey", ls="--", lw=0.8)
        ax.axhline(-np.log10(GWS), color="firebrick", lw=0.6, ls=":")
        ax.set_xlabel("Expected -log10(p)")
        ax.set_ylabel("Observed -log10(p)")
        l = lam.loc[lam["analysis"] == an, "lambda"]
        ax.set_title(f"{an}\n{len(l)} cell types; lambda median {l.median():.3f} (range {l.min():.3f}-{l.max():.3f})", fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(pdir, "qq.png"), dpi=200)
        plt.close(fig)

    # ---------------- Manhattan + hits
    hits = []
    if not thin.empty:
        for (an, ct), t in thin.groupby(["analysis", "phenotype_id"]):
            manhattan(t, f"{an}: {ct}", os.path.join(args.outdir, "plots", an, "manhattan", f"{ct}.png"))
            sig = t.loc[t["pval"] < GWS]
            if len(sig):
                leads = clump(sig, args.clump_kb)
                leads["n_sig_variants_in_celltype"] = len(sig)
                hits.append(leads)
    hits = pd.concat(hits, ignore_index=True) if hits else pd.DataFrame(
        columns=["analysis", "phenotype_id", "variant_id", "chrom", "pos", "pval", "b"])
    if len(hits):
        hits["study_wide_threshold"] = GWS / hits["analysis"].map(n_ct)
        hits["study_wide_sig"] = hits["pval"] < hits["study_wide_threshold"]
        hits = hits.sort_values("pval")
    hits.to_csv(os.path.join(args.outdir, "hits.tsv"), sep="\t", index=False)
    print(f"Lead variants at p < {GWS}: {len(hits)} "
          f"({int(hits['study_wide_sig'].sum()) if len(hits) else 0} also < 5e-8 / n cell types)")
    lookup.to_csv(os.path.join(args.outdir, "lookup.tsv"), sep="\t", index=False)

    # ---------------- Permutation calibration
    if args.perm_dir:
        ps = read_many(os.path.join(args.perm_dir, "*.perm_stats.tsv"))
        ph = read_many(os.path.join(args.perm_dir, "*.perm_hist.tsv"))
        if not ps.empty:
            perm_chroms = sorted(ps["chrom"].astype(str).unique(), key=chrom_key)
            obs_h = hist.loc[hist["chrom"].astype(str).isin(perm_chroms)]
            rows = []
            for (an, ct), d in ps.groupby(["analysis", "celltype"]):
                d = d.groupby("perm").agg(n_tests=("n_tests", "sum"), min_p=("min_p", "min"),
                                          n5=("n_p_lt_1e-5", "sum"), n8=("n_p_lt_5e-8", "sum"),
                                          lam=("lambda", "median"))
                oh = obs_h.loc[(obs_h["analysis"] == an) & (obs_h["celltype"] == ct)]
                obs_lam = hist_lambda(oh)[0] if len(oh) else np.nan
                rows.append({"analysis": an, "celltype": ct, "chromosomes": ",".join(perm_chroms),
                             "n_perm": len(d), "n_tests_per_perm": int(d["n_tests"].median()),
                             "observed_lambda": round(obs_lam, 4),
                             "perm_lambda_median": round(d["lam"].median(), 4),
                             "perm_lambda_2.5%": round(d["lam"].quantile(0.025), 4),
                             "perm_lambda_97.5%": round(d["lam"].quantile(0.975), 4),
                             "observed_lambda_percentile": round((d["lam"] <= obs_lam).mean(), 3),
                             "null_rate_p_lt_1e-5": d["n5"].sum() / d["n_tests"].sum(),
                             "perms_with_p_lt_5e-8": (d["n8"] > 0).mean()})
            cal = pd.DataFrame(rows)
            fw = ps.groupby(["analysis", "perm"])["min_p"].min().reset_index()
            fw = fw.groupby("analysis")["min_p"].apply(lambda x: (x < GWS).mean()).rename("perms_any_celltype_p_lt_5e-8")
            cal = cal.merge(fw, left_on="analysis", right_index=True, how="left")
            cal.to_csv(os.path.join(args.outdir, "permutation_calibration.tsv"), sep="\t", index=False)
            print(cal.groupby("analysis")[["observed_lambda", "perm_lambda_median", "null_rate_p_lt_1e-5",
                                           "perms_any_celltype_p_lt_5e-8"]].median())

            for an in cal["analysis"].unique():
                fig, ax = plt.subplots(figsize=(5, 5))
                e, o = hist_qq(ph.loc[ph["analysis"] == an])
                ax.plot(e, o, color="grey", lw=2, label="permuted (pooled)")
                for ct, h in obs_h.loc[obs_h["analysis"] == an].groupby("celltype"):
                    e, o = hist_qq(h)
                    ax.plot(e, o, lw=0.7, alpha=0.6)
                lim = max(ax.get_xlim()[1], 6)
                ax.plot([0, lim], [0, lim], color="black", ls="--", lw=0.8)
                ax.set_xlabel("Expected -log10(p)")
                ax.set_ylabel("Observed -log10(p)")
                c = cal.loc[cal["analysis"] == an]
                ax.set_title(f"{an}, chr {','.join(perm_chroms)}\nobserved (coloured, per cell type) vs "
                             f"{int(c['n_perm'].iloc[0])} permutations (grey)\nperm lambda median "
                             f"{c['perm_lambda_median'].median():.3f}", fontsize=8)
                ax.legend(fontsize=7, loc="upper left")
                fig.tight_layout()
                fig.savefig(os.path.join(args.outdir, "plots", an, "qq_permuted_chromosomes.png"), dpi=200)
                plt.close(fig)


if __name__ == "__main__":
    main()
