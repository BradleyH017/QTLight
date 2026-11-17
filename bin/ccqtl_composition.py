#!/usr/bin/env python
"""
Cell-type composition phenotypes and covariates for ccQTL mapping, following the
differential_abundance pipeline (Dart_BTNLs/differential_abundance/bin/001-compute_composition.py).

For one aggregation column (e.g. predicted_labels) and one 'within' setting:
  within = none : one analysis, composition across all aggregation values
  within = <col>: one analysis per value of <col>, composition across the aggregation
                  values assigned to it (--mapping, from ccqtl_within_map.py)

Per analysis:
  1. Count cells per (sample, aggregation value). Samples kept: genotyped (in the
     psam), with complete covariates.
  2. Drop samples with < --min-cells-per-sample cells (in this analysis), then
     aggregation values with > 0 cells in < --min-individuals-per-celltype samples.
     Analyses with < 2 cell types (CLR needs >= 2 parts) or < --min-samples samples
     left are skipped.
  3. Multiplicative zero replacement (frac x smallest non-zero proportion per cell type).
  4. CLR; optionally rank-inverse-normal transform each cell type across samples.
  5. Covariates: genotype PCs (plink2 eigenvec), --covariates obs columns (numeric as
     is, categorical dummy-coded against the most common level), and optionally
     log10 of the sample's total cells in this analysis. Constant covariates dropped.

Outputs (in --outdir), per analysis <name> = <aggregation>__<within>__<value>:
  <name>.phenotype.bed        tensorQTL-style bed (#chr start end phenotype_id samples...);
                              dummy positions on 'chrCC' (map_trans ignores them)
  <name>.covariates.tsv       covariates x samples
  <name>.composition.tsv.gz   long table: sample, celltype, count, proportion, CLR, phenotype
and <aggregation>__<within>.summary.tsv: one row per analysis (also the skipped ones)
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

from ccqtl_obs_utils import read_obs_columns, sanitize


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--aggregation-column", required=True)
    p.add_argument("--within", default="none", help="'none' or an obs column")
    p.add_argument("--mapping", default=None, help="ccqtl_within_map.py output (required unless --within none)")
    p.add_argument("--sample-column", required=True, help="obs column with the genotype sample ID")
    p.add_argument("--psam", required=True, help="Genotype samples (plink2 .psam)")
    p.add_argument("--eigenvec", required=True, help="Genotype PCs (plink2 .eigenvec)")
    p.add_argument("--n-genotype-pcs", type=int, default=5)
    p.add_argument("--covariates", default="", help="Comma-separated obs columns, constant within a sample")
    p.add_argument("--include-log10-n-cells", action="store_true")
    p.add_argument("--min-cells-per-sample", type=int, default=50)
    p.add_argument("--min-individuals-per-celltype", type=int, default=10)
    p.add_argument("--min-samples", type=int, default=50,
                   help="Skip an analysis with fewer samples than this after the filters")
    p.add_argument("--multiplicative-frac", type=float, default=0.65)
    p.add_argument("--rint", action="store_true", help="Rank-inverse-normal transform each cell type's CLR")
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def multiplicative_replacement(counts, frac):
    props = counts.div(counts.sum(axis=1), axis=0)
    delta = frac * props.where(props > 0).min(axis=0)
    is_zero = props == 0
    delta_mat = is_zero.mul(delta, axis=1)
    row_delta = delta_mat.sum(axis=1)
    if (row_delta >= 1).any():
        raise ValueError("Zero replacement would exceed the sample total; lower multiplicative_frac")
    return props, props.mul(1 - row_delta, axis=0).where(~is_zero, delta_mat)


def clr(mat):
    log_mat = np.log(mat)
    return log_mat.sub(log_mat.mean(axis=1), axis=0)


def rint(x):
    """Rank-based inverse normal transform (ties averaged), per column."""
    n = x.shape[0]
    return x.apply(lambda c: pd.Series(stats.norm.ppf((stats.rankdata(c) - 0.5) / n), index=c.index))


def read_psam_ids(psam):
    df = pd.read_csv(psam, sep="\t", dtype=str)
    col = "IID" if "IID" in df.columns else ("#IID" if "#IID" in df.columns else df.columns[0])
    return set(df[col].astype(str))


def read_eigenvec(path, n_pcs):
    df = pd.read_csv(path, sep=r"\s+", dtype={"#FID": str, "IID": str, "#IID": str})
    id_col = "IID" if "IID" in df.columns else "#IID"
    pcs = [c for c in df.columns if c.startswith("PC")][:n_pcs]
    if len(pcs) < n_pcs:
        raise ValueError(f"{path} has {len(pcs)} PCs, {n_pcs} requested")
    out = df.set_index(df[id_col].astype(str))[pcs]
    out.columns = [f"genotype_{c}" for c in pcs]
    return out


def sample_level(obs, sample_col, cols):
    """One value per sample; error if a covariate varies within a sample."""
    g = obs.groupby(sample_col, observed=True)
    varying = [c for c in cols if (g[c].nunique(dropna=True) > 1).any()]
    if varying:
        raise ValueError(f"Covariates vary within samples: {varying}")
    return g[cols].first()


def encode_covariates(df):
    """Numeric columns as is; categorical ones dummy-coded against the most common level."""
    out = []
    for c in df.columns:
        x = df[c]
        if pd.api.types.is_numeric_dtype(x) and not pd.api.types.is_bool_dtype(x):
            out.append(x.astype(float).rename(c))
        else:
            x = x.astype(str)
            ref = x.value_counts().index[0]
            for lvl in sorted(set(x) - {ref}):
                out.append((x == lvl).astype(float).rename(f"{c}_{sanitize(lvl)}_vs_{sanitize(ref)}"))
    return pd.concat(out, axis=1) if out else pd.DataFrame(index=df.index)


def run_analysis(name, obs, agg_col, sample_col, cov_base, args):
    summary = {"analysis": name, "aggregation_column": agg_col, "within": args.within}
    counts = obs.groupby([sample_col, agg_col], observed=True).size().unstack(fill_value=0)
    counts.index = counts.index.astype(str)
    counts.columns = counts.columns.astype(str)
    # samples: genotyped + complete covariates (cov_base already restricted to those)
    counts = counts.loc[counts.index.isin(cov_base.index)]
    summary["n_samples_input"] = counts.shape[0]
    summary["n_celltypes_input"] = counts.shape[1]

    n_cells = counts.sum(axis=1)
    counts = counts.loc[n_cells >= args.min_cells_per_sample]
    n_present = (counts > 0).sum(axis=0)
    counts = counts.loc[:, n_present >= max(args.min_individuals_per_celltype, 1)]
    counts = counts.loc[counts.sum(axis=1) > 0]
    summary["n_samples"] = counts.shape[0]
    summary["n_celltypes"] = counts.shape[1]
    summary["celltypes"] = ",".join(counts.columns)

    if counts.shape[1] < 2 or counts.shape[0] < args.min_samples:
        summary["status"] = f"skipped: {counts.shape[1]} cell types, {counts.shape[0]} samples after filters"
        print(f"  {name}: {summary['status']}")
        return summary

    props, replaced = multiplicative_replacement(counts, args.multiplicative_frac)
    clr_df = clr(replaced)
    pheno = rint(clr_df) if args.rint else clr_df

    cov = cov_base.loc[counts.index].copy()
    if args.include_log10_n_cells:
        cov["log10_n_cells"] = np.log10(counts.sum(axis=1))
    constant = [c for c in cov.columns if cov[c].nunique() < 2]
    if constant:
        print(f"  {name}: dropping constant covariates {constant}")
        cov = cov.drop(columns=constant)
    summary["covariates"] = ",".join(cov.columns)

    # tensorQTL-style outputs (phenotype IDs = sanitised cell-type names)
    ids = [sanitize(c) for c in pheno.columns]
    if len(set(ids)) != len(ids):
        raise ValueError(f"Cell-type names collide after sanitising in {name}")
    bed = pheno.T.copy()
    bed.index = ids
    bed.insert(0, "phenotype_id", ids)
    bed.insert(0, "end", np.arange(1, len(ids) + 1) + 1)
    bed.insert(0, "start", np.arange(1, len(ids) + 1))
    bed.insert(0, "#chr", "chrCC")
    bed.to_csv(os.path.join(args.outdir, f"{name}.phenotype.bed"), sep="\t", index=False, float_format="%.8g")
    cov.T.to_csv(os.path.join(args.outdir, f"{name}.covariates.tsv"), sep="\t", float_format="%.8g")

    long = pd.concat({
        "count": counts.stack(), "proportion": props.stack(), "proportion_zero_replaced": replaced.stack(),
        "clr": clr_df.stack(), "phenotype": pheno.stack(),
    }, axis=1).reset_index()
    long.columns = ["sample", "celltype"] + list(long.columns[2:])
    long["phenotype_id"] = long["celltype"].map(dict(zip(pheno.columns, ids)))
    long.to_csv(os.path.join(args.outdir, f"{name}.composition.tsv.gz"), sep="\t", index=False,
                float_format="%.8g", compression="gzip")

    summary["status"] = "ok"
    print(f"  {name}: {counts.shape[0]} samples x {counts.shape[1]} cell types; covariates: {summary['covariates']}")
    return summary


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    cov_cols = [c.strip() for c in args.covariates.split(",") if c.strip()]
    within = args.within if args.within and args.within.lower() != "none" else None
    obs = read_obs_columns(args.h5ad, [args.aggregation_column, args.sample_column] + cov_cols + ([within] if within else []))
    obs = obs.dropna(subset=[args.aggregation_column, args.sample_column])
    obs[args.sample_column] = obs[args.sample_column].astype(str)
    obs[args.aggregation_column] = obs[args.aggregation_column].astype(str)
    print(f"{len(obs):,} cells, {obs[args.sample_column].nunique()} samples, "
          f"{obs[args.aggregation_column].nunique()} {args.aggregation_column} values")

    # Sample-level covariates: genotype PCs + obs covariates, genotyped samples with complete data
    genotyped = read_psam_ids(args.psam)
    pcs = read_eigenvec(args.eigenvec, args.n_genotype_pcs)
    samples = sorted(set(obs[args.sample_column]))
    cov = pcs.reindex(samples)
    if cov_cols:
        cov = cov.join(encode_covariates(sample_level(obs, args.sample_column, cov_cols)).reindex(samples))
    keep = cov.notna().all(axis=1) & cov.index.isin(genotyped)
    print(f"Samples: {len(samples)} in h5ad; {int(cov.index.isin(genotyped).sum())} genotyped; "
          f"{int(keep.sum())} genotyped with complete covariates")
    cov = cov.loc[keep]

    agg = sanitize(args.aggregation_column)
    summaries = []
    if within is None:
        summaries.append(run_analysis(f"{agg}__none__all", obs, args.aggregation_column, args.sample_column, cov, args))
    else:
        if not args.mapping:
            raise ValueError("--mapping is required with --within")
        mapping = pd.read_csv(args.mapping, sep="\t", dtype=str)
        mapping = mapping.dropna(subset=["within_value"])
        for w, d in mapping.groupby("within_value"):
            sub = obs.loc[obs[args.aggregation_column].isin(set(d["aggregation_value"]))]
            name = f"{agg}__{sanitize(within)}__{sanitize(w)}"
            summaries.append(run_analysis(name, sub, args.aggregation_column, args.sample_column, cov, args))

    pd.DataFrame(summaries).to_csv(os.path.join(args.outdir, f"{agg}__{sanitize(within or 'none')}.summary.tsv"),
                                   sep="\t", index=False)


if __name__ == "__main__":
    main()
