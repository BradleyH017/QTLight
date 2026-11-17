#!/usr/bin/env python
"""
Work out how an aggregation column (e.g. predicted_labels) sits within a broader
'within' column (e.g. predicted_category), for ccQTL testing within each value of
the 'within' column.

For every aggregation value, the cells are cross-tabulated against the within
column. Each aggregation value is assigned to the within value holding most of
its cells, if that fraction is >= --min-frac; otherwise it is left unassigned
(and excluded from the within analyses). The whole aggregation value (all its
cells) then counts towards that within group, so a label's abundance is never
split across groups.

Output (--output, TSV), one row per aggregation value:
    aggregation_value, within_value (assigned, or NA), n_cells, n_cells_in_assigned,
    frac_in_assigned, n_within_values (with >0 cells), all_within (value:n;...), status
"""
import argparse
import sys

import pandas as pd

from ccqtl_obs_utils import read_obs_columns


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--aggregation-column", required=True)
    p.add_argument("--within-column", required=True)
    p.add_argument("--min-frac", type=float, default=0.9,
                   help="Min fraction of an aggregation value's cells in one within value to assign it [%(default)s]")
    p.add_argument("--output", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    obs = read_obs_columns(args.h5ad, [args.aggregation_column, args.within_column])
    obs = obs.dropna()
    agg = obs[args.aggregation_column].astype(str)
    within = obs[args.within_column].astype(str)

    ct = pd.crosstab(agg, within)
    n_cells = ct.sum(axis=1)
    best = ct.idxmax(axis=1)
    n_best = ct.max(axis=1)
    frac = n_best / n_cells
    assigned = frac >= args.min_frac

    out = pd.DataFrame({
        "aggregation_value": ct.index,
        "within_value": best.where(assigned).values,
        "n_cells": n_cells.values,
        "n_cells_in_assigned": n_best.values,
        "frac_in_assigned": frac.round(4).values,
        "n_within_values": (ct > 0).sum(axis=1).values,
        "all_within": [";".join(f"{w}:{n}" for w, n in row[row > 0].sort_values(ascending=False).items())
                       for _, row in ct.iterrows()],
    })
    out["status"] = "unassigned"
    out.loc[assigned.values & (out["n_within_values"] == 1).values, "status"] = "nested"
    out.loc[assigned.values & (out["n_within_values"] > 1).values, "status"] = "majority"
    out = out.sort_values(["within_value", "aggregation_value"], na_position="last")
    out.to_csv(args.output, sep="\t", index=False)

    n = out["status"].value_counts().to_dict()
    print(f"{args.aggregation_column} within {args.within_column}: {len(out)} values; "
          f"nested {n.get('nested', 0)}, majority (>= {args.min_frac}) {n.get('majority', 0)}, "
          f"unassigned {n.get('unassigned', 0)}")
    for w, d in out.dropna(subset=["within_value"]).groupby("within_value"):
        print(f"  {w}: {len(d)} values")
    if n.get("unassigned", 0):
        print("  Unassigned (excluded from the within analyses):", file=sys.stderr)
        print(out.loc[out["status"] == "unassigned", ["aggregation_value", "all_within"]].to_string(index=False),
              file=sys.stderr)


if __name__ == "__main__":
    main()
