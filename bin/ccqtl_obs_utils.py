"""Shared helpers for the ccQTL (cell-composition QTL) scripts."""
import re

import h5py
import numpy as np
import pandas as pd

# read_elem has moved between anndata versions
try:
    from anndata.io import read_elem  # anndata >= 0.11
except ImportError:
    try:
        from anndata.experimental import read_elem  # anndata 0.8-0.10
    except ImportError:
        from anndata._io.specs import read_elem  # older


def sanitize(name):
    """Make a value safe for file names and phenotype IDs."""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(name))


def read_obs_columns(h5ad, columns):
    """Read only the requested obs columns (and the index) of an h5ad, without
    loading X or the other obs columns."""
    columns = list(dict.fromkeys(c for c in columns if c))
    with h5py.File(h5ad, "r") as f:
        obs = f["obs"]
        missing = [c for c in columns if c not in obs]
        if missing:
            raise KeyError(f"Columns not in obs of {h5ad}: {missing}. Available: {sorted(obs.keys())}")
        index_key = obs.attrs.get("_index", "_index")
        index = read_elem(obs[index_key]) if index_key in obs else None
        data = {}
        for c in columns:
            x = read_elem(obs[c])
            data[c] = pd.Series(np.asarray(x) if not isinstance(x, pd.Categorical) else x)
    out = pd.DataFrame(data)
    if index is not None:
        out.index = pd.Index(np.asarray(index).astype(str))
    return out
