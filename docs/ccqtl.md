# ccQTL: cell-type composition QTL mapping

`params.mode = 'ccqtl'` runs `workflows/ccqtl.nf` instead of the eQTL workflow. Every variant is tested against the abundance of every cell type (genome-wide, like a GWAS), using the same genotype processing as the eQTL workflow.

## Steps

| Process | Script | What it does |
|---|---|---|
| (genotypes) | existing modules | `PREPROCESS_GENOTYPES` (bcftools filters) → `PGEN_CONVERT` → `GENOTYPE_PC_CALCULATION`, exactly as for eQTLs. You can skip the first two with `genotypes.preprocessed_pgen_file`, and the PCA with `covariates.genotype_pcs_file`. |
| `CCQTL_COMPOSITION` | `bin/ccqtl_within_map.py`, `bin/ccqtl_composition.py` | One task per `aggregation_columns` × `cc_within` value. Counts cells per sample, filters samples (`min_cells_per_sample`), cell types (`min_individuals_per_celltype`) and analyses (`min_samples`, ≥ 2 cell types). Applies multiplicative zero replacement and CLR, plus RINT if `RINT_composition`. Covariates are the genotype PCs, `ccQTL.covariates` (categorical ones dummy-coded) and `log10_n_cells`. Matches `Dart_BTNLs/differential_abundance` (identical CLR). |
| `CCQTL_MAP` | `bin/ccqtl_map.py` | One task per chromosome. Reads dosages from the pgen in chunks and runs tensorQTL `map_trans` for every variant (MAF ≥ `maf` in the analysis' samples) against every cell type. Output is **unfiltered**. |
| `CCQTL_PERMUTE` | `bin/ccqtl_map.py --n-perm` | Rescans `permutation_chromosomes` `n_permutations` times. Each permutation shuffles the samples of the covariate-residualised phenotypes, using one shuffle per permutation shared by all cell types of an analysis. |
| `CCQTL_SUMMARY` | `bin/ccqtl_summarise.py` | Produces λ, QQ, Manhattan, lead variants (p < 5e-8, distance-clumped), `lookup_variants`, and the permutation calibration. |

## `cc_within`

Set in the sample input, comma-separated:

- **`none`**: composition across all values of the aggregation column (e.g. all `predicted_labels`).
- **an obs column** (e.g. `predicted_category`): one analysis per value of that column. That analysis uses the composition of the aggregation values within it, i.e. relative proportions within the category only.

  `ccqtl_within_map.py` works out which within value each aggregation value belongs to. It's the value holding ≥ `within_min_frac` of the aggregation value's cells; otherwise the aggregation value is left unassigned and reported. The mapping is saved as `composition/<aggregation>__<within>.within_map.tsv`.

## Outputs (`${outdir}/ccQTL/`)

- `composition/`:
  - `<analysis>.{phenotype.bed,covariates.tsv,composition.tsv.gz}`, where `<analysis>` = `<aggregation>__<within>__<value>`;
  - `*.summary.tsv`, which also lists skipped analyses;
  - `*.within_map.tsv`.
- `assoc/<analysis>/<celltype>/<celltype>__chr<c>.tsv.gz`: every test, with columns `variant_id chrom pos ref alt af n_samples b b_se t pval`.
- `stats/`: per-chromosome `-log10(p)` histograms, thinned points and lookups (summary inputs).
- `perm/`: permutation statistics and histograms.
- `summary/`:
  - `lambda.tsv`
  - `hits.tsv`, with study-wide significance at 5e-8 / n cell types
  - `lookup.tsv`
  - `permutation_calibration.tsv`
  - `plots/<analysis>/{qq.png, qq_permuted_chromosomes.png, manhattan/<celltype>.png}`

## Running

`run_ccQTL.sh` uses `conf/extra_confs/sanger-cd_healthy.config` + `sample_input/input_cd_healthy.nf` + `conf/extra_confs/ccqtl_healthy.config`. All ccQTL settings and their defaults are in the `ccQTL` block of `nextflow.config`.
