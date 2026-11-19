params{
    method= 'single_cell' //or a [bulk | single_cell] (if single cell used the *phenotype_file* is a h5ad file)
    input_vcf ='/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/core_analysis_output/IBDverse_multi-tissue_eQTL_project/IBDverse_genotypes/2024_07_11-genotype_plate12345/imputed.vcf.gz' // Batch 5 genotypes
    genotype_phenotype_mapping_file = '' // annotation file containing Genotype | RNA and | Condiion
    annotation_file = '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/IBDVerse-sc-eQTL-code/data/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt' //assets file that has start and end positions for the genes, this one is using hg38
    phenotype_file = '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/core_analysis_output/IBDverse_multi-tissue_eQTL_project/IBDverse_scRNA/celltypist_0.5_ngene_ncount_mt_filt_nomiss.h5ad' // ALL
    aggregation_columns='predicted_labels,predicted_category' 
    cc_within='none,predicted_category' // ccQTL only: 'none' = composition across all aggregation_columns values; an obs column = composition within each of its values (comma-separated)
    gt_id_column='Genotyping_ID' //for the scrna h5ad defines which column has individual level id (corresponding to the VCF file)
    sample_column='sanger_sample_id' // for the scrna h5ad defines which column has the sample name (can be identical to gt_id_column if sample=individual)
    sample_covariates='' //covariates to be included in the model - LEAVE BLANK, DOES NOTHING AT PRESENT
}