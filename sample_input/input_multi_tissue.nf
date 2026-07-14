params{
    method= 'single_cell' //or a [bulk | single_cell] (if single cell used the *phenotype_file* is a h5ad file)
    input_vcf ='/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/core_analysis_output/IBDverse_multi-tissue_eQTL_project/IBDverse_genotypes/2024_07_11-genotype_plate12345/imputed.vcf.gz' // Batch 5 genotypes
    genotype_phenotype_mapping_file = '' // annotation file containing Genotype | RNA and | Condiion
    annotation_file = '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/IBDVerse-sc-eQTL-code/data/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt' //assets file that has start and end positions for the genes, this one is using hg38
    phenotype_file = '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/atlassing/results/objects/from_irods/celltypist_0.5_ngene_ncount_mt_filt_nomiss.h5ad' //this should point to h5ad file in a single cell experiments.
    aggregation_columns='unannotated_tissue,predicted_category_tissue,predicted_labels_tissue' // For interaction
    // aggregation_columns='unannotated,predicted_category,predicted_labels,unannotated_tissue,predicted_category_tissue,predicted_labels_tissue' // For base
    // aggregation_columns='unannotated,predicted_category' // For testing
    gt_id_column='Genotyping_ID' //for the scrna h5ad defines which column has individual level id (corresponding to the VCF file)
    sample_column='Genotyping_ID' // for the scrna h5ad defines which column has the sample name (can be identical to gt_id_column if sample=individual)
    sample_covariates='' //covariates to be included in the model - LEAVE BLANK, DOES NOTHING AT PRESENT
}