params{
    method= 'single_cell' //or a [bulk | single_cell] (if single cell used the *phenotype_file* is a h5ad file)
    input_vcf ='/lustre/scratch126/humgen/teams_v2/anderson/rh33/IBDverse_eqtl_fr002/data/2026_03_22-genotypes_plate123456/imputed-renamed.vcf.gz' // Batch 1-6 genotypes
    genotype_phenotype_mapping_file = '' // annotation file containing Genotype | RNA and | Condiion
    annotation_file = '/lustre/scratch126/humgen/teams_v2/anderson/rh33/IBDverse_eqtl_fr002/IBDVerse-sc-eQTL-code/data/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt' //assets file that has start and end positions for the genes, this one is using hg38
    //phenotype_file = '/lustre/scratch126/humgen/teams_v2/anderson/rh33/IBDverse_eqtl_fr002/celltypist_fig/adata_celltypist_annotated_filtered.h5ad' //this should point to h5ad file in a single cell experiments.
    phenotype_file = '/lustre/scratch126/humgen/teams_v2/anderson/rh33/IBDverse_eqtl_fr002/celltypist_fig/adata_celltypist_annotated_filtered_3-stripped_layers.h5ad' //this should point to h5ad file in a single cell experiments.
    aggregation_columns='leiden_tissue,leiden,Category_tissue,Category,unannotated_tissue,unannotated'
    gt_id_column='Batch123456_Genotyping_ID' //for the scrna h5ad defines which column has individual level id (corresponding to the VCF file)
    sample_column='Batch123456_Genotyping_ID' // for the scrna h5ad defines which column has the sample name (can be identical to gt_id_column if sample=individual)
    sample_covariates='' //covariates to be included in the model - LEAVE BLANK, DOES NOTHING AT PRESENT
}