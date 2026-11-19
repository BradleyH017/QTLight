// ccQTL: cell-type composition QTL mapping (see workflows/ccqtl.nf)

process CCQTL_COMPOSITION {
    // Composition phenotypes (CLR, optionally RINT) and covariates for one
    // aggregation column x 'within' setting; one analysis per within value.
    tag "${aggregation_column} within ${within}"
    label 'process_medium'

    publishDir  path: "${params.outdir}/ccQTL",
                mode: "${params.copy_mode}",
                overwrite: "true"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        tuple val(aggregation_column), val(within)
        path(h5ad)
        path(eigenvec)
        path(plink_dir)

    output:
        path("composition/*.phenotype.bed"), emit: beds, optional: true
        path("composition/*.covariates.tsv"), emit: covs, optional: true
        path("composition/*.composition.tsv.gz"), optional: true
        path("composition/*.within_map.tsv"), optional: true
        path("composition/*.summary.tsv"), emit: summary

    script:
        within_map = "composition/${aggregation_column}__${within}.within_map.tsv"
        within_cmd = within == 'none' ? '' : """
        ccqtl_within_map.py \\
            --h5ad ${h5ad} \\
            --aggregation-column "${aggregation_column}" \\
            --within-column "${within}" \\
            --min-frac ${params.ccQTL.within_min_frac} \\
            --output "${within_map}"
        """
        mapping_arg = within == 'none' ? '' : "--mapping \"${within_map}\""
        n_cells_arg = params.ccQTL.include_log10_n_cells ? '--include-log10-n-cells' : ''
        rint_arg = params.ccQTL.RINT_composition ? '--rint' : ''
        """
        mkdir -p composition
        ${within_cmd}
        ccqtl_composition.py \\
            --h5ad ${h5ad} \\
            --aggregation-column "${aggregation_column}" \\
            --within "${within}" ${mapping_arg} \\
            --sample-column "${params.gt_id_column}" \\
            --psam ${plink_dir}/plink_genotypes.psam \\
            --eigenvec ${eigenvec} \\
            --n-genotype-pcs ${params.covariates.nr_genotype_pcs} \\
            --covariates "${params.ccQTL.covariates}" ${n_cells_arg} \\
            --min-cells-per-sample ${params.ccQTL.min_cells_per_sample} \\
            --min-individuals-per-celltype ${params.ccQTL.min_individuals_per_celltype} \\
            --min-samples ${params.ccQTL.min_samples} \\
            --multiplicative-frac ${params.ccQTL.multiplicative_frac} ${rint_arg} \\
            --outdir composition
        """
}


process CCQTL_MAP {
    // All variants on one chromosome x every cell type of every analysis (unfiltered output)
    tag "chr${chrom}"
    label 'process_medium'

    publishDir  path: "${params.outdir}/ccQTL",
                mode: "${params.copy_mode}",
                overwrite: "true"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        path(analysis_files)
        path(plink_dir)
        each chrom

    output:
        path("assoc/**.tsv.gz"), emit: assoc
        path("stats/*"), emit: stats

    script:
        """
        ccqtl_map.py \\
            --pgen-prefix ${plink_dir}/plink_genotypes \\
            --chrom ${chrom} \\
            --analysis-dir . \\
            --maf ${params.maf} \\
            --lookup-variants "${params.ccQTL.lookup_variants}" \\
            --threads ${task.cpus} \\
            --outdir .
        """
}


process CCQTL_PERMUTE {
    // Calibration: rescan one chromosome with permuted (residualised) phenotypes.
    // One job per chromosome x block of permutations (ccQTL.permutation_jobs blocks)
    tag "chr${perm_job[0]}, permutations ${perm_job[1]}-${perm_job[1] + perm_job[2] - 1}"
    label 'process_medium'
    label 'process_long'

    publishDir  path: "${params.outdir}/ccQTL",
                mode: "${params.copy_mode}",
                overwrite: "true"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        path(analysis_files)
        path(plink_dir)
        each perm_job       // [chromosome, first permutation, number of permutations]

    output:
        path("perm/*"), emit: perm

    script:
        """
        ccqtl_map.py \\
            --pgen-prefix ${plink_dir}/plink_genotypes \\
            --chrom ${perm_job[0]} \\
            --analysis-dir . \\
            --maf ${params.maf} \\
            --n-perm ${perm_job[2]} \\
            --perm-start ${perm_job[1]} \\
            --seed ${params.ccQTL.permutation_seed} \\
            --threads ${task.cpus} \\
            --outdir .
        """
}


process CCQTL_SUMMARY {
    // lambda, QQ, Manhattan, lead variants, lookups and permutation calibration
    label 'process_low'

    publishDir  path: "${params.outdir}/ccQTL",
                mode: "${params.copy_mode}",
                overwrite: "true"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        path(stats_files, stageAs: "stats/*")
        path(perm_files, stageAs: "perm/*")

    output:
        path("summary"), emit: summary

    script:
        """
        ccqtl_summarise.py \\
            --stats-dir stats \\
            --perm-dir perm \\
            --clump-kb ${params.ccQTL.clump_kb} \\
            --outdir summary
        """
}
