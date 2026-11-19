/*
========================================================================================
    ccQTL: cell-type composition QTL mapping (params.mode = 'ccqtl')
========================================================================================
    Genotypes go through the same steps as the eQTL workflow (bcftools filters ->
    plink2 pgen -> genotype PCs), or are taken from genotypes.preprocessed_pgen_file /
    covariates.genotype_pcs_file. Then, for every aggregation column x cc_within value:
      composition per sample (CLR, optionally RINT) -> tensorQTL map_trans of every
      variant (MAF >= params.maf) against every cell type, one job per chromosome,
      unfiltered results -> permutation calibration on ccQTL.permutation_chromosomes ->
      summary.
    Results: ${params.outdir}/ccQTL/{composition,assoc,stats,perm,summary}
----------------------------------------------------------------------------------------
*/

include {PREPROCESS_GENOTYPES} from '../modules/nf-core/modules/preprocess_genotypes/main'
include {PGEN_CONVERT} from '../modules/nf-core/modules/plink_convert/main'
include {GENOTYPE_PC_CALCULATION} from '../modules/nf-core/modules/genotype_pc_calculation/main'
include {CCQTL_COMPOSITION; CCQTL_MAP; CCQTL_PERMUTE; CCQTL_SUMMARY} from '../modules/nf-core/modules/ccqtl/main'


def as_list(x) {
    if (x == null || x.toString().trim() == '') { return [] }
    return (x instanceof List ? x : x.toString().split(',')).collect { it.toString().trim() }.findAll { it }
}


workflow CCQTL {
    log.info '------ ccQTL: cell-type composition QTL mapping ------'

    if (!params.genotypes.use_gt_dosage) {
        error "ccQTL mode reads dosages from the plink2 pgen: set genotypes.use_gt_dosage = true"
    }
    if (params.genotypes.subset_genotypes_to_available) {
        log.warn "genotypes.subset_genotypes_to_available is ignored in ccQTL mode (samples are matched to the psam)"
    }

    aggregation_columns = as_list(params.aggregation_columns)
    within = as_list(params.cc_within ?: 'none')
    chromosomes = as_list(params.ccQTL.chromosomes)
    perm_chromosomes = (params.ccQTL.n_permutations as int) > 0 ? as_list(params.ccQTL.permutation_chromosomes) : []
    log.info "aggregation_columns : ${aggregation_columns.join(', ')}"
    log.info "cc_within           : ${within.join(', ')}"
    log.info "covariates          : genotype PCs (${params.covariates.nr_genotype_pcs}), ${params.ccQTL.covariates ?: '-'}${params.ccQTL.include_log10_n_cells ? ', log10_n_cells' : ''}"
    log.info "RINT_composition    : ${params.ccQTL.RINT_composition}"
    log.info "chromosomes         : ${chromosomes.join(',')}; permutations: ${params.ccQTL.n_permutations} on ${perm_chromosomes.join(',') ?: '-'} in ${params.ccQTL.permutation_jobs ?: 1} jobs per chromosome"

    // ---------------- Genotypes (as in the eQTL workflow)
    if (params.genotypes.preprocessed_pgen_file == '') {
        donorsvcf = Channel.from(params.input_vcf)
        if (params.genotypes.apply_bcftools_filters) {
            PREPROCESS_GENOTYPES(donorsvcf)
            plink_convert_input = PREPROCESS_GENOTYPES.out.filtered_vcf
        } else {
            plink_convert_input = donorsvcf
        }
        PGEN_CONVERT(plink_convert_input)
        plink_path = PGEN_CONVERT.out.plink_path.first()
    } else {
        plink_path = Channel.fromPath(params.genotypes.preprocessed_pgen_file, type: 'dir', checkIfExists: true).first()
    }

    if (params.covariates.genotype_pcs_file == '') {
        GENOTYPE_PC_CALCULATION(plink_path)
        genotype_pcs_file = GENOTYPE_PC_CALCULATION.out.gtpca_plink.first()
    } else {
        genotype_pcs_file = Channel.fromPath(params.covariates.genotype_pcs_file, checkIfExists: true).first()
    }

    // ---------------- Composition phenotypes: aggregation column x within
    pairs = Channel.fromList([aggregation_columns, within].combinations())
    h5ad = Channel.fromPath(params.phenotype_file, checkIfExists: true).first()
    CCQTL_COMPOSITION(pairs, h5ad, genotype_pcs_file, plink_path)

    analysis_files = CCQTL_COMPOSITION.out.beds.flatten()
        .mix(CCQTL_COMPOSITION.out.covs.flatten())
        .collect()

    // ---------------- Association (one job per chromosome) and permutations
    CCQTL_MAP(analysis_files, plink_path, chromosomes)
    if (perm_chromosomes) {
        // Split the permutations of each chromosome into ccQTL.permutation_jobs blocks of consecutive
        // permutation numbers (each permutation has its own seed, so the split does not change the results)
        n_perm = params.ccQTL.n_permutations as int
        n_jobs = Math.max(1, Math.min((params.ccQTL.permutation_jobs ?: 1) as int, n_perm))
        block = (int) Math.ceil(n_perm / (double) n_jobs)
        perm_jobs = []
        perm_chromosomes.each { c ->
            for (int s = 1; s <= n_perm; s += block) perm_jobs << [c, s, Math.min(block, n_perm - s + 1)]
        }
        CCQTL_PERMUTE(analysis_files, plink_path, perm_jobs)
        perm_files = CCQTL_PERMUTE.out.perm.flatten().collect()
    } else {
        perm_files = Channel.value([])
    }

    CCQTL_SUMMARY(CCQTL_MAP.out.stats.flatten().collect(), perm_files)
}
