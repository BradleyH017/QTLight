#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/eqtl
========================================================================================
    Github : https://github.com/nf-core/eqtl
    Website: https://nf-co.re/eqtl
    Slack  : https://nfcore.slack.com/channels/eqtl
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

// WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { EQTL } from './workflows/eqtl'
include { CCQTL } from './workflows/ccqtl'

//
// WORKFLOW: Run main nf-core/eqtl analysis pipeline
//
workflow NFCORE_EQTL {
    // params.mode: 'eqtl' (default) or 'ccqtl' (cell-type composition QTLs, workflows/ccqtl.nf)
    if (params.mode == 'ccqtl') {
        CCQTL ()
    } else {
        EQTL ()
    }
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_EQTL ()
}

/*
========================================================================================
    THE END
========================================================================================
*/


workflow.onComplete {

    log.info "Pipeline completed at: $workflow.complete"
    log.info "Command line: $workflow.commandLine"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"

    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}