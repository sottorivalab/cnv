#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sottorivalab/cnv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/sottorivalab/cnv
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CNV                     } from './workflows/cnv'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_cnv_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_cnv_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_cnv_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta     = getGenomeAttribute('fasta')
params.fasta_fai = getGenomeAttribute('fasta_fai')
params.gc_wiggle = getGenomeAttribute('gc_wiggle')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow SOTTORIVALAB_CNV {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:
    versions = Channel.empty()
    //
    // SUBWORKFLOW: Prepare genome if needed
    //
    // gather files or get them from params
    // Gather gc_wiggle file: either generate or use existing
    // TODO create warning if neither gc_wiggle is available nor create_gc_wiggle is supplied

    fasta_ch = params.fasta ?
        Channel.fromPath(params.fasta)
            .map{ it -> [ [id:it.baseName], it ] }
            .collect() :
        Channel.empty()

    PREPARE_GENOME (
        fasta_ch
    )

    gc_wiggle_ch = params.create_gc_wiggle
        ? PREPARE_GENOME.out.gc_wiggle
        : Channel.fromPath(params.gc_wiggle).map { it -> [ [id: it.baseName], it ] }.collect()

    fasta_fai_ch = params.fasta_fai
        ? Channel.fromPath(params.fasta_fai).map{ it -> [ [id:'fai'], it ] }.collect()
        : PREPARE_GENOME.out.fasta_fai

    bin_size_ch = params.bin_size
        ? Channel.value(params.bin_size)
        : Channel.value(50)

    purity_ch = params.purity
        ? Channel.value(params.purity)
        : Channel.value( [  20, 40, 60, 80 ] )

    ploidy_ch = params.purity
        ? Channel.value(params.ploidy)
        : Channel.value(7)

    gamma_ch = params.gamma
        ? Channel.value(params.gamma)
        : Channel.value(40)

    sex_ch = params.sex
        ? Channel.value(params.sex)
        : Channel.value("XX")

//
    // WORKFLOW: Run pipeline
    //
    CNV (
        samplesheet,
        fasta_ch,
        fasta_fai_ch,
        gc_wiggle_ch,
        bin_size_ch,
        purity_ch,
		ploidy_ch,
		gamma_ch,
		sex_ch
    )
    emit:
    multiqc_report = CNV.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    SOTTORIVALAB_CNV (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        SOTTORIVALAB_CNV.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
