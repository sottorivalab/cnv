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
params.gc_wiggle = getGenomeAttribute('wiggle')
params.fasta_gzi = getGenomeAttribute('fasta_gzi')

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
    // println(params.gc_wiggle)
    gc_wiggle_ch = params.create_gc_wiggle
        ? PREPARE_GENOME.out.gc_wiggle
        : Channel.fromPath(params.gc_wiggle).map { it -> [ [id: it.baseName], it ] }.collect()

    fasta_fai_ch = params.fasta_fai
        ? Channel.fromPath(params.fasta_fai).map{ it -> [ [id:'fai'], it ] }.collect()
        : PREPARE_GENOME.out.fasta_fai

    fasta_gzi_ch = params.fasta_gzi
        ? Channel.fromPath(params.fasta_gzi).map{ it -> [ [id:'gzi'], it ] }.collect()
        : PREPARE_GENOME.out.fasta_gzi // TODO make gzi if absent

    bin_size_ch = params.bin_size
        ? Channel.value(params.bin_size)
        : Channel.value(50)

    if (params.purity == "range" )
        { purity_ch =  Channel.value([20, 40, 60, 80, 100]) }
        else if (params.purity instanceof Number )
        { purity_ch = Channel.value(params.purity) }


    //
    // WORKFLOW: Run pipeline
    //

    CNV (
        samplesheet,
        fasta_ch,
        fasta_fai_ch,
		fasta_gzi_ch,
        gc_wiggle_ch,
        bin_size_ch,
        purity_ch
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
