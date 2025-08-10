#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sottorivalab/cnv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/sottorivalab/cnv
----------------------------------------------------------------------------------------
*/
nextflow.enable.dsl = 2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.ascat_alleles = getGenomeAttribute('ascat_alleles')
params.ascat_genome  = getGenomeAttribute('ascat_genome')
params.ascat_loci    = getGenomeAttribute('ascat_loci')
params.ascat_loci_gc = getGenomeAttribute('ascat_loci_gc')
params.ascat_loci_rt = getGenomeAttribute('ascat_loci_rt')
params.fasta         = getGenomeAttribute('fasta')
params.fasta_fai     = getGenomeAttribute('fasta_fai')
params.gc_wiggle     = getGenomeAttribute('gc_wiggle')
params.fasta_gzi     = getGenomeAttribute('fasta_gzi')
params.bed_file      = getGenomeAttribute('bed_file')

// Initialize value channels based on params, defined in the params.genomes[params.genome] scope
println "DEBUG params.ascat_genome=${params.ascat_genome}"
ascat_genome         = params.ascat_genome   ?:  Channel.empty()

println "DEBUG ascat_genome=${ascat_genome}"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CNV                     } from './workflows/cnv'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_cnv_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_cnv_pipeline'

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
        fasta_ch,
        params.ascat_alleles, 
        params.ascat_loci,    
        params.ascat_loci_gc, 
        params.ascat_loci_rt 
    )

    // For ASCAT, extracted from zip or tar.gz files
    allele_files           = PREPARE_GENOME.out.allele_files
    loci_files             = PREPARE_GENOME.out.loci_files
    gc_file                = PREPARE_GENOME.out.gc_file
    rt_file                = PREPARE_GENOME.out.rt_file

    gc_wiggle_ch = params.gc_wiggle
        ? Channel.fromPath(params.gc_wiggle).map { it -> [ [id:'gc_wiggle'], it ] }.collect()
        : PREPARE_GENOME.out.gc_wiggle

    bed_file_ch = params.bed_file
        ? Channel.fromPath(params.bed_file, checkIfExists: true)
        : PREPARE_GENOME.out.bed_file

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
        purity_ch,
        allele_files, 
        loci_files,    
        bed_file_ch,
        gc_file,
        rt_file       
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genome config file e.g. fasta
//

def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
