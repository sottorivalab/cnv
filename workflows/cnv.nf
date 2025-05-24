/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_cnv_pipeline'
include { SEQUENZAUTILS_BAM2SEQZ } from '../modules/nf-core/sequenzautils/bam2seqz/main'
include { TABIX_TABIX            } from '../modules/nf-core/tabix/tabix/main'
include { SEQUENZAUTILS_BIN      } from '../modules/local/sequenzautils/bin/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CNV {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta       // channel from reference folder
    ch_fasta_fai   // channel from reference folder
    ch_wig         // channel from reference folder
    ch_bin_size    // default sequenza-utils bin size

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // prepare chromosomes for bam2seqz this allows parallelization over chromosomes
    chromosome_list = ch_fasta_fai
                    .map{ it[1] }
                    .splitCsv(sep: "\t")
                    .map{ chr -> chr[0] }
                    .filter( ~/^chr\d+|^chr[X,Y]|^\d+|[X,Y]/ )
                    .collect()
    //
    // MODULE: BAM2SEQZ
    //

    SEQUENZAUTILS_BAM2SEQZ (
        ch_samplesheet,
        ch_fasta,
        ch_wig,
        chromosome_list
    )

    ch_versions = ch_versions.mix(SEQUENZAUTILS_BAM2SEQZ.out.versions)

	// collect bam2seqz output for binning
    SEQUENZAUTILS_BAM2SEQZ.out.seqz
                          .groupTuple()
                          .set{merge_seqz_input}

	//
    // MODULE: TABIX
    //
    // merge and index bam2seqz output
    TABIX_TABIX(merge_seqz_input)

    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    //
    // MODULE: BINNING
    //
    // binning of bam2seqz output default bin size is 50kb
    SEQUENZAUTILS_BIN (
        TABIX_TABIX.out.seqz,
        ch_bin_size
    )


	//
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'cnv_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
