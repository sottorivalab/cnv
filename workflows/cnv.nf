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
include { SEQUENZAUTILS_BAM2SEQZ } from '../modules/nf-core/sequenzautils/bam2seqz/'
include { TABIX_TABIX            } from '../modules/nf-core/tabix/tabix/'
include { SEQUENZAUTILS_BIN      } from '../modules/local/sequenzautils/bin/'
include { SEQUENZAUTILS_RSEQZ    } from '../modules/local/sequenza/rseqz/'
include { GAWK                   } from '../modules/nf-core/gawk/'
include { ASCAT                  } from '../modules/nf-core/ascat'
include { UNZIP                  } from '../modules/nf-core/unzip'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CNV {

    take:
    ch_samplesheet      // channel: samplesheet read in from --input
    ch_fasta            // channel from reference folder
    ch_fasta_fai        // channel from reference folder
    ch_fasta_gzi        // channel from reference folder
    ch_wig              // channel from reference folder
    ch_bin_size         // default or supplied  sequenza-utils bin size
    ch_purity           // purity default or supplied
    ch_ascat_alleles
    ch_ascat_loci
	ch_bed_file
    ch_ascat_loci_gc
    ch_ascat_loci_rt

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Branching logic based on presence of seqz
    ch_samplesheet
        .branch {
            seqz:  it.size() == 2 && it[1].toString().endsWith(".seqz.gz")
            bam:   it.size() == 3
        }
        .set { ch_input }

    // stub just run on three if it is a stub test
    if (!workflow.stubRun){
    chromosome_list = ch_fasta_fai
                    .map{ it[1] }
                    .splitCsv(sep: "\t")
                    .map{ chr -> chr[0] }
					.filter( ~/^chr\d+|^chr[X,Y]|^\d+|[X,Y]/ )
					.collect()
    } else {
        chromosome_list = Channel.value(['chr1', 'chr2', 'chr3'])
    }

    //
    // MODULE: BAM2SEQZ
    //
    //
	// split into two channels normal and tumour for bam2seqz and ascat
	//

    ch_input.bam
        .branch { meta, bam, bai ->
            normal: meta.status == 0
            tumour: meta.status == 1
        }
        .set { bam_ch }

    ch_input.bam
        .map { meta, bam, bai -> tuple(meta.patient, [meta, bam, bai]) }
        .groupTuple()
        .flatMap { patient_id, samples ->
            def normals = samples.findAll { it[0].status == 0 }
            def tumours = samples.findAll { it[0].status == 1 }

            if (normals.isEmpty()) return []

            def normal = normals[0]
            def meta_n = normal[0]
            def bam_n  = normal[1]
            def bai_n  = normal[2]

            return tumours.collect { tumour ->
                def meta_t = tumour[0]
                def bam_t  = tumour[1]
                def bai_t  = tumour[2]

                def new_meta = [ patient: patient_id, id: meta_t.id, sex: meta_t.sex, ploidy: meta_t.ploidy, gamma: meta_t.gamma ]

                tuple(new_meta, bam_n, bai_n, bam_t, bai_t)
            }
        }
       .set { ch_paired_tumour_normal }

    SEQUENZAUTILS_BAM2SEQZ (
        ch_paired_tumour_normal,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_gzi,
        ch_wig,
        chromosome_list
    )

    ch_versions = ch_versions.mix(SEQUENZAUTILS_BAM2SEQZ.out.versions)

	// collect bam2seqz output for binning
    SEQUENZAUTILS_BAM2SEQZ.out.seqz
        .groupTuple()
        .map { meta, files ->
            def chr_order = (1..22).collect { it.toString() } + ['X', 'Y']

            def chr_index = { file ->
                def match = (file.getName() =~ /chr(\d+|X|Y)/)
                def chr = match ? match[0][1] : null
                return chr && chr_order.contains(chr) ? chr_order.indexOf(chr) : 999
            }

            def canonical = files.findAll { it.getName() =~ /chr(\d+|X|Y)/ }
            def sorted = canonical.sort { a, b -> chr_index(a) <=> chr_index(b) }

            tuple(meta, sorted)
        }
        .set { merge_seqz_input }
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
        TABIX_TABIX.out.concat_seqz,
        ch_bin_size
    )

    ch_versions = ch_versions.mix(SEQUENZAUTILS_BIN.out.versions.first())

    ch_seqz_final = ch_input.seqz.mix(SEQUENZAUTILS_BIN.out.seqz_bin)

    SEQUENZAUTILS_RSEQZ(
        ch_seqz_final,
        ch_purity
        )


    ASCAT( ch_paired_tumour_normal,
           ch_ascat_alleles,
           ch_ascat_loci,
           ch_bed_file,
           ch_fasta,
           ch_ascat_loci_gc,
           ch_ascat_loci_rt
        )

    ch_versions = ch_versions.mix(ASCAT.out.versions.first())

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
