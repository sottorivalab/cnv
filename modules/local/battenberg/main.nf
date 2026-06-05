process BATTENBERG {
    tag "$meta.id"
    label 'process_high'

    conda (params.containsKey('enable_conda') && params.enable_conda ? "bioconda::battenberg" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/chelauk/battenberg-beagle:latest' :
        'ghcr.io/chelauk/battenberg-beagle:latest' }"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumour_bam), path(tumour_bai)
    path impute_info
    path g1000_loci
    path problem_loci
    path gc_correction
    path rt_correction
    path g1000_alleles
    path beagle_jar
    path beagle_ref
    path beagle_plink

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path "versions.yml"             , emit: versions, topic: 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args instanceof Map ? task.ext.args : [:]
    prefix = task.ext.prefix ?: "${meta.id}"
    tumour_name = task.ext.tumour_name ?: "${meta.id}"
    normal_name = task.ext.normal_name ?: (meta.normal_id ?: "${meta.patient}_normal")
    is_male = meta.sex in ['XY', 'M', 'male', 'Male', 'MALE'] ? 'TRUE' : 'FALSE'
    data_type = args.data_type ?: 'wgs'
    impute_exe = args.impute_exe ?: 'impute2'
    allelecounter_exe = args.allelecounter_exe ?: 'alleleCounter'
    platform_gamma = args.platform_gamma ?: 1
    phasing_gamma = args.phasing_gamma ?: 1
    segmentation_gamma = args.segmentation_gamma ?: 10
    segmentation_kmin = args.segmentation_kmin ?: 3
    phasing_kmin = args.phasing_kmin ?: 1
    clonality_dist_metric = args.clonality_dist_metric ?: 0
    ascat_dist_metric = args.ascat_dist_metric ?: 1
    min_ploidy = args.min_ploidy ?: 1.6
    max_ploidy = args.max_ploidy ?: 4.8
    min_rho = args.min_rho ?: 0.1
    min_goodness = args.min_goodness ?: 0.63
    uninformative_baf_threshold = args.uninformative_baf_threshold ?: 0.51
    min_normal_depth = args.min_normal_depth ?: 10
    min_base_qual = args.min_base_qual ?: 20
    min_map_qual = args.min_map_qual ?: 35
    calc_seg_baf_option = args.calc_seg_baf_option ?: 1
    skip_allele_counting = args.skip_allele_counting == null ? 'FALSE' : (args.skip_allele_counting in [true, 'true', 'TRUE', 'T', '1'] ? 'TRUE' : 'FALSE')
    skip_preprocessing = args.skip_preprocessing == null ? 'FALSE' : (args.skip_preprocessing in [true, 'true', 'TRUE', 'T', '1'] ? 'TRUE' : 'FALSE')
    skip_phasing = args.skip_phasing == null ? 'FALSE' : (args.skip_phasing in [true, 'true', 'TRUE', 'T', '1'] ? 'TRUE' : 'FALSE')
    use_beagle = args.use_beagle == null ? 'TRUE' : (args.use_beagle in [true, 'true', 'TRUE', 'T', '1'] ? 'TRUE' : 'FALSE')
    beagle_max_mem = args.beagle_max_mem ?: 8
    beagle_threads = args.beagle_threads ?: 1
    beagle_window = args.beagle_window ?: 40
    beagle_overlap = args.beagle_overlap ?: 4
    heterozygous_filter = args.heterozygous_filter ?: 'none'
    template 'battenberg.Rscript'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        battenberg: stub
    END_VERSIONS
    """
}
