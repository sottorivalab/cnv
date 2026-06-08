process INFERSEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.22.0--py39hcada746_1' :
        'quay.io/biocontainers/pysam:0.22.0--py39hcada746_1' }"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumour_bam), path(tumour_bai)

    output:
    tuple val(meta), path("${prefix}.sex.txt")       , emit: sex
    tuple val(meta), path("${prefix}.sex_depths.tsv"), emit: depths
    path "versions.yml"                              , emit: versions, topic: 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'infersex.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "unknown" > ${prefix}.sex.txt
    cat <<-END_DEPTHS > ${prefix}.sex_depths.tsv
    group	chrom	start	end	depth
    autosomes	chr1	10000000	10100000	30
    chrX	chrX	50000000	50100000	15
    chrY	chrY	7000000	7100000	3
    END_DEPTHS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: stub
    END_VERSIONS
    """
}
