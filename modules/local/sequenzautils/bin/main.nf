process SEQUENZAUTILS_BIN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py39h67e14b5_5' :
        'biocontainers/sequenza-utils:3.0.0--py39h67e14b5_5' }"

    input:
    tuple val(meta), path(concat_seqz), path(concat_seqz_tbi)
    val bin_size

    output:
    tuple val(meta), path("*bin${params.bin}.seqz.gz"), emit: seqz_bin
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sequenza-utils \\
        seqz_binning \\
        -w $bin_size \\
        --seqz $concat_seqz \\
        -o - |\\
        gzip > ${prefix}_bin${bin_size}.seqz.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "This is a stub run. ext.when = ${task.ext.when}"
    echo | gzip > ${prefix}.seqz.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: stub
    END_VERSIONS
    """

}
