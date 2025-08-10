process SEQUENZAUTILS_BAM2SEQZ {
    tag "${meta.id}_${chromosome}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py39he88f293_8' :
        'biocontainers/sequenza-utils:3.0.0--py39he88f293_8' }"

    input:
    tuple val(meta), path(normalbam), path(normalbam_bai), path(bam), path(bai)
    tuple val(meta_1), path(fasta)
    tuple val(meta_2), path(fasta_fai)
    tuple val(meta_3), path(fasta_gzi)
	tuple val(meta_4), path(wigfile)
    each chromosome

    output:
    tuple val(meta), path("*.gz"), emit: seqz
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: "-C $chromosome"
    def prefix = task.ext.prefix ?: "${meta.id}_${chromosome}"
    """
    sequenza-utils \\
        bam2seqz \\
        $args \\
        -n $normalbam \\
        -t $bam \\
        --fasta $fasta \\
        -gc $wigfile \\
        --het ${params.het} \\
        --hom ${params.hom} \\
        -o ${prefix}.gz 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: "-C $chromosome"
    def prefix = task.ext.prefix ?: "${meta.id}_${chromosome}"
    """
    echo "This is a stub run. ext.when = ${task.ext.when}"
        echo | gzip > ${prefix}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: stub
    END_VERSIONS
    """

}
