process SEQUENZAUTILS_GCWIGGLE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py38h6ed170a_2' :
        'biocontainers/sequenza-utils:3.0.0--py38h6ed170a_2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.wig.gz"), emit: gc_wiggle
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sequenza-utils \\
        gc_wiggle \\
        $args \\
        --fasta $fasta \\
        -o ${prefix}.wig.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "This is a stub run. ext.when = ${task.ext.when}"
	echo | gzip > ${prefix}.wig.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: stub 
    END_VERSIONS
    """
}
