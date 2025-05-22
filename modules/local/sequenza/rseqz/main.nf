process SEQUENZAUTILS_RSEQZ {
    tag "${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::r-sequenza=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-sequenza%3A3.0.0--r42h3342da4_5' :
        'biocontainers/r-sequenza%3A3.0.0--r42h3342da4_5' }"

    input:
    tuple val(meta), path(seqz_bin), val(tissue), val(purity)
    val  gender
    val  ploidy
    val  seq_gam

    output:
    tuple val(meta), val(tissue), val(purity), path("${tissue}"), emit: rseqz
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def seqz_in = "${seqz_bin.toString().minus(".gz")}"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${tissue}"
    """
    zcat ${seqz_bin} > $seqz_in
    analyse_cn_sequenza.R ${seqz_in} ${prefix} ${meta.id} ${gender} ${ploidy} ${purity} ${seq_gam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def seqz_in = "${seqz_bin.toString().minus(".gz")}"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${tissue}"
    """
    echo "analyse_cn_sequenza.R ${seqz_in} ${prefix} ${meta.id} ${gender} ${ploidy} ${purity} ${seq_gam}"
    mkdir ${tissue}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: 3.0.0
    END_VERSIONS
    """
}
