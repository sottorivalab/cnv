process SEQUENZAUTILS_RSEQZ {
    tag "${meta.id}_${meta.sex}_${purity}"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::r-sequenza=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-sequenza:2.1.2--r351h6115d3f_2' :
        'biocontainers/r-sequenza%2.1.2--r351h6115d3f_2' }"

    input:
    tuple val(meta), path(seqz_bin), path(seqz_tbi)
    each purity


    output:
    tuple val(meta), path("${meta.id}_${purity}"), emit: rseqz
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${purity}"
    template "cn_analysis.sh"

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${purity}"
    """
    echo "analyse_cn_sequenza.R ${seqz_bin} ${seqz_bin} ${prefix} ${meta.id} ${meta.gender} ${meta.ploidy} ${purity} ${meta.gamma}"
    mkdir ${prefix}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: 3.0.0
    END_VERSIONS
    """
}
