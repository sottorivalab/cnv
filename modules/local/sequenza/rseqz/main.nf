process SEQUENZAUTILS_RSEQZ {
    tag "${meta.id}_${meta.sex}_${purity}"

	conda (params.enable_conda ? "bioconda::r-sequenza=3.0.0" : null)
    container 'sottorivalab_sequenza.sif'

    input:
    tuple val(meta), path(seqz_bin)// , path(filtered_seqz_bin)
    each purity


    output:
    tuple val(meta), path("${meta.id}_${purity}"), emit: rseqz
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${purity}"
	def seqz = seqz_bin.toString().replaceAll(/\.gz/,'') 
    """
    zcat $seqz_bin > $seqz
	analyse_cn_sequenza.R ${seqz} ${prefix} ${meta.id} ${meta.sex} ${meta.ploidy} ${purity} ${meta.gamma}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """

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
