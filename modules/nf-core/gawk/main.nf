process GAWK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*filtered.seqz.gz"), emit: filtered_seqz
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    zcat ${input} | \
    awk 'BEGIN {OFS="\t"} 
    NR == 1 {
        print "chromosome", "position", "depth.normal", "depth.tumor", "depth.ratio", "Af", "GC.percent"; next
    }
    (\$1 ~ /^chr([0-9]+|X|Y)\$/ && \$4 >= 10 && \$5 >= 10 && \$7 > 0.1 && \$7 < 0.9) {
       print \$1, \$2, \$4, \$5, \$6, \$7, \$10
    }' | gzip > ${prefix}_filtered.seqz.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_filtered.seqz.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
