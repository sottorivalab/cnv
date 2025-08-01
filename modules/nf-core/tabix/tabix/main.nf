process TABIX_TABIX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92859404d861ae01afb87e2b789aebc71c0ab546397af890c7df74e4ee22c8dd/data' :
        'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa' }"

    input:
    tuple val(meta), path(seqz)

    output:
    tuple val(meta), path("*concat.seqz.gz"), path("*concat.seqz.gz.tbi"), emit: concat_seqz
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def all_seqz = seqz.collect { it.getName() }.join(' ')
    """
    echo "Merging these files:"
    for f in $all_seqz; do echo \$f; done
    head -n 1 <( zcat ${seqz[0]} ) > header
    cat header <( zcat $all_seqz | awk '{if (NR!=1 && \$1 != "chromosome") {print \$0}}' ) | bgzip > \
        ${prefix}_concat.seqz.gz
    tabix -f -s 1 -b 2 -e 2 -S 1 ${prefix}_concat.seqz.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(tabix -h 2>&1 | sed -n 's/^Version: //p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def all_seqz = seqz.collect { it.getName() }.join(' ')
    """
    touch ${prefix}_concat.seqz.gz
    touch ${prefix}_concat.seqz.gz.tbi
    echo $all_seqz 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: stub 
    END_VERSIONS
    """
}
