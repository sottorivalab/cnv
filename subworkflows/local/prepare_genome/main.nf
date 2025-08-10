//
// PREPARE GENOME
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments

include { SEQUENZAUTILS_GCWIGGLE } from '../../../modules/nf-core/sequenzautils/gcwiggle'
include { SAMTOOLS_FAIDX         } from '../../../modules/nf-core/samtools/faidx'
include { UNZIP as UNZIP_ALLELES } from '../../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_GC      } from '../../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_LOCI    } from '../../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_RT      } from '../../../modules/nf-core/unzip/main'

workflow PREPARE_GENOME {
    take:
    fasta                // channel: [mandatory] fasta
    ascat_alleles        // params.ascat_alleles
    ascat_loci           // params.ascat_loci
    ascat_loci_gc        // params.ascat_loci_gc
    ascat_loci_rt        // params.ascat_loci_rt

    main:
    versions = Channel.empty()
	SEQUENZAUTILS_GCWIGGLE(fasta)
    SAMTOOLS_FAIDX(fasta, [ [ id:'no_fai' ], [] ] )
    
    // prepare ascat and controlfreec reference files
    
    if (!ascat_alleles) allele_files = Channel.empty()
    else if (ascat_alleles.endsWith(".zip")) {
        UNZIP_ALLELES(Channel.fromPath(file(ascat_alleles)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
        allele_files = UNZIP_ALLELES.out.unzipped_archive.map{ it[1] }
        versions = versions.mix(UNZIP_ALLELES.out.versions)
    } else allele_files = Channel.fromPath(ascat_alleles).collect()

    if (!ascat_loci) loci_files = Channel.empty()
    else if (ascat_loci.endsWith(".zip")) {
        UNZIP_LOCI(Channel.fromPath(file(ascat_loci)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
        loci_files = UNZIP_LOCI.out.unzipped_archive.map{ it[1] }
        versions = versions.mix(UNZIP_LOCI.out.versions)
    } else loci_files = Channel.fromPath(ascat_loci).collect()

    if (!ascat_loci_gc) gc_file = Channel.value([])
    else if (ascat_loci_gc.endsWith(".zip")) {
        UNZIP_GC(Channel.fromPath(file(ascat_loci_gc)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
        gc_file = UNZIP_GC.out.unzipped_archive.map{ it[1] }
        versions = versions.mix(UNZIP_GC.out.versions)
    } else gc_file = Channel.fromPath(ascat_loci_gc).collect()

    if (!ascat_loci_rt) rt_file = Channel.value([])
    else if (ascat_loci_rt.endsWith(".zip")) {
        UNZIP_RT(Channel.fromPath(file(ascat_loci_rt)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
        rt_file = UNZIP_RT.out.unzipped_archive.map{ it[1] }
        versions = versions.mix(UNZIP_RT.out.versions)
    } else rt_file = Channel.fromPath(ascat_loci_rt).collect()

	versions = versions.mix(UNZIP_RT.out.versions)
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
	versions = versions.mix(SEQUENZAUTILS_GCWIGGLE.out.versions) // channel: versions.yml

    emit:
    gc_wiggle = SEQUENZAUTILS_GCWIGGLE.out.gc_wiggle.collect() // path: gcwiggle
    fasta_fai = SAMTOOLS_FAIDX.out.fai.collect()               // path: genome.fasta.fai
    allele_files   // path: allele_files
    loci_files      // path: loci_files
    gc_file         // path: gc_file
    rt_file         // path: rt_file
    versions        // channel: [ versions.yml ]
}
