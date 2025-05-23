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

workflow PREPARE_GENOME {
    take:
    fasta                // channel: [mandatory] fasta

    main:
    versions = Channel.empty()
	SEQUENZAUTILS_GCWIGGLE(fasta)
	SAMTOOLS_FAIDX(fasta, [ [ id:'no_fai' ], [] ] )
    
	versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
	versions = versions.mix(SEQUENZAUTILS_GCWIGGLE.out.versions) // channel: versions.yml

    emit:
    gc_wiggle = SEQUENZAUTILS_GCWIGGLE.out.gc_wiggle.collect() // path: gcwiggle
    fasta_fai = SAMTOOLS_FAIDX.out.fai.collect()               // path: genome.fasta.fai

    versions        // channel: [ versions.yml ]
}
