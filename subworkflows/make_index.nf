/*
* Make references
*/

include { MAKE_STAR } from '../modules/make_star.nf'
include { MAKE_BWA  } from '../modules/make_bwa.nf'


workflow MAKE_INDEX {
    take:
    genome_fa
    gtf
    aligner_indices

    main: 
    star_dir = Channel.empty()
    bwa_dir = Channel.empty()
    bowtie_dir = Channel.empty()
    bowtie2_dir = Channel.empty()
    salmon_dir = Channel.empty()

    // Init aligners
    def aligner_indices_list = ["bowtie", "bowtie2", "bwa", "salmon", "star"]

    if ((aligner_indices_list + aligner_indices).unique().size() != aligner_indices_list.size()) {
        exit 1, "Invalid aligner option: ${aligner_indices}. Valid options: ${aligner_indices_list.join(', ')}"
    }

    if ("star" in aligner_indices){
        MAKE_STAR(
            genome_fa,
            gtf
        )
        star_dir = MAKE_STAR.out.dir
    }

    if ("bwa" in aligner_indices){
        MAKE_BWA(
            genome_fa
        )
        bwa_dir = MAKE_BWA.out.dir
    }

    emit:
    star = star_dir
    bwa = bwa_dir


}
