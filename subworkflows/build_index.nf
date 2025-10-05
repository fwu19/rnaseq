/*
* Make references
*/

include { MAKE_STAR } from '../modules/make_star.nf'
include { MAKE_BWA  } from '../modules/make_bwa.nf'


workflow BUILD_INDEX {
    take:
    genome_fa
    gtf
    aligner_index

    main: 
    star_dir = Channel.empty()
    bwa_dir = Channel.empty()
    bowtie_dir = Channel.empty()
    bowtie2_dir = Channel.empty()
    salmon_dir = Channel.empty()

    // Init aligners
    def aligner_index_list = ["bowtie", "bowtie2", "bwa", "salmon", "star"]
    index_list = aligner_index.split(',').collect{ it.trim().toLowerCase()}
    if ((aligner_index_list + index_list).unique().size() != aligner_index_list.size()) {
        exit 1, "Invalid aligner option found in ${aligner_index}. Valid options: ${aligner_index_list.join(', ')}"
    } 

    if ("star" in index_list){
        MAKE_STAR(
            genome_fa.split(',').collect(),
            gtf
        )
        star_dir = MAKE_STAR.out.dir
    }

    if ("bwa" in index_list){
        MAKE_BWA(
            genome_fa.split(',').collect()
        )
        bwa_dir = MAKE_BWA.out.dir
    }

    emit:
    star = star_dir
    bwa = bwa_dir


}
