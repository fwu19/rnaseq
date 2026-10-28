/*
* Make references
*/

include { MAKE_STAR } from '../modules/make_star.nf'
include { MAKE_STAR as MAKE_STAR_HOST } from '../modules/make_star.nf'
include { MAKE_BWA  } from '../modules/make_bwa.nf'
include { MAKE_BWA as MAKE_BWA_HOST  } from '../modules/make_bwa.nf'
include { MAKE_XENGSORT } from '../modules/make_xengsort.nf'

workflow BUILD_INDEX {
    take:
    genome_fa
    gtf
    genome_fa_host
    gtf_host
    reference_index

    main: 
    star_dir = Channel.empty()
    star_host_dir = Channel.empty()
    bwa_dir = Channel.empty()
    bwa_host_dir = Channel.empty()
    bowtie_dir = Channel.empty()
    bowtie2_dir = Channel.empty()
    salmon_dir = Channel.empty()
    xengsort_dir = Channel.empty()
    ch_versions = Channel.empty()

    // Init aligners
    def reference_index_list = ["bowtie", "bowtie2", "bwa", "bwa-mem", "salmon", "star", "xengsort"]
    index_list = reference_index.split(',').collect{ it.trim().toLowerCase()}
    if ((reference_index_list + index_list).unique().size() != reference_index_list.size()) {
        exit 1, "Invalid aligner option found in ${reference_index}. Valid options: ${reference_index_list.join(', ')}"
    } 

    genome_fa_list = (genome_fa ?: '').split(',').collect{ it.trim() }.findAll{ it != '' }
    //println "genome_fa_list: ${genome_fa_list}"
    gtf_list = (gtf ?: '').split(',').collect{ it.trim() }.findAll{ it != '' }
    genome_fa_host_list = (genome_fa_host ?: '').split(',').collect{ it.trim() }.findAll{ it != '' }
    println "genome_fa_host_list: ${genome_fa_host_list}"
    gtf_host_list = (gtf_host ?: '').split(',').collect{ it.trim() }.findAll{ it != '' }

    if ("star" in index_list){
        MAKE_STAR(
            genome_fa_list ? Channel.value(genome_fa_list) : Channel.empty(),
            gtf_list ? Channel.value(gtf_list) : []
        )
        star_dir = MAKE_STAR.out.dir
        ch_versions = ch_versions.mix(MAKE_STAR.out.versions)

        MAKE_STAR_HOST(
            genome_fa_host_list ? Channel.value(genome_fa_host_list) : Channel.empty(),
            gtf_host_list ? Channel.value(gtf_host_list) : []
        )
        star_host_dir = MAKE_STAR_HOST.out.dir
        ch_versions = ch_versions.mix(MAKE_STAR_HOST.out.versions)
    }

    if ("bwa" in index_list || "bwa-mem" in index_list){
        MAKE_BWA(
            genome_fa_list ? Channel.value(genome_fa_list) : Channel.empty()
        )
        bwa_dir = MAKE_BWA.out.dir
        ch_versions = ch_versions.mix(MAKE_BWA.out.versions)

        MAKE_BWA_HOST(
            genome_fa_host_list ? Channel.value(genome_fa_host_list) : Channel.empty()
        )
        bwa_host_dir = MAKE_BWA_HOST.out.dir
        ch_versions = ch_versions.mix(MAKE_BWA_HOST.out.versions)
    }

    if ("xengsort" in index_list){
        MAKE_XENGSORT(
            genome_fa_list ? Channel.value(genome_fa_list) : Channel.empty(),
            genome_fa_host_list ? Channel.value(genome_fa_host_list) : Channel.empty()
        )
        xengsort_dir = MAKE_XENGSORT.out.dir
        ch_versions = ch_versions.mix(MAKE_XENGSORT.out.versions)
    }

    emit:
    star = star_dir
    star_host = star_host_dir
    bwa = bwa_dir
    bwa_host = bwa_host_dir
    xengsort = xengsort_dir
    versions = ch_versions

}
