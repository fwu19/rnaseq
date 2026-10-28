/*
* Align all fastq files
*/

include { BWA_MEM } from '../modules/bwa_mem.nf'
include { STAR } from '../modules/star.nf'
include { BWA_MEM as BWA_MEM_HOST } from '../modules/bwa_mem.nf'
include { STAR as STAR_HOST } from '../modules/star.nf'

workflow ALIGN_FASTQ {
    take:
    ch_reads
    index_dir
    gtf
    ch_reads_host
    index_host_dir
    gtf_host

    main: 
    ch_bam_bai = Channel.empty()
    ch_counts = Channel.empty()
    ch_tx_bam = Channel.empty() // STAR's transcript bam, used for downstream quantification with Salmon
    ch_star_log = Channel.empty()
    ch_bam_bai_host = Channel.empty()
    ch_counts_host = Channel.empty()
    ch_tx_bam_host = Channel.empty() // STAR's transcript bam, used for downstream quantification with Salmon
    ch_star_log_host = Channel.empty()
    ch_versions = Channel.empty()

    if (params.aligner == 'star'){            
        STAR(
            ch_reads,
            index_dir, 
            gtf
        )
        ch_bam_bai = STAR.out.bam_bai
        // [ [meta], val(out_prefix), path(bam), path(bai) ]
        ch_tx_bam = STAR.out.tx_bam
        // [ [meta], val(out_prefix), path(bam) ]
        ch_star_log = STAR.out.log
        // [ [meta], val(out_prefix), path(log) ]
        ch_counts = STAR.out.counts
        // [ [meta], val(out_prefix), path("ReadsPerGene.tab") ]
        ch_versions = ch_versions.mix(STAR.out.versions.first())

        if (params.workflow == 'pdx'){
            STAR_HOST(
                ch_reads_host,
                index_host_dir, 
                gtf_host
            )
            ch_bam_bai_host = STAR_HOST.out.bam_bai
            ch_tx_bam_host = STAR_HOST.out.tx_bam
            ch_star_log_host = STAR_HOST.out.log
            ch_counts_host = STAR_HOST.out.counts
            ch_versions = ch_versions.mix(STAR_HOST.out.versions.first())
        }
    }
    
    if (params.aligner == 'bwa-mem'){
        BWA_MEM(
            ch_reads, 
            index_dir
        )
        ch_bam_bai = BWA_MEM.out.bam_bai
        // [ [meta], val(out_prefix), path(bam), path(bai) ]
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

        if (params.workflow == 'pdx'){
            BWA_MEM_HOST(
                ch_reads_host, 
                index_host_dir
            )
            ch_bam_bai_host = BWA_MEM_HOST.out.bam_bai
            ch_versions = ch_versions.mix(BWA_MEM_HOST.out.versions.first())
        }
    }


    emit:
    bam_bai = ch_bam_bai
    tx_bam = ch_tx_bam
    counts = ch_counts
    star_log = ch_star_log
    bam_bai_host = ch_bam_bai_host
    tx_bam_host= ch_tx_bam_host
    counts_host = ch_counts_host
    star_log_host = ch_star_log_host
    versions = ch_versions

}
