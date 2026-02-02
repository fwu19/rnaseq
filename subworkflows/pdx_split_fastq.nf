/*
* Split fastq files into smaller chunks, run STAR and XenofilteR, and merge filtered bams
*/

include { SPLIT_FASTQ } from '../modules/split_fastq'
include { BWA_MEM } from '../modules/bwa_mem.nf'
include { BWA_MEM  as  BWA_MEM_HOST } from '../modules/bwa_mem.nf'
include { STAR } from '../modules/star.nf'
include { STAR  as  STAR_HOST } from '../modules/star.nf'
include { XENOFILTER } from '../modules/xenofilter.nf'
include { MERGE_BAM } from '../modules/merge_bam.nf'

workflow PDX_SPLIT_FASTQ {
    take:
    ch_reads
    split_size
    aligner_index
    aligner_index_host


    main:
    ch_bam_bai = Channel.empty()
    ch_bam_bai_host = Channel.empty()
    ch_bam_bai_xeno = Channel.empty()
    ch_versions = Channel.empty()

    SPLIT_FASTQ(
        ch_reads,
        split_size
    )
    ch_versions = ch_versions.mix(SPLIT_FASTQ.out.versions.first())

    ch_reads
        .map{ it -> [ it[1], it ] }
        .cross(
            SPLIT_FASTQ.out.csv
                .map { it -> it[2] }
                .splitCsv( header: false )
        )
        .map{ it -> [ it[0][1][0], it[1][1], it[1][2], it[1][3] ]}
        .set{ ch_reads_in}


    if (params.aligner == 'star'){
        STAR(
            ch_reads_in,
            params.genome, 
            aligner_index, 
            params.gtf
        )
        ch_bam_bai = STAR.out.bam_bai
        // [ [meta], val(out_prefix), path(bam), path(bai) ]
        ch_versions = ch_versions.mix(STAR.out.versions.first())

        /*
        * align to host genome for PDX samples
        */
        STAR_HOST(
                ch_reads_in, 
                params.genome_host, 
                aligner_index_host, 
                params.gtf_host
        )
        ch_bam_bai_host = STAR_HOST.out.bam_bai
        // [ [meta], val(out_prefix), path(bam), path(bai) ]
        ch_versions = ch_versions.mix(STAR_HOST.out.versions.first())
    }
    
    if (params.aligner == 'bwa-mem'){
        BWA_MEM(
        ch_reads_in, 
        params.genome, 
        aligner_index
        )
        ch_bam_bai = BWA_MEM.out.bam_bai
        // [ [meta], val(out_prefix), path(bam), path(bai) ]
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

        /*
        * align to host genome for PDX samples
        */
        BWA_MEM_HOST(
            ch_reads_in, 
            params.genome_host, 
            aligner_index_host
        )
        ch_bam_bai_host = BWA_MEM_HOST.out.bam_bai
        // [ [meta], val(out_prefix), path(bam), path(bai) ]
        ch_versions = ch_versions.mix(BWA_MEM_HOST.out.versions.first())
        
    }

    ch_bam_bai
        .map{ it -> [ [ it[0], it[1] ], it[2], it[3] ]}
        .join (
            ch_bam_bai_host
                .map{ it -> [ [ it[0], it[1] ], it[2], it[3] ]}
        )
        .map{ it -> [ it[0][0], it[0][1], it[1], it[2], it[3], it[4] ] }
        .set { ch_bam_bai_paired }
    // [ [meta], val(out_prefix), [path/to/graft.{bam,bai}], [path/to/host.{bam,bai}]]

    XENOFILTER(
        ch_bam_bai_paired, 
        params.genome, 
        params.mm_threshold
    )
    ch_bam_bai_xeno = XENOFILTER.out.bam_bai 
    // [ [meta], val(out_prefix), path(bam), path(bai) ]    
    ch_versions = ch_versions.mix(XENOFILTER.out.versions.first())

    MERGE_BAM(
        ch_bam_bai_xeno
        .map{ it -> [ it[0], it[2] ]}
        .groupTuple( by: [0] )
    )
    ch_bam_bai_xeno = MERGE_BAM.out.bam_bai    
    // ch_bam_xeno.view()
    // [ [meta], val(out_prefix), path("*.{bam,bai}") ]      
    ch_versions = ch_versions.mix(MERGE_BAM.out.versions.first())    

    emit:
    bam_bai = ch_bam_bai_xeno    
    versions = ch_versions
    
}
