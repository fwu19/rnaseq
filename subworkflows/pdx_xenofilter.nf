/*
* Align all fastq files, and run XenofilteR if needed
*/

include { ALIGN_FASTQ } from '../subworkflows/align_fastq.nf'
include { ALIGN_FASTQ as ALIGN_FASTQ_HOST } from '../subworkflows/align_fastq.nf'

include { SPLIT_FASTQ } from '../modules/split_fastq'
include { XENOFILTER } from '../modules/xenofilter.nf'
include { MERGE_BAM } from '../modules/merge_bam.nf'
include { BAM_TO_FASTQ } from '../modules/bam_to_fastq.nf'

workflow PDX_XENOFILTER {
    take:
    ch_reads
    reference_index
    reference_index_host
    ch_bam_bai
    ch_bam_bai_host

    main: 
    ch_bam_bai_xeno = Channel.empty()
    ch_graft_reads = Channel.empty()
    ch_versions = Channel.empty()

    if (!params.run_split_fastq){
        XENOFILTER(
            ch_bam_bai
                .map{ it -> [ [ it[0], it[1] ], it[2], it[3] ]}
                .join (
                    ch_bam_bai_host
                    .map{ it -> [ [ it[0], it[1] ], it[2], it[3] ]}
                )
                .map{ it -> [ it[0][0], it[0][1], it[1], it[2], it[3], it[4] ] }, 
            params.mm_threshold
        )
        ch_bam_bai_xeno = XENOFILTER.out.bam_bai 
        ch_versions = ch_versions.mix(XENOFILTER.out.versions)

    }else{        
        SPLIT_FASTQ(
            ch_reads,
            params.split_size
        )
        ch_versions = ch_versions.mix(SPLIT_FASTQ.out.versions)

        ch_reads
            .map{ it -> [ it[1], it ] }
            .cross(
                SPLIT_FASTQ.out.csv
                    .map { it -> it[2] }
                    .splitCsv( header: false )
            )
            .map{ it -> [ it[0][1][0], it[1][1], it[1][2], it[1][3] ]}
            .set{ ch_reads_in}

        /* align to graft genome */
        ALIGN_FASTQ(
            ch_reads_in,
            reference_index
        )
        
        ch_software_versions = ch_software_versions.mix(ALIGN_FASTQ.out.versions)

        /* align to host genome */
        ALIGN_FASTQ_HOST(
            ch_reads_in,
            reference_index_host
        )
        ch_software_versions = ch_software_versions.mix(ALIGN_FASTQ_HOST.out.versions)

        /* join graft and host bam files */
        ALIGN_FASTQ.out.bam_bai
            .map{ it -> [ [ it[0], it[1] ], it[2], it[3] ]}
            .join (
                ALIGN_FASTQ_HOST.out.bam_bai
                    .map{ it -> [ [ it[0], it[1] ], it[2], it[3] ]}
            )
            .map{ it -> [ it[0][0], it[0][1], it[1], it[2], it[3], it[4] ] }
            .set { ch_bam_bai_paired }
        // [ [meta], val(out_prefix), [path/to/graft.{bam,bai}], [path/to/host.{bam,bai}]]

        XENOFILTER(
            ch_bam_bai_paired, 
            params.mm_threshold
        )        
        // [ [meta], val(out_prefix), path(bam), path(bai) ]    
        ch_versions = ch_versions.mix(XENOFILTER.out.versions)

        MERGE_BAM(
            XENOFILTER.out.bam_bai
                .map{ it -> [ it[0], it[2] ]}
                .groupTuple( by: [0] )
        )
        ch_bam_bai_xeno = MERGE_BAM.out.bam_bai    
        // ch_bam_xeno.view()
        // [ [meta], val(out_prefix), path("*.{bam,bai}") ]      
        ch_versions = ch_versions.mix(MERGE_BAM.out.versions)    

    }

    /*
    * save graft-only reads
    */
    if(params.run_arriba || params.only_filter_fastq){
        BAM_TO_FASTQ(
            ch_bam_bai_xeno
        )
        ch_graft_reads = BAM_TO_FASTQ.out.fq
        ch_versions = ch_versions.mix(BAM_TO_FASTQ.out.versions)
    }


    emit:
    bam_bai_xeno = ch_bam_bai_xeno
    graft_reads = ch_graft_reads
    versions = ch_versions

}
