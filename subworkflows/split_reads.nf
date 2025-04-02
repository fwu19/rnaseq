/*
* Split fastq files into smaller chunks, run STAR and XenofilteR, and merge filtered bams
*/

include { SPLIT_FASTQ } from '../modules/split_fastq'
include { STAR } from '../modules/star.nf'
include { STAR  as  STAR_HOST } from '../modules/star.nf'
include { XENOFILTER } from '../modules/xenofilter.nf'
include { MERGE_BAM } from '../modules/merge_bam.nf'

workflow SPLIT_READS {
    take:
    ch_reads


    main:
    SPLIT_FASTQ(
        ch_reads
    )
    ch_reads = ch_reads
        .map{ it -> [ it[0].id, it ] }
        .cross(
            SPLIT_FASTQ.out.csv
                .map { it -> it[1] }
                .splitCsv( header: false )
        )
        .map{ it -> [ it[0][1][0], it[1][1], it[1][2], it[1][3] ]}

    /*
    * run STAR alignment 
    * for pdx workflow, align to both graft and host genomes 
    * for pdx workflow, run XenofilteR to remove reads with host origin
    */
    
    ch_bam = Channel.empty()
    ch_bam_host = Channel.empty()
    ch_bam_xeno = Channel.empty()
    STAR(
        ch_reads, 
        params.genome, 
        params.star, 
        params.gtf
    )
    ch_bam = STAR.out.bam
    // [ [meta], val(out_prefix), path(bam) ]

    STAR_HOST(
        ch_reads, 
        params.genome_host, 
        params.star_host, 
        params.gtf_host
    )
    ch_bam_host = STAR_HOST.out.bam
    // [ [meta], val(out_prefix), path(bam) ]

    ch_bam
        .map{ it -> [ [ it[0],it[1] ], it[2] ]}
        .join (
            ch_bam_host
                .map{ it -> [ [ it[0],it[1] ], it[2] ]}
        )
        .map{ it -> [ it[0][0], it[0][1], it[1], it[2] ] }
        .set { ch_bam_paired }
    // ch_bam_paired.view()
    // [ [meta], val(out_prefix), [path/to/graft.{bam,bai}], [path/to/host.{bam,bai}]]

    XENOFILTER(
        ch_bam_paired, 
        params.genome, 
        params.mm_threshold
    )
    ch_bam_xeno = XENOFILTER.out.bam 
    // [ [meta], path/to/filtered.bam ]

        
    MERGE_BAM(
        ch_bam_xeno
            .map{ it -> [ it[0], it[2] ]}
            .groupTuple( by: [0] )
    )
    ch_bam_xeno = MERGE_BAM.out.bam      
    // ch_bam_xeno.view()
    // [ [meta], path("*.{bam,bai}") ]      
    

    emit:
    //versions = ch_versions
    bam_xeno = ch_bam_xeno    

}
