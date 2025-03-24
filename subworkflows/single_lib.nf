/*
* Split fastq files into smaller chunks, run STAR and XenofilteR, and merge filtered bams
*/

include { XENOFILTER } from '../modules/xenofilter.nf'

workflow SINGLE_LIB {
    take:
    ch_bam
    ch_bam_host

    main:    
    ch_bam_xeno = Channel.empty()

    ch_bam
        .map{ it -> [ [ it[0], it[1] ], it[2] ]}
        .join (
            ch_bam_host
                .map{ it -> [ [ it[0], it[1] ], it[2] ]}
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
    // [ [meta], val(out_prefix), path/to/filtered.bam ]    

    emit:
    //versions = ch_versions]
    bam_host = ch_bam_host
    bam_xeno = ch_bam_xeno    

}
