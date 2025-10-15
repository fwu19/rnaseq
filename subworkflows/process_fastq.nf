/*
* Check fastq files and cat fastq if needed
*/

include { CAT_FASTQ } from '../modules/cat_fastq.nf'
include { CUTADAPT  } from '../modules/cutadapt.nf'

workflow PROCESS_FASTQ {
    take:
    samplesheet
    fq
    cat_fastq

    main: 
    ch_reads = Channel.empty()
    ch_reads_trimmed = Channel.empty()
    ch_cutadapt_js = Channel.empty()

    fq
        .splitCsv(header: true)
        .map {
            row -> [ row.id, row.fastq_1, row.fastq_2 ]
        }
        .groupTuple (by: [0])
        .branch {
            id, fastq_1, fastq_2 ->
                single  : fastq_1.size() == 1
                    return [ id, fastq_1, fastq_2 ]
                multiple: fastq_1.size() > 1
                    return [ id, fastq_1, fastq_2 ]
        }
        .set { ch_fastq }

    merged_fastq = Channel.empty()
    if ( cat_fastq ){
        CAT_FASTQ(
            ch_fastq.multiple
        )  
        merged_fastq = CAT_FASTQ.out.reads      
    }else{
        // add a suffix to meta.id to make unique
    }
    
    samplesheet
            .splitCsv( header: true )
            .map {
                row -> [ row.id, row ]
            }
            .join ( 
                merged_fastq 
                    .mix(ch_fastq.single)
            )
            .map { it -> [ it[1], it[1].id, it[2], it[3] ]}
            .set { ch_reads }

    /*
    * old codes
    if (cat_fastq){
        CAT_FASTQ(
            fq
                .splitCsv(header: true)
                .map {
                    row -> [ row.id, row.fastq_1, row.fastq_2 ]
                }
                .groupTuple (by: [0])
        )

        samplesheet
            .splitCsv( header: true )
            .map {
                row -> [ row.id, row ]
            }
            .join ( CAT_FASTQ.out.reads )
            .map { it -> [ it[1], it[1].id, it[2], it[3] ]}
            .set { ch_reads }
        
    }else {
        fq
            .splitCsv(header: true)
            .map {
                row -> [ row, row.id, row.fastq_1, row.fastq_2 ]
            }
            .set { ch_reads}
    }

    */


    /* consider to add splitting fastq here ? */

    /*
    * run cutadapt
    */
    if (params.run_cut_adapt){
        CUTADAPT(
            ch_reads
        )
        ch_reads_trimmed = CUTADAPT.out.fq
        ch_cutadapt_js = CUTADAPT.out.js

    }
    // ch_reads.view()
    // [ [meta], meta.id, path("R1.fastq.gz"), path("R2.fastq.gz") ]


    emit:
    reads = ch_reads
    reads_trimmed = ch_reads_trimmed
    cutadapt_js = ch_cutadapt_js

}
