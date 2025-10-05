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
