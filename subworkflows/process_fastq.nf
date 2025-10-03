/*
* Check fastq files and cat fastq if needed
*/

include { CAT_FASTQ } from '../modules/cat_fastq.nf'

workflow PROCESS_FASTQ {
    take:
    samplesheet
    fq
    cat_fastq

    main: 
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

    emit:
    reads = ch_reads

}
