/*
* QC fastq 
*/

include { FASTQC  } from '../modules/fastqc.nf'
include { FASTQC as FASTQC_TRIMMED  } from '../modules/fastqc.nf'


workflow QC_FASTQ {
    take:
    raw_reads
    trimmed_reads

    main: 
    ch_fastqc = Channel.empty()
    ch_fastq_trimmed = Channel.empty()
    

    if (params.run_fastqc){
        FASTQC(
            raw_reads
        )
        ch_fastqc = FASTQC.out.qc

        FASTQC_TRIMMED(
            trimmed_reads
        )
        ch_fastqc_trimmed = FASTQC_TRIMMED.out.qc
    }

    emit:
    fastqc = ch_fastqc
    fastqc_trimmed = ch_fastqc_trimmed


}
