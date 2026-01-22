/*
* Check fastq files and cat fastq if needed
*/

include { CAT_FASTQ } from '../modules/cat_fastq.nf'
include { CUTADAPT  } from '../modules/cutadapt.nf'
include { FASTP  } from '../modules/fastp.nf'
include { WRITE_CSV as WRITE_CSV_TRIM_FASTQ} from '../modules/write_csv.nf'

workflow PROCESS_FASTQ {
    take:
    samplesheet // one row per unique id
    fq // one row per fastq file
    cat_fastq
    trimmer

    main: 
    ch_reads = Channel.empty()
    ch_reads_trimmed = Channel.empty()
    ch_cutadapt_js = Channel.empty()
    ch_fastp_js = Channel.empty()
    ch_fastp_html = Channel.empty()
    ch_versions = Channel.empty()


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




    /* consider to add splitting fastq here ? 
    SPLIT_FASTQ(
    )
    */

    /*
    * run fastp
    */
    if (( trimmer == 'fastp' || !params.adapters ) && params.run_cut_adapt){
        FASTP(
            ch_reads
        )
        ch_reads_trimmed = FASTP.out.fq
        ch_fastp_js = FASTP.out.js
        ch_fastp_html = FASTP.out.html
        ch_versions = ch_versions.mix(FASTP.out.versions.first())

    } else if (trimmer == 'cutadapt' && params.run_cut_adapt){
        CUTADAPT(
            ch_reads
        )
        ch_reads_trimmed = CUTADAPT.out.fq
        ch_cutadapt_js = CUTADAPT.out.js
        ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())
    }
    // ch_reads.view()
    // [ [meta], meta.id, path("R1.fastq.gz"), path("R2.fastq.gz") ]

    /* Write trimmed fastq paths to csv */
    if (params.only_trim_fastq){
        def my_dir = new File("${params.outdir}")
        def outdir = my_dir.absolutePath
        WRITE_CSV_TRIM_FASTQ(
                ch_reads_trimmed
                    .map { 
                        it -> it[0] + [ trimmed_fastq_1: "${outdir}/trimmed_fastq/${it[2].name}" ] + [ trimmed_fastq_2: "${outdir}/trimmed_fastq/${it[3].name}" ] 
                    }
                    .collect(),
                "trim_fastq.csv"        
        )
        ch_versions = ch_versions.mix(WRITE_CSV_TRIM_FASTQ.out.versions.first())

    }

    emit:
    reads = ch_reads
    reads_trimmed = ch_reads_trimmed
    cutadapt_js = ch_cutadapt_js
    fastp_js = ch_fastp_js
    fastp_html = ch_fastp_html
    versions = ch_versions

}
