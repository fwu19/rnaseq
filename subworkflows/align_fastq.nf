/*
* Align all fastq files, and run XenofilteR if needed
*/

include { BWA_MEM } from '../modules/bwa_mem.nf'
include { BWA_MEM  as  BWA_MEM_HOST } from '../modules/bwa_mem.nf'
include { STAR } from '../modules/star.nf'
include { STAR  as  STAR_HOST } from '../modules/star.nf'

include { MAKE_INDEX } from './make_index.nf'
include { MAKE_INDEX as MAKE_INDEX_HOST } from './make_index.nf'

workflow ALIGN_FASTQ {
    take:
    ch_reads
    aligner

    main: 
    ch_bam = Channel.empty()
    ch_bai = Channel.empty()
    ch_counts = Channel.empty()
    ch_tx_bam = Channel.empty() // STAR's transcript bam, used for downstream quantification with Salmon
    ch_star_log = Channel.empty()
    ch_bam_host = Channel.empty()
    ch_bai_host = Channel.empty()
    ch_star_log_host = Channel.empty()


    if (params.aligner == 'star'){
        if (params.star == null){
            if (params.genome_fa == null || params.gtf == null){
                exit 1, "Need to specify valid paths to --genome_fa and --gtf"
            }else{
                MAKE_INDEX(
                    file(params.genome_fa, checkIfExists: true),
                    file(params.gtf, checkIfExists: true),
                    "star"
                )
                aligner_index = MAKE_INDEX.out.star
            }
        }else{
            aligner_index = file(params.star, checkIfExists: true)
        }

        STAR(
            ch_reads,
            params.genome, 
            aligner_index, 
            params.gtf
        )
        ch_bam = STAR.out.bam
        // [ [meta], val(out_prefix), path(bam) ]
        ch_bai = STAR.out.bai
        // [ [meta], val(out_prefix), path(bai) ]
        ch_tx_bam = STAR.out.tx_bam
        // [ [meta], val(out_prefix), path(bam) ]
        ch_star_log = STAR.out.log
        // [ [meta], val(out_prefix), path(log) ]
        ch_counts = STAR.out.counts
        // [ [meta], val(out_prefix), path("ReadsPerGene.tab") ]


        /*
        * align to host genome for PDX samples
        */
        if (params.workflow == 'pdx'){
            if (params.star_host == null){
                if (params.genome_fa_host == null || params.gtf_host == null){
                exit 1, "Need to specify valid paths to --genome_fa_host and --gtf_host"
            }else{
                MAKE_INDEX_HOST(
                    file(params.genome_fa_host, checkIfExists: true),
                    file(params.gtf_host, checkIfExists: true),
                    "star"
                )
                aligner_index_host = MAKE_INDEX_HOST.out.star
            }
        }else{
            aligner_index_host = file(params.star_host, checkIfExists: true)
        }

            STAR_HOST(
                ch_reads, 
                params.genome_host, 
                aligner_index_host,  
                params.gtf_host
            )
            ch_bam_host = STAR_HOST.out.bam
            // [ [meta], val(out_prefix), path(bam) ]
            ch_bai_host = STAR_HOST.out.bai
            // [ [meta], val(out_prefix), path(bam) ]
            ch_star_log_host = STAR_HOST.out.log
            // [ [meta], val(out_prefix), path(log) ]

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

        }
    }
    
    if (params.aligner == 'bwa-mem'){
        if (params.bwa == null){
            if (params.genome_fa == null){
                exit 1, "Need to specify valid paths to --genome_fa"
            }else{
                MAKE_INDEX(
                    file(params.genome_fa, checkIfExists: true),
                    '',
                    "bwa"
                )
                aligner_index = MAKE_INDEX.out.bwa
            }
        }else{
            aligner_index = file(params.bwa, checkIfExists: true)
        }

        BWA_MEM(
            ch_reads, 
            params.genome, 
            aligner_index
        )
        ch_bam = BWA_MEM.out.bam
        // [ [meta], val(out_prefix), path(bam) ]
        ch_bai = BWA_MEM.out.bai
        // [ [meta], val(out_prefix), path(bai) ]

        /*
        * align to host genome for PDX samples
        */
        if (params.workflow == 'pdx'){
            if (params.bwa_host == null){
                if (params.genome_fa_host == null){
                    exit 1, "Need to specify valid paths to --genome_fa_host"
                }else{
                    MAKE_INDEX_HOST(
                        file(params.genome_fa_host, checkIfExists: true),
                        '',
                        "bwa"
                    )
                    aligner_index_host = MAKE_INDEX_HOST.out.bwa
                }
            }else{
            aligner_index_host = file(params.bwa_host, checkIfExists: true)
            }
            BWA_MEM_HOST(
                ch_reads, 
                params.genome_host, 
                aligner_index_host
            )
            ch_bam_host = BWA_MEM_HOST.out.bam
            // [ [meta], val(out_prefix), path(bam) ]
            ch_bai_host = BWA_MEM_HOST.out.bai
            // [ [meta], val(out_prefix), path(bam) ]

        }
    }

    emit:
    //versions = ch_versions
    bam = ch_bam
    bai = ch_bai
    tx_bam = ch_tx_bam
    counts = ch_counts
    star_log = ch_star_log
    bam_host = ch_bam_host
    bai_host = ch_bai_host
    star_log_host = ch_star_log_host


}
