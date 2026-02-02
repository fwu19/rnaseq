/*
* Align all fastq files, and run XenofilteR if needed
*/

include { BUILD_INDEX as BUILD_INDEX_HOST } from './build_index.nf'
include { PDX_SPLIT_FASTQ } from './pdx_split_fastq.nf'

include { BWA_MEM } from '../modules/bwa_mem.nf'
include { BWA_MEM  as  BWA_MEM_HOST } from '../modules/bwa_mem.nf'
include { STAR } from '../modules/star.nf'
include { STAR  as  STAR_HOST } from '../modules/star.nf'
include { XENOFILTER } from '../modules/xenofilter.nf'
include { BAM_TO_FASTQ } from '../modules/bam_to_fastq.nf'
include { WRITE_CSV as WRITE_CSV_ALIGN_FASTQ} from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_STAR_COUNTS} from '../modules/write_csv.nf'

workflow ALIGN_FASTQ {
    take:
    ch_reads
    split_fastq // if true, don't run xenofilteR
    aligner_index

    main: 
    ch_bam_bai = Channel.empty()
    ch_counts = Channel.empty()
    ch_tx_bam = Channel.empty() // STAR's transcript bam, used for downstream quantification with Salmon
    ch_star_log = Channel.empty()
    ch_bam_bai_host = Channel.empty()
    ch_star_log_host = Channel.empty()
    ch_bam_bai_xeno = Channel.empty()
    ch_graft_reads = Channel.empty()
    ch_versions = Channel.empty()

    if (params.aligner == 'star'){
            
        STAR(
            ch_reads,
            params.genome, 
            aligner_index, 
            params.gtf
        )
        ch_bam_bai = STAR.out.bam_bai
        // [ [meta], val(out_prefix), path(bam), path(bai) ]
        ch_tx_bam = STAR.out.tx_bam
        // [ [meta], val(out_prefix), path(bam) ]
        ch_star_log = STAR.out.log
        // [ [meta], val(out_prefix), path(log) ]
        ch_counts = STAR.out.counts
        // [ [meta], val(out_prefix), path("ReadsPerGene.tab") ]
        ch_versions = ch_versions.mix(STAR.out.versions.first())

        /*
        * align to host genome for PDX samples
        */
        if (params.workflow == 'pdx'){
            /* check index */
            if (params.star_host == null){
                if (params.genome_fa_host == null){
                    exit 1, "Need to specify valid paths to --genome_fa_host"
                }else{
                    BUILD_INDEX_HOST(
                        file(params.genome_fa_host, checkIfExists: true),
                        file(params.gtf_host?:"$projectDir/assets/dummy_file.csv", checkIfExists: true),
                        "star"
                    )
                    aligner_index_host = BUILD_INDEX_HOST.out.star
                    ch_versions = ch_versions.mix(BUILD_INDEX_HOST.out.versions.first())
                }
            }else{
                aligner_index_host = file(params.star_host, checkIfExists: true)
            }

            /* align to host genome */
            STAR_HOST(
                ch_reads, 
                params.genome_host, 
                aligner_index_host,  
                params.gtf_host?:"$projectDir/assets/dummy_file.csv"
            )
            ch_bam_bai_host = STAR_HOST.out.bam_bai
            // [ [meta], val(out_prefix), path(bam), path(bai) ]
            ch_star_log_host = STAR_HOST.out.log
            // [ [meta], val(out_prefix), path(log) ]
            ch_versions = ch_versions.mix(STAR_HOST.out.versions.first())

        }
    }
    
    if (params.aligner == 'bwa-mem'){

        BWA_MEM(
            ch_reads, 
            params.genome, 
            aligner_index
        )
        ch_bam_bai = BWA_MEM.out.bam_bai
        // [ [meta], val(out_prefix), path(bam), path(bai) ]
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

        /*
        * align to host genome for PDX samples
        */
        if (params.workflow == 'pdx'){
            if (params.bwa_host == null){
                if (params.genome_fa_host == null){
                    exit 1, "Need to specify valid paths to --genome_fa_host"
                }else{
                    BUILD_INDEX_HOST(
                        file(params.genome_fa_host, checkIfExists: true),
                        "$projectDir/assets/dummy_file.csv",
                        "bwa"
                    )
                    aligner_index_host = BUILD_INDEX_HOST.out.bwa
                    ch_versions = ch_versions.mix(BUILD_INDEX_HOST.out.versions.first())
                }
            }else{
            aligner_index_host = file(params.bwa_host, checkIfExists: true)
            }
            BWA_MEM_HOST(
                ch_reads, 
                params.genome_host, 
                aligner_index_host
            )
            ch_bam_bai_host = BWA_MEM_HOST.out.bam_bai
            // [ [meta], val(out_prefix), path(bam), path(bai) ]
            ch_versions = ch_versions.mix(BWA_MEM_HOST.out.versions.first())
        }
    }

    /* filter host reads */
    if (params.workflow == 'pdx'){
        if (split_fastq){
            PDX_SPLIT_FASTQ(
                ch_reads,
                params.split_size,
                aligner_index,
                aligner_index_host
            )
            ch_bam_bai_xeno = PDX_SPLIT_FASTQ.out.bam_bai 
            ch_versions = ch_versions.mix(PDX_SPLIT_FASTQ.out.versions.first())
        }else{
            XENOFILTER(
                ch_bam_bai
                    .map{ it -> [ [ it[0], it[1] ], it[2], it[3] ]}
                    .join (
                        ch_bam_bai_host
                        .map{ it -> [ [ it[0], it[1] ], it[2], it[3] ]}
                    )
                    .map{ it -> [ it[0][0], it[0][1], it[1], it[2], it[3], it[4] ] }, 
                params.genome, 
                params.mm_threshold
            )
            ch_bam_bai_xeno = XENOFILTER.out.bam_bai 
            ch_versions = ch_versions.mix(XENOFILTER.out.versions.first())
        }

        // save graft-only reads
        if(params.run_arriba || params.only_filter_fastq){
            BAM_TO_FASTQ(
                ch_bam_bai_xeno
            )
            ch_graft_reads = BAM_TO_FASTQ.out.fq
            ch_software_versions = ch_software_versions.mix(BAM_TO_FASTQ.out.versions)
        }

    }


    emit:
    bam_bai = ch_bam_bai
    tx_bam = ch_tx_bam
    counts = ch_counts
    star_log = ch_star_log
    bam_bai_host = ch_bam_bai_host
    star_log_host = ch_star_log_host
    bam_bai_xeno = ch_bam_bai_xeno
    graft_reads = ch_graft_reads
    versions = ch_versions

}
