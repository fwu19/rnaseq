/*
* Align all fastq files, and run XenofilteR if needed
*/

include { BWA_MEM } from '../modules/bwa_mem.nf'
include { BWA_MEM  as  BWA_MEM_HOST } from '../modules/bwa_mem.nf'
include { STAR } from '../modules/star.nf'
include { STAR  as  STAR_HOST } from '../modules/star.nf'
include { XENOFILTER } from '../modules/xenofilter.nf'
include { WRITE_CSV as WRITE_CSV_ALIGN_FASTQ} from '../modules/write_csv.nf'

include { BUILD_INDEX as BUILD_INDEX_HOST } from './build_index.nf'
include { PDX_SPLIT_FASTQ } from './pdx_split_fastq.nf'

workflow ALIGN_FASTQ {
    take:
    ch_reads
    split_fastq // if true, don't run xenofilteR
    aligner_index
    workflow
    write_csv

    main: 
    ch_bam = Channel.empty()
    ch_bai = Channel.empty()
    ch_counts = Channel.empty()
    ch_tx_bam = Channel.empty() // STAR's transcript bam, used for downstream quantification with Salmon
    ch_star_log = Channel.empty()
    ch_bam_host = Channel.empty()
    ch_bai_host = Channel.empty()
    ch_star_log_host = Channel.empty()
    ch_bam_xeno = Channel.empty()
    ch_bai_xeno = Channel.empty()
    ch_versions = Channel.empty()

    if (params.aligner == 'star'){
            
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
        ch_versions = ch_versions.mix(STAR.out.versions)

        /*
        * align to host genome for PDX samples
        */
        if (workflow == 'pdx'){
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
                    ch_versions = ch_versions.mix(BUILD_INDEX_HOST.out.versions)
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
            ch_bam_host = STAR_HOST.out.bam
            // [ [meta], val(out_prefix), path(bam) ]
            ch_bai_host = STAR_HOST.out.bai
            // [ [meta], val(out_prefix), path(bam) ]
            ch_star_log_host = STAR_HOST.out.log
            // [ [meta], val(out_prefix), path(log) ]
            ch_versions = ch_versions.mix(STAR_HOST.out.versions)

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

        BWA_MEM(
            ch_reads, 
            params.genome, 
            aligner_index
        )
        ch_bam = BWA_MEM.out.bam
        // [ [meta], val(out_prefix), path(bam) ]
        ch_bai = BWA_MEM.out.bai
        // [ [meta], val(out_prefix), path(bai) ]
        ch_versions = ch_versions.mix(BWA_MEM.out.versions)

        /*
        * align to host genome for PDX samples
        */
        if (workflow == 'pdx'){
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
                    ch_versions = ch_versions.mix(BUILD_INDEX_HOST.out.versions)
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
            ch_versions = ch_versions.mix(BWA_MEM_HOST.out.versions)
        }
    }

    /* filter host reads */
    if (workflow == 'pdx'){
        if (split_fastq){
            PDX_SPLIT_FASTQ(
                ch_reads,
                params.split_size,
                aligner_index,
                aligner_index_host
            )
            ch_bam_xeno = PDX_SPLIT_FASTQ.out.bam 
            ch_bai_xeno = PDX_SPLIT_FASTQ.out.bai
            ch_versions = ch_versions.mix(PDX_SPLIT_FASTQ.out.versions)
        }else{
            XENOFILTER(
                ch_bam
                .map{ it -> [ [ it[0], it[1] ], it[2] ]}
                .join (
                            ch_bai
                                .map{ it -> [ [ it[0], it[1] ], it[2] ]}
                )
                .join (
                            ch_bam_host
                                .map{ it -> [ [ it[0], it[1] ], it[2] ]}
                )
                .join (
                            ch_bai_host
                                .map{ it -> [ [ it[0], it[1] ], it[2] ]}
                )
                .map{ it -> [ it[0][0], it[0][1], it[1], it[2], it[3], it[4] ] }, 
                params.genome, 
                params.mm_threshold
            )
            ch_bam_xeno = XENOFILTER.out.bam 
            ch_bai_xeno = XENOFILTER.out.bai
            ch_versions = ch_versions.mix(XENOFILTER.out.versions)
        }
    }

    /* write csv */
    if(workflow == 'pdx' && write_csv){
        subdir_aln_fq = params.aligner.toUpperCase()
        WRITE_CSV_ALIGN_FASTQ(
                ch_bam
                    .join(ch_bam_host)
                    .join(ch_bam_xeno)
                    .map { 
                        it -> it[0] + [graft_bam: "${subdir_aln_fq}/${params.genome}/_unfiltered/${it[0].id}.bam" ] + [host_bam: "${subdir_aln_fq}/${params.host_genome}/_unfiltered/${it[0].id}.bam" ]+ [filtered_graft_bam: "${subdir_aln_fq}/${params.genome}/${it[0].id}.bam" ]
                    }
                    .collect(),
                "mapped.csv"        
        )
    
    }else if (params.aligner == 'star' && write_csv){
        subdir_aln_fq = params.aligner.toUpperCase()
        WRITE_CSV_ALIGN_FASTQ(
                ch_bam
                    .map { 
                        it -> it[0] + [bam: "${subdir_aln_fq}/${params.genome}/${it[0].id}.bam"] + [tx_bam: "${subdir_aln_fq}/${params.genome}/_work/${it[0].id}/Aligned.toTranscriptome.out.bam" ] + [gene_count: "${subdir_aln_fq}/${params.genome}/_work/${it[0].id}/ReadsPerGene.out.tab" ] 
                    }
                    .collect(),
                "mapped.csv"        
        )
        
    }else if (params.aligner == 'bwa-mem' && write_csv){
        subdir_aln_fq = params.aligner.toUpperCase()
        WRITE_CSV_ALIGN_FASTQ(
                ch_bam
                    .map { 
                        it -> it[0] + [bam: "${subdir_aln_fq}/${params.genome}/${it[0].id}.bam"]
                    }
                    .collect(),
                "mapped.csv"        
        )
        
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
    bam_xeno = ch_bam_xeno
    bai_xeno = ch_bai_xeno
    versions = ch_versions

}
