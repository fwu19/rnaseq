#!/usr/bin/env nextflow 

include { GET_FASTQ_PATHS } from '../modules/get_fastq_paths'
include { CHECK_INPUT  } from '../modules/check_input.nf'
include { FASTQC  } from '../modules/fastqc.nf'
include { FASTQC as FASTQC_TRIMMED  } from '../modules/fastqc.nf'
include { CUTADAPT  } from '../modules/cutadapt.nf'
include { CAT_FASTQ  } from '../modules/cat_fastq.nf'
include { XENOFILTER } from '../modules/xenofilter.nf'
include { BAM_TO_FASTQ } from '../modules/bam_to_fastq.nf'
include { SALMON  } from '../modules/salmon.nf'
include { ARRIBA  } from '../modules/arriba.nf'
include { GTF2GENES  } from '../modules/gtf2genes.nf'
include { GTF2TRANSCRIPTS  } from '../modules/gtf2transcripts.nf'
include { GTF2BED  } from '../modules/gtf2bed.nf'
include { COLLAPSE_GTF  } from '../modules/collapse_gtf.nf'
include { GTF2FASTA  } from '../modules/gtf2fasta.nf'
include { RNASEQC  } from '../modules/rnaseqc.nf'
include { RSEQC  } from '../modules/rseqc.nf'
include { HS_METRICS  } from '../modules/hs_metrics'
include { SAMTOOLS_VIEW} from '../modules/samtools_view.nf'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_HOST } from '../modules/samtools_view.nf'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_XENO } from '../modules/samtools_view.nf'
include { MULTIQC  } from '../modules/multiqc.nf'
include { MULTIQC_PDX  } from '../modules/multiqc_pdx.nf'
include { FEATURECOUNTS } from '../modules/featureCounts.nf'
include { DIFFERENTIAL_EXPRESSION } from '../modules/differential_expression.nf'
include { DIFFERENTIAL_TRANSCRIPTS } from '../modules/differential_transcripts.nf'
include { GENERATE_GENE_COUNT_MATRIX } from '../modules/generate_gene_count_matrix.nf'
include { GENERATE_TRANSCRIPT_COUNT_MATRIX } from '../modules/generate_transcript_count_matrix.nf'
include { GENERATE_REPORT } from '../modules/generate_report.nf'
include { WRITE_CSV as WRITE_CSV_ALIGN_FASTQ} from '../modules/write_csv.nf'

//include { OUTPUT_PARAMS  } from '../modules/output_params.nf'
//include { TEST  } from './modules/test.nf'

include { MAKE_INDEX } from '../subworkflows/make_index.nf'
include { ALIGN_FASTQ } from '../subworkflows/align_fastq.nf'
include { PROCESS_FASTQ } from '../subworkflows/process_fastq.nf'

samplesheet = null
tools = params.tools.split(',').collect()
skip_tools = params.skip_tools.split(',').join()
if ('arriba' in skip_tools){
    params.run_arriba = false
}
if ('salmon' in skip_tools){
    params.run_salmon = false
    params.run_tx_count = false
    params.run_de = false
}

workflow RNASEQ {
    /*
    * Make references only
    */
    if (params.run_make_index && params.only_ref){
        if (params.genome_fa && params.gtf){
            /*
            gtf = params.gtf.split(',').collect()
            if (gtf.size() > 1){
                CAT(
                    gtf,
                    "genes.gtf"
                )
                gtf = CAT.out.file
            }
            */
            MAKE_INDEX(
                params.genome_fa.split(',').collect(),
                params.gtf,
                params.aligner_indices?:params.aligner.split(',').collect{ it.trim().toLowerCase() }
            )

        }else{
            exit 1, "Need to provide valid paths to genome fastq file and gene gtf file via --genome_fa and --gtf."
        }
    }

    /*
    * Run input check
    */
    ch_input = Channel.fromPath( params.input, checkIfExists: true )
    if (params.run_input_check){
        if ( params.input_dir =~ 'dummy' ){
            if ( params.input =~ 'dummy' ){
                exit 1, 'Neither --input nor --input_dir is specified!'
            }
        }else {
            GET_FASTQ_PATHS (
                Channel.fromPath("${params.input_dir}", checkIfExists: true)
            )
            ch_input = GET_FASTQ_PATHS.out.csv
        }

        CHECK_INPUT(
            ch_input,
            ch_metadata
        )
        samplesheet = CHECK_INPUT.out.csv
        // final samplesheet, one row per id
        //CHECK_INPUT.out.fq.view()
        // samplesheet, one row per pair of fastq
    }

    /*
    * cat fastq files if needed
    */
    ch_reads = Channel.empty()
    if (params.run_cat_fastq){
        CAT_FASTQ(
            CHECK_INPUT.out.fq
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
        
    }else if (!params.only_ref && !params.only_input){
        CHECK_INPUT.out.fq
                .splitCsv(header: true)
                .map {
                    row -> [ row, row.id, row.fastq_1, row.fastq_2 ]
                }
                .set { ch_reads}
    }
    ch_reads_raw = ch_reads
    // ch_reads.view()

    /*
    * run cutadapt if needed
    */
    ch_fastqc_trimmed = Channel.empty()
    ch_cutadapt_js = Channel.empty()
    if (params.run_cut_adapt){
        CUTADAPT(
            ch_reads_raw
        )
        ch_reads = CUTADAPT.out.fq
        ch_reads_trimmed = ch_reads
        ch_cutadapt_js = CUTADAPT.out.js

    }
    // ch_reads.view()
    // [ [meta], meta.id, path("R1.fastq.gz"), path("R2.fastq.gz") ]

    /*
    * run STAR alignment 
    * for pdx workflow, align to both graft and host genomes 
    * for pdx workflow, run XenofilteR to remove reads with host origin
    */
    ch_counts = Channel.empty()
    ch_bam = Channel.empty()
    ch_bai = Channel.empty()
    ch_tx_bam = Channel.empty() // STAR's transcript bam, used for downstream quantification with Salmon
    ch_star_log = Channel.empty()
    ch_bam_host = Channel.empty() 
    ch_bai_host = Channel.empty()
    ch_star_log_host = Channel.empty() 
    ch_bam_xeno = Channel.empty()
    ch_bai_xeno = Channel.empty()


    if (params.run_alignment){
        ALIGN_FASTQ(
            ch_reads,
            params.aligner
        )
        
        ch_bam = ALIGN_FASTQ.out.bam
        ch_bai = ALIGN_FASTQ.out.bai
        ch_counts = ALIGN_FASTQ.out.counts
        ch_tx_bam = ALIGN_FASTQ.out.tx_bam
        ch_star_log = ALIGN_FASTQ.out.star_log
        ch_bam_host = ALIGN_FASTQ.out.bam_host
        ch_bai_host = ALIGN_FASTQ.out.bai_host
        ch_star_log_host = ALIGN_FASTQ.out.star_log_host

        if(params.workflow == 'pdx'){
            if(params.run_split_fastq){
                ALIGN_SPLIT_FASTQ(
                    ch_reads,
                    params.split_size
                )
                ch_bam_xeno = ALIGN_SPLIT_FASTQ.out.bam_xeno
                ch_bai_xeno = ALIGN_SPLIT_FASTQ.out.bai_xeno

            } else {
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
                // [ [meta], val(out_prefix), path/to/filtered.bam ]    
            }

            subdir_aln_fq = params.aligner.toUpperCase()
            WRITE_CSV_ALIGN_FASTQ(
                ch_bam
                    .join(ch_bam_host)
                    .join(ch_bam_xeno)
                    .map { 
                        it -> it[0] + [graft_bam: "${params.outdir}/${subdir_aln_fq}/${params.genome}/_unfiltered/${it[0].id}.bam" ] + [host_bam: "${params.outdir}/${subdir_aln_fq}/${params.host_genome}/_unfiltered/${it[0].id}.bam" ]+ [filtered_graft_bam: "${params.outdir}/${subdir_aln_fq}/${params.genome}/${it[0].id}.bam" ]
                    }
                    .collect(),
                "align_fastq.csv"        
            )

        }else{
            subdir_aln_fq = params.aligner.toUpperCase()
            WRITE_CSV_ALIGN_FASTQ(
                ch_bam
                    .map { it -> it[0] + [bam: "${params.outdir}/${subdir_aln_fq}/${params.genome}/${it[0].id}.bam" ] }
                    .collect(),
                "align_fastq.csv"        
            )
        }
    
    }


    // save graft-only reads
    ch_graft_reads = Channel.empty()
    if(params.run_alignment && params.workflow == 'pdx'){
        if (params.run_arriba || params.only_filter_fastq){
            BAM_TO_FASTQ(
                ch_bam_xeno
            )
            ch_graft_reads = BAM_TO_FASTQ.out.fq
        }
    }


    ch_salmon = Channel.empty()
    if (params.run_alignment && params.run_salmon){
        if (params.tx_fa){
            tx_fa = file(params.tx_fa, checkIfExists:true)
        }else{
            GTF2FASTA(
                params.genome_fa,
                params.gtf
            )
            tx_fa = GTF2FASTA.out.fa
        }
        
        if (params.workflow == 'pdx'){
            //should filter Aligned.toTranscriptome.out.bam
        }else{
            SALMON(
                ch_tx_bam,
                tx_fa
            )
        }
        ch_salmon = SALMON.out.sf
        // [ [meta], val(out_prefix), path("${out_prefix}/") ]
    }

    if (params.run_alignment && params.run_arriba){
        ARRIBA(
            params.workflow == 'pdx' ? ch_graft_reads : ch_reads, 
            params.genome, 
            params.star, 
            params.gtf,
            params.genome_fa,
            params.blacklist,
            params.known_fusions,
            params.protein_domains
        )

    }

    /*
    * run QC
    */    
    ch_fastqc = Channel.empty()
    ch_rnaseqc = Channel.empty()
    ch_rseqc = Channel.empty()
    ch_hs_metrics = Channel.empty()
    ch_bam_stat = Channel.empty()
    ch_bam_stat_host = Channel.empty()
    ch_bam_stat_xeno = Channel.empty()

    ch_bam
        .map{ it -> [ [ it[0], it[1] ], it[2] ]}
        .join (
            ch_bai
                .map{ it -> [ [ it[0], it[1] ], it[2] ]}
        )
        .map { it -> [ it[0][0], it[0][1], it[1], it[2] ]}
        .set { ch_bam_bai }
                    
    
    if (params.workflow == 'pdx'){
        ch_bam_host
            .map{ it -> [ [ it[0], it[1] ], it[2] ]}
            .join (
                ch_bai_host
                    .map{ it -> [ [ it[0], it[1] ], it[2] ]}
            )
            .map { it -> [ it[0][0], it[0][1], it[1], it[2] ]}
            .set { ch_bam_bai_host }

        ch_bam_xeno
            .map{ it -> [ [ it[0], it[1] ], it[2] ]}
            .join (
                ch_bai_xeno
                    .map{ it -> [ [ it[0], it[1] ], it[2] ]}
            )
            .map { it -> [ it[0][0], it[0][1], it[1], it[2] ]}
            .set { ch_bam_bai_xeno }

        ch_bam_bai_qc = ch_bam_bai_xeno
        
    }else{
        ch_bam_bai_qc = ch_bam_bai
        
    }

    if (params.run_qc){

        /*
        * FastQC
        */
        if (params.run_fastqc){
            FASTQC(
                ch_reads_raw
            )
            ch_fastqc = FASTQC.out.qc

            FASTQC_TRIMMED(
                ch_reads_trimmed
            )
            ch_fastqc_trimmed = FASTQC_TRIMMED.out.qc

        }

        /*
        * RNASeQC
        */
        if (params.run_rnaseqc){
            if (params.rnaseqc_gtf){
                collapsed_gtf = file(params.rnaseqc_gtf, checkIfExists:true)
            }else{
                COLLAPSE_GTF(params.gtf)
                collapsed_gtf = COLLAPSE_GTF.out.gtf
            }
            RNASEQC(
                    ch_bam_bai_qc,
                    collapsed_gtf,
                    params.strand,
                    params.read_type
            )            
            ch_rnaseqc = RNASEQC.out.qc
            // [ [meta], path("*") ]
        }

        /*
        * RSeQC
        */
        if (params.run_rseqc){
            if (params.rseqc_bed){
                tx_bed = file(params.tx_bed, checkIfExists:true)
            }else{
                GTF2BED(params.gtf)
                tx_bed = GTF2BED.out.bed
            }
            RSEQC(
                ch_bam_bai_qc,
                    tx_bed
            )
            ch_rseqc = RSEQC.out.qc
            // [ [meta], path("*") ]
        }

        /*
        * GATK collect_hs_metrics
        */
        if (params.run_hs_metrics){
            HS_METRICS(
                ch_bam_bai_qc,
                params.genome_fa,
                params.target_region
            )
            ch_hs_metrics = HS_METRICS.out.qc
            // [ [meta], path("*") ]

        }
        

        /*
        * pdx
        */
        if (params.run_samtools){
            /*
            * samtools
            */
            SAMTOOLS_VIEW(
                ch_bam_bai
            )
            ch_bam_stat = SAMTOOLS_VIEW.out.data

            SAMTOOLS_VIEW_HOST(
                ch_bam_bai_host
            )
            ch_bam_stat_host = SAMTOOLS_VIEW_HOST.out.data

            SAMTOOLS_VIEW_XENO(
                ch_bam_bai_xeno
            )
            ch_bam_stat_xeno = SAMTOOLS_VIEW_XENO.out.data
        }


    }


    /*
    * MultiQC
    */
    ch_multiqc = Channel.empty()
    if (params.run_multiqc){
        if (params.workflow == 'pdx'){
            MULTIQC_PDX(
            ch_multiqc_config,
            ch_fastqc.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_cutadapt_js.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_fastqc_trimmed.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_rseqc.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_rnaseqc.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_hs_metrics.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_star_log.map{it[2]}.flatten().collect().ifEmpty([]), 
            ch_star_log_host.map{it[2]}.flatten().collect().ifEmpty([]), 
            ch_bam_stat.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_bam_stat_host.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_bam_stat_xeno.map{it[1]}.flatten().collect().ifEmpty([])
            )
            ch_multiqc = MULTIQC_PDX.out.data
            //ch_multiqc.view()

        }else{
            MULTIQC(
            ch_multiqc_config,
            ch_fastqc.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_cutadapt_js.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_fastqc_trimmed.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_rseqc.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_rnaseqc.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_hs_metrics.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_star_log.map{it[2]}.flatten().collect().ifEmpty([]), 
            ch_counts.map{it[2]}.flatten().collect().ifEmpty([]),
            ch_salmon.map{it[2]}.flatten().collect().ifEmpty([])
            )
            ch_multiqc = MULTIQC.out.data
            //ch_multiqc.view()
        }
    }

    /*
    * run featureCounts
    */
    if(params.workflow == 'pdx'){
        ch_bam_fc = ch_bam_xeno
    }else{
        ch_bam_fc = ch_bam
    }
    if (params.run_featurecounts){
        FEATURECOUNTS(
                ch_bam_fc, 
                params.gtf,
                params.read_type,
                params.strand
        )
        
        ch_counts = FEATURECOUNTS.out.counts
        // [ [meta], val(out_prefix), path("count.txt") ]

    }

    /*
    *   collect gene-level count matrix and run PCA
    */
    ch_gene_rds = Channel.empty()
    if (params.run_alignment && (params.run_gene_count || params.run_dt)){
        if (params.gene_txt){
            gene_txt = file(params.gene_txt, checkIfExists:true)
        }else if (params.gtf){
            GTF2GENES(params.gtf)
            gene_txt = GTF2GENES.out.txt
        }
    }

    if (params.run_alignment && params.run_gene_count){
        GENERATE_GENE_COUNT_MATRIX(
            samplesheet, 
            ch_counts.map{it[2]}.collect().ifEmpty([]), 
            gene_txt,
            params.length_col?:"gene_length",
            params.strand,
            params.workflow
        )
        ch_gene_rds = GENERATE_GENE_COUNT_MATRIX.out.rds
    }
    
    /*
    * differential genes
    */
    ch_de = Channel.empty()
    if (params.run_alignment && params.run_gene_count && params.run_de){
        DIFFERENTIAL_EXPRESSION(
            samplesheet, 
            Channel.fromPath(params.comparison, checkIfExists: true), 
            ch_gene_rds, 
            params.fdr,
            params.fc,
            params.fdr2,
            params.fc2,
            params.de_gene_txt?:gene_txt,
            params.de_gene_type?:'all'
        )
        ch_de = DIFFERENTIAL_EXPRESSION.out.rds
    }

    /*
    *   collect transcript-level count matrix and run PCA
    */
    ch_tx_rds = Channel.empty()
    if (params.run_alignment && params.run_salmon && params.run_tx_count){
        if (params.tx_txt){
            tx_txt = file(params.tx_txt, checkIfExists:true)
        }else if (params.gtf){
            GTF2TRANSCRIPTS(params.gtf)
            tx_txt = GTF2TRANSCRIPTS.out.txt
        }
        
        GENERATE_TRANSCRIPT_COUNT_MATRIX(
            samplesheet, 
            ch_salmon.map{it[2]}.collect().ifEmpty([]), 
            tx_txt,
            "EffectiveLength"
        )
        ch_tx_rds = GENERATE_TRANSCRIPT_COUNT_MATRIX.out.rds
    }


    /*
    * differential transcripts
    */
    ch_dt = Channel.empty()
    if (params.run_alignment && params.run_salmon && params.run_tx_count && params.run_dt){
        DIFFERENTIAL_TRANSCRIPTS(
            samplesheet, 
            Channel.fromPath(params.comparison, checkIfExists: true), 
            ch_tx_rds, 
            "EffectiveLength",
            params.fdr,
            params.fc,
            params.fdr2,
            params.fc2,
            params.de_gene_txt?:gene_txt,
            params.de_gene_type?:'all'
        )
        ch_dt = DIFFERENTIAL_TRANSCRIPTS.out.rds
    }


    /*
    * Generate a report
    */
    if (params.run_report){
        GENERATE_REPORT(
            params.workflow,
            samplesheet,
            ch_multiqc.ifEmpty([]),
            ch_gene_rds.ifEmpty([]),
            ch_de.ifEmpty([]),
            ch_tx_rds.ifEmpty([]),
            ch_dt.ifEmpty([]),
            ch_hs_metrics.map{it[1]}.flatten().collect().ifEmpty([]),
            params.report_dir ? Channel.fromPath(params.report_dir, type: 'dir', checkIfExists: true) : Channel.fromPath("$projectDir/assets/report/", type: 'dir', checkIfExists: true)
        )
    }
}
