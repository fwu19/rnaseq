#!/usr/bin/env nextflow 

include { GET_FASTQ_PATHS } from '../modules/get_fastq_paths'
include { CHECK_INPUT  } from '../modules/check_input.nf'
include { CUTADAPT  } from '../modules/cutadapt.nf'
include { BAM_TO_FASTQ } from '../modules/bam_to_fastq.nf'
include { ARRIBA  } from '../modules/arriba.nf'
include { GTF2GENES  } from '../modules/gtf2genes.nf'
include { MULTIQC  } from '../modules/multiqc.nf'
include { MULTIQC_PDX  } from '../modules/multiqc_pdx.nf'
include { GENERATE_REPORT } from '../modules/generate_report.nf'

include { BUILD_INDEX } from '../subworkflows/build_index.nf'
include { ALIGN_FASTQ } from '../subworkflows/align_fastq.nf'
include { PROCESS_FASTQ } from '../subworkflows/process_fastq.nf'
include { QC_FASTQ } from '../subworkflows/qc_fastq.nf'
include { QC_ALIGNMENT } from '../subworkflows/qc_alignment.nf'
include { QUANT_GENES } from '../subworkflows/quant_genes.nf'
include { MAP_TRANSCRIPTS } from '../subworkflows/map_transcripts.nf'
include { QUANT_TRANSCRIPTS } from '../subworkflows/quant_transcripts.nf'


ch_metadata = params.metadata ? Channel.fromPath( params.metadata, checkIfExists: true ) : Channel.fromPath("$projectDir/assets/dummy_file.csv", checkIfExists: true)

ch_multiqc_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true ) : Channel.fromPath("$projectDir/assets/multiqc_config.yml")

samplesheet = file("$projectDir/assets/dummy_file.csv", checkIfExists: true)
gene_txt = file("$projectDir/assets/dummy_file.csv", checkIfExists: true)

//tools = params.tools.split(',').collect()
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
    if (params.build_index) {
        if (params.genome_fa && params.gtf){
            BUILD_INDEX(
                params.genome_fa,
                params.gtf,
                params.aligner_index ?: params.aligner 
            )

        }else{
            exit 1, "Need to provide valid paths to genome fastq file and gene gtf file via --genome_fa and --gtf."
        }
    }

    /*
    * Run input check
    */
    ch_input = Channel.fromPath( params.input, checkIfExists: true )
    if (params.step == "mapping" && params.run_input_check){
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
        fq = CHECK_INPUT.out.fq
        // final samplesheet, one row per id
        //CHECK_INPUT.out.fq.view()
        // samplesheet, one row per pair of fastq
    }

    /*
    * cat fastq files
    */
    ch_reads = Channel.empty()
    if (params.step == "mapping" && params.run_process_fastq){
        PROCESS_FASTQ(
            samplesheet,
            fq,
            params.cat_fastq
        )
        ch_reads = PROCESS_FASTQ.out.reads
        ch_reads_raw = ch_reads
        // ch_reads.view()
    }

    /*
    * run cutadapt
    */
    ch_fastqc_trimmed = Channel.empty()
    ch_cutadapt_js = Channel.empty()
    if (params.step == "mapping" && params.run_cut_adapt){
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


    if (params.step == "mapping" && params.run_alignment){
        ALIGN_FASTQ(
            ch_reads,
            params.run_split_fastq // if true, dont run xenofilteR
        )
        
        ch_bam = ALIGN_FASTQ.out.bam
        ch_bai = ALIGN_FASTQ.out.bai
        ch_counts = ALIGN_FASTQ.out.counts
        ch_tx_bam = ALIGN_FASTQ.out.tx_bam
        ch_star_log = ALIGN_FASTQ.out.star_log
        ch_bam_host = ALIGN_FASTQ.out.bam_host
        ch_bai_host = ALIGN_FASTQ.out.bai_host
        ch_star_log_host = ALIGN_FASTQ.out.star_log_host
        ch_bam_xeno = ALIGN_FASTQ.out.bam_xeno 
        ch_bai_xeno = ALIGN_FASTQ.out.bai_xeno
    
    }


    // save graft-only reads
    ch_graft_reads = Channel.empty()
    if(params.step == "mapping" && params.run_alignment && params.workflow == 'pdx'){
        if (params.run_arriba || params.only_filter_fastq){
            BAM_TO_FASTQ(
                ch_bam_xeno
            )
            ch_graft_reads = BAM_TO_FASTQ.out.fq
        }
    }


    /*
    * get gene records from gtf
    */
    if (params.step in [ "mapping", "expression_quantification" ] && (params.run_gene_count || params.run_dt)){
        if (params.gene_txt){
            gene_txt = file(params.gene_txt, checkIfExists:true)
        }else if (params.gtf){
            GTF2GENES(params.gtf)
            gene_txt = GTF2GENES.out.txt
        }
    }

    /*
    *   collect gene-level count matrix and call differential expression
    */
    ch_gene_rds = Channel.empty()
    ch_de = Channel.empty()
    if (params.step in [ "mapping", "expression_quantification", "differential_expression" ] && ( params.run_gene_count || params.run_de)){
        QUANT_GENES(
            params.workflow == 'pdx' ? ch_bam_xeno.ifEmpty([]) : ch_bam.ifEmpty([]),
            ch_counts.ifEmpty([]),
            params.step == "mapping" ? samplesheet : file("${params.outdir}/csv/align_fastq.csv", checkIfExists: true), 
            gene_txt
        )
        
        ch_gene_rds = QUANT_GENES.out.gene_rds
        ch_de = QUANT_GENES.out.de
        // update ch_counts if necessary
        if (params.step in [ "expression_quantification", "differential_expression" ] || params.run_featurecounts){
            ch_counts = QUANT_GENES.out.counts
        }

    }


    // Transcript-level analysis
    ch_salmon = Channel.empty()
    if (params.step in [ "mapping", "expression_quanification" ] && params.run_alignment ){
        MAP_TRANSCRIPTS(
            ch_tx_bam
        )
        ch_salmon = MAP_TRANSCRIPTS.out.salmon
        // [ [meta], val(out_prefix), path("${out_prefix}/") ]
    }

    /*
    *   Collect transcript-level counts and call differential transcripts
    */
    ch_tx_rds = Channel.empty()
    ch_dt = Channel.empty()
    if (params.step in [ "mapping", "expression_quantification", "differential_expression" ] && (params.run_tx_count || params.run_dt)){
        QUANT_TRANSCRIPTS(
            ch_salmon,
            params.step == "mapping" ? samplesheet : file("${params.outdir}/csv/align_fastq.csv", checkIfExists: true),
            gene_txt
        )

        ch_tx_rds = QUANT_TRANSCRIPTS.out.tx_rds
        ch_dt = QUANT_TRANSCRIPTS.out.dt
    }


    /*
    * Identify gene fusions
    */
    if ((params.step == "mapping" || params.gene_fusion) && params.run_arriba){
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
    * QC fastq
    */
    ch_fastqc = Channel.empty()
    ch_fastqc_trimmed = Channel.empty()
    if ((params.step == "mapping" || params.qc_fastq ) && params.run_qc_fastq){
        QC_FASTQ(
            ch_reads_raw,
            ch_reads_trimmed
        )
        ch_fastqc = QC_FASTQ.out.fastqc
        ch_fastqc_trimmed = QC_FASTQ.out.fastqc_trimmed
    }

    /*
    * QC alignment
    */
    ch_rnaseqc = Channel.empty()
    ch_rseqc = Channel.empty()
    ch_hs_metrics = Channel.empty()
    ch_bam_stat = Channel.empty()
    ch_bam_stat_host = Channel.empty()
    ch_bam_stat_xeno = Channel.empty()

    if ((params.step == "mapping" || params.qc_alignment ) && params.run_qc_alignment){
        // need to add a parser to align_fastq.csv
        QC_ALIGNMENT(
            ch_bam,
            ch_bai,
            ch_bam_host,
            ch_bai_host,
            ch_bam_xeno,
            ch_bai_xeno,
        )
        ch_rnaseqc = QC_ALIGNMENT.out.rnaseqc
        ch_rseqc = QC_ALIGNMENT.out.rnaseqc
        ch_hs_metrics = QC_ALIGNMENT.out.hs_metrics
        ch_bam_stat = QC_ALIGNMENT.out.bam_stat
        ch_bam_stat_host = QC_ALIGNMENT.out.bam_stat_host  
        ch_bam_stat_xeno = QC_ALIGNMENT.out.bam_stat_xeno
        
    }


    /*
    * MultiQC
    */
    ch_multiqc = Channel.empty()
    if ((params.step == "mapping" || params.qc_fastq || params.qc_alignment) && params.run_multiqc){
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
    * Generate a report
    */
    if ((params.step in [ "mapping", "bam_qc", "expression_quantification", "differential_expression" ] || params.qc_reads || params.qc_alignment || params.transcript_expression) && params.run_report){
        GENERATE_REPORT(
            params.workflow,
            params.step == "mapping" ? samplesheet : file("${params.outdir}/csv/align_fastq.csv", checkIfExists: true), 
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
