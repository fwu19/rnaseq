#!/usr/bin/env nextflow 

include { GET_REFERENCE } from '../subworkflows/get_reference.nf'
include { GET_INPUT } from '../subworkflows/get_input.nf'
include { PROCESS_FASTQ } from '../subworkflows/process_fastq.nf'
include { ALIGN_FASTQ } from '../subworkflows/align_fastq.nf'
include { QUANT_GENES } from '../subworkflows/quant_genes.nf'
include { MAP_TRANSCRIPTS } from '../subworkflows/map_transcripts.nf'
include { QUANT_TRANSCRIPTS } from '../subworkflows/quant_transcripts.nf'
include { QC_FASTQ } from '../subworkflows/qc_fastq.nf'
include { QC_ALIGNMENT } from '../subworkflows/qc_alignment.nf'
include { GENERATE_REPORT } from '../subworkflows/generate_report.nf'

include { INFER_EXPERIMENT } from '../modules/infer_experiment.nf'
include { BAM_TO_FASTQ } from '../modules/bam_to_fastq.nf'
include { ARRIBA  } from '../modules/arriba.nf'
include { MULTIQC } from '../modules/multiqc.nf'


if ( params.genome_fa == null && params.run_build_index ){
    exit 1, "Need to provide a valid path with --genome_fa path/to/genome/fasta."
}

if ( params.gtf == null && (params.run_gene_count || params.run_tx_count)){
    exit 1, "Need to provide a valid path with --gtf path/to/genes/gtf."
}


workflow RNASEQ {
    ch_software_versions = Channel.empty()
    /*
    * Run reference check
    */
    gene_txt = Channel.empty()
    tx_bed = Channel.empty()
    tx_txt = Channel.empty()
    tx_fa = Channel.empty()
    collapsed_gtf = Channel.empty()

    if (params.run_reference_check) {
        GET_REFERENCE(
            params.genome_fa,
            params.gtf,
            params.aligner,
            params.aligner_index ?: params.aligner
        )
        index_dir = GET_REFERENCE.out.index_dir // match params.aligner
        gene_txt = GET_REFERENCE.out.gene_txt
        tx_txt = GET_REFERENCE.out.tx_txt
        tx_bed = GET_REFERENCE.out.tx_bed
        tx_fa = GET_REFERENCE.out.tx_fa
        collapsed_gtf = GET_REFERENCE.out.collapsed_gtf
        ch_software_versions = ch_software_versions.mix(GET_REFERENCE.out.versions)
    }


    /*
    * Run input check
    */
    ch_input = Channel.fromPath( params.input, checkIfExists: true )
    if (params.run_input_check){
        GET_INPUT(
            file(params.input, checkIfExists: true),
            file(params.input_dir, checkIfExists: true),
            file(params.metadata, checkIfExists: true)
        )
        samplesheet = GET_INPUT.out.samplesheet
        fq = GET_INPUT.out.fq

    }


    /*
    * cat/split fastq files
    */
    ch_reads = Channel.empty()
    ch_reads_trimmed = Channel.empty()
    if (params.run_process_fastq){
        PROCESS_FASTQ(
            samplesheet,
            fq,
            params.cat_fastq,
            params.trimmer
        )
        ch_reads = PROCESS_FASTQ.out.reads
        ch_reads_trimmed = PROCESS_FASTQ.out.reads_trimmed
        ch_cutadapt_js = PROCESS_FASTQ.out.cutadapt_js
        ch_fastp_js = PROCESS_FASTQ.out.fastp_js
        ch_software_versions = ch_software_versions.mix(PROCESS_FASTQ.out.versions)
    }


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

    if (params.gene_level){
        if (params.run_alignment){
        ALIGN_FASTQ(
            params.run_cut_adapt ? ch_reads_trimmed : ch_reads,
            params.run_split_fastq, // if true, dont run xenofilteR
            index_dir,
            params.workflow
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
        ch_software_versions = ch_software_versions.mix(ALIGN_FASTQ.out.versions)
    

        // save graft-only reads
        ch_graft_reads = Channel.empty()
        if(params.workflow == 'pdx' && (params.run_arriba || params.only_filter_fastq)){
            BAM_TO_FASTQ(
                ch_bam_xeno
            )
            ch_graft_reads = BAM_TO_FASTQ.out.fq
            ch_software_versions = ch_software_versions.mix(BAM_TO_FASTQ.out.versions)
        }

        } else {
        if (params.run_quant_genes && params.run_gene_count){
            ch_bam = file("${params.outdir}/csv/mapped.csv", checkIfExists: true)
                .splitCsv(header: true)
                .map { it -> [ [it.id, it.sample_group], it.id , path("${params.outdir}/${it.bam}") ] }

            ch_bai = file("${params.outdir}/csv/mapped.csv", checkIfExists: true)
                .splitCsv(header: true)
                .map { it -> [ [ it.id, it.sample_group ], it.id, path("${params.outdir}/${it.bai}") ] }

            if (params.workflow == 'pdx'){
                ch_bam_xeno = ch_bam
                ch_bai_xeno = ch_bai
            }

            ch_counts = file("${params.outdir}/csv/gene_counts.${params.aligner.toUpperCase()}.csv")
                    .splitCsv(header: true)
                    .map { it -> [ [ it.id, it.sample_group ], it.id, path("${params.outdir}/${it.gene_count}") ] }

        }

        
        if (params.run_quant_transcripts && params.run_tx_count ){
            ch_tx_bam = file("${params.outdir}/csv/mapped.csv", checkIfExists: true)
                .splitCsv(header: true)
                .map { it -> [ [ it.id, it.sample_group ], it.id, path("${params.outdir}/${it.tx_bam}") ] }
        }
        }
    }
    /*
    * Infer strandedness and read type
    */
    infer_experiment = Channel.empty()
    if (params.run_infer_experiment ){
        INFER_EXPERIMENT(
                ch_bam
                .map{ it -> [ [ it[0], it[1] ], it[2] ]}
                .join (
                    ch_bai
                    .map{ it -> [ [ it[0], it[1] ], it[2] ]}
                )
                .map { it -> [ it[0][0], it[0][1], it[1], it[2] ]}
                .first(),
                tx_bed
        )
        infer_experiment = INFER_EXPERIMENT.out.csv
        ch_software_versions = ch_software_versions.mix(INFER_EXPERIMENT.out.versions)
        
    }else {
        if (params.run_quant_genes && params.run_gene_count){
            infer_experiment = Channel.fromPath("${params.outdir}/csv/infer_experiment.csv")
        }
    }    
    

    /*
    *   collect gene-level count matrix and call differential expression
    */
    ch_gene_expr = Channel.empty()
    ch_de = Channel.empty()
    if (params.gene_level){
        if ( params.run_quant_genes ){
        QUANT_GENES(
            samplesheet,
            params.workflow == 'pdx' ? ch_bam_xeno : ch_bam,
            params.workflow == 'pdx' ? ch_bai_xeno : ch_bai,
            ch_counts.ifEmpty([]),
            gene_txt,
            tx_bed,
            infer_experiment
        )
        
        ch_counts = QUANT_GENES.out.counts // if params.run_featurecounts, ch_counts is updated
        ch_gene_expr = QUANT_GENES.out.expr 
        ch_de = QUANT_GENES.out.de
        ch_software_versions = ch_software_versions.mix(QUANT_GENES.out.versions)

        } 
    }


    /*
    * Transcript-level analysis
    */
    ch_salmon = Channel.empty()
    ch_tx_expr = Channel.empty()
    ch_dt = Channel.empty()
    if (params.transcript_expression ){
        /* Map to transcripts */
        if (params.run_alignment){
            MAP_TRANSCRIPTS(
                ch_tx_bam,
                tx_fa
            )
            ch_salmon = MAP_TRANSCRIPTS.out.salmon
            // [ [meta], val(out_prefix), path("${out_prefix}/") ]
            ch_software_versions = ch_software_versions.mix(MAP_TRANSCRIPTS.out.versions)
        }else{
            // generate ch_salmon by parse csv/salmon.csv
        }

        /*  Collect transcript-level counts and call differential transcripts */
        if ((params.run_quant_transcripts)){
            QUANT_TRANSCRIPTS(
                ch_salmon,
                samplesheet,
                gene_txt,
                tx_txt
            )

            ch_tx_expr = QUANT_TRANSCRIPTS.out.tx_expr
            ch_dt = QUANT_TRANSCRIPTS.out.dt
            ch_software_versions = ch_software_versions.mix(QUANT_TRANSCRIPTS.out.versions)
        }
    }


    /*
    * Identify gene fusions
    */
    if (params.gene_fusion){
        // can build a subworkflow for gene_fusion analysis
        if (params.run_arriba){
            ARRIBA(
            params.workflow == 'pdx' ? ch_graft_reads : (params.run_cut_adapt ? ch_reads_trimmed : ch_reads), 
            params.genome, 
            index_dir, 
            params.gtf,
            params.genome_fa,
            params.blacklist,
            params.known_fusions,
            params.protein_domains
            )
            ch_software_versions = ch_software_versions.mix(ARRIBA.out.versions)
        }
    }


    /*
    * QC fastq
    */
    ch_fastqc = Channel.empty()
    ch_fastqc_trimmed = Channel.empty()
    if (params.run_qc_fastq){
        QC_FASTQ(
            ch_reads,
            ch_reads_trimmed
        )
        ch_fastqc = QC_FASTQ.out.fastqc
        ch_fastqc_trimmed = QC_FASTQ.out.fastqc_trimmed
        ch_software_versions = ch_software_versions.mix(QC_FASTQ.out.versions)
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

    if (params.run_qc_alignment){
        // need to add a parser to align_fastq.csv
        QC_ALIGNMENT(
            ch_bam,
            ch_bai,
            ch_bam_host,
            ch_bai_host,
            ch_bam_xeno,
            ch_bai_xeno,
            collapsed_gtf,
            tx_bed,
            infer_experiment
        )
        ch_rnaseqc = QC_ALIGNMENT.out.rnaseqc
        ch_rseqc = QC_ALIGNMENT.out.rseqc
        ch_hs_metrics = QC_ALIGNMENT.out.hs_metrics
        ch_bam_stat = QC_ALIGNMENT.out.bam_stat
        ch_bam_stat_host = QC_ALIGNMENT.out.bam_stat_host  
        ch_bam_stat_xeno = QC_ALIGNMENT.out.bam_stat_xeno
        ch_software_versions = ch_software_versions.mix(QC_ALIGNMENT.out.versions)
        
    }


    /*
    * MultiQC
    */
    ch_multiqc = Channel.empty()
    if (params.run_multiqc){
        MULTIQC(
            params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true ) : Channel.fromPath("$projectDir/assets/multiqc_config.yml"),
            ch_fastqc.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_cutadapt_js.map{it[2]}.flatten().collect().ifEmpty([]),  
            ch_fastp_js.map{it[2]}.flatten().collect().ifEmpty([]),  
            ch_fastqc_trimmed.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_rseqc.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_rnaseqc.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_hs_metrics.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_star_log.map{it[2]}.flatten().collect().ifEmpty([]), 
            ch_counts.map{it[2]}.flatten().collect().ifEmpty([]),
            ch_salmon.map{it[2]}.flatten().collect().ifEmpty([])
        )
 
        ch_multiqc = MULTIQC.out.data // multiqc_data/
        //ch_multiqc.view()
        ch_software_versions = ch_software_versions.mix(MULTIQC.out.versions)
    
    }

    /*
    * Generate a report
    */
    if (params.run_report){
        GENERATE_REPORT(
            samplesheet, 
            ch_multiqc.ifEmpty([]), 
            ch_gene_expr.ifEmpty([]), 
            ch_de.ifEmpty([]), 
            ch_tx_expr.ifEmpty([]), 
            ch_dt.ifEmpty([]), 
            ch_hs_metrics.map{it[1]}.collect().ifEmpty([]), 
            params.report_dir ? Channel.fromPath(params.report_dir, type: 'dir', checkIfExists: true) : Channel.fromPath("$projectDir/assets/report/", type: 'dir', checkIfExists: true)
        )
        ch_software_versions = ch_software_versions.mix(GENERATE_REPORT.out.versions)
        
    }

    /*
    * Collect software versions
    */
    if (!params.update){
        ch_software_versions
            .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'software_versions.yml', sort: true, newLine: true)
    }


}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

import groovy.json.JsonOutput
workflow.onComplete {
    def jsonStr = JsonOutput.toJson(params)
    def pretty  = JsonOutput.prettyPrint(jsonStr)
    file("${params.outdir ?: '.'}/pipeline_info/params.json").text = pretty

}
