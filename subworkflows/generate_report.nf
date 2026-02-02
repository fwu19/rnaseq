/*
* Generate analysis report 
*/

include { MULTIQC } from '../modules/multiqc.nf'
include { RENDER_REPORT } from '../modules/render_report.nf'

workflow GENERATE_REPORT {
    take:
    samplesheet
    ch_star_counts
    ch_fc_counts
    ch_salmon
    ch_fastqc
    ch_cutadapt_js
    ch_fastp_js
    ch_fastqc_trimmed
    ch_star_log
    ch_rnaseqc
    ch_rseqc
    ch_bam_stat
    ch_bam_stat_host
    ch_bam_stat_xeno
    ch_hs_metrics
    ch_gene_expr
    ch_de
    ch_tx_expr
    ch_dt
    srcdir

    main: 
    ch_versions  = Channel.empty()
    ch_multiqc = Channel.empty()

    if (params.run_multiqc){

        if (!params.run_gene_count){
            ch_fc_counts = Channel.fromPath( "${srcdir}/csv/gene_counts.featurecounts.csv" )
                .splitCsv(header: true)
                .map { it -> [ "${srcdir}/${it.gene_count}" ] }
                .collect()
        }

        if (!params.run_alignment && params.aligner == 'star'){
            ch_star_counts = Channel.fromPath( "${srcdir}/csv/gene_counts.star.csv" )
                .splitCsv(header: true)
                .map { it -> [ "${srcdir}/${it.gene_count}" ] }
                .collect()
        }

        if (!params.run_salmon){
            ch_salmon = Channel.fromPath( "${srcdir}/csv/salmon.csv" )
                .splitCsv(header: true)
                .map { it -> [ "${srcdir}/${it.salmon}" ] }
                .collect()
        }

        if (!params.run_samtools){
            if (params.workflow == 'pdx'){
                ch_bam_stat = Channel.fromPath("${srcdir}/QC/samtools/${params.genome}/_unfiltered/*").collect()
                ch_bam_stat_host = Channel.fromPath("${srcdir}/QC/samtools/${params.genome_host}/_unfiltered/*").collect()
                ch_bam_stat_xeno = Channel.fromPath("${srcdir}/QC/samtools/${params.genome}/_filtered/*").collect()
            } else {
                ch_bam_stat = Channel.fromPath("${srcdir}/QC/samtools/*").collect()
            }
        }

        MULTIQC(
            params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true ) : Channel.fromPath("$projectDir/assets/multiqc_config.yml"),
            params.run_fastqc ? ch_fastqc.ifEmpty([]) : Channel.fromPath("${srcdir}/QC/fastqc/*").flatten().collect().ifEmpty([]),
            params.run_cut_adapt ? ch_cutadapt_js.ifEmpty([]) : Channel.fromPath("${srcdir}/QC/cutadapt/*").flatten().collect().ifEmpty([]),  
            params.run_cut_adapt ? ch_fastp_js.ifEmpty([]) : Channel.fromPath("${srcdir}/QC/fastp/*").flatten().collect().ifEmpty([]),  
            params.run_fastqc ? ch_fastqc_trimmed.ifEmpty([]) : Channel.fromPath("${srcdir}/QC/fastqc_trimmed/*").flatten().collect().ifEmpty([]),  
            params.run_alignment ? ch_star_log.ifEmpty([]) : Channel.fromPath("${srcdir}/qc/star/*.Log.final.out").collect().ifEmpty([]), 
            params.run_rnaseqc ? ch_rnaseqc.ifEmpty([]) : Channel.fromPath("${srcdir}/QC/rnaseqc/*").collect().ifEmpty([]),
            params.run_rseqc ? ch_rseqc.ifEmpty([]) : Channel.fromPath("${srcdir}/QC/rseqc/*").collect().ifEmpty([]),  
            ch_bam_stat.ifEmpty([]),
            ch_bam_stat_host.ifEmpty([]),
            ch_bam_stat_xeno.ifEmpty([]),
            params.run_hs_metrics ? ch_hs_metrics : Channel.fromPath("${srcdir}/QC/gatk/*").collect().ifEmpty([]),
            ch_star_counts.ifEmpty([]),
            ch_fc_counts.ifEmpty([]),
            ch_salmon.ifEmpty([])
        )
        ch_multiqc = MULTIQC.out.data // multiqc_data/
        //ch_multiqc.view()
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    
    }

    if (params.run_report){
        RENDER_REPORT(
        samplesheet, 
        params.run_multiqc ? ch_multiqc.ifEmpty([]) : Channel.fromPath("${srcdir}/multiqc/multiqc_data/", type: 'dir').ifEmpty([]), 
        params.run_gene_count ? ch_gene_expr : Channel.fromPath("${srcdir}/expression_quantification/all_samples.gene_raw_counts.txt").ifEmpty([]), 
        params.run_de ? ch_de : Channel.fromPath("${srcdir}/differential_genes/", type: 'dir').ifEmpty([]), 
        params.run_tx_count ? ch_tx_expr : Channel.fromPath("${srcdir}/expression_quantification/all_samples.transcript_raw_counts.txt").ifEmpty([]),  
        params.run_dt ? ch_dt : Channel.fromPath("${srcdir}/differential_transcripts/", type: 'dir').ifEmpty([]), 
        params.run_hs_metrics ? ch_hs_metrics : Channel.fromPath("${srcdir}/QC/gatk/*", type: 'dir').collect().ifEmpty([]), 
        params.report_dir ? Channel.fromPath(params.report_dir, type: 'dir', checkIfExists: true) : Channel.fromPath("$projectDir/assets/report/", type: 'dir', checkIfExists: true)
        )
        ch_versions = ch_versions.mix(RENDER_REPORT.out.versions)
    }

    emit:
    versions = ch_versions

}
