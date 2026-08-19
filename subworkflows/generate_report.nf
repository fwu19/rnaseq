/*
* Generate analysis report 
*/

include { MULTIQC } from '../modules/multiqc.nf'
include { RENDER_REPORT } from '../modules/render_report.nf'

workflow GENERATE_REPORT {
    take:
    samplesheet
    ch_star_counts
    ch_fc_summary
    ch_gene_expr
    ch_de
    ch_salmon
    ch_tx_expr
    ch_dt
    ch_arriba
    ch_fastqc
    ch_cutadapt_js
    ch_fastp_js
    ch_trimgalore_js
    ch_fastqc_trimmed
    ch_star_log
    ch_star_log_host
    ch_rnaseqc
    ch_rseqc
    ch_bam_stat
    ch_bam_stat_host
    ch_bam_stat_xeno
    ch_hs_metrics
    srcdir

    main: 
    ch_versions  = Channel.empty()
    ch_multiqc = Channel.empty()
    rmd = Channel.empty()

    if (params.run_multiqc){
        File star_csv = new File("${srcdir}/csv/gene_counts.star.csv")
        star_counts = star_csv.exists() ? Channel.fromPath( star_csv )
                .splitCsv(header: true)
                .map { it -> [ "${srcdir}/${it.gene_count}" ] }
                .collect()
                .ifEmpty([]) : []
        
        File fc_csv = new File("${srcdir}/csv/gene_counts.featurecounts.csv")
        fc_summary = fc_csv.exists() ? Channel.fromPath( fc_csv )
                .splitCsv(header: true)
                .map { it -> [ "${srcdir}/${it.gene_count_summary}" ] }
                .collect()
                .ifEmpty([]) : []

        File salmon_csv = new File("${srcdir}/csv/salmon.csv")
        salmon = salmon_csv.exists() ? Channel.fromPath( salmon_csv )
                .splitCsv(header: true)
                .map { it -> [ "${srcdir}/${it.salmon}" ] }
                .collect()
                .ifEmpty([]) : []
        
        File fastqc_dir = new File("${srcdir}/${params.fastqc_dir}")
        fastqc = fastqc_dir.exists() ? Channel.fromPath("${fastqc_dir}/*").flatten().collect().ifEmpty([]) : []
        
        File fastqc_trimmed_dir = new File("${srcdir}/${params.fastqc_trimmed_dir}")
        fastqc_trimmed = fastqc_trimmed_dir.exists() ? Channel.fromPath("${fastqc_trimmed_dir}/*").flatten().collect().ifEmpty([]) : []

        File cutadapt_dir = new File("${srcdir}/${params.cutadapt_dir}")
        cutadapt_js = cutadapt_dir.exists() ? Channel.fromPath("${cutadapt_dir}/*").flatten().collect().ifEmpty([]) : []

        File fastp_dir = new File("${srcdir}/${params.fastp_dir}")
        fastp_js = fastp_dir.exists() ? Channel.fromPath("${fastp_dir}/*").flatten().collect().ifEmpty([]) : []

        File trimgalore_dir = new File("${srcdir}/${params.trimgalore_dir}")
        trimgalore_js = trimgalore_dir.exists() ? Channel.fromPath("${trimgalore_dir}/*").flatten().collect().ifEmpty([]) : []

        File stat_dir = new File("${srcdir}/${params.samtools_dir}")       
        bam_stat = stat_dir.exists() ? Channel.fromPath("${stat_dir}/*").collect().ifEmpty([]) : []
        
        File stat_host_dir = new File("${srcdir}/${params.samtools_host_dir}")
        bam_stat_host = stat_host_dir.exists() ? Channel.fromPath("${stat_host_dir}/*").collect().ifEmpty([]) : []
        
        File stat_xeno_dir = new File("${srcdir}/${params.samtools_xeno_dir}")  
        bam_stat_xeno = stat_xeno_dir.exists() ? Channel.fromPath("${stat_xeno_dir}/*").collect().ifEmpty([]) : []
            
        File star_log_dir = new File("${srcdir}/${params.star_log_dir}")
        star_log = star_log_dir.exists() ?  Channel.fromPath("${star_log_dir}/**.final.out").collect().ifEmpty([]) : []

        File star_log_host_dir = new File("${srcdir}/${params.star_log_host_dir}")
        star_log_host = star_log_host_dir.exists() ?  Channel.fromPath("${star_log_host_dir}/*").collect().ifEmpty([]) : []

        File rnaseqc_dir = new File("${srcdir}/${params.rnaseqc_dir}")
        rnaseqc = rnaseqc_dir.exists() ? Channel.fromPath("${rnaseqc_dir}/*").collect().ifEmpty([]) : []

        File rseqc_dir = new File("${srcdir}/${params.rseqc_dir}")
        rseqc = rseqc_dir.exists() ? Channel.fromPath("${rseqc_dir}/*").collect().ifEmpty([]) : []  

        File gatk_dir = new File("${srcdir}/${params.gatk_dir}")
        hs_metrics = gatk_dir.exists() ? Channel.fromPath("${gatk_dir}/*").collect().ifEmpty([]) : []   

        MULTIQC(
            params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true ) : Channel.fromPath("$projectDir/assets/multiqc_config.yml"),
            params.run_fastqc ? ch_fastqc : fastqc,
            params.run_cut_adapt && params.trimmer == 'cutadapt' ? ch_cutadapt_js : cutadapt_js,  
            params.run_cut_adapt && params.trimmer == 'fastp' ? ch_fastp_js : fastp_js,  
            params.run_cut_adapt && params.trimmer == 'trimgalore' ? ch_trimgalore_js : trimgalore_js,  
            params.run_cut_adapt && params.run_fastqc ? ch_fastqc_trimmed : fastqc_trimmed,  
            params.run_alignment && params.aligner == 'star' ? ch_star_log : star_log, 
            params.run_alignment && params.aligner == 'star' ? ch_star_log_host : star_log_host, 
            params.run_rnaseqc ? ch_rnaseqc : rnaseqc,
            params.run_rseqc ? ch_rseqc : rseqc,  
            params.run_samtools ? ch_bam_stat : bam_stat,
            params.run_samtools ? ch_bam_stat_host : bam_stat_host,
            params.run_samtools ? ch_bam_stat_xeno : bam_stat_xeno,
            params.run_hs_metrics ? ch_hs_metrics : hs_metrics,
            params.workflow == 'pdx' ? [] : (params.run_alignment && params.aligner == 'star' ? ch_star_counts : star_counts),
            params.run_featurecounts ? ch_fc_summary : fc_summary,
            params.run_salmon ? ch_salmon : salmon
        )
        ch_multiqc = MULTIQC.out.data // multiqc_data/
        //multiqc.view()
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    
    }

    if (params.run_report){
        File multiqc_dir = new File("${srcdir}/${params.multiqc_dir}")
        multiqc = multiqc_dir.exists() ? Channel.fromPath("${multiqc_dir}/multiqc_data/", type: 'dir').ifEmpty([]) : []   
        //println "multiqc: ${multiqc.collect().flatten()}"
        
        File arriba_csv = new File("${srcdir}/csv/arriba.csv")
        arriba = arriba_csv.exists() ? Channel.fromPath( arriba_csv )
                .splitCsv(header: true)
                .map { it -> [ "${srcdir}/${it.arriba}" ] }
                .collect()
                .ifEmpty([]) : []

        File gene_expr_dir = new File("${srcdir}/${params.gene_expr_dir}/")
        gene_expr = gene_expr_dir.exists() ? Channel.fromPath("${gene_expr_dir}/all_samples.gene_raw_counts.txt").ifEmpty([]) : []
            
        File de_dir = new File("${srcdir}/${params.de_dir}/")
        de = de_dir.exists() ? Channel.fromPath("${de_dir}", type: 'dir').ifEmpty([]) : []

        File tx_expr_dir = new File("${srcdir}/${params.tx_expr_dir}/")
        tx_expr = tx_expr_dir.exists() ? Channel.fromPath("${srcdir}/${params.tx_expr_dir}/all_samples.transcript_raw_counts.txt").ifEmpty([]) : []

        File dt_dir = new File("${srcdir}/${params.dt_dir}/")
        dt = dt_dir.exists() ? Channel.fromPath("${srcdir}/${params.dt_dir}/", type: 'dir').ifEmpty([]) : []

        File gatk_dir = new File("${srcdir}/${params.gatk_dir}")
        gatk_hs_metrics = gatk_dir.exists() ? Channel.fromPath("${gatk_dir}/*").collect().ifEmpty([]) : []

        RENDER_REPORT(
            samplesheet, 
            params.run_multiqc ? ch_multiqc.ifEmpty([]) : multiqc, 
            params.run_gene_count ? ch_gene_expr : gene_expr, 
            params.run_de ? ch_de : de, 
            params.run_tx_count ? ch_tx_expr : tx_expr,  
            params.run_dt ? ch_dt : dt, 
            params.run_hs_metrics ? ch_hs_metrics : gatk_hs_metrics, 
            params.run_arriba ? ch_arriba : arriba,
            params.report_rmd ? Channel.fromPath(params.report_rmd, checkIfExists: true) : (params.report_dir ? Channel.fromPath("${params.report_dir}/rnaseq_${params.workflow}.Rmd", checkIfExists: true) : Channel.fromPath("$projectDir/assets/report/rnaseq_${params.workflow}.Rmd", checkIfExists: true))
        )
        rmd = RENDER_REPORT.out.rmd
        ch_versions = ch_versions.mix(RENDER_REPORT.out.versions)
    }

    emit:
    rmd
    versions = ch_versions

}
