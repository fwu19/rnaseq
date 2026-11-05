/*
* Generate analysis report 
*/

include { MULTIQC } from '../modules/multiqc.nf'
include { RENDER_REPORT } from '../modules/render_report.nf'

workflow GENERATE_REPORT {
    take:
    samplesheet
    ch_fastqc
    ch_cutadapt_js
    ch_fastp_js
    ch_trimgalore_js
    ch_fastqc_trimmed
    ch_xengsort_log
    ch_star_log
    ch_star_log_host
    ch_rnaseqc
    ch_rnaseqc_host
    ch_rseqc
    ch_rseqc_host
    ch_bam_stat
    ch_bam_stat_host
    ch_bam_stat_xeno
    ch_hs_metrics
    ch_star_counts
    ch_star_counts_host
    ch_fc_summary
    ch_fc_summary_host
    ch_gene_expr
    ch_gene_expr_host
    ch_de
    ch_de_host
    ch_salmon
    ch_salmon_host
    ch_tx_expr
    ch_tx_expr_host
    ch_dt
    ch_dt_host
    ch_arriba

    srcdir

    main: 
    ch_versions  = Channel.empty()
    ch_multiqc = Channel.empty()
    rmd = Channel.empty()

    if (params.run_multiqc){
        
        File fastqc_dir = new File("${srcdir}/${params.fastqc_dir}")
        fastqc = fastqc_dir.exists() ? Channel.fromPath("${fastqc_dir}/*") : Channel.empty()
        
        File fastqc_trimmed_dir = new File("${srcdir}/${params.fastqc_trimmed_dir}")
        fastqc_trimmed = fastqc_trimmed_dir.exists() ? Channel.fromPath("${fastqc_trimmed_dir}/*") : Channel.empty()

        File cutadapt_dir = new File("${srcdir}/${params.cutadapt_dir}")
        cutadapt_js = cutadapt_dir.exists() ? Channel.fromPath("${cutadapt_dir}/*") : Channel.empty()

        File fastp_dir = new File("${srcdir}/${params.fastp_dir}")
        fastp_js = fastp_dir.exists() ? Channel.fromPath("${fastp_dir}/*") : Channel.empty()

        File trimgalore_dir = new File("${srcdir}/${params.trimgalore_dir}")
        trimgalore_js = trimgalore_dir.exists() ? Channel.fromPath("${trimgalore_dir}/*") : Channel.empty()

        File xengsort_dir = new File("${srcdir}/${params.xengsort_dir}")
        xengsort_log = xengsort_dir.exists() ? Channel.fromPath("${xengsort_dir}/*") : Channel.empty()
        println "xengsort_log: ${xengsort_log.collect().flatten()}" 

        File stat_dir = new File("${srcdir}/${params.samtools_dir}")       
        bam_stat = stat_dir.exists() ? Channel.fromPath("${stat_dir}/*") : Channel.empty()
        
        File stat_host_dir = new File("${srcdir}/${params.samtools_host_dir}")
        bam_stat_host = stat_host_dir.exists() ? Channel.fromPath("${stat_host_dir}/*") : Channel.empty()
        
        File stat_xeno_dir = new File("${srcdir}/${params.samtools_xeno_dir}")  
        bam_stat_xeno = stat_xeno_dir.exists() ? Channel.fromPath("${stat_xeno_dir}/*") : Channel.empty()
            
        File star_log_dir = new File("${srcdir}/${params.star_log_dir}")
        star_log = star_log_dir.exists() ?  Channel.fromPath("${star_log_dir}/**.final.out") : Channel.empty()

        File star_log_host_dir = new File("${srcdir}/${params.star_log_host_dir}")
        star_log_host = star_log_host_dir.exists() ?  Channel.fromPath("${star_log_host_dir}/*") : Channel.empty()

        File rnaseqc_dir = new File("${srcdir}/${params.rnaseqc_dir}")
        rnaseqc = rnaseqc_dir.exists() ? Channel.fromPath("${rnaseqc_dir}/*") : Channel.empty()

        File rnaseqc_host_dir = new File("${srcdir}/${params.rnaseqc_host_dir}")
        rnaseqc_host = rnaseqc_host_dir.exists() ? Channel.fromPath("${rnaseqc_host_dir}/*") : Channel.empty()

        File rseqc_dir = new File("${srcdir}/${params.rseqc_dir}")
        rseqc = rseqc_dir.exists() ? Channel.fromPath("${rseqc_dir}/*") : Channel.empty()  

        File rseqc_host_dir = new File("${srcdir}/${params.rseqc_host_dir}")
        rseqc_host = rseqc_host_dir.exists() ? Channel.fromPath("${rseqc_host_dir}/*") : Channel.empty()  

        File gatk_dir = new File("${srcdir}/${params.gatk_dir}")
        hs_metrics = gatk_dir.exists() ? Channel.fromPath("${gatk_dir}/*") : Channel.empty()   

        File gatk_host_dir = new File("${srcdir}/${params.gatk_host_dir}")
        hs_metrics_host = gatk_host_dir.exists() ? Channel.fromPath("${gatk_host_dir}/*") : Channel.empty()   

        File star_csv = new File("${srcdir}/csv/gene_counts.star.${params.genome}.csv")
        star_counts = star_csv.exists() ? Channel.fromPath( star_csv )
                .splitCsv(header: true)
                .map { 
                    it -> [ "${srcdir}/${it.gene_count}" ] 
                } : Channel.empty()
        
        File star_host_csv = new File("${srcdir}/csv/gene_counts.star.${params.genome_host}.csv")
        star_counts_host = star_host_csv.exists() ? Channel.fromPath( star_host_csv )
                .splitCsv(header: true)
                .map { 
                    it -> [ "${srcdir}/${it.gene_count}" ] 
                } : Channel.empty()

        File fc_csv = new File("${srcdir}/csv/gene_counts.featurecounts.${params.genome}.csv")
        fc_summary = fc_csv.exists() ? Channel.fromPath( fc_csv )
                .splitCsv(header: true)
                .map { 
                    it -> ["${srcdir}/${it.gene_count}.summary" ] 
                } : Channel.empty()

        File fc_host_csv = new File("${srcdir}/csv/gene_counts.featurecounts.${params.genome_host}.csv")
        fc_summary_host = fc_host_csv.exists() ? Channel.fromPath( fc_host_csv )
                .splitCsv(header: true)
                .map { 
                    it -> [ "${srcdir}/${it.host_gene_count}.summary" ] 
                } : Channel.empty()

        File salmon_csv = new File("${srcdir}/csv/salmon.${params.genome}.csv")
        salmon = salmon_csv.exists() ? Channel.fromPath( salmon_csv )
                .splitCsv(header: true)
                .map { 
                    it -> ["${srcdir}/${it.salmon}" ]
                } : Channel.empty()

        File salmon_host_csv = new File("${srcdir}/csv/salmon.${params.genome_host}.csv")
        salmon_host = salmon_host_csv.exists() ? Channel.fromPath( salmon_host_csv )
                .splitCsv(header: true)
                .map { 
                    it -> [ "${srcdir}/${it.salmon}" ] 
                } : Channel.empty()


        MULTIQC(
            params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true ) : Channel.fromPath("$projectDir/assets/multiqc_config.yml"),
            (params.run_fastqc ? ch_fastqc.flatten().collect() : fastqc.collect()).ifEmpty([]),
            (params.run_cut_adapt && params.trimmer == 'cutadapt' ? ch_cutadapt_js.flatten().collect() : cutadapt_js.flatten().collect()).ifEmpty([]),  
            (params.run_cut_adapt && params.trimmer == 'fastp' ? ch_fastp_js.flatten().collect() : fastp_js.flatten().collect()).ifEmpty([]),  
            (params.run_cut_adapt && params.trimmer == 'trimgalore' ? ch_trimgalore_js.flatten().collect() : trimgalore_js.flatten().collect()).ifEmpty([]),  
            (params.run_cut_adapt && params.run_fastqc ? ch_fastqc_trimmed.flatten().collect() : fastqc_trimmed.flatten().collect()).ifEmpty([]),  
            (params.workflow == 'pdx' && params.run_xengsort ? ch_xengsort_log.collect() : xengsort_log.collect()).ifEmpty([]),
            (params.run_alignment && params.aligner == 'star' ? ch_star_log.flatten().collect() : star_log.flatten().collect()).ifEmpty([]), 
            (params.run_alignment && params.aligner == 'star' ? ch_star_log_host.flatten().collect() : star_log_host.flatten().collect()).ifEmpty([]), 
            (params.run_rnaseqc ? ch_rnaseqc.flatten().collect() : rnaseqc.flatten().collect()).ifEmpty([]),
            (params.run_rnaseqc ? ch_rnaseqc_host.flatten().collect() : rnaseqc_host.flatten().collect()).ifEmpty([]),
            (params.run_rseqc ? ch_rseqc.flatten().collect() : rseqc.flatten().collect()).ifEmpty([]),  
            (params.run_rseqc ? ch_rseqc_host.flatten().collect() : rseqc_host.flatten().collect()).ifEmpty([]),  
            (params.run_samtools ? ch_bam_stat.flatten().collect() : bam_stat.flatten().collect()).ifEmpty([]),
            (params.run_samtools ? ch_bam_stat_host.flatten().collect() : bam_stat_host.flatten().collect()).ifEmpty([]),
            (params.run_samtools ? ch_bam_stat_xeno.flatten().collect() : bam_stat_xeno.flatten().collect()).ifEmpty([]),
            (params.run_hs_metrics ? ch_hs_metrics.flatten().collect() : hs_metrics.flatten().collect()).ifEmpty([]),
            params.workflow == 'pdx' && params.run_xenofilter ? [] : (params.run_alignment && params.aligner == 'star' ? ch_star_counts.map{it[2]}.flatten().collect() : star_counts.collect()).ifEmpty([]),
            params.workflow == 'pdx' && params.run_xenofilter ? [] : (params.run_alignment && params.aligner == 'star' ? ch_star_counts_host.map{it[2]}.flatten().collect() : star_counts_host.collect()).ifEmpty([]),
            (params.run_featurecounts ? ch_fc_summary.map{it[2]}.flatten().collect() : fc_summary.collect()).ifEmpty([]),
            (params.run_featurecounts_host ? ch_fc_summary_host.map{it[2]}.flatten().collect() : fc_summary_host.collect()).ifEmpty([]),
            (params.run_salmon ? ch_salmon.map{it[2]}.flatten().collect() : salmon.collect()).ifEmpty([]),
            (params.run_salmon_host ? ch_salmon_host.map{it[2]}.flatten().collect() : salmon_host.collect()).ifEmpty([])
        )
        ch_multiqc = MULTIQC.out.data // multiqc_data/
        //multiqc.view()
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    
    }

    if (params.run_report){
        File multiqc_dir = new File("${srcdir}/${params.multiqc_dir}")
        multiqc = multiqc_dir.exists() ? Channel.fromPath("${multiqc_dir}/multiqc_data/", type: 'dir') : Channel.empty()
        //println "multiqc: ${multiqc.collect().flatten()}"
        
        File arriba_csv = new File("${srcdir}/csv/arriba.${params.genome}.csv")
        arriba = arriba_csv.exists() ? Channel.fromPath( arriba_csv )
                .splitCsv(header: true)
                .map { it -> [ "${srcdir}/${it.arriba}" ] } : Channel.empty()

        File gene_expr_file = new File("${srcdir}/${params.gene_expr_dir}/all_samples.gene_raw_counts.txt")
        gene_expr = gene_expr_file.exists() ? Channel.fromPath("${srcdir}/${params.gene_expr_dir}/all_samples.gene_raw_counts.txt") : Channel.empty()

        File gene_expr_host_file = new File("${srcdir}/${params.gene_expr_host_dir}/all_samples.gene_raw_counts.txt")
        gene_expr_host = gene_expr_host_file.exists() ? Channel.fromPath("${srcdir}/${params.gene_expr_host_dir}/all_samples.gene_raw_counts.txt") : Channel.empty()
        
        File de_dir = new File("${srcdir}/${params.de_dir}/")
        de = de_dir.exists() ? Channel.fromPath("${de_dir}", type: 'dir') : Channel.empty()

        File de_host_dir = new File("${srcdir}/${params.de_host_dir}/")
        de_host = de_host_dir.exists() ? Channel.fromPath("${de_host_dir}", type: 'dir') : Channel.empty()

        File tx_expr_file = new File("${srcdir}/${params.tx_expr_dir}/all_samples.transcript_raw_counts.txt")
        tx_expr = tx_expr_file.exists() ? Channel.fromPath("${srcdir}/${params.tx_expr_dir}/all_samples.transcript_raw_counts.txt") : Channel.empty()

        File tx_expr_host_file = new File("${srcdir}/${params.tx_expr_host_dir}/all_samples.transcript_raw_counts.txt")
        tx_expr_host = tx_expr_host_file.exists() ? Channel.fromPath("${srcdir}/${params.tx_expr_host_dir}/all_samples.transcript_raw_counts.txt") : Channel.empty()

        File dt_dir = new File("${srcdir}/${params.dt_dir}/")
        dt = dt_dir.exists() ? Channel.fromPath("${srcdir}/${params.dt_dir}/", type: 'dir') : Channel.empty()

        File dt_host_dir = new File("${srcdir}/${params.dt_host_dir}/")
        dt_host = dt_host_dir.exists() ? Channel.fromPath("${srcdir}/${params.dt_host_dir}/", type: 'dir') : Channel.empty()

        File gatk_dir = new File("${srcdir}/${params.gatk_dir}")
        hs_metrics = gatk_dir.exists() ? Channel.fromPath("${gatk_dir}/*") : Channel.empty()

        def report_rmd_path 
        if (params.report_rmd) {
            report_rmd_path = params.report_rmd
        } else if (params.report_dir) {
            report_rmd_path = "${params.report_dir}/rnaseq_${params.workflow}.Rmd"
        } else {
            report_rmd_path = "${projectDir}/assets/report/rnaseq_${params.workflow}.Rmd"
        }

        RENDER_REPORT(
            samplesheet, 
            (params.run_multiqc ? ch_multiqc.collect() : multiqc.collect()).ifEmpty([]), 
            (params.run_gene_count ? ch_gene_expr.collect() : gene_expr.collect()).ifEmpty([]), 
            (params.run_gene_count_host ? ch_gene_expr_host.collect() : gene_expr_host.collect()).ifEmpty([]), 
            (params.run_de ? ch_de.collect() : de.collect()).ifEmpty([]), 
            (params.run_de_host ? ch_de_host.collect() : de_host.collect()).ifEmpty([]), 
            (params.run_tx_count ? ch_tx_expr.collect() : tx_expr.collect()).ifEmpty([]), 
            (params.run_tx_count_host ? ch_tx_expr_host.collect() : tx_expr_host.collect()).ifEmpty([]),   
            (params.run_dt ? ch_dt.collect() : dt.collect()).ifEmpty([]), 
            (params.run_dt_host ? ch_dt_host.collect() : dt_host.collect()).ifEmpty([]), 
            (params.run_hs_metrics ? ch_hs_metrics.flatten().collect() : hs_metrics.collect()).ifEmpty([]), 
            (params.run_arriba ? ch_arriba.map{it[2]}.flatten().collect() : arriba.collect()).ifEmpty([]),
            channel.fromPath(report_rmd_path, checkIfExists: true)
        )
        rmd = RENDER_REPORT.out.rmd
        ch_versions = ch_versions.mix(RENDER_REPORT.out.versions)
    }

    emit:
    rmd
    versions = ch_versions

}
