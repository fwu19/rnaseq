/*
* QC alignment
*/

include { RNASEQC  } from '../modules/rnaseqc.nf'
include { RSEQC  } from '../modules/rseqc.nf'
include { HS_METRICS  } from '../modules/hs_metrics'
include { SAMTOOLS_VIEW} from '../modules/samtools_view.nf'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_XENO } from '../modules/samtools_view.nf'
include { RNASEQC as RNASEQC_HOST  } from '../modules/rnaseqc.nf'
include { RSEQC as RSEQC_HOST } from '../modules/rseqc.nf'
include { HS_METRICS as HS_METRICS_HOST  } from '../modules/hs_metrics'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_HOST } from '../modules/samtools_view.nf'


workflow QC_ALIGNMENT {
    take:
    ch_bam_bai
    ch_bam_bai_host
    ch_bam_bai_xeno
    collapsed_gtf
    tx_bed
    collapsed_gtf_host
    tx_bed_host
    infer_experiment
    srcdir
    
    main: 
    ch_rnaseqc = Channel.empty()
    ch_rseqc = Channel.empty()
    ch_hs_metrics = Channel.empty()
    ch_rnaseqc_host = Channel.empty()
    ch_rseqc_host = Channel.empty()
    ch_hs_metrics_host = Channel.empty()
    ch_bam_stat = Channel.empty()
    ch_bam_stat_host = Channel.empty()
    ch_bam_stat_xeno = Channel.empty()
    ch_versions = Channel.empty()
    
    if (!params.run_alignment){
        ch_bam_bai = Channel.fromPath("${srcdir}/csv/map2genome.${params.aligner}.${params.genome}.csv")
            .splitCsv(header: true)
            .map { it -> [ [ it ], it.id , "${srcdir}/${it.bam}", "${srcdir}/${it.bai}" ] }

        if (params.workflow == 'pdx'){
            ch_bam_bai = Channel.fromPath("${srcdir}/csv/map2genome.${params.aligner}.${params.genome}.csv")
                .splitCsv(header: true)
                .map { it -> [ [ it ], it.id , "${srcdir}/${it.graft_bam}", "${srcdir}/${it.graft_bai}" ] }

            ch_bam_bai_host = Channel.fromPath("${srcdir}/csv/map2genome.${params.aligner}.${params.genome_host}.csv")
                .splitCsv(header: true)
                .map { it -> [ [ it ], it.id , "${srcdir}/${it.host_bam}", "${srcdir}/${it.host_qbai}" ] }

            ch_bam_bai_xeno = Channel.fromPath("${srcdir}/csv/map2genome.${params.aligner}.${params.genome}.xenofilteR_filtered.csv")
                .splitCsv(header: true)
                .map { it -> [ [ it ], it.id , "${srcdir}/${it.bam}", "${srcdir}/${it.bai}" ] }
        }
    }

    /*
    * RNASeQC
    */
    if (params.run_rnaseqc){
        RNASEQC(
            params.run_xenofilter ? ch_bam_bai_xeno : ch_bam_bai,
            collapsed_gtf,
            params.run_infer_experiment ? infer_experiment : file("${srcdir}/csv/infer_experiment.csv")
        )            
        ch_rnaseqc = RNASEQC.out.qc
        ch_versions = ch_versions.mix(RNASEQC.out.versions.first())

        RNASEQC_HOST(
            params.run_xenofilter ? Channel.empty() : ch_bam_bai_host,
            collapsed_gtf_host,
            params.run_infer_experiment ? infer_experiment : file("${srcdir}/csv/infer_experiment.csv")
        )            
        ch_rnaseqc_host = RNASEQC_HOST.out.qc
        ch_versions = ch_versions.mix(RNASEQC_HOST.out.versions.first())
        
    }

    /*
    * RSeQC
    */
    if (params.run_rseqc){
        RSEQC(
            params.run_xenofilter ? ch_bam_bai_xeno : ch_bam_bai,
            tx_bed
        )
        ch_rseqc = RSEQC.out.qc
        ch_versions = ch_versions.mix(RSEQC.out.versions.first())

        RSEQC_HOST(
            params.run_xenofilter ? Channel.empty() : ch_bam_bai_host,
            tx_bed_host
        )
        ch_rseqc_host = RSEQC_HOST.out.qc
        ch_versions = ch_versions.mix(RSEQC_HOST.out.versions.first())

    }

    /*
    * GATK collect_hs_metrics
    */
    if (params.run_hs_metrics){
            HS_METRICS(
                params.run_xenofilter ? ch_bam_bai_xeno : ch_bam_bai,
                params.genome_fa,
                params.target_region
            )
            ch_hs_metrics = HS_METRICS.out.qc
            ch_versions = ch_versions.mix(HS_METRICS.out.versions.first())

            HS_METRICS_HOST(
                params.run_xenofilter ? Channel.empty() : ch_bam_bai_host,
                params.genome_fa,
                params.target_region
            )
            ch_hs_metrics_host = HS_METRICS_HOST.out.qc
            ch_versions = ch_versions.mix(HS_METRICS_HOST.out.versions.first())
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
            ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

            SAMTOOLS_VIEW_HOST(
                ch_bam_bai_host
            )
            ch_bam_stat_host = SAMTOOLS_VIEW_HOST.out.data
            ch_versions = ch_versions.mix(SAMTOOLS_VIEW_HOST.out.versions.first())

            SAMTOOLS_VIEW_XENO(
                ch_bam_bai_xeno
            )
            ch_bam_stat_xeno = SAMTOOLS_VIEW_XENO.out.data
            ch_versions = ch_versions.mix(SAMTOOLS_VIEW_XENO.out.versions.first())
    }


    emit:
    rnaseqc = ch_rnaseqc
    rseqc = ch_rseqc
    hs_metrics = ch_hs_metrics
    rnaseqc_host = ch_rnaseqc_host
    rseqc_host = ch_rseqc_host
    hs_metrics_host = ch_hs_metrics_host
    bam_stat = ch_bam_stat
    bam_stat_host  = ch_bam_stat_host
    bam_stat_xeno = ch_bam_stat_xeno
    versions = ch_versions

}
