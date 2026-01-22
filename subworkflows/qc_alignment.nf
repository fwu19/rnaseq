/*
* QC alignment
*/

include { RNASEQC  } from '../modules/rnaseqc.nf'
include { RSEQC  } from '../modules/rseqc.nf'
include { HS_METRICS  } from '../modules/hs_metrics'
include { SAMTOOLS_VIEW} from '../modules/samtools_view.nf'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_HOST } from '../modules/samtools_view.nf'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_XENO } from '../modules/samtools_view.nf'


workflow QC_ALIGNMENT {
    take:
    bam
    bai
    bam_host
    bai_host
    bam_xeno
    bai_xeno
    collapsed_gtf
    tx_bed
    experiment
    
    main: 
    ch_rnaseqc = Channel.empty()
    ch_rseqc = Channel.empty()
    ch_hs_metrics = Channel.empty()
    ch_bam_stat = Channel.empty()
    ch_bam_stat_host = Channel.empty()
    ch_bam_stat_xeno = Channel.empty()
    ch_versions = Channel.empty()
    
    /*
    * Build bam/bai channels
    */
    bam
        .map{ it -> [ [ it[0], it[1] ], it[2] ]}
        .join (
            bai
                .map{ it -> [ [ it[0], it[1] ], it[2] ]}
        )
        .map { it -> [ it[0][0], it[0][1], it[1], it[2] ]}
        .set { ch_bam_bai }
                    
    
    if (params.workflow == 'pdx'){
        bam_host
            .map{ it -> [ [ it[0], it[1] ], it[2] ]}
            .join (
                bai_host
                    .map{ it -> [ [ it[0], it[1] ], it[2] ]}
            )
            .map { it -> [ it[0][0], it[0][1], it[1], it[2] ]}
            .set { ch_bam_bai_host }

        bam_xeno
            .map{ it -> [ [ it[0], it[1] ], it[2] ]}
            .join (
                bai_xeno
                    .map{ it -> [ [ it[0], it[1] ], it[2] ]}
            )
            .map { it -> [ it[0][0], it[0][1], it[1], it[2] ]}
            .set { ch_bam_bai_xeno }
        
    }

    /*
    * RNASeQC
    */
    if (params.run_rnaseqc){
        RNASEQC(
            params.workflow == 'pdx' ? ch_bam_bai_xeno : ch_bam_bai,
            collapsed_gtf,
            experiment
        )            
        ch_rnaseqc = RNASEQC.out.qc
        // [ [meta], path("*") ]
        ch_versions = ch_versions.mix(RNASEQC.out.versions.first())
    }

    /*
    * RSeQC
    */
    if (params.run_rseqc){
        RSEQC(
            params.workflow == 'pdx' ? ch_bam_bai_xeno : ch_bam_bai,
            tx_bed
        )
        ch_rseqc = RSEQC.out.qc
        // [ [meta], path("*") ]
        ch_versions = ch_versions.mix(RSEQC.out.versions.first())
    }

    /*
    * GATK collect_hs_metrics
    */
    if (params.run_hs_metrics){
            HS_METRICS(
                params.workflow == 'pdx' ? ch_bam_bai_xeno : ch_bam_bai,
                params.genome_fa,
                params.target_region
            )
            ch_hs_metrics = HS_METRICS.out.qc
                .map { it[1] }
            // [ [meta], path("*") ]
            ch_versions = ch_versions.mix(HS_METRICS.out.versions.first())

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
    bam_stat = ch_bam_stat
    bam_stat_host  = ch_bam_stat_host
    bam_stat_xeno = ch_bam_stat_xeno
    versions = ch_versions

}
