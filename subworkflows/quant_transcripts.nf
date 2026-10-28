/*
* Analyze transcripts
*/

include { SALMON  } from '../modules/salmon.nf'
include { GENERATE_TRANSCRIPT_COUNT_MATRIX } from '../modules/generate_transcript_count_matrix.nf'
include { DIFFERENTIAL_TRANSCRIPTS } from '../modules/differential_transcripts.nf'
include { SALMON as SALMON_HOST  } from '../modules/salmon.nf'
include { GENERATE_TRANSCRIPT_COUNT_MATRIX as GENERATE_TRANSCRIPT_COUNT_MATRIX_HOST } from '../modules/generate_transcript_count_matrix.nf'
include { DIFFERENTIAL_TRANSCRIPTS as DIFFERENTIAL_TRANSCRIPTS_HOST } from '../modules/differential_transcripts.nf'


workflow QUANT_TRANSCRIPTS {
    take:
    samplesheet
    ch_tx_bam
    tx_fa
    tx_txt
    gene_txt
    ch_tx_bam_host
    tx_fa_host
    tx_txt_host
    gene_txt_host
    srcdir

    main: 
    ch_salmon = Channel.empty()
    ch_tx_expr = Channel.empty()
    ch_dt = Channel.empty()
    dt_csv = Channel.empty()
    ch_salmon_host = Channel.empty()
    ch_tx_expr_host = Channel.empty()
    ch_dt_host = Channel.empty()
    dt_csv_host = Channel.empty()
    ch_versions = Channel.empty()

    /* generate count matrix */
    if (params.run_tx_count){
        if (!params.run_alignment ){
            ch_tx_bam = Channel.fromPath("${srcdir}/csv/map2transcriptome.${params.aligner}.${params.genome}.csv")
                .splitCsv(header: true)
                .map { it -> [ [ it ], it.id, "${srcdir}/${it.tx_bam}" ] }
        }

        if (params.run_salmon){        
            SALMON(
                    ch_tx_bam,
                    tx_fa
            )
            ch_salmon = SALMON.out.sf
            ch_versions = ch_versions.mix(SALMON.out.versions.first())

        }

        GENERATE_TRANSCRIPT_COUNT_MATRIX(
            samplesheet, 
            ch_salmon.map{it[2]}.collect(), 
            tx_txt,
            "EffectiveLength"
        )
        ch_tx_expr = GENERATE_TRANSCRIPT_COUNT_MATRIX.out.rds
        ch_versions = ch_versions.mix(GENERATE_TRANSCRIPT_COUNT_MATRIX.out.versions.first())

    }

    /* differential transcripts */
    if (params.run_dt){
        DIFFERENTIAL_TRANSCRIPTS(
            samplesheet, 
            params.comparison_transcripts ? file(params.comparison_transcripts) : (params.comparison ? file(params.comparison) : file("$projectDir/assets/dummy_file.csv")),
            params.run_tx_count ? ch_tx_expr : Channel.fromPath("${srcdir}/${params.tx_expr_dir}/all_samples.transcript_raw_counts.txt", checkIfExists: true), 
            "EffectiveLength",
            gene_txt,
            params.de_method
        )
        dt_csv = DIFFERENTIAL_TRANSCRIPTS.out.csv
        ch_dt = DIFFERENTIAL_TRANSCRIPTS.out.rds
        ch_versions = ch_versions.mix(DIFFERENTIAL_TRANSCRIPTS.out.versions)
    }

    /* generate count matrix for host genome */
    if (params.workflow == 'pdx' && params.run_tx_count_host){
        if (!params.run_alignment ){
            ch_tx_bam_host = Channel.fromPath("${srcdir}/csv/map2transcriptome.${params.aligner}.${params.genome_host}.csv")
                .splitCsv(header: true)
                .map { it -> [ [ it ], it.id, "${srcdir}/${it.tx_bam}" ] }
        }

        if (params.run_salmon_host){        
            SALMON_HOST(
                ch_tx_bam_host,
                tx_fa_host
            )
            ch_salmon_host = SALMON_HOST.out.sf
            ch_versions = ch_versions.mix(SALMON_HOST.out.versions.first())

        }

        GENERATE_TRANSCRIPT_COUNT_MATRIX_HOST(
            samplesheet, 
            ch_salmon_host.map{it[2]}.collect(), 
            tx_txt_host,
            "EffectiveLength"
        )
        ch_tx_expr_host = GENERATE_TRANSCRIPT_COUNT_MATRIX_HOST.out.rds
        ch_versions = ch_versions.mix(GENERATE_TRANSCRIPT_COUNT_MATRIX_HOST.out.versions.first())

    }

    /* differential transcripts */
    if (params.workflow == 'pdx' && params.run_dt_host){
        DIFFERENTIAL_TRANSCRIPTS_HOST(
            samplesheet, 
            params.comparison_transcripts ? file(params.comparison_transcripts) : (params.comparison ? file(params.comparison) : file("$projectDir/assets/dummy_file.csv")),
            params.run_tx_count_host ? ch_tx_expr_host : Channel.fromPath("${srcdir}/${params.tx_expr_host_dir}/all_samples.transcript_raw_counts.txt", checkIfExists: true), 
            "EffectiveLength",
            gene_txt_host,
            params.de_method
        )
        dt_csv_host = DIFFERENTIAL_TRANSCRIPTS_HOST.out.csv
        ch_dt_host = DIFFERENTIAL_TRANSCRIPTS_HOST.out.rds
        ch_versions = ch_versions.mix(DIFFERENTIAL_TRANSCRIPTS_HOST.out.versions)
    }

    emit:
    salmon = ch_salmon
    tx_expr = ch_tx_expr
    dt = ch_dt
    dt_csv
    salmon_host = ch_salmon_host
    tx_expr_host = ch_tx_expr_host
    dt_host = ch_dt_host
    dt_csv_host
    versions = ch_versions

}
