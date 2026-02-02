/*
* Analyze transcripts
*/

include { SALMON  } from '../modules/salmon.nf'
include { GENERATE_TRANSCRIPT_COUNT_MATRIX } from '../modules/generate_transcript_count_matrix.nf'
include { DIFFERENTIAL_TRANSCRIPTS } from '../modules/differential_transcripts.nf'


workflow QUANT_TRANSCRIPTS {
    take:
    samplesheet
    ch_tx_bam
    tx_fa
    tx_txt
    gene_txt
    srcdir

    main: 
    ch_salmon = Channel.empty()
    ch_tx_expr = Channel.empty()
    ch_dt = Channel.empty()
    ch_versions = Channel.empty()

    /* generate count matrix */
    if (params.run_tx_count){
        if (!params.run_alignment ){
            def my_dir = new File("${srcdir}")
            def srcdir = my_dir.absolutePath
            ch_tx_bam = Channel.fromPath("${srcdir}/csv/mapped.${params.aligner}.csv")
                .splitCsv(header: true)
                .map { it -> [ [ it ], it.id, "${srcdir}/${it.tx_bam}" ] }
        }

        if (params.run_salmon){        
            if (params.workflow == 'pdx'){
            //should filter Aligned.toTranscriptome.out.bam
            }else{
                SALMON(
                    ch_tx_bam,
                    tx_fa
                )
            }
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
            file(params.comparison),
            params.run_tx_count ? ch_tx_expr : Channel.fromPath("${srcdir}/expression_quantification/all_samples.transcript_raw_counts.txt", checkIfExists: true), 
            "EffectiveLength",
            gene_txt
        )
        ch_dt = DIFFERENTIAL_TRANSCRIPTS.out.rds
        ch_versions = ch_versions.mix(DIFFERENTIAL_TRANSCRIPTS.out.versions)
    }


    emit:
    salmon = ch_salmon
    tx_expr = ch_tx_expr
    dt = ch_dt
    versions = ch_versions

}
