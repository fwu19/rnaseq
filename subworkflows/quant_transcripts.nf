/*
* Analyze transcripts
*/

include { GENERATE_TRANSCRIPT_COUNT_MATRIX } from '../modules/generate_transcript_count_matrix.nf'
include { DIFFERENTIAL_TRANSCRIPTS } from '../modules/differential_transcripts.nf'


workflow QUANT_TRANSCRIPTS {
    take:
    ch_salmon
    samplesheet
    gene_txt
    tx_txt

    main: 
    ch_tx_expr = Channel.empty()
    ch_dt = Channel.empty()
    ch_versions = Channel.empty()

    /* generate count matrix */
    if (params.run_tx_count){
        
        GENERATE_TRANSCRIPT_COUNT_MATRIX(
            samplesheet, 
            ch_salmon.map{it[2]}.collect(), 
            tx_txt,
            "EffectiveLength"
        )
        ch_tx_expr = GENERATE_TRANSCRIPT_COUNT_MATRIX.out.rds
        ch_versions = ch_versions.mix(GENERATE_TRANSCRIPT_COUNT_MATRIX.out.versions.first())

    }else{
        if (params.run_dt){
            ch_tx_expr = Channel.fromPath("${params.outdir}/expression_quantification/all_samples.transcript_raw_counts.txt")
        }
        
    }

    /* differential transcripts */
    if (params.run_dt){
        DIFFERENTIAL_TRANSCRIPTS(
            samplesheet, 
            file(params.comparison),
            ch_tx_expr, 
            "EffectiveLength",
            gene_txt
        )
        ch_dt = DIFFERENTIAL_TRANSCRIPTS.out.rds
        ch_versions = ch_versions.mix(DIFFERENTIAL_TRANSCRIPTS.out.versions)
    }


    emit:
    tx_expr = ch_tx_expr
    dt = ch_dt
    versions = ch_versions

}
