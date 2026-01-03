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
    ch_tx_rds = Channel.empty()
    ch_dt = Channel.empty()
    ch_versions = Channel.empty()

    /* generate count matrix */
    if (params.run_tx_count){
        
        GENERATE_TRANSCRIPT_COUNT_MATRIX(
            samplesheet, 
            ch_salmon.map{it[2]}.collect().ifEmpty([]), 
            tx_txt,
            "EffectiveLength"
        )
        ch_tx_rds = GENERATE_TRANSCRIPT_COUNT_MATRIX.out.rds
        ch_versions = ch_versions.mix(GENERATE_TRANSCRIPT_COUNT_MATRIX.out.versions)

    }

    /* differential transcripts */
    if (params.run_dt && params.comparison){

        DIFFERENTIAL_TRANSCRIPTS(
            samplesheet, 
            Channel.fromPath(params.comparison, checkIfExists: true), 
            ch_tx_rds, 
            "EffectiveLength",
            params.fdr,
            params.fc,
            params.fdr2,
            params.fc2,
            gene_txt,
            params.de_gene_type
        )
        ch_dt = DIFFERENTIAL_TRANSCRIPTS.out.rds
        ch_versions = ch_versions.mix(DIFFERENTIAL_TRANSCRIPTS.out.versions)
    }


    emit:
    tx_rds = ch_tx_rds
    dt = ch_dt
    versions = ch_versions

}
