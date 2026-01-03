/*
* Map to transcripts
*/

include { SALMON  } from '../modules/salmon.nf'


workflow MAP_TRANSCRIPTS {
    take:
    ch_tx_bam
    tx_fa

    main: 
    ch_salmon = Channel.empty()
    ch_versions = Channel.empty()

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
        ch_versions = ch_versions.mix(SALMON.out.versions)
    }

    emit:
    salmon = ch_salmon
    versions = ch_versions

}
