/*
* Map to transcripts
*/

include { SALMON  } from '../modules/salmon.nf'
include { WRITE_CSV as WRITE_CSV_SALMON} from '../modules/write_csv.nf'


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
        ch_versions = ch_versions.mix(SALMON.out.versions.first())

        WRITE_CSV_SALMON(
            ch_salmon
                .map { 
                    it -> it[0] + [salmon: "Salmon/${params.genome}/${it[0].id}/"]
                }
                .collect(),
            "salmon.csv"        
        )
    }

    emit:
    salmon = ch_salmon
    versions = ch_versions

}
