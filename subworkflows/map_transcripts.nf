/*
* Map to transcripts
*/

include { SALMON  } from '../modules/salmon.nf'
include { GTF2FASTA  } from '../modules/gtf2fasta.nf'


workflow MAP_TRANSCRIPTS {
    take:
    ch_tx_bam

    main: 
    ch_salmon = Channel.empty()

    if (params.run_salmon){
        // generate tx.fa
        File ref = new File("${params.outdir}/references/${params.genome}/transcripts.fa")
        if (params.tx_fa){
            tx_fa = file(params.tx_fa, checkIfExists:true)
        } else if (ref.exists()){
            tx_fa = file(ref, checkIfExists: true)
        }else{
            GTF2FASTA(
                params.genome_fa,
                params.gtf
            )
            tx_fa = GTF2FASTA.out.fa
        }
        
        if (params.workflow == 'pdx'){
            //should filter Aligned.toTranscriptome.out.bam
        }else{
            SALMON(
                ch_tx_bam,
                tx_fa
            )
        }
        ch_salmon = SALMON.out.sf
    }

    emit:
    salmon = ch_salmon

}
