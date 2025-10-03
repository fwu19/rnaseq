/*
* Analyze transcripts
*/

include { GTF2TRANSCRIPTS  } from '../modules/gtf2transcripts.nf'
include { GTF2GENES  } from '../modules/gtf2genes.nf'
include { GENERATE_TRANSCRIPT_COUNT_MATRIX } from '../modules/generate_transcript_count_matrix.nf'
include { DIFFERENTIAL_TRANSCRIPTS } from '../modules/differential_transcripts.nf'


workflow QUANT_TRANSCRIPTS {
    take:
    ch_salmon
    samplesheet
    gene_txt

    main: 
    ch_tx_rds = Channel.empty()
    ch_dt = Channel.empty()

    // generate transcripts.txt
    File ref = new File("${params.outdir}/references/${params.genome}/transcripts.txt")
    if (params.tx_txt){
        tx_txt = file(params.tx_txt, checkIfExists:true)
    } else if (ref.exists()){
        tx_txt = file(ref, checkIfExists: true)
    } else if (params.gtf){
        GTF2TRANSCRIPTS(params.gtf)
        tx_txt = GTF2TRANSCRIPTS.out.txt
    }
        
    GENERATE_TRANSCRIPT_COUNT_MATRIX(
        samplesheet, 
        ch_salmon.map{it[2]}.collect().ifEmpty([]), 
        tx_txt,
        "EffectiveLength"
    )
    ch_tx_rds = GENERATE_TRANSCRIPT_COUNT_MATRIX.out.rds


    /*
    * differential transcripts
    */
    if (params.run_dt && params.comparison){
        if (! gene_txt.name =~ 'dummy'){
            File refgene = new File("${params.outdir}/references/${params.genome}/genes.txt")
            if (params.gene_txt){
                gene_txt = file(params.gene_txt, checkIfExists: true)
            } else if (refgene.exists()){
                gene_txt = file(refgene, checkIfExists: true)
            } else if (params.gtf){
                GTF2GENES(params.gtf)
                gene_txt = GTF2GENES.out.txt
            }
        }

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
            params.de_gene_type?:'all'
        )
        ch_dt = DIFFERENTIAL_TRANSCRIPTS.out.rds
    }


    emit:
    tx_rds = ch_tx_rds
    dt = ch_dt

}
