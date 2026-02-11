/*
* Check fastq files and cat fastq if needed
*/

include { ARRIBA } from '../modules/arriba.nf'

workflow DETECT_FUSIONS {
    take:
    ch_reads
    genome
    index_dir
    gtf
    genome_fa
    blacklist
    known_fusions
    protein_domains

    main: 
    ch_arriba = Channel.empty()
    ch_versions = Channel.empty()

    if (params.run_arriba){
        ARRIBA(
            ch_reads, 
            genome, 
            index_dir, 
            gtf,
            genome_fa,
            blacklist,
            known_fusions,
            protein_domains
        )
        ch_arriba = ARRIBA.out.tsv
        ch_versions = ARRIBA.out.versions
    }

    emit:
    versions = ch_versions
    arriba = ch_arriba
}
