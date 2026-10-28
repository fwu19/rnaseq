/*
* Make references
*/

include { CAT as CAT_GENOME_FA } from '../modules/cat.nf'
include { CAT as CAT_GTF } from '../modules/cat.nf'
include { GTF2GENES  } from '../modules/gtf2genes.nf'
include { GTF2TRANSCRIPTS  } from '../modules/gtf2transcripts.nf'
include { COLLAPSE_GTF  } from '../modules/collapse_gtf.nf'
include { GTF2BED  } from '../modules/gtf2bed.nf'
include { GTF2FASTA  } from '../modules/gtf2fasta.nf'
include { MAKE_STAR } from '../modules/make_star.nf'
include { MAKE_BWA  } from '../modules/make_bwa.nf'
include { MAKE_XENGSORT } from '../modules/make_xengsort.nf'


workflow CHECK_REFERENCE {
    take:
    genome
    genome_fa
    gtf_path
    aligner
    star_path
    bwa_path
    rnaseqc_gtf_path
    rseqc_bed_path
    gene_txt_path
    tx_bed_path
    tx_fa_path
    tx_txt_path
    srcdir

    main: 
    index_dir = Channel.empty()
    gene_txt = Channel.empty()
    tx_fa = Channel.empty()
    tx_txt = Channel.empty()
    tx_bed = Channel.empty()
    collapsed_gtf = Channel.empty()
    xengsort_dir = Channel.empty()
    ch_versions = Channel.empty()

    /*
    *   Process GTF
    */
    /* generate genes.gtf */
    if (params.run_gtf2genes){
        File ref = new File("${srcdir}/references/${genome}/genes.txt")
        if (gene_txt_path){
            gene_txt = file(gene_txt_path, checkIfExists: true)
        } else if (ref.exists()){
            gene_txt = file(ref, checkIfExists: true)
        } else {
            GTF2GENES(file(gtf_path, checkIfExists: true))
            gene_txt = GTF2GENES.out.txt
            ch_versions = ch_versions.mix(GTF2GENES.out.versions)
        }
    }

    /* generate transcripts.fa */
    if (params.run_gtf2fasta){
        File ref = new File("${srcdir}/references/${genome}/transcripts.fa")
        if (tx_fa_path){
            tx_fa = file(tx_fa_path, checkIfExists:true)
        } else if (ref.exists()){
            tx_fa = file(ref, checkIfExists: true)
        }else{
            GTF2FASTA(
                file(genome_fa, checkIfExists: true),
                file(gtf_path, checkIfExists: true)
            )
            tx_fa = GTF2FASTA.out.fa
            ch_versions = ch_versions.mix(GTF2FASTA.out.versions)
        }

    }

    /* generate transcripts.txt */
    if (params.run_gtf2tx){
        File ref = new File("${srcdir}/references/${genome}/transcripts.txt")
        if (tx_txt_path){
            tx_txt = file(tx_txt_path, checkIfExists:true)
        } else if (ref.exists()){
            tx_txt = file(ref, checkIfExists: true)
        } else {
            GTF2TRANSCRIPTS(file(gtf_path, checkIfExists: true))
            tx_txt = GTF2TRANSCRIPTS.out.txt
            ch_versions = ch_versions.mix(GTF2TRANSCRIPTS.out.versions)
        }
    }

    /* collapse gtf for RNA-SeQC */
    if (params.run_collapse_gtf){
        File ref = new File("${srcdir}/references/${genome}/genes.collapsed.gtf")
        if (rnaseqc_gtf_path){
            collapsed_gtf = file(rnaseqc_gtf_path, checkIfExists:true)
        } else if (ref.exists()){
            collapsed_gtf = file(ref, checkIfExists: true)
        } else {
            COLLAPSE_GTF(file(gtf_path, checkIfExists: true))
            collapsed_gtf = COLLAPSE_GTF.out.gtf
            ch_versions = ch_versions.mix(COLLAPSE_GTF.out.versions)
        }
    }

    /* transcripts.bed for RSeQC */
    if (params.run_gtf2bed){
        File ref = new File("${srcdir}/references/${genome}/transcripts.bed")
        if (rseqc_bed_path){
            tx_bed = file(rseqc_bed_path, checkIfExists:true)
        } else if (ref.exists()){
            tx_bed = file(ref, checkIfExists: true)
        }else{
            GTF2BED(file(gtf_path, checkIfExists: true))
            tx_bed = GTF2BED.out.bed
            ch_versions = ch_versions.mix(GTF2BED.out.versions)
        }
    }

    /* 
    * Build index if needed
    */
    genome_fa_list = (genome_fa ?: '').split(',').collect{ it.trim() }.findAll{ it != '' }
    gtf_list = (gtf_path ?: '').split(',').collect{ it.trim() }.findAll{ it != '' }

    if (params.run_alignment && aligner == 'star'){
        if (star_path){
            index_dir = file(star_path, checkIfExists: true)
        }else{
            MAKE_STAR(
            genome_fa_list ? Channel.value(genome_fa_list) : Channel.empty(),
            file(gtf_list, checkIfExists: true)
            )
            index_dir = MAKE_STAR.out.dir
            ch_versions = ch_versions.mix(MAKE_STAR.out.versions)

        }

    } else if (params.run_alignment && aligner == 'bwa-mem'){
        if (bwa_path){
            index_dir = file(bwa_path, checkIfExists: true)
        }else{
            MAKE_BWA(
                genome_fa_list ? Channel.value(genome_fa_list) : Channel.empty()
            )
            index_dir = MAKE_BWA.out.dir
            ch_versions = ch_versions.mix(MAKE_BWA.out.versions)

        }
    }

    if (params.run_classify_fastq && params.run_xengsort){
        if (params.xengsort){
            xengsort_dir = file(params.xengsort, checkIfExists: true)
        }else{
            genome_fa_list = (params.genome_fa ?: '').split(',').collect{ it.trim() }.findAll{ it != '' }
            genome_fa_host_list = (params.genome_fa_host ?: '').split(',').collect{ it.trim() }.findAll{ it != '' }

            MAKE_XENGSORT(
                genome_fa_list ? Channel.value(genome_fa_list) : Channel.empty(),
                genome_fa_host_list ? Channel.value(genome_fa_host_list) : Channel.empty()
            )
            xengsort_dir = MAKE_XENGSORT.out.dir
            ch_versions = ch_versions.mix(MAKE_XENGSORT.out.versions)
        }
    }

    emit:
    index_dir
    gene_txt
    tx_fa
    tx_txt
    tx_bed
    collapsed_gtf
    xengsort_dir
    versions = ch_versions
}
