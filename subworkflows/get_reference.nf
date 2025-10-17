/*
* Make references
*/

include { BUILD_INDEX } from './build_index.nf'

include { CAT as CAT_GENOME_FA } from '../modules/cat.nf'
include { CAT as CAT_GTF } from '../modules/cat.nf'
include { GTF2GENES  } from '../modules/gtf2genes.nf'
include { GTF2TRANSCRIPTS  } from '../modules/gtf2transcripts.nf'
include { COLLAPSE_GTF  } from '../modules/collapse_gtf.nf'
include { GTF2BED  } from '../modules/gtf2bed.nf'
include { GTF2FASTA  } from '../modules/gtf2fasta.nf'


workflow GET_REFERENCE {
    take:
    fa_str
    gtf_str
    aligner
    aligner_index // could be null

    main: 
    index_dir = Channel.empty()
    gene_txt = Channel.empty()
    tx_fa = Channel.empty()
    tx_txt = Channel.empty()
    tx_bed = Channel.empty()
    collapsed_gtf = Channel.empty()

    genome_fa = fa_str
    gtf = gtf_str
    
    /*
    * Concatenate multiple fasta and gtf files
    */
    /*
    def fa_list = fa_str.split(',')
    def fa_no = fa_list.size()
    def gtf_list = gtf_str.split(',')
    def gtf_no = gtf_list.size()

    genome_fa = "$projectDir/assets/dummy_file.csv"
    gtf = "$projectDir/assets/dummy_file.csv"
    
    if ( fa_no > 1){
        CAT_GENOME_FA(
            Channel.fromList(fa_list).collect(),
            "genome.fa",
            "fa"
        )
        genome_fa = CAT_GENOME_FA.out.file
    }else{
        genome_fa = file(fa_str, checkIfExists: true)
    }

    if ( gtf_no > 1 ){
        CAT_GTF(
            Channel.fromList(gtf_list).collect(),
            "genes.gtf",
            "gtf"
        )
        gtf = CAT_GTF.out.file
    }else{
        gtf = gtf_str
    }
    */

    /*
    *   Process GTF
    */
    /* generate genes.gtf */
    if (params.run_gene_count || params.run_dt){
        File ref = new File("${params.outdir}/references/${params.genome}/genes.txt")
        if (params.gene_txt){
            gene_txt = file(params.gene_txt, checkIfExists: true)
        } else if (ref.exists()){
            gene_txt = file(ref, checkIfExists: true)
        } else {
            GTF2GENES(gtf)
            gene_txt = GTF2GENES.out.txt
        }
    }

    /* generate transcripts.fa */
    if (params.run_salmon){
        File ref = new File("${params.outdir}/references/${params.genome}/transcripts.fa")
        if (params.tx_fa){
            tx_fa = file(params.tx_fa, checkIfExists:true)
        } else if (ref.exists()){
            tx_fa = file(ref, checkIfExists: true)
        }else{
            GTF2FASTA(
                genome_fa,
                gtf
            )
            tx_fa = GTF2FASTA.out.fa
        }

    }

    /* generate transcripts.txt */
    if (params.run_tx_count){
        File ref = new File("${params.outdir}/references/${params.genome}/transcripts.txt")
        if (params.tx_txt){
            tx_txt = file(params.tx_txt, checkIfExists:true)
        } else if (ref.exists()){
            tx_txt = file(ref, checkIfExists: true)
        } else {
            GTF2TRANSCRIPTS(gtf)
            tx_txt = GTF2TRANSCRIPTS.out.txt
        }
    }

    /* collapse gtf for RNA-SeQC */
    if (params.run_rnaseqc){
        if (params.rnaseqc_gtf){
                collapsed_gtf = file(params.rnaseqc_gtf, checkIfExists:true)
        }else{
                COLLAPSE_GTF(gtf)
                collapsed_gtf = COLLAPSE_GTF.out.gtf
        }
    }

    /* transcripts.bed for RSeQC */
    if (params.run_rseqc){
        if (params.rseqc_bed){
                tx_bed = file(params.tx_bed, checkIfExists:true)
        }else{
                GTF2BED(gtf)
                tx_bed = GTF2BED.out.bed
        }
    }

    /* 
    * Build index 
    */
    if (params.run_build_index){
        BUILD_INDEX(
            genome_fa,
            gtf == null ? "$projectDir/assets/dummy_file.csv" : gtf,
            aligner_index
        )

        if (aligner == 'star'){
            index_dir = BUILD_INDEX.out.star
        }else if (aligner == 'bwa'){
            index_dir = BUILD_INDEX.out.bwa
        }

    }else if (aligner == 'star'){
        if (params.star){
            index_dir = file(params.star, checkIfExists: true)
        }else{
            BUILD_INDEX(
                genome_fa,
                gtf == null ? "$projectDir/assets/dummy_file.csv" : gtf,
                "star"
            )
            index_dir = BUILD_INDEX.out.star
        }
    }else if (aligner == 'bwa-mem'){
        if (params.bwa){
            index_dir = file(params.bwa, checkIfExists: true)
        }else{
            BUILD_INDEX(
                genome_fa,
                gtf == null ? "$projectDir/assets/dummy_file.csv" : gtf,
                "bwa"
            )
            index_dir = BUILD_INDEX.out.bwa
        }
    }else{
        exit 1, "Use --aligner to provide an aligner. Possible values are star, bwa."
    }



    emit:
    index_dir = index_dir
    gene_txt = gene_txt
    tx_fa = tx_fa
    tx_txt = tx_txt
    tx_bed = tx_bed
    collapsed_gtf = collapsed_gtf

}
