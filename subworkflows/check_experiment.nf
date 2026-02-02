/*
* Make references
*/

include { SAMPLE_FASTQ } from '../modules/sample_fastq.nf'
include { STAR } from '../modules/star.nf'
include { BWA_MEM } from '../modules/bwa_mem.nf'
include { INFER_EXPERIMENT } from '../modules/infer_experiment.nf'


workflow CHECK_EXPERIMENT {
    take:
    fq
    aligner_index
    tx_bed
    srcdir

    main: 
    ch_reads= Channel.empty()
    ch_bam = Channel.empty()
    ch_bai = Channel.empty()
    ch_versions = Channel.empty()
    csv = Channel.empty()

    File saved = new File("${srcdir}/csv/infer_experiment.csv")
    if (saved.exists()){
        csv = file(saved)
    }else{

        ch_reads = fq
        .splitCsv(header: true)
        .map {
            row -> [ row, row.id, row.fastq_1, row.fastq_2 ]
        }

        SAMPLE_FASTQ(
        ch_reads.first(),
        1000000
        )
        sub_reads = SAMPLE_FASTQ.out.fq

        if (params.aligner == 'star'){            
        STAR(
            sub_reads,
            params.genome, 
            aligner_index, 
            params.gtf
        )
        ch_bam_bai = STAR.out.bam_bai
        // [ [meta], val(out_prefix), path(bam), path(bai) ]
        ch_versions = ch_versions.mix(STAR.out.versions.first())

        }
    
        if (params.aligner == 'bwa-mem'){

        BWA_MEM(
            sub_reads, 
            params.genome, 
            aligner_index
        )
        ch_bam = BWA_MEM.out.bam_bai
        // [ [meta], val(out_prefix), path(bam), path(bai) ]
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

        }

        INFER_EXPERIMENT(
        ch_bam_bai,
        tx_bed
        )
        csv = INFER_EXPERIMENT.out.csv
        ch_versions = ch_versions.mix(INFER_EXPERIMENT.out.versions)

    }

    emit:
    csv
    versions = ch_versions

}
