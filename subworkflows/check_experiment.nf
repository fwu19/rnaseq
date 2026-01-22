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

    main: 
    ch_reads= Channel.empty()
    ch_bam = Channel.empty()
    ch_bai = Channel.empty()
    ch_versions = Channel.empty()
    csv = Channel.empty()

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
        ch_bam = STAR.out.bam
        // [ [meta], val(out_prefix), path(bam) ]
        ch_bai = STAR.out.bai
        // [ [meta], val(out_prefix), path(bai) ]
        ch_versions = ch_versions.mix(STAR.out.versions.first())

    }
    
    if (params.aligner == 'bwa-mem'){

        BWA_MEM(
            sub_reads, 
            params.genome, 
            aligner_index
        )
        ch_bam = BWA_MEM.out.bam
        // [ [meta], val(out_prefix), path(bam) ]
        ch_bai = BWA_MEM.out.bai
        // [ [meta], val(out_prefix), path(bai) ]
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    }

    INFER_EXPERIMENT(
        ch_bam
            .map{ it -> [ [ it[0], it[1] ], it[2] ]}
            .join (
                ch_bai
                    .map{ it -> [ [ it[0], it[1] ], it[2] ]}
            ) 
            .map { it -> [ it[0][0], it[0][1], it[1], it[2] ]},
        tx_bed
    )
    csv = INFER_EXPERIMENT.out.csv
    ch_versions = ch_versions.mix(INFER_EXPERIMENT.out.versions.first())

    emit:
    csv
    versions = ch_versions

}
