process RNASEQC {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['RNA-SeQC/2.3.4-foss-2019b', 'SAMtools/1.11-GCC-10.2.0']


    tag "RNA-SeQC on ${sample_id}"

    publishDir "${params.outdir}/QC/rnaseqc/", pattern: "*.{tsv,gct}", mode: 'copy'

    input:
    tuple val(sample_id), path("*.{bam,bai}", stageAs: "input/*")

    output:
    tuple path("*"), emit: qc

    script:
    """
    rnaseqc.sh $sample_id ${params.rnaseqc_gtf} ${params.strand} ${params.read_type} input/*.bam 

    """
}


	