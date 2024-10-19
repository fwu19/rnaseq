process RSEQC {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['R/4.1.2-foss-2020b', 'SAMtools/1.11-GCC-10.2.0']


    tag "RSeQC on ${sample_id}"

    publishDir "${params.outdir}/RSeQC/", pattern: "${sample_id}.*.{txt,pdf,r,bed,log,xls,xlsx}", mode: 'copy'

    input:
    tuple val(sample_id), path("*.{bam,bai}", stageAs: "input/*")

    output:
    tuple path("${sample_id}.*"), emit: qc

    script:
    """
    rseqc.sh $sample_id ${params.rseqcBed} ${params.txBed} ${params.geneBed} input/*.bam 

    """
}

