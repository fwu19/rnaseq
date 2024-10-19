process FASTQC {
    time = '1d'
    cpus = 1
    memory = '12G'
    module = 'FastQC/0.11.9-Java-11'


    tag "FASTQC on $sample_id"

    publishDir "${params.outdir}/FastQC/", pattern: '*.{html,zip}', mode: 'copy'

    input:
    tuple val(sample_id), path("*.fastq.gz", stageAs: "input/*")
    

    output:
    tuple path("*.{html,zip}"), emit: qc

    script:
    """
    fastqc.sh ${sample_id} input/
    """
}
