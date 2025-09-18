process FASTQC {
    time = '1d'
    cpus = 1
    memory = '12G'
    module = 'FastQC/0.12.1-Java-11'


    tag "FASTQC on ${meta.id}"

    publishDir "${params.outdir}/QC/fastqc/", pattern: '*.{html,zip}', mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path("${meta.id}_R1.fastq.gz"), path("${meta.id}_R2.fastq.gz")
    

    output:
    tuple val(meta), path("*.{html,zip}"), emit: qc

    script:
    """
    fastqc -o ./ --casava ${meta.id}_R1.fastq.gz ${meta.id}_R2.fastq.gz
    """
}
