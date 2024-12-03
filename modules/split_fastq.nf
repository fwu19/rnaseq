/*
* Split fastq files into smaller chunks
*/


process SPLIT_FASTQ {
    time = '1d'
    cpus = 6
    memory = '36G'

    tag "Split fastq on ${meta.id}"

    publishDir "${params.outdir}/fastq/", mode: 'copy'

    input:
    tuple val(meta), path(fq1), path(fq2)

    output:
    tuple val(meta), path( "${meta.id}.csv" ), emit: csv

    script:
    def args = task.ext.args ?: ""
    """
    split_fastq.sh ${meta.id} $fq1 $fq2 "$args" 
    """
    

}
