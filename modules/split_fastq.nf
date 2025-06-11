/*
* Split fastq files into smaller chunks
*/


process SPLIT_FASTQ {
    time = '1d'
    cpus = 6
    memory = '36G'

    tag "Split fastq on ${out_prefix}"

    publishDir "${params.outdir}/fastq/", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(fq1), path(fq2)
    val(size)

    output:
    tuple val(meta), val(out_prefix), path( "${out_prefix}.csv" ), emit: csv

    script:
    def args = task.ext.args ?: ""
    """
    split_fastq.sh ${out_prefix} $fq1 $fq2 $size "$args" 
    """
    

}
