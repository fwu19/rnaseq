/*
* Sample reads from a fastq file
*/


process SAMPLE_FASTQ {
    label "process_low"

    tag "Sample reads from ${out_prefix}"

    input:
    tuple val(meta), val(out_prefix), path(fq1), path(fq2)
    val(n_reads)

    output:
    tuple val(meta), val(out_prefix), path("${out_prefix}_R1.head_${n_reads}.fastq.gz"), path("${out_prefix}_R2.head_${n_reads}.fastq.gz"), emit: fq

    script:
    def n_lines = n_reads * 4
    """
    sample_fastq.sh $fq1 $fq2 $out_prefix $n_reads

    """
}

