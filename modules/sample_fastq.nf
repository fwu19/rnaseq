/*
* Sample reads from a fastq file
*/


process SAMPLE_FASTQ {
    label "process_single"

    tag "Sample reads from ${out_prefix}"

    input:
    tuple val(meta), val(out_prefix), path(fq1), path(fq2)
    val(n_reads)

    output:
    tuple val(meta), val(out_prefix), path("${out_prefix}_R1.head_${n_reads}.fastq.gz"), path("${out_prefix}_R2.head_${n_reads}.fastq.gz"), emit: fq

    script:
    def n_lines = n_reads * 4
    """
    zcat ${fq1} | head -n ${n_lines} | gzip -c > ${out_prefix}_R1.head_${n_reads}.fastq.gz
    zcat ${fq2} | head -n ${n_lines} | gzip -c > ${out_prefix}_R2.head_${n_reads}.fastq.gz
    """
}

