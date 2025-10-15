
process BAM_TO_FASTQ {
    label "process_medium"
    module = ['BEDTools/2.30.0-GCC-12.2.0', 'SAMtools/1.17-GCC-12.2.0']


    tag "convert bam to fastq on ${out_prefix}"

    publishDir "${params.outdir}/fastq/${out_prefix}", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(bam)

    output:
    tuple val(meta), val(out_prefix), path( "${out_prefix}_R1.fastq.gz" ), path( "${out_prefix}_R2.fastq.gz" ), emit: fq
    path("*.fastq.gz")

    script:
    def args = task.ext.args ?: ""
    """
    samtools sort -n $bam | bedtools bamtofastq -i - -fq ${out_prefix}_R1.fastq -fq2 ${out_prefix}_R2.fastq
    gzip ${out_prefix}_R1.fastq
    gzip ${out_prefix}_R2.fastq
    """
}
