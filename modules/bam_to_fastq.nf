
process BAM_TO_FASTQ {
    time = '1d'
    cpus = 6
    memory = '24G'
    //module = ['BEDTools/2.30.0-GCC-12.2.0', 'SAMtools/1.17-GCC-12.2.0']
    container "quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:96081f93c4cf2811a53d5112d138ad07ebe1b815-0" // bedtools 2.27.1 samtools 1.9

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
