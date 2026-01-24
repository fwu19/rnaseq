
process BAM_TO_FASTQ {
    label "process_medium"

    container = 'quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:ac2171da3c40903e835de564d81f41fb9eb475d2-2' 
//'community.wave.seqera.io/library/bamtools_salmon_samtools:22d508928f8d86c7'
    //module = ['BEDTools/2.30.0-GCC-12.2.0', 'SAMtools/1.17-GCC-12.2.0']

    tag "convert bam to fastq on ${out_prefix}"

    publishDir "${params.outdir}/fastq/${out_prefix}", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(bam)

    output:
    tuple val(meta), val(out_prefix), path( "${out_prefix}_R1.fastq.gz" ), path( "${out_prefix}_R2.fastq.gz" ), emit: fq
    path("*.fastq.gz")
    path ('versions.yml'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    """
    samtools sort -n $bam | bedtools bamtofastq -i - -fq ${out_prefix}_R1.fastq -fq2 ${out_prefix}_R2.fastq
    gzip ${out_prefix}_R1.fastq
    gzip ${out_prefix}_R2.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$( bedtools --version | head -n 1 | sed -e "s/.* //g" )
        samtools: \$( samtools --version | head -n 1 | sed -e "s/.* //g" )
    END_VERSIONS
    """
}
