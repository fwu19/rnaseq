
process CAT_FASTQ {
    time = '1d'
    cpus = 1
    memory = '12G'


    tag "CAT_FASTQ on $sample_id"

    //publishDir "${params.outdir}/merged_fastq/", pattern: '*.fastq.gz', mode: 'copy'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)
    
    output:
    tuple val(sample_id), path("${sample_id}_merged_R1.fastq.gz"), path("${sample_id}_merged_R2.fastq.gz"), emit: reads

    script:
    """
    cat_fastq.sh $sample_id

    """

}
