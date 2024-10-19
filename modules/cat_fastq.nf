
process CAT_FASTQ {
    time = '1d'
    cpus = 1
    memory = '12G'


    tag "CAT_FASTQ on $sample_id"

    //publishDir "${params.outdir}/merged_fastq/", pattern: '*.fastq.gz', mode: 'copy'

    input:
    tuple val(sample_id), path("*.fastq.gz", stageAs: "input/*")
    
    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: reads

    script:
    """
    cat_fastq.sh $sample_id input/
    """
}
