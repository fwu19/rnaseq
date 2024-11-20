
process CAT_FASTQ {
    time = '1d'
    cpus = 1
    memory = '12G'


    tag "CAT_FASTQ on $id"

    //publishDir "${params.outdir}/merged_fastq/", pattern: '*.fastq.gz', mode: 'copy'

    input:
    tuple val(id), path( "*" ), path( "*" )
    
    output:
    tuple val(id), path("${id}_merged_R1.fastq.gz"), path("${id}_merged_R2.fastq.gz"), emit: reads

    script:
    """
    cat_fastq.sh $id

    """

}
