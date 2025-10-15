
process CAT_FASTQ {
    label "process_single"


    tag "CAT_FASTQ on $id"

    //publishDir "${params.outdir}/merged_fastq/", pattern: '*.fastq.gz', mode: 'copy'

    input:
    tuple val(id), path( "read1/*" ), path( "read2/*" )
    
    output:
    tuple val(id), path("${id}_merged_R1.fastq.gz"), path("${id}_merged_R2.fastq.gz"), emit: reads

    script:
    """
    cat_fastq.sh $id read1/ read2/

    """

}
