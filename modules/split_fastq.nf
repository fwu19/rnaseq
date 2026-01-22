/*
* Split fastq files into smaller chunks
*/


process SPLIT_FASTQ {
    label "process_medium"

    tag "Split fastq on ${out_prefix}"

    publishDir "${params.outdir}/fastq/", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(fq1), path(fq2)
    val(size)

    output:
    tuple val(meta), val(out_prefix), path( "${out_prefix}.csv" ), emit: csv
    path ('versions.yml'), emit: versions

    script:
    def args = task.ext.args ?: ""
    """
    split_fastq.sh ${out_prefix} $fq1 $fq2 $size "$args" 

        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit -h 2>&1 | sed -n 3,3p | sed -e "s/.* //g" )
    END_VERSIONS

    """
    

}
