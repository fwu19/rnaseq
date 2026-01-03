process SAMTOOLS_VIEW {
    label "process_medium"
    module = ['SAMtools/1.17-GCC-12.2.0']


    tag "samtools on ${out_prefix}"

    publishDir "${params.outdir}/BAM/", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(bam), path(bai)
    
    output:
    tuple val(meta), path("${out_prefix}.{stat,flagstat}"), emit: data
    path("${out_prefix}.{stat,flagstat}")
    path ('versions.yml'), emit: versions

    script:
    def suffix = task.ext.suffix ?: ""
    def args = task.ext.args ?: ""
    """
    samtools stats -@ ${task.cpus} $bam >${out_prefix}.stat
    samtools flagstat -@ ${task.cpus} -O tsv $bam >${out_prefix}.flagstat
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version | head -n 1 | sed -e "s/.* //g" )
    END_VERSIONS

    """
}
