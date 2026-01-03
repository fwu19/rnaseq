process RNASEQC {
    label "process_high"
    module = ['RNA-SeQC/2.4.2-foss-2021b', 'SAMtools/1.17-GCC-12.2.0']


    tag "RNA-SeQC on ${meta.id}"

    publishDir "${params.outdir}/QC/rnaseqc/", pattern: "*.{tsv,gct}", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(bam), path(bai)
    path(gtf)
    val(strand)
    val(read_type)

    output:
    tuple val(meta), path("${out_prefix}*.{tsv,gct}"), emit: qc
    path("${out_prefix}*.{tsv,gct}")
    path ('versions.yml'), emit: versions

    script:
    """
    rnaseqc.sh ${out_prefix} ${gtf} ${strand} ${read_type} ${bam} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version | head -n 1 | sed -e "s/.* //g" )
        rnaseqc: \$( rnaseqc -v 2>&1 | sed -n 3,3p | sed -e "s/.* //g" )
    END_VERSIONS

    """
}


	