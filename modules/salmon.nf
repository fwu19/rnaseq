
process SALMON {
    label "process_high"

    container = 'community.wave.seqera.io/library/bamtools_salmon_samtools:22d508928f8d86c7'
    //module = ['Salmon/1.10.1-GCC-12.2.0', 'SAMtools/1.17-GCC-12.2.0']

    tag "Salmon on ${out_prefix}"

    publishDir "${params.outdir}/Salmon/${params.genome}/", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(bam)
    path (tx_fa)
    
    output:
    tuple val(meta), val(out_prefix), path( "${out_prefix}/" ), emit: sf
    path( "${out_prefix}/" )
    path ('versions.yml'), emit: versions

    script:

    """
    salmon.sh ${task.cpus} ${bam} ${tx_fa} ${out_prefix}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version 2>&1 | head -n 1 | sed -e "s/.* //g" )
        salmon: \$( salmon --version 2>&1 | head -n 1 | sed -e "s/.* //g" )
    END_VERSIONS

    """
}
