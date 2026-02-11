process HS_METRICS {
    label "process_medium"

    container = 'community.wave.seqera.io/library/bcftools_gatk4_samtools:460fbebd289c03b4'
    //module = ['GATK/4.4.0.0-GCCcore-12.2.0-Java-17', 'SAMtools/1.17-GCC-12.2.0']

    tag "collect hs_metrics on ${meta.id}"

    publishDir "${params.outdir}/QC/gatk/", pattern: "*.{tsv,gct}", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(bam), path(bai)
    path("genome.fa")
    path(target_region)

    output:
    path("${out_prefix}*.hs_metrics.txt"), emit: qc
    path ('versions.yml'), emit: versions

    script:
    """
    collect_hs_metrics.sh ${out_prefix} $bam $target_region genome.fa 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version 2>&1 | head -n 1 | sed -e "s/.* //g" )
        gatk: \$( gatk --version 2>&1 | sed -n 4,4p | sed -e "s/.* //g" )
    END_VERSIONS
    """
}


	