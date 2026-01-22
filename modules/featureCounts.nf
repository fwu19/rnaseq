process FEATURECOUNTS {
    label "process_high"

    container = 'quay.io/biocontainers/subread:2.0.0--hed695b0_0'
    //module = ['Subread/2.0.0-GCC-8.3.0']

    tag "featureCounts on ${meta.id}"

    publishDir "${params.outdir}/featureCounts/", pattern: "${meta.id}*.{txt,summary}", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path( bam )
    path (gtf)
    val (experiment)
    
    output:
    tuple val(meta), val(out_prefix), path("${out_prefix}*.txt"), emit: counts
    path ('versions.yml'), emit: versions
    path("${out_prefix}*.{txt,summary}")

    script:
    """
    featureCounts.sh ${out_prefix} ${gtf} ${bam} ${experiment.read_type} ${experiment.strand} ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        featureCounts: \$( featureCounts -v 2>&1 | sed -n 2,2p | sed -e "s/.* //g" )
    END_VERSIONS

    """
}