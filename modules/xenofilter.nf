
process XENOFILTER {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['fhR/4.1.2-foss-2021b', 'SAMtools/1.11-GCC-10.2.0']


    tag "Differential expression"

    publishDir "${params.outdir}/STAR/${genome}_filtered", mode: 'copy'
    
    input:
    tuple val(sample_id), path(graft_bam, stageAs: "graft/*"), path(host_bam, stageAs: "host/*")
    val (genome)
    val (mm_threshold)


    output:
    tuple val (sample_id), path ("*"), emit: bam
    
    script:
    """
    xenofilteR.r graft/*.bam host/*.bam $sample_id $mm_threshold ${task.cpus}
    rm -r graft/ host/
    mv Filtered_bams/XenofilteR.log ${sample_id}.XenofilteR.log
    mv Filtered_bams/* .
    rm -r Filtered_bams
    """
}
