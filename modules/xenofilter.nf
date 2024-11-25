
process XENOFILTER {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['fhR/4.1.2-foss-2021b', 'SAMtools/1.11-GCC-10.2.0']


    tag "Filter bam on ${meta.id}"

    publishDir "${params.outdir}/STAR/${genome}_filtered", mode: 'copy'
    
    input:
    tuple val(meta), path(graft_bam, stageAs: "graft/*"), path(host_bam, stageAs: "host/*")
    val (genome)
    val (mm_threshold)


    output:
    tuple val(meta), path ("*.{bam,bai,log}"), emit: bam
    path("*.{bam,bai,log}")
    
    script:
    """
    xenofilteR.r graft/*.bam host/*.bam ${meta.id} $mm_threshold ${task.cpus}
    mv Filtered_bams/XenofilteR.log ${meta.id}.XenofilteR.log
    mv Filtered_bams/* .
    
    """
}
