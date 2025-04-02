
process XENOFILTER {
    time = '1d'
    cpus = 8
    memory = '48G'
    module = ['fhR/4.1.2-foss-2021b', 'SAMtools/1.17-GCC-12.2.0']


    tag "Filter bam on ${out_prefix}"

    publishDir "${params.outdir}/STAR/${genome}_filtered", mode: 'copy'
    
    input:
    tuple val(meta), val(out_prefix), path(graft_bam, stageAs: "graft/*"), path(host_bam, stageAs: "host/*")
    val (genome)
    val (mm_threshold)


    output:
    tuple val(meta), val(out_prefix), path ("*.bam"), emit: bam
    tuple val(meta), val(out_prefix), path ("*.bai"), emit: bai
    path("*.{bam,bai,log}")
    
    script:
    """
    xenofilteR.r graft/*.bam host/*.bam ${out_prefix} $mm_threshold ${task.cpus}
    mv Filtered_bams/XenofilteR.log ${out_prefix}.XenofilteR.log
    mv Filtered_bams/*.bam* .
    
    """
}
