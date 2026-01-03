
process XENOFILTER {
    label "process_high"
    container "docker://fwu19/r-libs:4.1.2" 


    tag "Filter bam on ${out_prefix}"

    publishDir "${params.outdir}/STAR/${genome}_filtered", mode: 'copy'
    
    input:
    tuple val(meta), val(out_prefix), path(graft_bam, stageAs: "graft/*"), path(graft_bai, stageAs: "graft/*"), path(host_bam, stageAs: "host/*"), path(host_bai, stageAs: "host/*")
    val (genome)
    val (mm_threshold)


    output:
    tuple val(meta), val(out_prefix), path ("*.bam"), emit: bam
    tuple val(meta), val(out_prefix), path ("*.bai"), emit: bai
    path("*.{bam,bai,log}")
    path ('versions.yml'), emit: versions
    
    script:
    """
    xenofilteR.r graft/*.bam host/*.bam ${out_prefix} $mm_threshold ${task.cpus}
    mv Filtered_bams/XenofilteR.log ${out_prefix}.XenofilteR.log
    mv Filtered_bams/*.bam* .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS

    
    """
}
