
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
    tuple val(meta), val(out_prefix), path ("*.bam"), path ("*.bai"), emit: bam_bai
    path("*.{bam,bai,log}")
    path ('versions.yml'), emit: versions
    
    script:
    """
    xenofilteR.r graft_dir="graft" host_dir="host" out_prefix=${out_prefix} mm_threshold=$mm_threshold nworkers=${task.cpus}
    mv Filtered_bams/XenofilteR.log ${out_prefix}.XenofilteR.log
    mv Filtered_bams/*.bam* .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1 | head -n 1 | sed -e "s/R version //g; s/ .*//g" )
    END_VERSIONS

    
    """
}
