
process GENERATE_TRANSCRIPT_COUNT_MATRIX {
    time = '1d'
    cpus = 1
    memory = '24G'
    container "docker://fwu19/r-libs:4.1.2" 

    tag "Differential transcripts"

    publishDir "${params.outdir}/transcript_expression/", mode: 'copy'
    
    input:
    path(samplesheet)
    path(comparison)
    path("*", stageAs: "counts/*")
    path(gene_txt)
    val(length_col)
    
    output:
    path ("transcript.y0.rds"), emit: rds
    path ("*")
    
    script:
    """
    generate_transcript_count_matrix.r input=$samplesheet comparison=$comparison count.dir=counts gene.txt=$gene_txt length.col=${length_col}
    mv y0.rds transcript.y0.rds
    """
}
