
process GENERATE_GENE_COUNT_MATRIX {
    label "process_single"
    container "docker://fwu19/r-libs:4.1.2" 


    tag "generate count matrix"

    publishDir "saved_data/", mode: 'copy'
    
    input:
    path(samplesheet)
    //tuple val(meta) val(id) path("counts/*")
    path("counts/*")
    path(gene_txt)
    val(length_col)
    val(strand)
    val(workflow)

    output:
    path ("gene.y0.rds"), emit: rds
    path ("*")
    
    script:
    """
    generate_gene_count_matrix.r input=$samplesheet count.dir=counts gene.txt=$gene_txt length.col=${length_col} strand=${strand} workflow=${workflow}
    mv y0.rds gene.y0.rds
    
    """
}
