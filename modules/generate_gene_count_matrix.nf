
process GENERATE_GENE_COUNT_MATRIX {
    time = '1d'
    cpus = 1
    memory = '24G'
    module = ['fhR/4.1.2-foss-2021b']


    tag "generate count matrix"

    publishDir "saved_data/", mode: 'copy'
    
    input:
    path(samplesheet)
    path("*", stageAs: "counts/*")
    path(gene_txt)
    val(length_col)
    val(strand)
    val(workflow)

    output:
    path ("y0.rds"), emit: rds, optional: true
    path ("*"), optional: true
    
    script:
    """
    generate_gene_count_matrix.r input=$samplesheet count.dir=counts gene.txt=$gene_txt length.col=${length_col} strand=${strand} workflow=${workflow}
    rm -r counts/
    ln -s y0.rds gene.y0.rds
    ln -s pca.rds gene.pca.rds
    """
}
