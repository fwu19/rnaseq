process GENERATE_GENE_COUNT_MATRIX {
    label "process_single"
    container "docker://fwu19/r-libs:4.1.2" 


    tag "generate count matrix"

    publishDir "saved_data/", mode: 'copy'
    
    input:
    path(samplesheet)
    path("counts/*")
    path(gene_txt)
    val(length_col)
    path(experiment)

    output:
    path ("gene.y0.rds"), emit: rds
    path ('versions.yml'), emit: versions
    path ("*")
    
    script:
    """
    generate_gene_count_matrix.r input=$samplesheet count.dir=counts gene.txt=$gene_txt length.col=$length_col experiment=$experiment workflow=${params.workflow}
    mv y0.rds gene.y0.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1 | head -n 1 | sed -e "s/R version //g; s/ .*//g" )
    END_VERSIONS
    
    """
}
