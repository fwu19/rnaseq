
process GENERATE_TRANSCRIPT_COUNT_MATRIX {
    time = '1d'
    cpus = 1
    memory = '24G'
    module = ['fhR/4.1.2-foss-2021b']


    tag "Differential transcripts"

    publishDir "${params.outdir}/transcript_expression/", mode: 'copy'
    
    input:
    path(samplesheet)
    path(comparison)
    path("*", stageAs: "counts/*")
    path(gene_txt)
    val(length_col)
    

    output:
    path ("y0.rds"), emit: rds, optional: true
    path ("*"), optional: true
    
    script:
    """
    generate_transcript_count_matrix.r input=$samplesheet comparison=$comparison count.dir=counts gene.txt=$gene_txt length.col=${length_col}
    rm -r counts/
    ln -s y0.rds transcript.y0.rds
    ln -s pca.rds transcript.pca.rds
    """
}
