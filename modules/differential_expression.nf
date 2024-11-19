
process DIFFERENTIAL_EXPRESSION {
    time = '1d'
    cpus = 1
    memory = '24G'
    module = ['fhR/4.1.2-foss-2021b']


    tag "Differential expression"

    publishDir "${params.outdir}/differential_expression/", mode: 'copy'
    
    input:
    path (samplesheet)
    path (comparison)
    path ("*.txt", stageAs: "counts/*")
    path (gene_txt)

    output:
    tuple path ("*")
    
    script:
    """
    differential_expression.r input=$input comparison=$comparison count.dir=counts gene.txt=$gene_txt length.col=${params.length_col} strand=${params.strand} fdr=${params.fdr} fc=${params.fc} fdr2=${params.fdr2} fc2=${params.fc2}
    rm -r counts/
    """
}
