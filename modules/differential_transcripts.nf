
process DIFFERENTIAL_TRANSCRIPTS {
    time = '1d'
    cpus = 1
    memory = '24G'
    container "docker://fwu19/r-libs:4.1.2" 


    tag "Differential transcripts"

    publishDir "${params.outdir}/differential_transcripts/", mode: 'copy'
    
    input:
    path(samplesheet)
    path(comparison)
    path(rds)
    val(length_col)
    val(fdr)
    val(fc)
    val(fdr2)
    val(fc2)
    path(gene_txt)
    val(gene_type)

    output:
    path ("differential_transcripts.rds"), emit: rds, optional: true
    path ("*"), optional: true
    
    script:
    """
    differential_transcripts.r input=$samplesheet comparison=$comparison rds=${rds} length.col=${length_col} fdr=${fdr} fc=${fc} fdr2=${fdr2} fc2=${fc2} gene_txt=$gene_txt gene_type=$gene_type
    """
}
