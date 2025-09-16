process MULTIQC_PDX {
    time = '1d'
    cpus = 6
    memory = '36G'
    //module = ['MultiQC/1.21-foss-2023a']
    container "quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0" // multiqc 1.21

    tag "MultiQC on all samples"

    publishDir "${params.outdir}/MultiQC/", mode: 'copy'

    input:
    path ( multiqc_config )
    path ( 'fastqc/*' )
    path ( 'cutadapt/*' )
    path ( 'fastqc_trimmed/*' )
    path ( 'rseqc/*' )
    path ( 'rnaseqc/*' )
    path ( 'gatk/*' )

    path ( 'star_log_graft/*' )
    path ( "star_log_host/*" )

    path ( "samtools_graft/*" )
    path ( "samtools_host/*" )
    path ( "samtools_xeno/*" )

    output:
    path ('multiqc_data/*'), emit: data
    path ('multiqc_data/', type: 'dir' )
    path ('multiqc_report.html')

    script:
    def args = task.ext.args ?: ''
    """
    multiqc -f --ignore _STARpass1/ --config $multiqc_config .
    """
}