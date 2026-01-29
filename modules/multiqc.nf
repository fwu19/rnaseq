process MULTIQC {
    label "process_medium"

    container = 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'
    //module = ['MultiQC/1.21-foss-2023a']

    tag "MultiQC on all samples"

    publishDir "${params.outdir}/MultiQC/", mode: 'copy'

    input:
    path ( multiqc_config)
    path ( 'fastqc/*' )
    path ( 'cutadapt/*' )
    path ( 'fastp/*' )
    path ( 'fastqc_trimmed/*' )
    path ( 'star_log/*' )
    path ( 'rnaseqc/*' )
    path ( 'rseqc/*' )
    path ( 'samtools/*' )
    path ( 'samtools_host/*' )
    path ( 'samtools_xeno/*' )
    path ( 'gatk/*' )
    path ( 'star_count/*' )
    path ( 'featureCounts/*' )
    path ( 'salmon/*' )

    output:
    path ('multiqc_data/'), emit: data
    path ('multiqc_data/', type: 'dir' )
    path ('multiqc_report.html')
    path ('versions.yml'), emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    [[ ${params.workflow} == 'pdx' ]] && mv samtools/ samtools_graft/
    multiqc -f --ignore _STARpass1/ --config $multiqc_config .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version 2>&1 | sed "s/.* //g" )
    END_VERSIONS
    """
}
