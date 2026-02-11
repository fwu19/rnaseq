process MULTIQC {
    label "process_medium"

    container = 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'
    //module = ['MultiQC/1.21-foss-2023a']

    tag "MultiQC on all samples"

    publishDir "${params.outdir}/MultiQC/", mode: 'copy'

    input:
    path ( multiqc_config)
    path fastqc_files,       stageAs: 'fastqc/*'
    path cutadapt_files,     stageAs: 'cutadapt/*'
    path fastp_files,        stageAs: 'fastp/*'
    path fastqc_trimmed_files, stageAs: 'fastqc_trimmed/*'
    path star_log_files,     stageAs: 'star_log/*/*'
    path star_log_host_files,     stageAs: 'star_log_host/*/*'
    path rnaseqc_files,      stageAs: 'rnaseqc/*'
    path rseqc_files,        stageAs: 'rseqc/*'
    path bam_stat_files,     stageAs: 'samtools/*'
    path bam_stat_host_files, stageAs: 'samtools_host/*'
    path bam_stat_xeno_files, stageAs: 'samtools_xeno/*'
    path hs_metrics_files,   stageAs: 'gatk/*'
    path star_count_files,   stageAs: 'star_count/*/*'
    path featureCounts_files, stageAs: 'featureCounts/*'
    path salmon_files,       stageAs: 'salmon/*'

    output:
    path ('multiqc_data/'), emit: data
    path ('multiqc_data/', type: 'dir' )
    path ('multiqc_report.html')
    path ('versions.yml'), emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    [[ ${params.workflow} == 'pdx' ]] && mv star_log/ star_log_graft/
    [[ ${params.workflow} == 'pdx' ]] && mv samtools/ samtools_graft/
    multiqc -f --ignore _STARpass1/ --config $multiqc_config .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version 2>&1 | sed "s/.* //g" )
    END_VERSIONS
    """
}
