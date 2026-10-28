process MULTIQC {
    label "process_medium"

    //container = 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'
    //module = ['MultiQC/1.21-foss-2023a']
    container 'multiqc/multiqc:v1.29'

    tag "MultiQC on all samples"

    publishDir "${params.outdir}/MultiQC/", mode: 'copy'

    input:
    path ( multiqc_config)
    path fastqc_files,       stageAs: 'fastqc/*'
    path cutadapt_files,     stageAs: 'cutadapt/*'
    path fastp_files,        stageAs: 'fastp/*'
    path trimgalore_files,    stageAs: 'trimgalore/*'
    path fastqc_trimmed_files, stageAs: 'fastqc_trimmed/*'
    path xengsort_log_files, stageAs: 'xengsort/*'
    path star_log_files,     stageAs: 'star_log/*/*'
    path star_log_host_files,     stageAs: 'star_log_host/*/*'
    path rnaseqc_files,      stageAs: 'rnaseqc/*'
    path rnaseqc_host_files,      stageAs: 'rnaseqc_host/*'
    path rseqc_files,        stageAs: 'rseqc/*'
    path rseqc_host_files,        stageAs: 'rseqc_host/*'
    path bam_stat_files,     stageAs: 'samtools/*'
    path bam_stat_host_files, stageAs: 'samtools_host/*'
    path bam_stat_xeno_files, stageAs: 'samtools_xeno/*'
    path hs_metrics_files,   stageAs: 'gatk/*'
    path star_count_files,   stageAs: 'star_count/*/*'
    path star_count_host_files,   stageAs: 'star_count_host/*/*'
    path featureCounts_files, stageAs: 'featureCounts/*'
    path featureCounts_host_files, stageAs: 'featureCounts_host/*'
    path salmon_files,       stageAs: 'salmon/*'
    path salmon_host_files,       stageAs: 'salmon_host/*'

    output:
    path ('multiqc_data/'), emit: data
    path ('multiqc_data/', type: 'dir' )
    path ('multiqc_report.html')
    path ('versions.yml'), emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    [[ ${params.workflow} == 'pdx' && -d star_log ]] && mv star_log star_log_graft
    [[ ${params.workflow} == 'pdx' && -d star_count ]] && mv star_count star_count_graft
    [[ ${params.workflow} == 'pdx' && -d samtools ]] && mv samtools samtools_graft
    [[ ${params.workflow} == 'pdx' && -d rnaseqc ]] && mv rnaseqc rnaseqc_graft
    [[ ${params.workflow} == 'pdx' && -d rseqc ]] && mv rseqc rseqc_graft
    [[ ${params.workflow} == 'pdx' && -d featureCounts ]] && mv featureCounts featureCounts_graft
    [[ ${params.workflow} == 'pdx' && -d salmon ]] && mv salmon salmon_graft

    TMPDIR=./ multiqc -f --ignore _STARpass1/ --config $multiqc_config .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version 2>&1 | sed "s/.* //g" )
    END_VERSIONS
    """
}
