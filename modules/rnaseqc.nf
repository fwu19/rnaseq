params.outdir = 'Analysis/'
params.gtfQC = '/fh/fast/_SR/Genomics/proj/fwu/reference/hg38/genes/gencode.v38.annotation.collapsed.gtf'

process RNASEQC {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['RNA-SeQC/2.3.4-foss-2019b', 'SAMtools/1.11-GCC-10.2.0']


    tag "RNA-SeQC on ${bam_path.baseName}"

    publishDir "${outdir}/RNA-SeQC/", mode: 'copy'

   input:
    tuple val(sample_id), path(bam_path)
    path(gtfQC)
    val(strand)
    val(readType)
    path(outdir)

    output:
    tuple path("*")

    script:
    """
    rnaseqc.sh $sample_id $bam_path $gtfQC $strand $readType

    """
}


	