params.outdir = 'Analysis/'

process RSEQC {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['R/4.1.2-foss-2020b', 'SAMtools/1.11-GCC-10.2.0']


    tag "RSeQC on ${bam_path.baseName}"

    publishDir "${outdir}/RSeQC/", mode: 'copy'

   input:
    tuple val(sample_id), path(bam_path) 
    path(rseqcBed) 
    path(txBed) 
    path(geneBed) 
    path(outdir) 

    output:
    tuple path("*")

    script:
    """
    rseqc.sh $sample_id $bam_path $rseqcBed $txBed $geneBed 

    """
}

