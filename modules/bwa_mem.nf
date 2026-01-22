
process BWA_MEM {
    label "process_high"

    container = 'community.wave.seqera.io/library/bwa_htslib_samtools:56c9f8d5201889a4'
    //module = ['BWA/0.7.17-GCCcore-12.2.0', 'SAMtools/1.17-GCC-12.2.0']

    tag "BWA_MEM on ${out_prefix}"

    publishDir "${params.outdir}/BWA/$genome/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(read1), path(read2)
    val (genome)
    path (bwa_index)
    
    output:
    tuple val(meta), val(out_prefix), path( "${out_prefix}.bam" ), emit: bam 
    tuple val(meta), val(out_prefix), path( "${out_prefix}.bam.bai" ), emit: bai
    path("${out_prefix}.{bam,bam.bai}")
    path ('versions.yml'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    """
    bwa_mem.sh ${task.cpus} ${bwa_index} ${read1} ${read2} ${out_prefix}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$( bwa 2>&1 | sed -n 3,3p | sed -e "s/.* //g" )
        samtools: \$( samtools --version | head -n 1 | sed -e "s/.* //g" )
    END_VERSIONS
    """
}
