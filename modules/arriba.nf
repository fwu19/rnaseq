
process ARRIBA {
    label "process_high"
    module = [ 'Arriba/2.4.0-GCC-12.2.0', 'SAMtools/1.17-GCC-12.2.0' ]


    tag "Arriba on ${out_prefix}"

    publishDir "${params.outdir}/arriba/${params.genome}/", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(read1), path(read2)
    val (genome)
    path (star_index)
    path (gtf)
    path (genome_fa)
    path (blacklist)
    path (known_fusions)
    path (protein_domains)

    
    output:
    path("${out_prefix}/", type: 'dir')
    path("${out_prefix}.{bam,bam.bai,fusions.tsv}")
    path ('versions.yml'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    """
    mkdir ${out_prefix}
    cd ${out_prefix}
    run_arriba.sh ../${star_index} ../${gtf} ../${genome_fa} ../${blacklist} ../${known_fusions} ../${protein_domains} ${task.cpus} ../${read1} ../${read2}
    mv Aligned.sortedByCoord.out.bam ../${out_prefix}.bam
    mv Aligned.sortedByCoord.out.bam.bai ../${out_prefix}.bam.bai
    mv fusions.tsv ../${out_prefix}.fusions.tsv 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba: \$( arriba -h | head -n 1 | sed -e "s/.* //g" )
        samtools: \$( samtools --version | head -n 1 | sed -e "s/.* //g" )
    END_VERSIONS
    """
}
