
process STAR {
    label "process_medium"

    container = 'community.wave.seqera.io/library/bedtools_samtools_star:4c1d7f700be70377'
    //container = 'biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0'
    //module = ['STAR/2.7.10b-GCC-12.2.0', 'SAMtools/1.17-GCC-12.2.0']

    tag "STAR on ${out_prefix}"

    publishDir "${params.outdir}/STAR/$genome/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(read1), path(read2)
    val (genome)
    path (star_index)
    path (gtf)
    
    output:
    tuple val(meta), val(out_prefix), path( "${out_prefix}.bam" ), path( "${out_prefix}.bam.bai" ), emit: bam_bai
    tuple val(meta), val(out_prefix), path( "${out_prefix}/${out_prefix}.Aligned.toTranscriptome.out.bam" ),emit: tx_bam 
    tuple val(meta), val(out_prefix), path( "${out_prefix}/${out_prefix}.ReadsPerGene.out.tab" ), emit: counts
    path( "${out_prefix}.Log.final.out" ), emit: log
    path("${out_prefix}/", type: 'dir')
    path ('versions.yml'), emit: versions

    script:
    def args = task.ext.args ?: ""
    def overhang = params.sjdb_overhang ? "--sjdbOverhang ${params.sjdb_overhang}" : ""
    def outSAMattributes = params.out_sam_attributes ? "--outSAMattributes ${params.out_sam_attributes}" : ""
    def rg = task.ext.rg ?: "none"
    """
    star.sh ${out_prefix} ${star_index} ${gtf} ${task.cpus} $read1 $read2 "$rg" $args $overhang $outSAMattributes

    mv ${out_prefix}/Aligned.sortedByCoord.out.bam ${out_prefix}.bam
    mv ${out_prefix}/Aligned.sortedByCoord.out.bam.bai ${out_prefix}.bam.bai
    mv ${out_prefix}/ReadsPerGene.out.tab ${out_prefix}/${out_prefix}.ReadsPerGene.out.tab
    mv ${out_prefix}/Aligned.toTranscriptome.out.bam ${out_prefix}/${out_prefix}.Aligned.toTranscriptome.out.bam
    cp ${out_prefix}/Log.final.out ${out_prefix}.Log.final.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR: \$(STAR --version | head -n 1)
        samtools: \$( samtools --version | head -n 1 | sed -e "s/.* //g" )
    END_VERSIONS

    """
}
