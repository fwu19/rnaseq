
process MAKE_STAR {
    label "process_medium"
    module = ['STAR/2.7.10b-GCC-12.2.0', 'SAMtools/1.17-GCC-12.2.0']


    tag "Make STAR index"

    publishDir "${params.outdir}/references/$genome/STAR/", mode: 'copy'

    input:
    path ("genome/*") 
    path ("genes/*")
    
    output:
    path ( "STAR2Index/", emit: dir)
    path  ("versions.yml", emit: versions)

    when:
    task.ext.when == null || task.ext.when

    script:
    def overhang = params.sjdb_overhang ? "--sjdbOverhang ${params.sjdb_overhang}" : ""
    def base = params.genome_sa_index_bases ? "--genomeSAindexNbases ${params.genome_sa_index_bases}" : ""
    def args = task.ext.args ?: ""
    """
    make_star.sh ${task.cpus} genome/ genes/ $overhang $base $args 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR: \$(STAR --version | head -n 1)
        samtools: \$( samtools --version | head -n 1 | sed -e "s/.* //g" )
    END_VERSIONS
    
    """

    stub:
    """
    mkdir STAR2Index
    touch STAR2Index/Genome
    touch versions.yml
    """
}
