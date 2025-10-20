
process MAKE_STAR {
    label "process_high"
    module = ['STAR/2.7.10b-GCC-12.2.0', 'SAMtools/1.17-GCC-12.2.0']


    tag "Make STAR index"

    publishDir "${params.outdir}/references/$genome/STAR/", mode: 'copy'

    input:
    path ("genome/*") 
    path ("genes/*")
    
    output:
    path ( "STAR2Index/", emit: dir)
    path  ("versions.yml", emit: versions)
    path ( "STAR2Index/", type: 'dir')

    when:
    task.ext.when == null || task.ext.when

    script:
    def overhang = params.sjdbOverhang ? "--sjdbOverhang ${params.sjdbOverhang}" : ""
    def base = params.genomeSAindexNbases ? "--genomeSAindexNbases ${params.genomeSAindexNbases}" : ""
    def args = task.ext.args ?: ""
    """
    make_star.sh ${task.cpus} genome/ genes/ $overhang $base $args 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR: \$(STAR --version | head -n 1)
    END_VERSIONS
    """

    stub:
    """
    mkdir STAR2Index
    touch STAR2Index/Genome
    touch versions.yml
    """
}
