process TRIMGALORE {
    label "process_high"

    container = 'ghcr.io/felixkrueger/trimgalore:v2.3.0'

    tag "TRIMGALORE on ${meta.id}"

    publishDir "${params.outdir}/trim_galore/", pattern: '*.gz', mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(read1), path(read2)
    
    output:
    tuple val(meta), val(out_prefix), path("${out_prefix}_val_1.fq.gz"), path("${out_prefix}_val_2.fq.gz"), emit: fq
    path( "*.{json,txt}" ), emit: js
    path ('versions.yml'), emit: versions

    script:
    def args = task.ext.args ?: ""
    if (params.adapters){
        def adapter_list = params.adapters.split(',').collect()
        if ( adapter_list.size == 1){
        adapter1 = adapter_list[0]
        adapter2 = adapter_list[0]
        } else {
        adapter1 = adapter_list[0]
        adapter2 = adapter_list[1]
        }
        """
        trim_galore --cores ${task.cpus} \
        $args \
        -a $adapter1 -a2 $adapter2 \
        --output_dir ./ --basename ${out_prefix} \
        --paired $read1 $read2

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trim_galore: \$( trim_galore -v 2>&1 | sed -n 4,4p | sed -e "s/.* //g" )
        END_VERSIONS
        """

    }else{
        """
        trim_galore --cores ${task.cpus} \
        $args --output_dir ./ --basename ${out_prefix} \
        --paired $read1 $read2

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trim_galore: \$( trim_galore -V | sed -e "s/.* //g" )
        END_VERSIONS
        """
    }
}
