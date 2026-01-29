/*
* Read input and generate samplesheet
*/

include { GET_FASTQ_PATHS } from '../modules/get_fastq_paths'
include { GET_INPUT  } from '../modules/get_input.nf'

workflow CHECK_INPUT {
    take:
    input
    input_dir
    metadata

    main: 
    fq = "$projectDir/assets/dummy_file.csv"
    samplesheet = "$projectDir/assets/dummy_file.csv"
    ch_versions = Channel.empty()

    if (params.run_input_check){
        if ( input_dir =~ 'dummy' ){
            if ( input =~ 'dummy' ){
                exit 1, 'Need to provide --input or --input_dir!'
            }
        }else {
            GET_FASTQ_PATHS (
                file(input_dir, checkIfExists: true)
            )
            input = GET_FASTQ_PATHS.out.csv
            ch_versions = ch_versions.mix(GET_FASTQ_PATHS.out.versions.first())
        }
        
        GET_INPUT(
            input,
            metadata
        )
        samplesheet = GET_INPUT.out.csv
        fq = GET_INPUT.out.fq
        ch_versions = ch_versions.mix(GET_INPUT.out.versions.first())
        
    }



    emit:
    samplesheet = samplesheet
    fq = fq
    versions = ch_versions

}
