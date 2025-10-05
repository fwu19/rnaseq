/*
* Read input and generate samplesheet
*/

include { GET_FASTQ_PATHS } from '../modules/get_fastq_paths'
include { CHECK_INPUT  } from '../modules/check_input.nf'

workflow GET_INPUT {
    take:
    input
    input_dir
    metadata

    main: 
    fq = "$projectDir/assets/dummy_file.csv"
    samplesheet = "$projectDir/assets/dummy_file.csv"

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
        }

        CHECK_INPUT(
            file( input, checkIfExists: true ),
            file( metadata, checkIfExists: true )
        )
        samplesheet = CHECK_INPUT.out.csv
        fq = CHECK_INPUT.out.fq
        
    }



    emit:
    samplesheet = samplesheet
    fq = fq

}
