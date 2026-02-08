nextflow.enable.dsl=2

        /*
         * Helper to normalize config values
         */
        
        def parseConfigValue(String raw) {
            if( raw == null ) {return null}

            def v = raw.trim()

            if( (v.startsWith('"') && v.endsWith('"')) || (v.startsWith("'") && v.endsWith("'")) )
                {v = v.substring(1, v.length()-1)}

            if( v == 'true' )  {return true}
            if( v == 'false' ) {return false}

            if( v.isInteger() )    {return v.toInteger()}
            if( v.isBigDecimal() ) {return v.toBigDecimal()}
            if( v == 'null' )      {return null}

            return v
        }

include { WRITE_JSON } from '../modules/write_json.nf'
include { WRITE_JSON as WRITE_JSON_BASE } from '../modules/write_json.nf'

workflow WRITE_PARAMS {

    take:
    ch_software_versions

    main:

    // Collect software versions
    ch_software_versions
            .collectFile(storeDir: "${params.info_dir}/", name: 'software_versions.yml', sort: false, newLine: true)


        /*
         * 1) Build BASE_PARAMS from config files
         */
        def BASE_PARAMS = [:]

        [ 'nextflow.config', 'conf/flowswitch.config' ]
        .each { relPath ->

            def cfgFile = new File("${projectDir}/${relPath}")
            if( !cfgFile.exists() ) return

            def lines    = cfgFile.readLines()
            def inParams = false

            lines.each { line ->
                def t = line.trim()
                if( t.startsWith('params {') ) {
                    inParams = true
                    return
                }
                if( inParams && t.startsWith('}') ) {
                    inParams = false
                    return
                }
                if( inParams && t && !t.startsWith('//') ) {
                    def parts = t.split('=', 2)
                    if( parts.size() == 2 ) {
                        def key      = parts[0].trim()
                        def rawValue = parts[1].replaceFirst(/\/\/.*/, '').trim()
                        //BASE_PARAMS[key] = rawValue
                        def parsed   = parseConfigValue(rawValue)
                        BASE_PARAMS[key] = parsed
                    }
                }
            }
        }

        /*
         * 2) Compare with final params
         */
        def FINAL_PARAMS      = params
        def OVERRIDDEN_PARAMS = [:]
        def SKIP_KEYS        = [ 'info_dir', 'input', 'input_dir', 'comparison', 'comparison_transcripts', 'report_dir', 'report_rmd', 'max_cpus', 'max_memory', 'max_time' ]

        BASE_PARAMS.each { key, defaultVal ->
            def finalVal = FINAL_PARAMS.containsKey(key) && FINAL_PARAMS[key] != null \
                ? FINAL_PARAMS[key]
                : defaultVal

            if( finalVal != defaultVal && !(key in SKIP_KEYS)) {
                OVERRIDDEN_PARAMS[key] = finalVal
            }
        }

        FINAL_PARAMS.each { key, v ->
            if( !BASE_PARAMS.containsKey(key) && !(key in SKIP_KEYS) ) {
                OVERRIDDEN_PARAMS[key] = v
            }
        }

        if (params.run_input_check && params.input){
            OVERRIDDEN_PARAMS['input'] = "${params.outdir}/csv/samplesheet.valid.csv"
        }        

        if (params.run_de && params.comparison){
            def de = de_csv.first()
            OVERRIDDEN_PARAMS['comparison'] = "${params.outdir}/csv/comparisons.differential_genes.csv"
        }

        if (params.run_dt && (params.comparison_transcripts || params.comparison)){
            def dt = dt_csv.first()
            OVERRIDDEN_PARAMS['comparison_transcripts'] = "${params.outdir}/csv/comparisons.differential_transcripts.csv"
        }

        if( params.run_report ) {
            def report = params.workflow == 'regular' ? "00_RNAseq_analysis_report.Rmd" : (params.workflow == 'exome' ? "00_RNAexome_analysis_report.Rmd" : (params.workflow == 'pdx' ? "00_PDX_RNAseq_analysis_report.Rmd" : null))
            OVERRIDDEN_PARAMS['report_rmd'] = "${params.outdir}/pipeline_info/${report}"
        }
        

        /*
         * 3) Emit JSON
         */
        Channel
            .value(OVERRIDDEN_PARAMS)
            .set { overridden_params_ch }

        WRITE_JSON(
            overridden_params_ch,
            "params.json"
        )

        /*write base json for reference
        Channel
            .value(BASE_PARAMS)
            .set { base_params_ch } 
        WRITE_JSON_BASE(
            base_params_ch,
            "base_params.json"
        )
        */
}
