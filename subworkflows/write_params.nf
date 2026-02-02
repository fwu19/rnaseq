/*
* Write user-defined params different from default params to a json file
*/

include { WRITE_JSON } from '../modules/write_json.nf'

workflow WRITE_PARAMS {

    main:

    /*
     * Build BASE_PARAMS from two config files:
     */

    def BASE_PARAMS = [:]

    [ 'conf/params_default.config',
      'conf/params_flowswitch_default.config'  // second file overrides first
    ].each { relPath ->

        def cfgFile = new File("${projectDir}/${relPath}")
        if( !cfgFile.exists() )
            return

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
                    def key   = parts[0].trim()
                    def value = parts[1].trim()
                    BASE_PARAMS[key] = value    // later files overwrite earlier
                }
            }
        }
    }

    // 2) Compare with final params
    def FINAL_PARAMS     = params
    def OVERRIDDEN_PARAMS = [:]

    // same keys, changed values
    BASE_PARAMS.each { key, defaultVal ->
        def finalVal = FINAL_PARAMS.containsKey(key) && FINAL_PARAMS[key] != null \
            ? FINAL_PARAMS[key]
            : defaultVal
        if( finalVal != defaultVal ) {
            OVERRIDDEN_PARAMS[key] = finalVal
        }
    }

    // extra keys only in FINAL_PARAMS
    FINAL_PARAMS.each { key, v ->
        if( !BASE_PARAMS.containsKey(key) ) {
            OVERRIDDEN_PARAMS[key] = v
        }
    }

    Channel
        .value(OVERRIDDEN_PARAMS)
        .set { overridden_params_ch }

    WRITE_JSON(
        overridden_params_ch,
        "overridden_params.json"
    )
}