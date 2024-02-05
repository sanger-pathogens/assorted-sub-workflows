//
// This file holds several functions used to perform JSON parameter validation, help and summary rendering.
//
// Modified from NF-Core's template: https://github.com/nf-core/tools

import org.yaml.snakeyaml.Yaml
import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.json.JsonBuilder
import nextflow.extension.FilesEx


class NextflowTool {
    //
    // Dump pipeline parameters in a json file
    //
    public static void dump_parameters(workflow, params) {
        def timestamp  = new java.util.Date().format( 'yyyy-MM-dd_HH-mm-ss')
        def filename   = "params_${timestamp}.json"
        def temp_pf    = new File(workflow.launchDir.toString(), ".${filename}")
        def jsonStr    = JsonOutput.toJson(params)
        temp_pf.text   = JsonOutput.prettyPrint(jsonStr)

        FilesEx.copyTo(temp_pf.toPath(), "${params.results_dir}/pipeline_info/params_${timestamp}.json")
        temp_pf.delete()
    }

    public static void help_message(pipeline_schema, schema_path_list, monochrome_logs, log) {
        Map colors = logColours(monochrome_logs)
        def indent = "      "
        //master schema is pipeline specific it only contains what needs to be there other stuff can be imported from the path_list
        //master schema contains a overwrite section which is not printed and instead used to overwrite defaults or whatever in the
        //general imported schema
        def master_in = new File(pipeline_schema).text
        def master_schema = new JsonSlurper().parseText(master_in)

        for (schema_path in schema_path_list) {
            def schema_in = new File(schema_path).text
            def json = new JsonSlurper().parseText(schema_in)

            //create a list of keys you want to supress in the json
            def banned_keylist = master_schema.overwrite_param.keySet() as Set

            json.params.each{
                if (!banned_keylist.contains(it.key)) {
                    log.info "${colors.purple} ${it.key} ${colors.reset}"
                    it.value.each {
                        if (it.key in master_schema.overwrite_param) {
                            // A groovy json path cannot be queried using a variable unless you force it using eval
                            // set up the query and then use eval.x to append the search term onto the end of the json object used to overwrite
                            // this then replaces the default value with the input

                            def overwrite_path = 'overwrite_param.' + it.key

                            def overwrite_param = Eval.x( master_schema, 'x.' + overwrite_path)

                            if (overwrite_param.help_text != "") {
                                log.info indent + "--" + it.key
                                log.info indent + indent + "default: " + overwrite_param.default
                                log.info indent + indent + overwrite_param.help_text
                            }

                        } else {
                            //if nothing needs to be overwritten just print what is there
                            log.info indent + "--" + it.key
                            log.info indent + indent + "default: " + it.value.default
                            log.info indent + indent + it.value.help_text
                        }
                    }
                    //put a line to seperate
                    log.info dashedLine(monochrome_logs)
                }
            }
        }
        //finally print the params in the master manifest
        master_schema.params.each {
            log.info "${colors.purple} ${it.key} ${colors.reset}"
            it.value.each {
                if (it.key.toString().contains('header')) {
                    if (it.value.title){
                        log.info indent
                        log.info "${colors.red} ${it.value.title} ${colors.reset}"
                    }
                    log.info indent + it.value.subtext
                    log.info indent

                } else {
                log.info indent + "--" + it.key
                log.info indent + indent + "default: " + it.value.default
                log.info indent + indent + it.value.help_text
                log.info indent
                }
            }
        //put a line to seperate
        log.info dashedLine(monochrome_logs)
        }
    }

    //
    // Print pipeline summary on completion
    //
    public static void summary(workflow, params, log) {
        Map colors = logColours(params.monochrome_logs)
        if (workflow.success) {
            if (workflow.stats.ignoredCount == 0) {
                log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} Pipeline completed successfully${colors.reset}-"
            } else {
                log.info "-${colors.purple}[$workflow.manifest.name]${colors.yellow} Pipeline completed successfully, but with errored process(es) ${colors.reset}-"
            }
        } else {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} Pipeline completed with errors${colors.reset}-"
        }
    }

    //
    // ANSII Colours used for terminal logging
    //
    public static Map logColours(Boolean monochrome_logs) {
        Map colorcodes = [:]

        // Reset / Meta
        colorcodes['reset']      = monochrome_logs ? '' : "\033[0m"
        colorcodes['bold']       = monochrome_logs ? '' : "\033[1m"
        colorcodes['dim']        = monochrome_logs ? '' : "\033[2m"
        colorcodes['underlined'] = monochrome_logs ? '' : "\033[4m"
        colorcodes['blink']      = monochrome_logs ? '' : "\033[5m"
        colorcodes['reverse']    = monochrome_logs ? '' : "\033[7m"
        colorcodes['hidden']     = monochrome_logs ? '' : "\033[8m"

        // Regular Colors
        colorcodes['black']      = monochrome_logs ? '' : "\033[0;30m"
        colorcodes['red']        = monochrome_logs ? '' : "\033[0;31m"
        colorcodes['green']      = monochrome_logs ? '' : "\033[0;32m"
        colorcodes['yellow']     = monochrome_logs ? '' : "\033[0;33m"
        colorcodes['blue']       = monochrome_logs ? '' : "\033[0;34m"
        colorcodes['purple']     = monochrome_logs ? '' : "\033[0;35m"
        colorcodes['cyan']       = monochrome_logs ? '' : "\033[0;36m"
        colorcodes['white']      = monochrome_logs ? '' : "\033[0;37m"

        // Bold
        colorcodes['bblack']     = monochrome_logs ? '' : "\033[1;30m"
        colorcodes['bred']       = monochrome_logs ? '' : "\033[1;31m"
        colorcodes['bgreen']     = monochrome_logs ? '' : "\033[1;32m"
        colorcodes['byellow']    = monochrome_logs ? '' : "\033[1;33m"
        colorcodes['bblue']      = monochrome_logs ? '' : "\033[1;34m"
        colorcodes['bpurple']    = monochrome_logs ? '' : "\033[1;35m"
        colorcodes['bcyan']      = monochrome_logs ? '' : "\033[1;36m"
        colorcodes['bwhite']     = monochrome_logs ? '' : "\033[1;37m"

        // Underline
        colorcodes['ublack']     = monochrome_logs ? '' : "\033[4;30m"
        colorcodes['ured']       = monochrome_logs ? '' : "\033[4;31m"
        colorcodes['ugreen']     = monochrome_logs ? '' : "\033[4;32m"
        colorcodes['uyellow']    = monochrome_logs ? '' : "\033[4;33m"
        colorcodes['ublue']      = monochrome_logs ? '' : "\033[4;34m"
        colorcodes['upurple']    = monochrome_logs ? '' : "\033[4;35m"
        colorcodes['ucyan']      = monochrome_logs ? '' : "\033[4;36m"
        colorcodes['uwhite']     = monochrome_logs ? '' : "\033[4;37m"

        // High Intensity
        colorcodes['iblack']     = monochrome_logs ? '' : "\033[0;90m"
        colorcodes['ired']       = monochrome_logs ? '' : "\033[0;91m"
        colorcodes['igreen']     = monochrome_logs ? '' : "\033[0;92m"
        colorcodes['iyellow']    = monochrome_logs ? '' : "\033[0;93m"
        colorcodes['iblue']      = monochrome_logs ? '' : "\033[0;94m"
        colorcodes['ipurple']    = monochrome_logs ? '' : "\033[0;95m"
        colorcodes['icyan']      = monochrome_logs ? '' : "\033[0;96m"
        colorcodes['iwhite']     = monochrome_logs ? '' : "\033[0;97m"

        // Bold High Intensity
        colorcodes['biblack']    = monochrome_logs ? '' : "\033[1;90m"
        colorcodes['bired']      = monochrome_logs ? '' : "\033[1;91m"
        colorcodes['bigreen']    = monochrome_logs ? '' : "\033[1;92m"
        colorcodes['biyellow']   = monochrome_logs ? '' : "\033[1;93m"
        colorcodes['biblue']     = monochrome_logs ? '' : "\033[1;94m"
        colorcodes['bipurple']   = monochrome_logs ? '' : "\033[1;95m"
        colorcodes['bicyan']     = monochrome_logs ? '' : "\033[1;96m"
        colorcodes['biwhite']    = monochrome_logs ? '' : "\033[1;97m"

        return colorcodes
    }

    //
    // Does what is says on the tin
    //
    public static String dashedLine(monochrome_logs) {
        Map colors = logColours(monochrome_logs)
        return "-${colors.dim}----------------------------------------------------${colors.reset}-"
    }
    
    //
    // pam-info logo
    //
    public static String logo(workflow, monochrome_logs) {
        Map colors = logColours(monochrome_logs)
        String.format(
            """\n
            ${dashedLine(monochrome_logs)}
            ${colors.blue} _____________________  ___     ____________   _________________ ${colors.reset}
            ${colors.blue} ___  __ \\__    |__   |/  /     ____  _/__  | / /__  ____/_  __ \\ ${colors.reset}
            ${colors.blue} __  /_/ /_  /| |_  /|_/ /_________  / __   |/ /__  /_   _  / / / ${colors.reset}
            ${colors.blue} _  ____/_  ___ |  /  / /_/_____/_/ /  _  /|  / _  __/   / /_/ / ${colors.reset}
            ${colors.blue} /_/     /_/  |_/_/  /_/        /___/  /_/ |_/  /_/      \\____/ ${colors.reset}                                         
            \n
            ${colors.purple}  ${workflow.manifest.name} ${workflow.manifest.version} ${colors.reset}
            ${dashedLine(monochrome_logs)}
            """.stripIndent()
        )
    }
}