//
// Parameter validation, run summary and software-version collation.
// Written without plugins so that the pipeline runs offline.
//

def validateParams(String schema_path) {
    def errors = []

    // Unknown parameters and enum / type checks against nextflow_schema.json.
    def schema = new groovy.json.JsonSlurper().parse(new File(schema_path))
    def props = [:]
    schema['$defs'].each { _group, definition -> props.putAll(definition.properties) }
    params.each { name, _value ->
        if (!props.containsKey(name) && !(name.toString().contains('-'))) {
            log.warn("Parameter --${name} is not defined in nextflow_schema.json and will be ignored")
        }
    }
    props.each { name, spec ->
        def value = params[name]
        if (value == null) {
            return
        }
        def allowed = spec['enum']
        if (allowed && !(value in allowed)) {
            errors.add("--${name} must be one of ${allowed.join(', ')} (got '${value}')")
        }
        if (spec.type == 'integer' && !(value instanceof Integer)) {
            errors.add("--${name} must be an integer (got '${value}')")
        }
        if (spec.type == 'boolean' && !(value instanceof Boolean)) {
            errors.add("--${name} must be true or false (got '${value}')")
        }
    }

    // Required inputs and conditional requirements.
    ['input', 'outdir', 'fasta', 'dbsnp'].each { p ->
        if (!params[p]) {
            errors.add("--${p} is required")
        }
    }
    if (params.scatter_count instanceof Integer && params.scatter_count < 1) {
        errors.add("--scatter_count must be >= 1")
    }
    if (params.run_germline && params.germline_filter == 'vqsr') {
        ['vqsr_hapmap', 'vqsr_omni', 'vqsr_1000g', 'vqsr_mills'].each { p ->
            if (!params[p]) {
                errors.add("--${p} is required with --germline_filter vqsr")
            }
        }
    }
    if (params.run_somatic && !params.germline_resource) {
        errors.add("--germline_resource (e.g. af-only-gnomad.hg38.vcf.gz) is required when --run_somatic is true")
    }
    if (params.run_vep && !params.vep_cache) {
        errors.add("--vep_cache is required when --run_vep is true")
    }
    if (!params.run_germline && !params.run_somatic) {
        errors.add("at least one of --run_germline and --run_somatic must be true")
    }
    if (errors) {
        error("Invalid parameters:\n  - " + errors.join("\n  - "))
    }

    if (params.run_somatic && !params.pon) {
        log.warn("No --pon given: Mutect2 will run without a panel of normals, which increases recurrent artefact calls")
    }
    if (params.run_somatic && !params.contamination_variants) {
        log.warn("No --contamination_variants given: cross-sample contamination will not be estimated before FilterMutectCalls")
    }
}

def pipelineSummary() {
    def keys = ['input', 'outdir', 'fasta', 'intervals', 'aligner', 'scatter_count', 'run_germline',
                'joint_caller', 'germline_filter', 'run_somatic', 'pon', 'germline_resource', 'run_vep']
    def width = keys.collect { k -> k.size() }.max()
    def lines = keys.collect { k -> "  ${k.padRight(width)} : ${params[k]}" }
    def rule = '-' * 60
    return (["${workflow.manifest.name} v${workflow.manifest.version}", rule] + lines + [rule]).join('\n')
}

def completionMessage() {
    def status = workflow.success ? 'completed successfully' : 'failed'
    return "Pipeline ${status}. Duration: ${workflow.duration}. Results: ${params.outdir}"
}
