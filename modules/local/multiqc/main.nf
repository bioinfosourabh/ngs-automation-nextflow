process MULTIQC {
    label 'process_single'
    container 'quay.io/biocontainers/multiqc:1.30--pyhdfd78af_1'

    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    path multiqc_files, stageAs: "inputs/?/*"
    path multiqc_config

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data"       , emit: data
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def config_cmd = multiqc_config ? "--config ${multiqc_config}" : ''
    """
    multiqc --force ${config_cmd} ${args} inputs

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/^multiqc, version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_report.html

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: 1.30
    END_VERSIONS
    """
}
