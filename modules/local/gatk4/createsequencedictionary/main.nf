process GATK4_CREATESEQUENCEDICTIONARY {
    tag "${fasta.name}"
    label 'process_single'
    container 'broadinstitute/gatk:4.6.1.0'

    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' || !params.save_reference ? null : fn }

    input:
    path fasta

    output:
    path "*.dict"       , emit: dict
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        CreateSequenceDictionary \\
        --REFERENCE ${fasta} \\
        --OUTPUT ${fasta.baseName}.dict \\
        --TMP_DIR .

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta.baseName}.dict

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
