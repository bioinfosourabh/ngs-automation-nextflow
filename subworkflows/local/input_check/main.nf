//
// Samplesheet parsing and validation.
//
// Expected columns (header required, order free):
//   patient, sample, status (0 = normal, 1 = tumour; optional, default 0),
//   lane (optional), fastq_1, fastq_2
//
// Relative FASTQ paths are resolved against the directory of the samplesheet.
//

def resolveFastq(String value, csv, int line) {
    if (!value) {
        error("Samplesheet line ${line}: FASTQ path is empty")
    }
    if (!(value ==~ /.+\.f(ast)?q\.gz$/)) {
        error("Samplesheet line ${line}: '${value}' must end in .fastq.gz or .fq.gz")
    }
    def f = (value ==~ /^[a-zA-Z][a-zA-Z0-9+.-]*:\/\/.*/ || value.startsWith('/')) ? file(value) : csv.parent.resolve(value)
    if (!f.exists()) {
        error("Samplesheet line ${line}: file not found: ${f}")
    }
    return f
}

def readSamplesheet(String path, String platform) {
    def csv = file(path, checkIfExists: true)
    def rows = csv.splitCsv(header: true, strip: true)
    if (!rows) {
        error("Samplesheet ${path} has no data rows")
    }
    def required = ['patient', 'sample', 'fastq_1', 'fastq_2']
    def missing = required.findAll { c -> !rows[0].containsKey(c) }
    if (missing) {
        error("Samplesheet ${path} is missing required column(s): ${missing.join(', ')}")
    }

    def records = []
    def seen = [:]
    rows.eachWithIndex { row, i ->
        def line = i + 2
        def sample  = row.sample
        def patient = row.patient
        if (!(sample ==~ /^[A-Za-z0-9._-]+$/)) {
            error("Samplesheet line ${line}: sample '${sample}' must be non-empty and contain only letters, digits, '.', '_' or '-'")
        }
        if (!(patient ==~ /^[A-Za-z0-9._-]+$/)) {
            error("Samplesheet line ${line}: patient '${patient}' must be non-empty and contain only letters, digits, '.', '_' or '-'")
        }
        def status_str = row.status ?: '0'
        if (!(status_str in ['0', '1'])) {
            error("Samplesheet line ${line}: status must be 0 (normal) or 1 (tumour), found '${status_str}'")
        }
        def n_prev = records.count { r -> r.sample == sample }
        def lane = row.lane ?: "${n_prev + 1}"
        def key = "${sample}|${lane}"
        if (seen.containsKey(key)) {
            error("Samplesheet line ${line}: sample '${sample}' lane '${lane}' already listed on line ${seen[key]}")
        }
        seen[key] = line
        def r1 = resolveFastq(row.fastq_1, csv, line)
        def r2 = resolveFastq(row.fastq_2, csv, line)
        records.add([sample: sample, patient: patient, status: status_str.toInteger(), lane: lane, r1: r1, r2: r2, line: line])
    }

    // Consistency checks across rows of the same sample and patient.
    records.groupBy { r -> r.sample }.each { sample, rs ->
        if (rs*.patient.unique().size() > 1 || rs*.status.unique().size() > 1) {
            error("Samplesheet: sample '${sample}' has rows with different patient or status values")
        }
    }
    def normals_by_patient = [:]
    records.findAll { r -> r.status == 0 }.each { r ->
        def known = normals_by_patient[r.patient]
        if (known && known != r.sample) {
            error("Samplesheet: patient '${r.patient}' has more than one normal sample (${known}, ${r.sample}); use one normal per patient")
        }
        normals_by_patient[r.patient] = r.sample
    }
    def tumour_patients = records.findAll { r -> r.status == 1 }*.patient.unique()

    return records.collect { r ->
        def n_lanes = records.count { x -> x.sample == r.sample }
        def meta = [
            id        : "${r.sample}_${r.lane}",
            sample    : r.sample,
            patient   : r.patient,
            status    : r.status,
            lane      : r.lane,
            num_lanes : n_lanes,
            normal_id : r.status == 1 ? normals_by_patient[r.patient] : null,
            has_tumour: r.status == 0 && (r.patient in tumour_patients),
            read_group: "@RG\\tID:${r.sample}.${r.lane}\\tSM:${r.sample}\\tLB:${r.sample}\\tPL:${platform}\\tPU:${r.sample}.${r.lane}"
        ]
        [meta, [r.r1, r.r2]]
    }
}
