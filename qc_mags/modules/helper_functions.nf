def ensureList(Object maybeListLike) {
    if (maybeListLike instanceof Collection) {
        return maybeListLike as List
    } else {
        return [maybeListLike]
    }
}

// params.fasta_ext without any leading '.'
def fastaExt() {
    return params.fasta_ext.replaceAll(/^\./, '')
}

// Genome name = fasta file name without .<fasta_ext>[.gz]; this is what the reports call preqc_genome_name
def genomeName(fasta) {
    String name = fasta.name.endsWith('.gz') ? fasta.name[0..<-3] : fasta.name
    String ext = ".${fastaExt()}"
    return name.endsWith(ext) ? name[0..<-ext.length()] : name.replaceAll(/\.[^.]*$/, '')
}

// Genome type per genome (preqc name -> 'mag'|'isolate') for MIXED runs, from --genome_types (curated
// provenance); genomes not listed use --genome_type. Validated as in bin/report.py read_genome_types.
def genomeTypeMap() {
    if (!params.genome_types) { return [:] }
    def rows = file(params.genome_types, checkIfExists: true).readLines()
    def header = rows[0].split('\t').collect { it.trim() }
    def gi = header.indexOf('genome')
    def ti = header.indexOf('genome_type')
    if (gi < 0 || ti < 0) { error "--genome_types ${params.genome_types}: needs columns 'genome' and 'genome_type'" }
    def types = rows.drop(1).findAll { it.trim() }.collectEntries { l ->
        def c = l.split('\t', -1)
        [(c[gi].trim()): c[ti].trim().toLowerCase()]
    }
    def bad = (types.values() as Set) - ['mag', 'isolate']
    if (bad) { error "--genome_types ${params.genome_types}: genome_type must be mag or isolate, found ${bad.sort()}" }
    return types
}

// Genome type of a single genome: --genome_types entry if listed, otherwise --genome_type
def genomeType(Map types, String name) {
    String run_type = (params.genome_type == 'isolate') ? 'isolate' : 'mag'
    return types.getOrDefault(name, run_type)
}

// Re-batch [meta, [fastas]] into [batch meta, [fastas]] of params.batch_size genomes each
def batchFastas(fasta_lists) {
    // Newer Nextflow (strict syntax) does not let a closure update a local variable (++i gives batch_0 for every batch)
    // and rejects the CLI param as buffer size without the cast
    def batchIndex = new java.util.concurrent.atomic.AtomicInteger()
    return fasta_lists
        .map { meta, fasta_list -> fasta_list }
        .flatten()
        .buffer( size: params.batch_size as int, remainder: true )
        .map { batch -> [ [ID: "batch_${batchIndex.incrementAndGet()}"], batch ] }
}
