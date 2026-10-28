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

// Genome type per genome (preqc name -> 'mag'|'isolate') for MIXED runs, read from the --genome_types TSV
// (header: genome<TAB>genome_type). Genomes not listed use --genome_type. bin/report.py relies on this validation.
def genomeTypeMap() {
    def genome_types = ['mag', 'isolate']
    if (!(params.genome_type in genome_types)) { error "--genome_type must be mag or isolate, found '${params.genome_type}'" }
    if (!params.genome_types) { return [:] }

    def lines = file(params.genome_types, checkIfExists: true).readLines()
    if (lines[0].trim() != "genome\tgenome_type") {
        error "--genome_types ${params.genome_types}: header must be 'genome<TAB>genome_type'"
    }
    def types = [:]
    lines.eachWithIndex { line, i ->
        if (i == 0 || !line.trim()) { return }
        def fields = line.split('\t', -1)*.trim()
        def where = "--genome_types ${params.genome_types} line ${i + 1}"
        if (fields.size() != 2) { error "${where}: expected 2 tab-separated columns, found ${fields.size()}" }
        def genome = fields[0]
        def type = fields[1].toLowerCase()
        if (!(type in genome_types)) { error "${where}: genome_type must be mag or isolate, found '${fields[1]}'" }
        if (types.containsKey(genome)) { error "${where}: genome '${genome}' is listed twice" }
        types[genome] = type
    }
    return types
}

// Genome type of a single genome: --genome_types entry if listed, otherwise --genome_type
def genomeType(Map types, String name) {
    return types.getOrDefault(name, params.genome_type)
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
