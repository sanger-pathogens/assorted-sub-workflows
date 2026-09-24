// Where each process publishes, so every output lands under results/<species>/.
//
// Species-wide items carry meta [ID: <species>]; per-group items carry
// meta [ID: <group>, species: <species>, ...].

// results/<species>
def species_outdir(meta) {
    return "${params.outdir}/${meta.species ?: meta.ID}"
}

// results/<species>/index/<tool_step>/species, or .../groups/<group> for a group's
// candidate index, e.g. index_outdir(meta, 'sbwt/build').
def index_outdir(meta, tool_step) {
    return "${species_outdir(meta)}/index/${tool_step}/${meta.species ? "groups/${meta.ID}" : 'species'}"
}
