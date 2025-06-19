def translateKey(in_key) {
    switch(in_key){ 
        case 'studyid':
        return 'study_id'

        case 'runid':
        return 'id_run'

        case 'laneid':
        return 'lane'

        case 'plexid':
        return 'tag_index'

        //if none of the cases return in_key
        default:
        return in_key
    }
}

def avuIdQuery(meta_query) {
    /**
    * Constructs a condensed JSON AVU query array string from a single metadata map or a list of metadata maps.
    *
    * Groups input metadata by key and performs value consolidation:
    * - If a key maps to a single unique value, emits a simple equality AVU clause: {"a": "<attr>", "v": "<value>"}.
    * - If a key maps to multiple unique values, emits a set-membership clause using the "in" operator:
    *   {"a": "<attr>", "v": ["<val1>", "<val2>", ...], "o": "in"}.
    *
    * Value filtering criteria:
    * - Retains values that are either:
    *   - Java Strings (`value.class.simpleName == 'String'`), or
    *   - Numeric values >= 0 (e.g., integers, floats).
    * - Nulls and negative numeric values are excluded.
    *
    * Input normalization:
    * - If `meta_query` is a single Map, it is wrapped into a singleton list.
    * - Assumes each metadata map uses abstract attribute keys (e.g., 'sample', 'plexid').
    * - Uses `translateKey(key)` to map abstract keys to iRODS AVU attribute names.
    *
    * Returns:
    * - A JSON-compatible Groovy string representing an AVU query array, suitable for embedding into
    *   a metaquery `target` block in iRODS/Baton-compatible workflows.
    *
    * Example output:
    * [
    *   {"a": "id_run", "v": "518005"},
    *   {"a": "tag_index", "v": ["1", "2", "3"], "o": "in"}
    * ]
    *
    * param meta_query A Map or List<Map> of metadata key-value pairs
    * return A JSON array string representing the AVU query clauses
    */
    def grouped = [:].withDefault { [] } //needs something to add to on line 63

    // Normalize input: if its a single map put into list
    def metas = (meta_query instanceof Map) ? [meta_query] : meta_query

    // Group values that pass the condition
    metas.each { m ->
        m.each { key, value ->
            if (value != null && (value.class.simpleName == 'String' || value >= 0)) {
                grouped[key] << value
            }
        }
    }

    // Build condensed AVU query list
    def query_list = []
    grouped.each { key, values ->
        def irods_key = translateKey(key)
        def unique_values = values.unique()

        if (unique_values.size() == 1) {
            def val = unique_values[0]
            query_list.add("""{"a": "${irods_key}", "v": "${val}"}""")
        } else {
            def val_array = unique_values.collect { "\"${it}\"" }.join(", ")
            query_list.add("""{"a": "${irods_key}", "v": [${val_array}], "o": "in"}""")
        }
    }

    return "[${query_list.join(", ")}]"
}


def objectOrCollection() {
    switch(params.read_type) {
        case 'ont':
            return '"collection": true'
        case 'all':
            return '"object": true, "collection": true'
        case 'illumina':
            return '"object": true'
        default:
            log.error("read_type '${params.read_type}' was not one of illumina|ont|all")
    }
}

process JSON_PREP {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    container 'quay.io/sangerpathogens/jq:1.6'

    input:
    val(meta)

    output:
    path(json_file), emit: path_channel

    script:
    json_file="input.json"
    """
    jq -n '{op: "metaquery", args: {${objectOrCollection()}, "avu": true}, target: {avus: ${avuIdQuery(meta)}}}' > ${json_file}
    """
}
