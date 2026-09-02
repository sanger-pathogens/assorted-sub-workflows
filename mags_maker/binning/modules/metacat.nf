process METACAT {
    tag "${meta.ID}"
    label 'cpu_8'
    label 'mem_4'
    label 'time_12'

    // No official biocontainer exists for MetaCAT as of writing (pip-only GitHub release,
    // no bioconda recipe) - installing the wheel at runtime into a generic python image.
    // Use the full (non-slim) image - Nextflow's own task wrapper needs `ps` for resource
    // tracing, which -slim doesn't ship; its absence silently aborts the whole task wrapper.
    // Replace `container` with a proper pinned image if/when one is published.
    container 'python:3.12'

    input:
    tuple val(meta), path(bam), path(bai), path(assembly)

    output:
    tuple val(meta), path(bins_out), emit: bins

    script:
    bins_out = "metacat_out"
    // pip's default --user install lands in the HOST's $HOME/.local (auto-bind-mounted into
    // the container by singularity), not anywhere on this container's PATH. Force it into an
    // ephemeral, work-dir-local prefix instead. That alone isn't enough though: Python's
    // site module silently disables user-site packages when running as UID 0 (root, the
    // default inside this container), so pip installs the files but nothing ever adds them
    // to sys.path - add the site-packages dir to PYTHONPATH explicitly to bypass that.
    """
    export PYTHONUSERBASE=\$PWD/.local_metacat
    pip install --quiet --user ${params.metacat_wheel_url}
    export PATH=\$PYTHONUSERBASE/bin:\$PATH
    export PYTHONPATH="\$(echo \$PYTHONUSERBASE/lib/python*/site-packages):\${PYTHONPATH:-}"

    mkdir -p ${bins_out}
    MetaCAT coverage -b ${bam} -o ${bins_out}/coverage -tc ${task.cpus}
    MetaCAT seed -f ${assembly} -o ${bins_out}/seed -t ${task.cpus}
    MetaCAT cluster -f ${assembly} -c ${bins_out}/coverage -s ${bins_out}/seed -o ${bins_out}/metacat -t ${task.cpus}
    """
}
