process COMEBIN {
    tag "${meta.ID}"
    label 'gpu'
    label 'cpu_5'
    label 'mem_8'
    label 'time_12'

    // No `container` here, deliberately: the published biocontainers image
    // (quay.io/biocontainers/comebin:1.1.0--hdfd78af_2) ships a CPU-only PyTorch build
    // (torch.version.cuda == None, confirmed directly) - not an env/scheduling issue, a
    // real limitation of that image, so no amount of GPU config fixes it. Use the
    // self-installed, GPU-enabled env from the ensemble-binning benchmark instead
    // (see that skill for how it was built) until/unless a GPU-capable image is published.
    beforeScript "export PATH=/data/pam/team162/shared_scratch/software_dbs/envs/comebin_1_1_0/bin:\$PATH"

    input:
    tuple val(meta), path(bam), path(bai), path(assembly)

    output:
    tuple val(meta), path(bins_out), emit: bins

    script:
    bins_out = "comebin/comebin_res/comebin_res_bins"
    """
    mkdir -p bamfiles
    ln -s "\$PWD/${bam}" "bamfiles/${bam}"
    ln -s "\$PWD/${bai}" "bamfiles/${bai}"

    # LSF's -gpu resource string (see the `gpu` label in binning.config) allocates a GPU for
    # scheduling but does not itself export CUDA_VISIBLE_DEVICES - set it explicitly, same as
    # the validated ensemble-binning skill recipe.
    export CUDA_VISIBLE_DEVICES=0
    run_comebin.sh -a ${assembly} -p bamfiles -o comebin -t ${task.cpus} -b 256
    """
}
