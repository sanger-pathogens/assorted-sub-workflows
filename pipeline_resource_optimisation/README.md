# GET_LABEL Process

## Overview

`GET_LABEL` is a Nextflow process that uses a pre-trained machine learning model to
predict a resource label (a memory value, in MB) for another downstream process —
in this instance `METASPADES`. Rather than hard-coding a memory allocation for
metaSPAdes, this process calls out to a `modelling` tool that predicts an
appropriate memory requirement based on characteristics of the sample being
assembled. The process itself is lightweight and only runs the prediction step.

## What the process does

1. Runs the `modelling` CLI tool inside the
   `quay.io/sangerpathogens/pipeline_resource_optimisation:test_6` container.
2. Calls it in `get_label` mode, using a specific pre-trained model
   (`METASPADES_GradientBoosting`, tagged `0.0.0`) to make a prediction
   (`--predict`).
3. Points the tool at a JSON file (`-d`) describing the sample/dataset the
   prediction should be based on.
4. Applies an adjustment to the raw prediction: `--adjust-value mae` adjusts the
   prediction using the model's mean absolute error, scaled by a factor of
   `1.5` (`--adjust-scale 1.5`). This padding helps avoid under-predicting the
   memory required.
5. Leaves the resulting predicted value on stdout.

The JSON is written by the process itself, from the feature values handed to it
on the input channel. An example:
```
{
    "Richness":371,
    "Shannon":3.64520578535183,
    "Inverse_Simpson":10.2939546318467,
    "total_reads":75243196,
    "q30_rate":0.927577
}
```

## Inputs

| Input | Description |
|---|---|
| `tuple val(meta), val(features_json)` | `meta` is the sample's meta map; `features_json` is a single JSON object holding one row of the model x table (every column except `sample_ID`) |

One task runs per sample, and each task writes its own `features.json` into its
work directory before calling `modelling`.

## Outputs

| Output | Emit name | Description |
|---|---|---|
| `tuple val(meta), stdout` | `label` | The predicted memory value as a string (a plain number, no units) |

Because the prediction is captured from stdout, anything else `modelling` has to
say must go to stderr, or it ends up in the label.

## How it fits into the subworkflow

```groovy
if (params.metaspades) {

    // one JSON object per row of the model x table, keyed by sample
    channel
    .fromPath ( params.model_x_table, checkIfExists: true )
    .splitCsv ( header:true, sep:',', strip:true )
    .map { row ->
        def features = row.findAll { column, value -> column != 'sample_ID' }
            .collectEntries { column, value -> [column, value?.isNumber() ? value as BigDecimal : value] }
        [ row.sample_ID, JsonOutput.toJson(features) ]
    }
    .set { model_x_ch }

    reads_ch
    .map { meta, read_1, read_2 -> [ meta.ID, meta, read_1, read_2 ] }
    .join ( model_x_ch, failOnMismatch: true, failOnDuplicate: true )
    .set { reads_with_features }

    GET_LABEL( reads_with_features.map { id, meta, read_1, read_2, features -> [ meta, features ] } )

    GET_LABEL.out.label
    .map { meta, label -> [ meta.ID, meta + [ memory_label: "${label.trim()} MB" ] ] }
    .join ( reads_with_features.map { id, meta, read_1, read_2, features -> [ id, read_1, read_2 ] } )
    .map { id, meta, read_1, read_2 -> [ meta, read_1, read_2 ] }
    .set { reads_with_label }

    METASPADES(reads_with_label)
}
```

This block only runs when `params.metaspades` is enabled, and works as follows:

1. **`model_x_ch`** turns each row of `params.model_x_table` into
   `[sample_ID, features_json]`. The columns are never named in the code, so
   retraining the model on different features needs no pipeline change — only
   `sample_ID` is special, and it is dropped from the JSON because it is an
   identifier rather than a feature.
2. **`reads_with_features`** joins those features onto the reads by sample ID.
   `failOnMismatch` means a sample missing from the table stops the run rather
   than being silently dropped from the assembly.
3. **`GET_LABEL`** runs once per sample and emits the prediction as a value.
4. **`reads_with_label`** merges that prediction into the sample's meta map as
   `meta.memory_label` (e.g. `4500 MB`) and re-attaches the reads. The `" MB"`
   suffix becomes unnecessary if `modelling` starts emitting units itself.
5. **`METASPADES(reads_with_label)`** reads the prediction straight off the meta
   map in its own dynamic memory directive, so no second channel has to be kept
   in step with the reads.

In short: **`GET_LABEL` predicts how much memory metaSPAdes will need for this
sample, and that prediction is fed straight into the `METASPADES` process's
resource allocation, instead of using a static memory setting.**

**METASPADES PROCESS**
in assorted-sub-workflows

```groovy
process METASPADES {
    tag "${meta.ID}"
    label 'cpu_8'
    memory { meta.memory_label }
    label 'time_12'

    container 'quay.io/biocontainers/spades:3.15.5--h95f258a_1'

    input:
    tuple val(meta), path(first_read), path(second_read)
```