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
5. Writes the resulting predicted value to `label.txt`.

An example json:
```
{
    "ID":"1761247.0",
    "Richness":371,
    "Shannon":3.64520578535183,
    "Inverse_Simpson":10.2939546318467,
    "total_reads":75243196,
    "q30_rate":0.927577
}
```

## Inputs

None — the process takes no channel inputs. The path to the input JSON
(`-d /lustre/scratch124/.../DRR671405.json`) is currently hard-coded in the
script block. This will change as the project is developed


## Outputs

| Output | Emit name | Description |
|---|---|---|
| `label.txt` | `label_file` | A text file containing the predicted memory value (a plain number, no units) |

## How it fits into the subworkflow

```groovy
if (params.metaspades) {
    GET_LABEL()

    label_ch = GET_LABEL.out.label_file
        .map { file -> file.text.trim() + " MB" } // the "+ MB" part may be unnecessary if the output of the GET_LABEL process also outputs the units. This can be changed by changing the code in the pipeline resource optimisation codebase itseld

    METASPADES(reads_ch, label_ch)
}
```

This block only runs when `params.metaspades` is enabled, and works as follows:

1. **`GET_LABEL()`** is called with no arguments, producing `label.txt`
   containing the predicted memory value.
2. **`label_ch`** is built by mapping over `GET_LABEL.out.label_file`: the file's
   contents are read (`file.text`), trimmed of whitespace/newlines
   (`.trim()`), and have `" MB"` appended. This turns the raw numeric
   prediction (e.g. `4500`) into a Nextflow-friendly memory string
   (e.g. `4500 MB`). This will likely be updated an made unnecessary.
3. **`METASPADES(reads_ch, label_ch)`** runs the metaSPAdes assembly process,
   passing in the reads channel (`reads_ch`) as normal, and the predicted
   memory string (`label_ch`) so the process can use it — typically to set its
   own dynamic memory directive (e.g. `memory { label }`) rather than a fixed
   value.

In short: **`GET_LABEL` predicts how much memory metaSPAdes will need for this
sample, and that prediction is fed straight into the `METASPADES` process's
resource allocation, instead of using a static memory setting.**

**METASPADES PROCESS**
in assorted-sub-workflows

```groovy
process METASPADES {
    tag "${meta.ID}"
    label 'cpu_8'
    memory "${MEMORY_LABEL}"
    label 'time_12'

    container 'quay.io/biocontainers/spades:3.15.5--h95f258a_1'

    input:
    tuple val(meta), path(first_read), path(second_read)
    val(MEMORY_LABEL)
```