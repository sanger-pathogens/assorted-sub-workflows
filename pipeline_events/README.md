# Pipeline Events Subworkflow

This workflow interacts with the Pipeline Events Database API to record the opening and closing of "batches", entities intended to embody pipeline runs, or similar units of analysis. Each batch has a number of parameters/metadata associated with it, which are logged in the form of a JSON file and can be tracked. Additionally, files created/published in the run can be associated with the batch and tracked.

[[_TOC_]]

## Workflows

### pipeline_events.nf

This file contains two workflows:

#### `PIPELINE_EVENTS_INIT`

- **Purpose**: Open a batch and generate a JSON file to record parameters supplied to nextflow

- **Inputs**: None

- **Outputs**:
    - `batch_id` (value): ID of the newly opened batch (here an automatically generated unique id (UUID)
    - `batch_manifest_info` (path|empty): Path if a it is intended that the batch metadata info JSON file is to be associated with the batch (currently disabled, enable with param `associate_batch_metadata`). Empty channel otherwise.

#### `PIPELINE_EVENTS_END`

- **Purpose**: Close the given batch and log the number of files that should now be tracked by the database.

- **Inputs**:
    - `batch_id` (value): ID of the opened batch to be closed
    - `batch_manifest_info` (path|empty): Path of the batch metadata info JSON file, or empty channel if not being tracked.
    - `created_file_infos` ([value, path, value]): Tuple of runid, path of the published file, and a filetype string that have been in the pipeline events database.

- **Outputs**:
    - `all_created_file_paths` (path): Filepath component (second element) of `all_created_file_infos`.
    - `all_created_file_infos` ([value, path, value]): Similar to `created_file_infos` except now also includes a tuple for the batch metadata info JSON file (if tracked).
    - `created_file_count` (value): Number of files that have been created/tracked in the pipeline events database.
    - `created_file_count_per_type` ([value, value]): Tuple of filetype and the number of files created/tracked for that specific type.

## Modules

### pipeline_events.nf

#### `PIPELINE_GET_METHOD`

- **Purpose**:  
  Retrieve pipeline metadata and construct consistent identifiers (name, URL, short name) for the running pipeline.  
  This process captures key manifest and parameter details for downstream tracking and reporting.

- **Inputs**:  
  - None

- **Outputs**:  
  - `method_url` (value): URL to the pipeline’s repository or version-specific branch/tag.  
  - `method_name` (value): Full pipeline name combined with version string.  
  - `method_short` (value): Short identifier for the pipeline, typically derived from its repository name.  
  - `pipeline_manifest_params` (value): Map containing the pipeline’s manifest, parameters, and environment information.

---

#### `PIPELINE_EVENTS_OPEN_BATCH`

- **Purpose**:  
  Open a new batch for the pipeline run, generate a JSON manifest capturing runtime parameters, and send an “open” event to the tracking system.

- **Inputs**:  
  - `method_url` (value): URL identifying the pipeline source.  
  - `methodname` (value): Full pipeline name with version.  
  - `methodshort` (value): Short name identifier for the pipeline.  
  - `batch_mani_params` (value): Map containing pipeline and parameter metadata to associate with the batch.

- **Outputs**:  
  - `batch_manifest_params` ([value, path]): Tuple containing the updated manifest map and the path to the generated batch JSON file.  
  - `batch_id` (value): ID assigned to the newly opened batch.

---

#### `PIPELINE_EVENTS_CLOSE_BATCH`

- **Purpose**:  
  Close a batch that was previously opened and record the number of files created within it. Sends a “close” event to the tracking system.

- **Inputs**:  
  - `batch_id` (value): ID of the batch to close.  
  - `filecreatedcount` (value): Total number of files created/tracked for this batch.

- **Outputs**:  
  - None (side effect: sends batch closure event to tracking system)

---

#### `GATHER_RESULTFILE_INFO`

- **Purpose**:  
  Compute and provide the published directory paths of result files to be tracked.

- **Inputs**:  
  - `[meta, resultfileWorkPath]` ([value, path]): Metadata map and the working path of the result file.  
  - `outputfoldertag` (value): Name of the output subdirectory where files are stored.  
  - `file_type` (value): File type string used for tracking purposes.  
  - `batch_id` (value): ID of the active batch associated with these result files.

- **Outputs**:  
  - `file_info` ([value, path, value, value, value]): Tuple containing metadata map, result file path, published directory absolute path, file type, and batch UUID.

---

#### `PIPELINE_EVENTS_CREATE_FILE`

- **Purpose**:  
  Record the creation of a result file in the tracking system, associating it with the correct batch and optionally a run ID.  
  Generates a `file` pipeline event containing metadata such as MD5 checksum, file type, and user/group ownership.

- **Inputs**:  
  - `[meta, resultfileWorkPath, resultfilePublishedDir, file_type, batch_id]` ([value, path, value, value, value]):  
    Tuple including metadata, working file path, published directory path (as a value), file type, and associated batch ID.

- **Outputs**:  
  - `created_file_info` ([value, value, value]): Tuple containing the output/run ID, the full published file path, and the file type.

## Parameters

For a description of available parameters and their uses, see [schema.json](./schema.json) and [pipeline_events.config](./pipeline_events.config).

## Pipelines/Workflows using this subworkflow
- [irods_extractor](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/pipelines/irods_extractor)