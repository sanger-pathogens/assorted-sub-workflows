# Taxo Profile subworkflow
## Subworkflow overview
Includes modules for sketching paired end reads with Sylph and querying/profiling those sketches against a pre-sketched Sylph database. Subworkflow produces reports including in Bracken style.

See the schema.json for parameters and help messages.
 
 ## Input
 - Paired end reads (meta, read_1, read_2)

 ## Output
- Sylph report
- Sylph-tax report
- Kraken2 Bracken style reports
- Metaphlan Abundance reports

<!---
## Subworkflow-specific notes
(none)
--->

## Associated pipelines and subworkflows
<!---
This subworkflow uses the following subworkflows:
(none)
--->

This subworkflow is used in the following pipelines/subworkflows:
- [QC subworkflow](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/pipelines/assorted-sub-workflows/-/tree/main/qc)
- [sylph_refset subworkflow](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/pipelines/assorted-sub-workflows/-/tree/main/sylph_refset)
- [QC Short Read pipeline](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/pipelines/qc-short-read)
- [Irods Extractor pipeline](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/pipelines/irods_extractor)

Please ensure any changes are carried through to those, including testing and documentation updates.
 
## Software versions
| Software      | Version | Image URL                                            |
| ------------- | ------- | ---------------------------------------------------- |
| sylph         | 0.8.1   | quay.io/biocontainers/sylph:0.8.1--ha6fb395_0                    |
| sylph-tax         | 1.7.0   | quay.io/biocontainers/sylph-tax:1.7.0--pyhdfd78af_0                    |

## Dependencies
- Kraken2 database
- Sylph database
- Sylph-tax metadata TSV