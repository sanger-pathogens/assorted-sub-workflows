# Sylph Refset subworkflow
## Subworkflow overview
Profiles a set of reads and creates a reference set of all genomes in the supplied database matching the detected taxa. Gets references at a defined taxonomic rank.

See the schema.json for parameters and help messages.
 
 ## Input
 - Paired end reads (meta, read_1, read_2)

 ## Output
- Collected references (meta, reference)
- Combined sylph report

## Subworkflow-specific notes
### Pooling latin taxa
Certain genus/species in GTDB are further divided by appended alphabet suffixes, for example Escherichia coli in GTDB has 3 species rank taxonomic groups; Escherichia_coli, Escherichia_coli_E and Escherichia_coli_F. Further explanation is available in the [GTDB documentation](https://gtdb.ecogenomic.org/faq#why-do-some-genus-and-species-names-end-with-an-alphabetic-suffix).

If you wanted to consider these as one group you can use the advanced option to pool latin taxa. Note that:

a) generated groups are no longer compliant with GTDB taxonomic definitions, consider if this affects downstream

b) the size of the produced group may be considerably larger, for example in GTDB release 232 at the genus level g__Clostridium has 1607 genomes but all 34 GTDB genuses in g__Clostridium* total at 2931 genomes.

Note that not all taxa belonging to a "traditional" species might be pooled this way due to certain GTDB species being named differently; for instance in GTDB r232, a new species called `ECMA0423 sp047199055` has been created out of genomes previously classified as `Escherichia_coli`.

## Associated pipelines and subworkflows
This subworkflow uses the following subworkflows:

- [taxo_profile](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/pipelines/assorted-sub-workflows/-/tree/main/taxo_profile)

This subworkflow is used in the following pipelines/subworkflows:

- [gemsweep](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/pipelines/gemsweep)

Please ensure any changes are carried through to those, including testing and documentation updates.
 
## Software versions
| Software      | Version | Image URL                                            |
| ------------- | ------- | ---------------------------------------------------- |
| pandas         | 2.2.1   | quay.io/sangerpathogens/pandas:2.2.1                    |
 
## Dependencies
- GTDB database
- Requires a file mapping genome IDs to their path