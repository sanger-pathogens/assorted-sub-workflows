# Combined Input subworkflow

## Mixed_inputs.nf

This workflow determines which processes should be called for the incoming parameters.

`validate_parameters` takes no input but uses the params map from the Nextflow runtime environment and checks for either IRODS or READ type parameters, as defined below:

There are two ways of specifying input for the pipeline:

- Manifests (batch input files):
- --manifest_ena: ENA
- --manifest_of_lanes: IRODS
- --manifest_of_reads or --manifest: READS

CLI arguments:
- --studyid: IRODS
- --runid: IRODS
- --laneid: IRODS
- --plexid: IRODS

`validate_parameters` returns a list of active workflows based on the runtime parameters

This list is then used to trigger the correct workflows for the data specified a input to the pipeline: either IRODS (merging of CLI + Manifest input) or parsing the read manifest.

The subworkflows can be executed separately through direct reference.  
For example, to specify the use of data from IRODS, one can call either:
```
COMBINE_IRODS
```
or directly
```
IRODS_EXTRACTOR
```

Or, in the scenario of using only reads files available on disk, one can call directly:
```
INPUT_CHECK
```

When data needs to be downloaded from the ENA:
```
ENA_DOWNLOAD
```

Mixed input should be used in cases that you expect manifests of reads, IRODS extraction, and/or download from ENA. 

## Scripts

### filter_metadata.py

This script facilitates filtering an input TSV file using a number of per-column filters (while ensuring that columns are correctly interpreted with a given type).

The column, filter and datatype are specified using the path to an input TSV manifest using the `-f` option:

- **column**: name of the column the filter applies to\; there should only be one row per column
- **filter**: filter to be applied to the column, which is a string that could be passed to [`pandas.DataFrame.query()`](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html)
- **datatype**: expected datatype for values in the column, either `int`, `float`, `datetime`, `bool` or `str` (use other types at your peril i.e. anything you can use as dtype in [`pandas.DataFrame.astype()`](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.astype.html)) 

Example manifest:

```
column	filter	datatype
center_name	center_name.str.contains("Wellcome Sanger Institute", na=False)	str
read_count	"read_count > 2500000"	int
collection_date	"2012 < collection_date"	datetime
```

By default, it will automatically remove any rows where a column value cannot be converted to the given type. For stricter validation, errors can be raised upon an invaid value (`-error_on_invalid_type` option). Missing values will be inferred, but additional values to be considered missing can be supplied with `--missing_values [MISSING_VALUES ...]`.

For the output TSV, you can select which columns you would like to output with the `--select` option, e.g. `--select sample_accession collection_date`. The header can optionally be removed from output using `--remove_header`.

Example command:
```
filter_metadata.py -f filter_manifest.tsv -i metadata.tsv -m missing not_provided not_available -s sample_accession -o test_output.tsv -r
```

You can find more help on any options using the `-h` option:
```
usage: filter_metadata.py [-h] --input INPUT --filter_manifest FILTER_MANIFEST
                          [--select SELECT [SELECT ...]]
                          [--missing_values MISSING_VALUES [MISSING_VALUES ...]]
                          [--remove_header] [--error_on_invalid_type]
                          [--drop_duplicates_from DROP_DUPLICATES_FROM [DROP_DUPLICATES_FROM ...]]
                          [--output OUTPUT] [--logfile LOGFILE]

Filter rows of a TSV file based on conditions.

options:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        Path to the input TSV file. (default: None)
  --filter_manifest FILTER_MANIFEST, -f FILTER_MANIFEST
                        Path to a TSV manifest specifying the columns in the input
                        TSV to which filter and datatype conversion should be
                        applied. The manifest should contain 3 columns: column,
                        filter, datatype. Entries in the 'column' column should
                        match a column in the input TSV and be unique. The filter
                        should match a string that could be supplied to
                        `pd.DataFrame.query()`, e.g. 'age > 30'. The datatype can be
                        int, float, datetime, bool and str. (default: None)
  --select SELECT [SELECT ...], -s SELECT [SELECT ...]
                        Specify columns to select in the output DataFrame. By
                        default, all columns will be selected. (default: None)
  --missing_values MISSING_VALUES [MISSING_VALUES ...], -m MISSING_VALUES [MISSING_VALUES ...]
                        Specify values that should be interpreted as NA (or NaT) by
                        pandas - see https://pandas.pydata.org/pandas-
                        docs/version/2.2/reference/missing_value.html). These values
                        are added to the list of values usually considered as NA
                        (see `na_values` argument of pandas.read_csv -
                        https://pandas.pydata.org/pandas-
                        docs/stable/reference/api/pandas.read_csv.html). (default:
                        None)
  --remove_header, -r   Remove the header from the output TSV. (default: False)
  --error_on_invalid_type, -e
                        During type conversion, upon encountering a value in the
                        column that cannot be converted to the given datatype, raise
                        an error. Default behaviour is to remove the row that
                        contains an invalid value. (default: False)
  --drop_duplicates_from DROP_DUPLICATES_FROM [DROP_DUPLICATES_FROM ...], -d DROP_DUPLICATES_FROM [DROP_DUPLICATES_FROM ...]
                        Specify columns that should not have duplicate values. Only
                        the first row with the value will be kept. (default: None)
  --output OUTPUT, -o OUTPUT
                        Path to the output file to save the filtered TSV. If not
                        specified, the metadata DataFrame will be printed to stdout.
                        (default: None)
  --logfile LOGFILE, -l LOGFILE
                        Path to the log file, defaults to
                        filter_metadata-<timestamp>.log (default:
                        filter_metadata-2024-28-31_12-28-32.log)
```