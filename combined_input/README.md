# Combined Input subworkflow

## Scripts

### filter_metadata.py

This script facilitates filtering an input TSV file using a number of per-column filters (while ensuring that columns are correctly interpreted with a given type).

The column, filter and datatype are specified using the path to an input TSV manifest using the `-f` option. Valid filters are any string that can be provided to [`pandas.DataFrame.query()`](https://pandas.pydata.org/pandas-docs/version/2.2/reference/api/pandas.DataFrame.query.html) and datatypes are `int`, `float`, `datetime`, `bool` and `str` (use other types at your peril). Example manifest:

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