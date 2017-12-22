# statistical_tools
Universal tools for group statistical analysis

## groupdata2statistics
This tool counts statistics from linker file containing two tab-delimited columns without a header: file path and group name.
```
python3 groupdata2statistics.py -g groupdata_processed.tsv -i reference_id -v id_mapped_reads_per_million_sample_total_reads -a 0.05 -o some_dir
```


## suffixes2statistics
Given two tab-delimited tables: the first table contains aliases, path prefixes and suffixes, the second contains abstract sample names and corresponding group IDs. So the ending name of **existing** file would be composed as `<path prefix><sample name><path suffix>`

```
python3 suffixes2statistics.py -s suffixes.tsv -g groupdata_raw.tsv -f bad_words.tsv -t 10 -c control_group_1 control_group_2 -m u-test -r fdr_bh -i reference_id -v id_mapped_reads_per_million_sample_total_reads -a 0.05 -o some_dir
```
