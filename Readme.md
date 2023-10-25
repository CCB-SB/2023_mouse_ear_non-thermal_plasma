# Skin treatment with non-thermal plasma modulates the immune system through miR-223-3p and its target genes

This repository contains the code necessary to produce the publication figures via a snakemake workflow.

## Necessary files
Provide the following files in the appropriate directories.
| name  | source | target path |
| ------------- | ------------- | ------------- |
| raw counts  | GEO Series GSE236788 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236788)  | `data/expression/raw/GSE236788_miRNA_expression_raw.tsv` |
| KEGG  | miRPathDB v2.0 (https://mpd.bioinf.uni-sb.de/download/version_2/miRPathDB2_mmu_p_values.tar.gz)  | `data/external_data/mirpathd/KEGG_prediction_union.csv` |
| human cellular microRNAome database  | Supplementary tables (https://doi.org/10.1093/gigascience/giac083)  | `data/external_data/giac083_supplemental_files/Supplementary_Table_S2.xlsx` `data/external_data/giac083_supplemental_files/Supplementary_Table_S7.xlsx`|

## Running the pipeline
```bash
# run snakemake in this folder
snakemake
# use conda a backend if mamba is not available
snakemake --use-conda --conda-frontend conda
```

## Results
If the pipeline suceeded, results are availabe in the `publication` folder. The names are according to the figures in the publication whereas `fig_2b_1` corresponds to the top left plot of Fig. 2b.

## License
[MIT © CCB-SB](/LICENSE)