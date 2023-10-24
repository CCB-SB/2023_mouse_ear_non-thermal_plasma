suppressPackageStartupMessages(library(data.table))

tbl = fread(snakemake@input[[1]], sep='\t', header=T)
tbl[, (2:ncol(tbl)):=(log2(tbl[,2:ncol(tbl)] + 1))]

fwrite(tbl, "data/expression/norm/miRNA_quantification_rpmmm_norm_log2.csv", sep='\t')


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
