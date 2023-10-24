library(data.table)

tbl = fread(snakemake@input[[1]], sep='\t', header=T)

sample_cols = colnames(tbl)[2:ncol(tbl)]
tbl[, (sample_cols):=lapply(.SD, function(x) as.numeric(x >= snakemake@params$min_detection)), .SDcols = sample_cols]

fwrite(tbl, "data/expression/raw/raw_detection_miRNA.csv", sep='\t', na="NA")


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
