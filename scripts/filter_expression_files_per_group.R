suppressPackageStartupMessages(library(data.table))

expr = fread(snakemake@input$expr, sep='\t', header=T)
expr_log2 = fread(snakemake@input$expr_log2, sep='\t', header=T)
expr_log10 = fread(snakemake@input$expr_log10, sep='\t', header=T)

detect = fread(snakemake@input$detect, sep='\t', header=T)

detection_rate = snakemake@params$detection_rate

feature_col <- "miRNA"

if(colnames(detect)[1] == "V1"){
  setnames(detect, "V1", feature_col)
}

annot = fread(snakemake@input$annot, sep='\t', colClasses=c(ID="character"))
grouping <- snakemake@params$filter_variable

to_keep = rep(F, nrow(detect))
for(g in unique(annot[,get(grouping)])){
  if(g != ""){
    sub = annot[get(grouping) == g, ID]
    to_keep = to_keep | rowMeans(detect[,..sub]) >= snakemake@params$detection_rate
  }
}

stopifnot(assertthat::are_equal(expr[,1], detect[,1]))
stopifnot(assertthat::are_equal(expr_log2[,1], detect[,1]))
stopifnot(assertthat::are_equal(expr_log10[,1], detect[,1]))

# fix filename for expr
fwrite(expr[to_keep], sprintf("data/expression/norm/miRNA_filtered_quantification_rpmmm_norm.detection_rate_%sp_per_group.csv", (detection_rate * 100)), sep='\t')
fwrite(expr_log2[to_keep], sprintf("data/expression/norm/miRNA_filtered_quantification_rpmmm_norm_log2.detection_rate_%sp_per_group.csv", (detection_rate * 100)), sep='\t')
fwrite(expr_log10[to_keep], sprintf("data/expression/norm/miRNA_filtered_quantification_rpmmm_norm_log10.detection_rate_%sp_per_group.csv", (detection_rate * 100)), sep='\t')


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
