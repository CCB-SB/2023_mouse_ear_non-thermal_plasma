suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(lme4))


snakemake@source("./custom_pvca.R")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_porps$set_seed)

print("pvca")

results_folder = snakemake@params$results_folder
output_folder = sprintf("%s/pvca", results_folder)
dir.create(output_folder, recursive=TRUE)

annot = fread(snakemake@input$annot, sep='\t', colClasses=c(ID="character"))
props = snakemake@params$props

prop_mapping = make.names(props)
names(prop_mapping) = props

tbl = fread(snakemake@input$expr_log10, sep='\t', header=T)

# keep only expression
expr = tbl[, sapply(tbl, is.numeric), with=F]

# remove constants
expr = expr[apply(expr, 1, var) != 0,]

# check if this is logarithmized data
if((is.null(snakemake@params$force) || !snakemake@params.force) && max(expr) >= 100) {
  stop("It seems that your data is not logarithmized. PVCA expects logarithmized data! You can force execution by setting params.force=True")
}

# Make factors
annot <- annot[,names(prop_mapping) := lapply(.SD, as.factor), .SDcols = names(prop_mapping)]

# Prepare for PVCA
pheno_dat <- data.frame(annot)
rownames(pheno_dat) <- pheno_dat$ID
pheno_dat <- pheno_dat[match(colnames(expr), rownames(pheno_dat)),]

# Run PVCA
pvcaObj <- pvcaBatchAssess(expr, pheno_dat, batch.factors=unname(prop_mapping), threshold=snakemake@params$parameters$min_var, threads=1, skip.unique=T)

# Prepare visualization
pvcaObj.df <- data.frame(variance=round(pvcaObj$dat, 3)*100, variable=gsub(".", " ", sub("resid", "Residual", pvcaObj$label, fixed=T), fixed=T), variance_label=paste0(as.character(round(pvcaObj$dat , 3)*100), "%"))
if(!snakemake@params$parameters$keep_zeros){
    pvcaObj.df <- pvcaObj.df[pvcaObj.df$variance != 0,]
}
other_var <- sum(pvcaObj.df[pvcaObj.df$variance < 100*snakemake@params$parameters$min_cutoff,]$variance)
if (snakemake@params$parameters$min_cutoff != 0) {
  pvcaObj.df <- rbind(pvcaObj.df, data.frame(variance=other_var, variable="Other", variance_label=sub(",", ".", paste0(as.character(other_var), "%"))), stringsAsFactors=FALSE)
}
pvcaObj.df <- pvcaObj.df[pvcaObj.df$variance >= 100*snakemake@params$parameters$min_cutoff,]
pvcaObj.df <- pvcaObj.df[order(-pvcaObj.df$variance),]
pvcaObj.df$variance <- pmax(pvcaObj.df$variance, 1)
pvcaObj.df$variable <- factor(pvcaObj.df$variable, levels=rev(pvcaObj.df$variable))

# Create plot
pvca_plot <- ggplot(pvcaObj.df, aes(x = variance, y = variable)) + 
  geom_bar(stat = "identity", width = 0.75, fill = "#D5A021") +
  scale_x_log10(limits = c(1, 900), expand = c(0, 0)) + 
  geom_text(aes(label=variance_label), nudge_x=log10(3.6), size=snakemake@params$plots_props$font_size/11.04*3.88, angle = 0, nudge_y=0) + 
  labs(x="Observed variance (%)", y="") +
  theme_classic() +
  theme(aspect.ratio=1, 
        axis.text.x = element_text(size=snakemake@params$plots_props$font_size), 
        legend.position="none", 
        text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
        axis.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
        axis.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
        plot.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size_header),
        legend.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
        legend.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size)
  )

# Save results to storage
ggsave(sprintf("%s/pvca.png", output_folder), pvca_plot, dpi=snakemake@params$plots_props$dpi, width=snakemake@params$plots_props$image_width, height=snakemake@params$plots_props$image_height, units = snakemake@params$plots_props$image_units)
ggsave(sprintf("%s/pvca.svg", output_folder), pvca_plot, width=snakemake@params$plots_props$image_width, height=snakemake@params$plots_props$image_height, units = snakemake@params$plots_props$image_units, limitsize = FALSE, scale = 1)

pvca_plot_ply = ggplotly(pvca_plot) %>% toWebGL()

saveRDS(pvca_plot_ply, file=sprintf("%s/pvca.rds", output_folder))


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
