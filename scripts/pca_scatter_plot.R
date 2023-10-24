suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(cowplot))


source("./scripts/iwanthue.R")

set.seed(snakemake@params$parameters_porps$set_seed)

print("pca") 


#---------------------------------- Load data ---------------------------------- 
expr = fread(snakemake@input$expr, sep='\t', header=T)
annot = fread(snakemake@input$annot, sep='\t', colClasses=c(ID="character"))

properties = snakemake@params$properties
plots_props = snakemake@params$plots_props
snakemake_colors = snakemake@params$colors
results_folder = snakemake@params$results_folder


# ------------------------------- Parameters -----------------------------------
proper_names = c("Control_left ear" = "left ear of the untreated mice", "Control_right ear" = "right ear of the untreated mice", "Treatment_left ear" = "left ear of the treated mice", "Treatment_right ear" = "right ear of the treated mice")


# --------------------------------- Script -------------------------------------

for (props in properties) {
  output_folder = sprintf("%s/pca/%s_%s", results_folder, props[1], props[2])
  dir.create(output_folder, recursive=TRUE)
  
  # set output filename
  output_filename = sprintf("%s/pca_cropped", output_folder)
    
  # remove constants
  expr_filtered = expr[apply(expr[,-"miRNA"], 1, var) != 0,]
  
  comp = prcomp(t(expr_filtered[,-"miRNA"]), center = T, scale=T)
  imp = summary(comp)$importance
  comp_df = as.data.frame(comp$x)
  plot_df = data.frame(D1=comp_df$PC1, D2=comp_df$PC2, ID=colnames(expr_filtered[,-"miRNA"]), stringsAsFactors = F)
  
  prop_mapping = make.names(props)
  names(prop_mapping) = props
  
  if(is.vector(props)){
    for(a in props){
      plot_df[[prop_mapping[a]]] = as.factor(unlist(lapply(plot_df$ID, function(n) annot[annot$ID == n, a, with=F][1])))
    }
    major_prop = props[1]
  } else {
    plot_df[[prop_mapping[props]]] = as.factor(unlist(lapply(plot_df$ID, function(n) annot[annot$ID == n, props, with=F][1])))
    major_prop = props
  }
  
  if (props[2] == "Time") {
    legend_title = "Timepoint"
  } else {
    legend_title = props[2]
  }
  
  if(!is.null(snakemake_colors)){
    colors = unlist(snakemake_colors[[major_prop]])
    if(length(colors) == 1 && colors == "auto") {
      colors = unname(iwanthue(nlevels(plot_df[[prop_mapping[major_prop]]])))
    } else if(length(colors) == 1 && colors == "viridis_c"){
      # make continuous instead of discrete
      plot_df[[prop_mapping[major_prop]]] = unlist(lapply(plot_df$ID, function(n) annot[annot$ID == n, props, with=F][1]))
    } else if(length(colors) > 1 && length(colors) < length(unique(prop_mapping[major_prop]))) {
      stop(sprintf("Not enough colors provided for the property %s! Expected %d, got %d (%s)", major_prop, length(unique(prop_mapping[major_prop])), length(color), paste(names(colors), sep=', ')))
    }
  } else {
    colors = unname(iwanthue(nlevels(plot_df[[prop_mapping[major_prop]]])))
  }
  
  if(is.vector(props) && length(props) == 2){
    p = ggplot(plot_df, aes_string(x="D1",y="D2", color=prop_mapping[major_prop], shape=prop_mapping[props[2]], text="ID"))
  } else {
    p = ggplot(plot_df, aes_string(x="D1",y="D2", color=prop_mapping[major_prop], text="ID"))
  }
  
  p = p + geom_point(size = 2) + xlim(-15, 18) + ylim(-10, 15)
  
  if(length(colors) == 1 && colors == "viridis_c"){
    p = p + scale_colour_viridis_c(name=major_prop)
  } else if(length(colors) == 1 && colors == "viridis_d") {
    p = p + scale_colour_viridis_d(name=major_prop)
  } else {
    p = p + scale_color_manual(name=major_prop, labels = proper_names, values=colors) +
      scale_shape_manual(labels = proper_names)
  }
  
  if(is.vector(props) && length(props) == 2){
    p = p + scale_shape(name=props[2])
  }
  
  if(nrow(plot_df) < 25){
    p = p + geom_text_repel(aes(label = plot_df$ID), show.legend = FALSE)
  }
  
  if(!is.numeric(annot[[major_prop]])) {
    p = p + theme_classic() + theme(legend.position="bottom")
  }
  
  if(TRUE){
    p = p + xlab(sprintf("PC1 (%.2f%%)", imp["Proportion of Variance","PC1"]*100)) +
      ylab(sprintf("PC2 (%.2f%%)", imp["Proportion of Variance", "PC2"]*100))
  } else {
    p = p + theme(legend.position=ifelse(is.numeric(annot[[major_prop]]), "right", "bottom"), 
                  axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                  axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.background = element_blank())
  }
  p = p + theme(text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size), 
                axis.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size), 
                axis.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size), 
                plot.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size_header), 
                legend.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
                legend.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size, face = "bold"),
                legend.margin=margin(t = -8), legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.4, "cm"))

  needed_rows = 1
  if(is.vector(props) && length(props) == 2){
    total_len = sum(nchar(unique(as.character(plot_df[[prop_mapping[major_prop]]])))) + sum(nchar(unique(as.character(plot_df[[prop_mapping[props[2]]]]))))
    if(nlevels(plot_df[[prop_mapping[major_prop]]]) + nlevels(plot_df[[prop_mapping[props[2]]]]) > 4 || total_len > 50){
      needed_rows = ceiling(total_len / 50)
      p = p + theme(legend.box = "vertical")
      p = p + guides(col=FALSE, fill=guide_legend(nrow=needed_rows, title=legend_title), shape=guide_legend(nrow=needed_rows, title=legend_title))
    }
  } else {
    if(length(colors) != 1 || colors == "viridis_d") {
      total_len = sum(nchar(unique(as.character(plot_df[[prop_mapping[major_prop]]]))), na.rm=T)
      if(nlevels(plot_df[[prop_mapping[major_prop]]]) > 4 || total_len > 50){
        needed_rows = ceiling(total_len / 50)
        p = p + guides(col=FALSE, fill=guide_legend(nrow=needed_rows, title=legend_title), shape=guide_legend(nrow=needed_rows, title=legend_title))
      }
    }
  }
  ggsave(file = sprintf("%s.svg",output_filename), plot = p, width = plots_props$image_width, height = plots_props$image_height, unit = plots_props$image_units, dpi = plots_props$dpi)
  ggsave(file = sprintf("%s.png",output_filename), plot = p, width = plots_props$image_width, height = plots_props$image_height, unit = plots_props$image_units, dpi = plots_props$dpi)
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])



