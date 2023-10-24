suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))

set.seed(snakemake@params$parameters_porps$set_seed)

print("violin")


#---------------------------------- Functions ---------------------------------- 
median_fun <- function(x){
  return(data.frame(y = median(x), label = "median"))
}


#---------------------------------- Load data ---------------------------------- 
expr = fread(snakemake@input$expr, sep='\t', header=T)
annot = fread(snakemake@input$annot, sep='\t', colClasses=c(ID="character"))
interesting_miRNA = fread(snakemake@input$miRNA_list)

comps = snakemake@params$comparisons

results_folder = snakemake@params$results_folder


#------------------------------------ Script ----------------------------------- 
for (miR in interesting_miRNA$miRNA) {

  # extract the corresponding miRNA row from expression matrix
  expr_filtered = expr[expr$miRNA == miR,]

  for(i in 1:length(comps)) {
    
    # which comp do we investigate and which groups are considered
    prop = names(comps)[i]
    c = comps[[i]]
    
    # generate plot_title and output folder which the specific information of comp and possible subset
    if (!is.null(c$subset)) {
      plot_title = sprintf("%s (%s=%s)", miR, c$subset[1], c$subset[2])
      output_folder = sprintf("%s/violin_plots/expr/%s_%s=%s", results_folder, prop, c$subset[1], c$subset[2])
    } else {
      plot_title = sprintf("%s", miR)
      output_folder = sprintf("%s/violin_plots/expr/%s", results_folder, prop)
    }

    plot_df = c()
    for (group_comp in c$g) {
      
      # select the IDs from the annot which matches with the one given in this comp and possible subset
      if (!is.null(c$subset)) {
        group_comp_IDs = annot[(annot[[prop]] == group_comp & annot[[c$subset[1]]] == c$subset[2]),]$ID
      } else {
        group_comp_IDs = annot[annot[[prop]] == group_comp,]$ID
      }

      # filter these IDs (samples) from the expression row and transpose it for plotting purposes
      tmp_group_comp = colnames(expr_filtered) %in% group_comp_IDs
      expr_filtered_group_comp = t(expr_filtered[, ..tmp_group_comp])
      
      # generate a vector as long as the filtered expression values containing only the considered group of the comp
      if(prop == "Status_Time") {
        tmp = strsplit(group_comp, "_")[[1]][1]
        expr_filtered_group_comp_prop = data.frame(rep(tmp, length(group_comp_IDs)))
      } else if (prop == "Tissue_Time"){
        tmp = strsplit(group_comp, "_")[[1]][2]
        expr_filtered_group_comp_prop = data.frame(rep(tmp, length(group_comp_IDs)))
      } else {
        expr_filtered_group_comp_prop = data.frame(rep(group_comp, length(group_comp_IDs)))
      }
 
      # concatenate the expression and the group names vector 
      tmp = data.table(expr_filtered_group_comp, expr_filtered_group_comp_prop)
      
      # hang the resulting tables among each other the create the plotting table
      plot_df = data.table(rbind(plot_df, tmp))
    }
    # set the colnames of the plotting table
    colnames(plot_df) = c("expr","comp")
    
    # to rename the x ticks labels in a proper way 
    if (prop == "Status") {
      plot_x_ticks = c("Treatment" = "Treated mice", "Control" = "Untreated mice")
    } else if (prop == "Group") {
      plot_x_ticks = c("Treatment_left ear" = "TL", "Treatment_right ear" = "TR", 
                     "Control_left ear" = "UL", "Control_right ear" = "UR")
    } else if (prop == "Tissue") {
      plot_x_ticks = c("left ear" = "Left ear of\n treated mice", "right ear" = "Right ear of\n untreated mice")
    } else if (prop == "Time") {
      plot_x_ticks = c("120m" = "Timepoint 120min", "0m" = "Timepoint 0min")
    } else if (prop == "Tissue_Time") {
      plot_x_ticks = c("120m" = "Left ear of\n timepoint 120min", "0m" = "Right ear of\n timepoint 0min")
    } else if (prop == "Status_Time") {
        plot_x_ticks = c("Treatment" = "Treated mice\n for timepoint 120min", "Control" = "Untreated mice\n for timepoint 0min")
    } else {
      plot_x_ticks = c("Treatment" = "Treated mice\n for timepoint 120min", "Control" = "Untreated mice\n for timepoint 120min")
    }
    
    if (prop == "Tissue") {
      if (c$subset[2] == "Treatment" || c$subset[2] == "Treatment_120m" || c$subset[2] == "Treatment_0m") {
        tmp = str_subset(names(snakemake@params$colors$Group), "Treatment")
        colour_palette = snakemake@params$colors$Group[tmp]
        names(colour_palette) = str_split_fixed(names(colour_palette), "_", 2)[, 2]
      } else {
        tmp = str_subset(names(snakemake@params$colors$Group), "Control")
        colour_palette = snakemake@params$colors$Group[tmp]
        names(colour_palette) = str_split_fixed(names(colour_palette), "_", 2)[, 2]
      }
    } else if (prop == "Tissue_Time") {
        colour_palette = snakemake@params$colors$Time
    } else {
      colour_palette = snakemake@params$colors[[strsplit(prop, "_")[[1]][1]]]
    }
    
    # create plot 
    violin_plot <- ggplot(plot_df, aes(x = comp, y = expr, colour = comp, fill = comp)) +
      geom_violin(trim = TRUE) +
      scale_fill_manual(values = alpha(colour_palette, 0.3)) +
      scale_colour_manual(values = alpha(colour_palette, 1)) +
      geom_jitter(size = 0.75, position = position_jitter(0.2)) + 
      stat_summary(fun=median, geom="crossbar", size=.2, aes(color = comp)) + 
      scale_x_discrete(labels = plot_x_ticks, limits = names(plot_x_ticks)) +
      ggtitle(plot_title) +
      xlab("") +
      ylab("Expression value (rpmmm)") +
      theme_classic() +
      theme(legend.position="none", text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            axis.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            axis.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            plot.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size_header),
            legend.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            legend.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size))   
   
    # save plot as png and svg
    dir.create(output_folder, recursive=TRUE)
    ggsave(sprintf("%s/%s.png", output_folder, miR), plot = violin_plot, dpi = snakemake@params$plots_props$dpi, width = snakemake@params$plots_props$image_width, height = snakemake@params$plots_props$image_height, units = snakemake@params$plots_props$image_units)
    ggsave(sprintf("%s/%s.svg", output_folder, miR), plot = violin_plot, width = snakemake@params$plots_props$image_width, height = snakemake@params$plots_props$image_height, units = snakemake@params$plots_props$image_units, limitsize = FALSE, scale = 1)
   
    # create plot expr in log10
    violin_plot <- ggplot(plot_df, aes(x = comp, y = expr, colour = comp, fill = comp)) +
      geom_violin(trim = TRUE) +
      scale_fill_manual(values = alpha(colour_palette, 0.3)) +
      scale_colour_manual(values = alpha(colour_palette, 1)) +
      geom_jitter(size = 0.75, position = position_jitter(0.2)) + 
      stat_summary(fun=median, geom="crossbar", size=.2, aes(color = comp)) + 
      scale_x_discrete(labels = plot_x_ticks, limits = names(plot_x_ticks)) +
      ggtitle(plot_title) +
      xlab("") +
      ylab("Expression value (rpmmm, log10)") +
      theme_classic() +
      theme(legend.position="none", text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            axis.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            axis.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            plot.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size_header),
            legend.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            legend.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size))   
    
    # save plot as png and svg
    dir.create(output_folder, recursive=TRUE)
    ggsave(sprintf("%s/%s_log10.png", output_folder, miR), plot = violin_plot, dpi = snakemake@params$plots_props$dpi, width = snakemake@params$plots_props$image_width, height = snakemake@params$plots_props$image_height, units = snakemake@params$plots_props$image_units)
    ggsave(sprintf("%s/%s_log10.svg", output_folder, miR), plot = violin_plot, width = snakemake@params$plots_props$image_width, height = snakemake@params$plots_props$image_height, units = snakemake@params$plots_props$image_units, limitsize = FALSE, scale = 1)
    
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
