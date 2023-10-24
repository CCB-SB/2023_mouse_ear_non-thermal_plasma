suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))


set.seed(snakemake@params$parameters_props$set_seed)

print("box_plots")


#---------------------------------- Load data ---------------------------------- 
expr = fread(snakemake@input$expr_log10, sep='\t', header=T)
annot = fread(snakemake@input$annot, sep='\t', colClasses=c(ID="character"))
interesting_miRNA = fread(snakemake@input$miRNA_list)

polynom_degree = snakemake@params$polynom_degree

plots_props = snakemake@params$plots_props

results_folder = snakemake@params$results_folder
output_folder = sprintf("%s/box_plot_l=%s", results_folder, polynom_degree)
dir.create(output_folder, recursive=TRUE)

expr_rownames = expr$miRNA
expr$miRNA = c()
timepoints = gsub("m","", annot$Time)
annot$Timepoint = as.numeric(timepoints)


#------------------------------------ Script ----------------------------------- 
for (miR in interesting_miRNA$miRNA) {
  # create the table to include all the sample points into the plot
  plot_table_2_y = c() 
  plot_table_2_x = c()
  plot_table_2_groups = c()
  for (treated in unique(annot$Status)){
    for (tissue in unique(annot$Tissue)){
      for (tp in unique(sort(annot$Timepoint))){
        group = annot[(annot$Status == treated) & (annot$Tissue == tissue) & (annot$Timepoint == tp)]$ID
        plot_table_2_y = append(plot_table_2_y, as.numeric(expr[expr_rownames == miR,..group]))
        plot_table_2_x = append(plot_table_2_x, rep(list(tp),length(expr[expr_rownames == miR,..group])))
        tmp = rep(sprintf("%s_%s", treated, tissue), length(expr[expr_rownames == miR,..group]))
        plot_table_2_groups = append(plot_table_2_groups, tmp)
      }
    }
  }
  
  plot_table_2 = data.frame(Time = unlist(plot_table_2_x), expr = unlist(plot_table_2_y), Groups = unlist(plot_table_2_groups))
  
  # make Time and Groups factors and give Groups an order
  group_ticks = c("Treatment_left ear" = "Left ear\n of treated\n mice", "Treatment_right ear" = "Right ear\n of treated\n mice", 
                   "Control_left ear" = "Left ear\n of untreated\n mice", "Control_right ear" = "Right ear\n of untreated\n mice")
  
  plot_table_2$Time = factor(plot_table_2$Time, ordered = TRUE, levels = seq(0, 120, 10))
  
  # create line plot
  line_plot = ggplot(plot_table_2, aes(x=Time, y=expr, fill=factor(Groups, levels = names(group_ticks), ordered = TRUE))) +  # group=Groups, colour=Groups
    geom_boxplot(outlier.size = 0.25, lwd=0.2, fatten=3) +
    scale_fill_manual(values=unlist(snakemake@params$colors$Group)) +
    scale_color_manual(values=unlist(snakemake@params$colors$Group)) +
    scale_x_discrete(drop=F, breaks=c("0", "10", "30", "60", "120")) +
    geom_smooth(method = "lm", formula = y ~ poly(x, polynom_degree), aes(x=as.numeric(Time), y=expr, color=factor(plot_table_2$Groups, levels = names(group_ticks), ordered = TRUE)), se=FALSE) +
    ggtitle(sprintf("%s", miR)) +
    xlab("Timepoint (min)") +
    ylab("Expression value (rpmmm, log10)") +
    theme_classic() +
    theme(legend.position = "none",
          text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size)
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) 
  
  ggsave(sprintf("%s/%s.png", output_folder, miR), plot = line_plot, dpi = plots_props$dpi, width = plots_props$image_width, height = plots_props$image_height, units = plots_props$image_units)
  ggsave(sprintf("%s/%s.svg", output_folder, miR), plot = line_plot, width = plots_props$image_width, height = plots_props$image_height, units = plots_props$image_units, limitsize = FALSE, scale = 1)
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])