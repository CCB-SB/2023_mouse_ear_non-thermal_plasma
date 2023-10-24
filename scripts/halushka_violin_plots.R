suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library("readxl"))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(gridExtra))

set.seed(snakemake@params$parameters_porps$set_seed)

print("halushka")

# functions
median_fun <- function(x){
  return(data.frame(y = median(x), label = "median"))
}

miRNA_list = c("hsa-miR-144-3p", "hsa-miR-451a", "hsa-miR-223-3p", "hsa-miR-142-5p")

colour_list_cell_classes = c("Immune" = "#42B93E", "Fat" = "#9BC9D5", "Stem" = "#6194A3", "RBC" = "#870005", "Plasma" = "#E61549", "Platelet" = "#F47024", "Epithelial" = "#E89CFF", "Endothelial" = "#B60095", "Fibroblast" = "#710093", "Brain" = "#B5B4B4", "Muscle" = "#6F7A85", "Sperm" = "#525457", "Other" = "#000000")

# load data 
expr_deseq2_norm = read_excel("data/external_data/giac083_supplemental_files/Supplementary_Table_S7.xlsx")
cell_type_mapping = read_excel("data/external_data/giac083_supplemental_files/Supplementary_Table_S2.xlsx")

# create save folder 
results_folder = snakemake@params$results_folder
output_folder = sprintf("%s/halushka_plots", results_folder)
dir.create(output_folder, recursive=TRUE)

# script
# filter for relevant miRNAs
tmp = expr_deseq2_norm$miRNA %in% miRNA_list
expr_deseq2_norm_filtered = data.frame(expr_deseq2_norm[tmp,])

expr_deseq2_norm_filtered_melt = data.frame(melt(expr_deseq2_norm_filtered))

cell_type_list = c()
for (sample in unique(expr_deseq2_norm_filtered_melt$variable)) {
  cell_type_mapping_filtered = cell_type_mapping[cell_type_mapping$Sample == sample,]
  if (dim(cell_type_mapping_filtered)[1] > 1) {
    if (length(unique(cell_type_mapping_filtered$CellType)) != 1) {
      print(sprintf("For %s we have different cell types given", sample))
    } else {
      cell_type_list[[sample]] = cell_type_mapping_filtered$CellType[1]
    } 
  } else {
    cell_type_list[[sample]] = cell_type_mapping_filtered$CellType
  }  
}  

cell_type_list_melt = data.frame(melt(cell_type_list))
colnames(cell_type_list_melt) = c("cell_type", "variable")

plot_table = merge(expr_deseq2_norm_filtered_melt, cell_type_list_melt, by = 'variable')

# To colour the celltypes by their classes
celltypes = unique(plot_table$cell_type)
cellclasses = c()
for (sample in plot_table$variable){
  tmp = cell_type_mapping[cell_type_mapping$Sample == sample, ]$Class
  if (length(tmp) == 1) {
    cellclasses = append(cellclasses, tmp)
  } else if (length(unique(tmp)) != 1) {
    print(sprintf("For %s we have different cell classes given", sample))
  } else {
    cellclasses = append(cellclasses, tmp[1])
  }
}

plot_table$cell_class = cellclasses
colour_palette = colour_list_cell_classes

legend_title = "Cell class"
for (num_cell_type in c(0, 50)) {
  p <- list()
  i = 1
  for (miR in unique(plot_table$miRNA)) {
    plot_table_filtered = plot_table[plot_table$miRNA == miR,]
  
    plot_table_filtered$cell_type_without_spacing = gsub("_", " ", plot_table_filtered$cell_type)
    
    # here we calculate the same what reorder does in the ggplot call
    # the "with" is needed but I do not know why
    # we then use levels to get the first 50 (or all) levels which are now ordered
    # this vector is then used as xlim
    tmp = with(plot_table_filtered, reorder(cell_type_without_spacing, -value, FUN = median))
    if (num_cell_type == 0) {
      x_show = levels(tmp)
    }
    else {
      x_show = levels(tmp)[1:num_cell_type]
    }

    p[[i]] <- ggplot(plot_table_filtered, aes(x = reorder(cell_type_without_spacing, -value, FUN = median), y = value, colour=cell_class)) +
      geom_violin()  +
      geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.5) + 
      scale_colour_manual(legend_title, values = colour_palette) +
      stat_summary(fun=median, geom="crossbar", size=.2, show.legend = FALSE) + # color="black"
      xlim(x_show) +
      xlab("") + 
      ylab("Expression values\n(normalizedCounts_DESeq)") +
      ggtitle(sprintf("%s", miR)) +
      theme_classic() +
      theme(legend.position="none",
            text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            axis.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            axis.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            plot.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size_header),
      )
    if (i == 1) {
      p[[i]] = p[[i]] + theme(plot.margin = unit(c(5.5, 5.5, -2, 5.5), "pt")) #t, , b, 
    } else if (i == length(unique(plot_table$miRNA))){
      p[[i]] = p[[i]] + theme(legend.position="none", 
                              legend.key=element_rect(fill=colour_palette), 
                              legend.spacing.y = unit(0.2, "cm"), legend.spacing.x = unit(0.5, "cm"), legend.key.size = unit(0.3, "cm"), legend.margin = margin(t = -18, l = -25),
                              plot.margin = unit(c(5.5, 5.5, -2, 5.5), "pt"), #t, , b, 
                              legend.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
                              legend.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size, face = "bold")
      )
    } else
      p[[i]] = p[[i]] + theme(plot.margin = unit(c(5.5, 5.5, -2, 5.5), "pt")) #t, , b, 
    
    if (num_cell_type == 0) {
      image_width = snakemake@params$plots_props$image_width*2
      p[[i]] = p[[i]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
    } else {
      image_width = snakemake@params$plots_props$image_width
      p[[i]] = p[[i]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
    }
    
    ggsave(sprintf("%s/%s_sorted_by_median_num_cell_type=%s.png", output_folder, miR, num_cell_type), p[[i]], dpi=snakemake@params$plots_props$dpi, width=image_width, height=snakemake@params$plots_props$image_height, units = snakemake@params$plots_props$image_units)
    ggsave(sprintf("%s/%s_sorted_by_median_num_cell_type=%s.svg", output_folder, miR, num_cell_type), p[[i]], width=image_width, height=snakemake@params$plots_props$image_height, units = snakemake@params$plots_props$image_units, limitsize = FALSE, scale = 1)
    i = i + 1
  }
  
  # creating concatinated plot
  if (num_cell_type == 0) {
    all_p = arrangeGrob(grobs=p, ncol=1)  
    image_height_all = snakemake@params$plots_props$image_height * length(p)
  } else {
    all_p = arrangeGrob(grobs=p, ncol=2, nrow = 2)
    image_height_all = snakemake@params$plots_props$image_height * 2
  }
  
  # Save results to storage
  ggsave(sprintf("%s/all_miRNAs_sorted_by_median_num_cell_type=%s.png", output_folder, num_cell_type), all_p, dpi=snakemake@params$plots_props$dpi, width=snakemake@params$plots_props$image_width*2, height=image_height_all, units = snakemake@params$plots_props$image_units)
  ggsave(sprintf("%s/all_miRNAs_sorted_by_median_num_cell_type=%s.svg", output_folder, num_cell_type), all_p, width=snakemake@params$plots_props$image_width*2, height=image_height_all, units = snakemake@params$plots_props$image_units, limitsize = FALSE, scale = 1)
  ggsave(sprintf("%s/all_miRNAs_sorted_by_median_height=18_num_cell_type=%s.png", output_folder, num_cell_type), all_p, dpi=snakemake@params$plots_props$dpi, width=snakemake@params$plots_props$image_width*2, height=snakemake@params$plots_props$image_height * 3, units = snakemake@params$plots_props$image_units)
  ggsave(sprintf("%s/all_miRNAs_sorted_by_median_height=18_num_cell_type=%s.svg", output_folder, num_cell_type), all_p, width=snakemake@params$plots_props$image_width*2, height=snakemake@params$plots_props$image_height * 3, units = snakemake@params$plots_props$image_units, limitsize = FALSE, scale = 1)
  }


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

