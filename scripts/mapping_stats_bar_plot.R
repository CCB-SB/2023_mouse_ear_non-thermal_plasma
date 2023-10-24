suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))


source("scripts/iwanthue.R")

set.seed(snakemake@params$parameters_props$set_seed)

print("bar_mapping")

colors = c("miRNA" = "#788A38", "other" = "grey")
colors_group_time = snakemake@params$colors
class_list = c("miRNA")

results_folder = snakemake@params$results_folder
output_folder = sprintf("%s/bar_plot/mapping", results_folder)
dir.create(output_folder, recursive=TRUE)

annot_complete = fread(snakemake@input$annot, sep='\t') 

# mapped
rna_composition_ori = fread(snakemake@input$rna_composition_mapped, sep='\t')

rna_composition_ori$ID = rna_composition_ori$Sample
rna_composition = merge(rna_composition_ori, annot_complete, by = "ID")

timepoints = unique(rna_composition$Time)
status = unique(rna_composition$Status)
ear_sides = unique(rna_composition$Tissue)

group_name = c()
group_color = c()
time_color = c()
reads_mapped = c()
for (timepoint in timepoints){
  for (treated in status){
    for (ear_side in ear_sides){
      ids = rna_composition[(rna_composition$Time == timepoint) & (rna_composition$Status == treated) & (rna_composition$Tissue == ear_side),]$ID
      if (treated == "Control") {
        treated_label = "U"
      } else {
        treated_label = "T"
      }
      if (ear_side == "left ear") {
        ear_side_label = "L"
      } else {
        ear_side_label = "R"
      }
      group_name = append(group_name, paste(timepoint, paste(treated_label, ear_side_label, sep = "")))
      name_for_color = paste(treated, ear_side, sep = "_")
      group_color = append(group_color, name_for_color)
      time_color = append(time_color, timepoint)
      tmp_mapped = c()
      for (RNA_class in class_list){
        tmp_mapped = append(tmp_mapped, mean(rna_composition[ids,][[RNA_class]]))
      }
      names(tmp_mapped) = class_list
      reads_mapped = rbind(reads_mapped, t(data.frame(tmp_mapped)))
    }
  }
}
reads_mapped_df = data.frame(reads_mapped)
rownames(reads_mapped_df) = c()
reads_mapped_df$other = 100 - rowSums(reads_mapped_df)
reads_mapped_df$Group = group_name
reads_mapped_df$group_color = group_color
reads_mapped_df$time_color = time_color

# melt for ggplot
plot_df = melt(reads_mapped_df, variable.name="mapping", value.name="percentage")
plot_df$mapping = factor(plot_df$mapping, levels=rev(c(class_list, "other")))

# order decreasing by column "percentage" (only values where column mapping == miRNA)
# get the numeric order
ord = order(plot_df[plot_df$mapping == "miRNA",]$percentage, decreasing = TRUE)
# get the names
x_order = plot_df[ord,]$Group

# merge all color mappings
all_colors = c(unlist(colors_group_time[["Group"]]), unlist(colors_group_time[["Time"]]), colors)
p = ggplot(plot_df, aes(x=factor(Group, levels=x_order), y=percentage, fill=mapping)) + 
  geom_bar(stat="identity") +
  labs(x="", y="Mapped reads (%)") +
  scale_fill_manual(breaks=c("miRNA", "other"), values=all_colors) +
  geom_tile(aes(x=factor(Group, levels=x_order), y=-3, fill=group_color), height=4) +
  geom_tile(aes(x=factor(Group, levels=x_order), y=-8, fill=time_color), height=4) +
  theme_classic() +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
        axis.text = element_text(angle = 0, family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size * 2/3),
        axis.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
        plot.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size_header),
        legend.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
        # make x invisible
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.margin=margin(t = -8), legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.4, "cm")
  )

# Save results to storage
ggsave(sprintf("%s/rna_composition_mapped.png", output_folder), p, dpi=snakemake@params$plots_props$dpi, width=snakemake@params$plots_props$image_width, height=snakemake@params$plots_props$image_height, units = snakemake@params$plots_props$image_units)
ggsave(sprintf("%s/rna_composition_mapped.svg", output_folder), p, width=snakemake@params$plots_props$image_width, height=snakemake@params$plots_props$image_height, units = snakemake@params$plots_props$image_units, limitsize = FALSE, scale = 1)


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
