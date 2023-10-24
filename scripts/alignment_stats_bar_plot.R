suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(svglite))

set.seed(snakemake@params$parameters_props$set_seed)

colors = c("aligned"="#D5A021", "not aligned"="grey")
colors_group_time = snakemake@params$colors

print("bar_alignment")

results_folder = snakemake@params$results_folder
output_folder = sprintf("%s/bar_plot/alignment", results_folder)
dir.create(output_folder, recursive=TRUE)

# load data
annot = fread(snakemake@input$annot, sep='\t') 
mapping_info = fread(snakemake@input$mapping_info, sep='\t')

mapping_info_1 = mapping_info[mapping_info$Mismatches == 1, ]

mapping_info_samples = merge(mapping_info_1, annot , by = "fastq_name")

timepoints = unique(mapping_info_samples$Time)
status = unique(mapping_info_samples$Status)
ear_sides = unique(mapping_info_samples$Tissue)

group_name = c()
group_color = c()
time_color = c()
reads_aligned = c()
reads_processed = c()
for (timepoint in timepoints){
  for (treated in status){
    for (ear_side in ear_sides){
      ids = mapping_info_samples[(mapping_info_samples$Time == timepoint) & (mapping_info_samples$Status == treated) & (mapping_info_samples$Tissue == ear_side)]$fastq_name
      if (treated == "Control") {
        treated_label = "U"
      }
      else {
        treated_label = "T"
      }
      if (ear_side == "left ear") {
        ear_side_label = "L"
      }
      else {
        ear_side_label = "R"
      }
      group_name = append(group_name, paste(timepoint, paste(treated_label, ear_side_label, sep = "")))
      name_for_color = paste(treated, ear_side, sep = "_")
      group_color = append(group_color, name_for_color)
      time_color = append(time_color, timepoint)
      tmp_total = mean(mapping_info_samples[ids,]$reads_processed) / 1000000
      tmp_aligned = mean(mapping_info_samples[ids,]$reads_aligned) / 1000000
      reads_processed = append(reads_processed, tmp_total)
      reads_aligned = append(reads_aligned, tmp_aligned)
      
    }
  }
}

plot_df = data.frame(group_name=group_name, group_color=group_color, time_color=time_color, aligned=reads_aligned, processed=reads_processed, check.names = F)
plot_df["not aligned"] = plot_df$processed - plot_df$aligned
# remove the total column (will show as the sum in the bar plot)
plot_df$processed = c()
# melt for ggplot
plot_df = melt(plot_df, variable.name="alignment", value.name="counts")
plot_df$alignment = factor(plot_df$alignment, levels=c("not aligned", "aligned"))

# order decreasing by column "counts" (only values where column alignment == aligned)
# get the numeric order
ord = order(plot_df[plot_df$alignment == "aligned",]$counts, decreasing = TRUE)
# get the names
x_order = plot_df[ord,]$group_name

# merge all color mappings
all_colors = c(unlist(colors_group_time[["Group"]]), unlist(colors_group_time[["Time"]]), colors)
p = ggplot(plot_df, aes(x=factor(group_name, levels=x_order), y=counts, fill=alignment)) + 
            geom_bar(stat="identity") +
            labs(x="", y="mean processed reads\nin millions") +
            scale_fill_manual(breaks = c("aligned", "not aligned"), values=all_colors) +
            geom_tile(aes(x=factor(group_name, levels=x_order), y=-0.75, fill=group_color)) +
            geom_tile(aes(x=factor(group_name, levels=x_order), y=-2.0, fill=time_color)) +
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
ggsave(sprintf("%s/mapping_stats.png", output_folder), p, dpi=snakemake@params$plots_props$dpi, width=snakemake@params$plots_props$image_width, height=snakemake@params$plots_props$image_height, units = snakemake@params$plots_props$image_units)
ggsave(sprintf("%s/mapping_stats.svg", output_folder), p, width=snakemake@params$plots_props$image_width, height=snakemake@params$plots_props$image_height, units = snakemake@params$plots_props$image_units, limitsize = FALSE, scale = 1)


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])



