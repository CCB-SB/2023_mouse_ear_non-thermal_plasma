suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))


set.seed(snakemake@params$parameters_porps$set_seed)

print("heatmap")

#---------------------------------- Functions ----------------------------------
plot_heatmap = function(df, row_names, plots_props, filename, test_value) {

  if(test_value == "fc") {
    color_bar_name = "Fold change (log2)" 
    legend_ticks = seq(from = -max(t(df)), to = max(t(df)), length.out = 5)
    mini = min(df)
    maxi = max(df)
    col_fun = colorRamp2(c(-maxi,0,maxi), hcl_palette = "Green-Brown", reverse = TRUE)
    col_fun_legend = col_fun
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (t(df)[i, j] >= log2(5)){
        grid.text("***", x, y, gp = gpar(fontsize = 6, col = "white"))
      } else if (t(df)[i, j] <= log2(1/5)){
        grid.text("***", x, y, gp = gpar(fontsize = 6, col = "white"))
      } else if (t(df)[i, j] >= log2(2)){
        grid.text("**", x, y, gp = gpar(fontsize = 6, col = "white"))
      } else if (t(df)[i, j] <= log2(1/2)){
        grid.text("**", x, y, gp = gpar(fontsize = 6, col = "white"))
      } else if (t(df)[i, j] >= log2(1.5)){
        grid.text("*", x, y, gp = gpar(fontsize = 6))
      } else if (t(df)[i, j] <= log2(1/1.5)) {
        grid.text("*", x, y, gp = gpar(fontsize = 6))
      }
    }
  } else if (test_value == "ttest_adjp") {
    color_bar_name = "Adj. p-value\nfrom t-test (-log10)"
    legend_ticks = seq(from = 0, to = max(t(df)), length.out = 5)
    maxi = max(df)
    col_fun = colorRamp2(c(0,maxi), hcl_palette = "Reds 3", reverse = TRUE)
    col_fun_legend = col_fun
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (t(df)[i, j] >= -log10(0.001)){
        grid.text("***", x, y, gp = gpar(fontsize = 6, col = "white"))
      } else if (t(df)[i, j] >= -log10(0.01)){ 
        grid.text("**", x, y, gp = gpar(fontsize = 6))
      } else if (t(df)[i, j] >= -log10(0.05)){ 
        grid.text("*", x, y, gp = gpar(fontsize = 6))
      }
    }
  } else {
    color_bar_name = "Area under curve"
    legend_ticks = seq(from = 0, to = 1, length.out = 5)
    col_fun = colorRamp2(c(0,0.5,1), hcl_palette = "Green-Brown", reverse = TRUE)
    col_fun_legend = col_fun
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (t(df)[i, j] >= 0.9){ 
        grid.text("***", x, y, gp = gpar(fontsize = 6, col = "white"))
      } else if (t(df)[i, j] <= 0.1){
        grid.text("***", x, y, gp = gpar(fontsize = 6, col = "white"))
      } else if (t(df)[i, j] >= 0.8){ 
        grid.text("**", x, y, gp = gpar(fontsize = 6))
      } else if (t(df)[i, j] <= 0.2){
        grid.text("**", x, y, gp = gpar(fontsize = 6))
      } else if (t(df)[i, j] >= 0.7){ 
        grid.text("*", x, y, gp = gpar(fontsize = 6))
      } else if (t(df)[i, j] <= 0.3){
        grid.text("*", x, y, gp = gpar(fontsize = 6))
      }
    }
  }
  
  legend_tick_labels = round(legend_ticks, digits = 2)

  rownames(df) = row_names

  p = ComplexHeatmap::Heatmap(t(df), col = col_fun,
                              rect_gp = gpar(col = "white", lwd = .5),
                              cell_fun = cell_fun,
                              height = unit(5.5, "cm"), width = unit(4.2, "cm"),
                              row_names_side = "left", row_names_gp = gpar(fontsize = 4, fontfamily=plots_props$font_family), column_names_side = "bottom", column_names_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                              column_dend_reorder = FALSE, cluster_columns = FALSE,
                              row_dend_reorder = FALSE, cluster_rows = FALSE,
                              show_heatmap_legend = FALSE
                              )
  legend = ComplexHeatmap::Legend(title = color_bar_name, col_fun = col_fun_legend, at = legend_tick_labels, direction="horizontal", #col_fun = circlize::colorRamp2(legend_ticks, viridis(5))
                                  title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                                  legend_width=unit(4, "cm"))
  
  png(sprintf("%s.png",filename), width = plots_props$image_width, height = 9, units = plots_props$image_units, res = plots_props$dpi)
  draw(p)
  draw(legend, x = unit(.5, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
  dev.off()
  svglite(sprintf("%s.svg",filename), width = plots_props$image_width / 2.54, height = 9 / 2.54)
  draw(p)
  draw(legend, x = unit(.5, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
  dev.off()
}


#---------------------------------- Load data ---------------------------------- 
diff_exp_logs = snakemake@params$diff_exp_logs
test_values = snakemake@params$test_values
comps = snakemake@params$comparisons
plots_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder

diff_exp = fread(sprintf("%s/matrices/diff_exp/diff_exp%s.csv", results_folder, diff_exp_logs[1]), sep='\t', header=T)
interesting_miRNA = fread(snakemake@input$miRNA_list)


#------------------------------------ Script ----------------------------------- 
for (test_value in test_values) {
  
  comp_col_names = c()
  all_test_value_mat_df_filtered = data.frame(RNA=interesting_miRNA$miRNA)
  for (i in 1:length(comps)) {
    
    prop = names(comps)[i]
    c = comps[[i]]
    
    c$g1 = c$g[1]
    c$g2 = c$g[2]
        
    prop = names(comps)[i]
    
    paired_string = ifelse(!is.null(c$paired), sprintf("_paired_%s", c$paired), "")
    sub_string = ifelse(!is.null(c$subset), sprintf("_%s=%s", c$subset[1], c$subset[2]), "")
    
    plot_df = data.table(RNA=diff_exp[["RNA"]],
                         comp_values=diff_exp[[sprintf("%s_%s__%s_vs_%s%s%s", test_value, prop, c$g1, c$g2, paired_string, sub_string)]])
    
    if (test_value == "ttest_adjp") {
      plot_df[,comp_values_log := -log10(comp_values)]
    } else if (test_value == "fc") {
      plot_df[,comp_values_log := log2(comp_values)]
    } else {
      plot_df[,comp_values_log := comp_values]
    }
    
    tmp = plot_df$RNA %in% interesting_miRNA$miRNA
    test_value_mat_df_filtered = plot_df[tmp,]
    rownames(test_value_mat_df_filtered) = c()
    
    all_test_value_mat_df_filtered = merge(all_test_value_mat_df_filtered, test_value_mat_df_filtered[, -"comp_values"], by="RNA")
    comp_col_names = append(comp_col_names, sprintf("%s_vs_%s%s", c$g1, c$g2, ifelse(!is.null(c$subset[2]), sprintf("_%s",c$subset[2]), "")) )
  }
  
  output_folder = sprintf("%s/heatmap_complex", results_folder)
  dir.create(output_folder, recursive=TRUE)
  
  colnames(all_test_value_mat_df_filtered) = c("RNA", comp_col_names)
  
  plot_heatmap(all_test_value_mat_df_filtered[,-1], all_test_value_mat_df_filtered$RNA, plots_props, sprintf("%s/%s", output_folder, test_value), test_value)
  
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
