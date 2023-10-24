suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))

set.seed(snakemake@params$parameters_props$set_seed)

print("scatter expr")


# ------------------------------- Input data -----------------------------------
diff_exp_logs = snakemake@params$diff_exp_logs
test_values_th = snakemake@params$test_values_th
comps = snakemake@params$comparisons
colors = snakemake@params$colors
results_folder = snakemake@params$results_folder

# read data
diff_exp = fread(sprintf("%s/matrices/diff_exp/diff_exp%s.csv", results_folder, diff_exp_logs[1]), sep='\t', header=T)
diff_exp_log = fread(sprintf("%s/matrices/diff_exp/diff_exp%s.csv", results_folder, diff_exp_logs[2]), sep='\t', header=T)

annot = fread(snakemake@input$annot, sep='\t', colClasses=c(ID="character"))

# save folder
output_folder = sprintf("%s/scatter_plot", results_folder)
dir.create(output_folder, recursive=TRUE)


# ----------------------------------- 4 ears -----------------------------------
# ------------------------- plot fc versus expr median -------------------------
  
plot_df_all = c()
for(i in 1:length(comps)) {
  # which comp do we investigate and which groups are considered
  prop = names(comps)[i]
  c = comps[[i]]

  c$g1 = c$g[1]
  c$g2 = c$g[2]

  paired_string = ifelse(!is.null(c$paired), sprintf("_paired_%s", c$paired), "")
  sub_string = ifelse(!is.null(c$subset), sprintf("_%s=%s", c$subset[1], c$subset[2]), "")

  plot_df = data.table(RNA=diff_exp[["RNA"]],
                       expr=diff_exp[[sprintf("median_%s__%s%s", prop, c$g1, sub_string)]], 
                       expr_log=diff_exp_log[[sprintf("median_%s__%s%s", prop, c$g1, sub_string)]],
                       ttest_adj_pval=diff_exp[[sprintf("ttest_adjp_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                       fc=diff_exp[[sprintf("fc_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                       auc=diff_exp[[sprintf("auc_value_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                       group=rep(sprintf("%s", c$subset[2]), length(diff_exp[["RNA"]])),
                       time=rep(sprintf("%s", c$g1), length(diff_exp[["RNA"]]))
                       )
  plot_df[,expr := expr_log]
  plot_df[,ttest_adj_pval_log10 := -log10(ttest_adj_pval)]
  plot_df[,logfc := log2(fc)]
  
  plot_df_all = rbind(plot_df_all, plot_df)
}

for (test_value in c("ttest_adj_pval_log10", "logfc", "auc")) {
  xlab_name = sprintf("Expression value (rpmmm, %s)", diff_exp_logs[2])
  if (test_value == "ttest_adj_pval_log10") {
    ylab_name = "Adj. p-value from ttest (-log10)"
  } else if (test_value == "logfc") {
    ylab_name = "Fold change (log2)"
  } else {
    ylab_name = "Area under curve"
  }

  # coloured by Group
  scatter_plot <- ggplot(plot_df_all, aes(x=expr, y=.data[[test_value]], color=group, shape=time, text="ID")) +
    geom_point(size = 2, alpha = 0.75) + 
    scale_colour_manual(values = colors[["Status"]]) +
    scale_shape() +
    xlab(xlab_name) +
    ylab(ylab_name) + 
    theme_classic() +
    theme(legend.position="bottom", 
          text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
          axis.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
          axis.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
          plot.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size_header),
          legend.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
          legend.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size, face = "bold"),
          legend.margin=margin(t = -8), legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.4, "cm")) +
    guides(col="none", shape=guide_legend(nrow=2, title="Timepoint"))
  
  if (test_value == "ttest_adj_pval_log10") {
    th = -log10(test_values_th$adj_p_value)
    scatter_plot = scatter_plot + geom_hline(yintercept=th, linetype="dashed", color="lightgrey", size=0.5)
  } else if (test_value == "logfc") {
    th_1 = log2(test_values_th$fc)
    th_2 = log2(1/test_values_th$fc)
    scatter_plot = scatter_plot + 
      geom_hline(yintercept=th_1, linetype="dashed", color="lightgrey", size=0.5) + 
      geom_hline(yintercept=th_2, linetype="dashed", color="lightgrey", size=0.5) +
      geom_vline(xintercept=log10(15), linetype="dashed", color="lightgrey", size=0.5) 
  }
  
  ggsave(sprintf("%s/%s_vs_expr_coloured=Group.png", output_folder, test_value), plot = scatter_plot, dpi = snakemake@params$plots_props$dpi, width = snakemake@params$plots_props$image_width, height = snakemake@params$plots_props$image_height, units = snakemake@params$plots_props$image_units)
  ggsave(sprintf("%s/%s_vs_expr_coloured=Group.svg", output_folder, test_value), plot = scatter_plot, width = snakemake@params$plots_props$image_width, height = snakemake@params$plots_props$image_height, units = snakemake@params$plots_props$image_units, limitsize = FALSE, scale = 1)

  # coloures by Time
  scatter_plot <- ggplot(plot_df_all, aes(x=expr, y=.data[[test_value]], color=time, shape=group, text="ID")) +
    geom_point(size = 2, alpha = 0.8) + 
    scale_colour_manual(values = colors[["Time"]]) +
    scale_shape(breaks = c("Treatment", "Control")) +
    xlab(xlab_name) +
    ylab(ylab_name) + 
    theme_classic() +
    theme(legend.position="bottom",
          text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
          axis.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
          axis.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
          plot.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size_header),
          legend.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
          legend.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size)) +
    guides(col="none", shape=guide_legend(nrow=1, title="Status"))  
  
  if (test_value == "ttest_adj_pval_log10") {
    th = -log10(test_values_th$adj_p_value)
    scatter_plot = scatter_plot + geom_hline(yintercept=th, linetype="dashed", color="lightgrey", size=0.5) 
  } else if (test_value == "logfc") {
    th_1 = log2(test_values_th$fc)
    th_2 = log2(1/test_values_th$fc)
    scatter_plot = scatter_plot + 
      geom_hline(yintercept=th_1, linetype="dashed", color="lightgrey", size=0.5) + 
      geom_hline(yintercept=th_2, linetype="dashed", color="lightgrey", size=0.5) +
      geom_vline(xintercept=log10(150), linetype="dashed", color="lightgrey", size=0.5) 
  }
  
  ggsave(sprintf("%s/%s_vs_expr_coloured=Time.png", output_folder, test_value), plot = scatter_plot, dpi = snakemake@params$plots_props$dpi, width = snakemake@params$plots_props$image_width, height = snakemake@params$plots_props$image_height, units = snakemake@params$plots_props$image_units)
  ggsave(sprintf("%s/%s_vs_expr_coloured=Time.svg", output_folder, test_value), plot = scatter_plot, width = snakemake@params$plots_props$image_width, height = snakemake@params$plots_props$image_height, units = snakemake@params$plots_props$image_units, limitsize = FALSE, scale = 1) 
}

  
#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

  
