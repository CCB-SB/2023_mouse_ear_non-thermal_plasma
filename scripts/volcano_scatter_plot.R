suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ggrepel))

set.seed(snakemake@params$parameters_porps$set_seed)

print("volcano")

diff_exp_logs = snakemake@params$diff_exp_logs
plots_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder

tbl = fread(sprintf("%s/matrices/diff_exp/diff_exp%s.csv", results_folder, diff_exp_logs[1]), sep='\t', header=T)

repel_max_time =  10 

output_folder = sprintf("%s/volcano", results_folder)
dir.create(output_folder)

font_size = 6
font_size_header = 9

comps = snakemake@params$comparisons
log_fc_thres = log2(snakemake@params$parameters$updownregulated)
sig_lvl = snakemake@params$parameters$significance
effect_size_thres = snakemake@params$parameters$mineffectsize

s_string = sprintf("Significant (P<%s)", sig_lvl)
ns_string = sprintf("Not significant (P\u2265%s)", sig_lvl)
up_string = sprintf("Sig. upregulated (FC\u2265%s)", snakemake@params$parameters$updownregulated)
down_string = sprintf("Sig. downregulated (FC\u22641/%s)", snakemake@params$parameters$updownregulated)
e_string = sprintf("Considerable effect (d\u2265%s)", effect_size_thres)
ne_string = sprintf("Neglectable effect (d<%s)", effect_size_thres)
e_up_string = sprintf("Considerably upregulated (FC\u2265%s)", snakemake@params$parameters$updownregulated)
e_down_string = sprintf("Considerably downregulated (FC\u22641/%s)", snakemake@params$parameters$updownregulated)
colors = c()
colors[ns_string] = "lightgrey"
colors[up_string] = "#31363B"
colors[down_string] = "#31363B"
colors[s_string] = "lightgrey"
colors[e_string] = "lightgrey"
colors[ne_string] = "lightgrey"
colors[e_up_string] = "#31363B"
colors[e_down_string] = "#31363B"

# set font-family for all plots
theme_set(theme_cowplot(font_family=snakemake@params$plots_props$font_family))

groups = list()
group_y_max = list()
group_x_max = list()
group_x_min = list()
group_plots = list()

for(i in 1:length(comps)){
  c = comps[[i]]
  
  c$g1 = c$g[1]
  c$g2 = c$g[2]
  
  group <- NULL
  if("volcano" %in% names(c)) {
    group <- c$volcano
    if(!(group %in% names(groups))) {
      groups[[group]] <- TRUE
      group_y_max[[group]] <- list("ttest_raw_pval" = -Inf,
                                    "ttest_adj_pval" = -Inf,
                                    "wilcox_raw_pval" = -Inf,
                                    "wilcox_adj_pval" = -Inf,
                                    "cohend_estimate" = -Inf)
      group_x_max[[group]] <- list("logfc" = -Inf)
      group_x_min[[group]] <- list("logfc" = Inf)
      group_plots[[group]] <- list()
    }
  }
  
  prop = names(comps)[i]
 
  paired_string = ifelse(!is.null(c$paired), sprintf("_paired_%s", c$paired), "")
  sub_string = ifelse(!is.null(c$subset), sprintf("_%s=%s", c$subset[1], c$subset[2]), "")
 
  plot_df = data.table(RNA=tbl[["RNA"]],
                       ttest_raw_pval=tbl[[sprintf("ttest_rawp_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                       ttest_adj_pval=tbl[[sprintf("ttest_adjp_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                       wilcox_raw_pval=tbl[[sprintf("wilcox_rawp_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                       wilcox_adj_pval=tbl[[sprintf("wilcox_adjp_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                       fc=tbl[[sprintf("fc_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                       cohend=tbl[[sprintf("cohend_estimate_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]])
  plot_df[,ttest_raw_pval_log10 := -log10(ttest_raw_pval)]
  plot_df[,ttest_adj_pval_log10 := -log10(ttest_adj_pval)]
  plot_df[,wilcox_raw_pval_log10 := -log10(wilcox_raw_pval)]
  plot_df[,wilcox_adj_pval_log10 := -log10(wilcox_adj_pval)]
  plot_df[,logfc := log2(fc)]
  plot_df[,cohend_abs := abs(cohend)]
  
  
  if(!is.null(group)) {
    group_y_max[[group]] <- list("ttest_raw_pval" = max(plot_df[,ttest_raw_pval_log10], group_y_max[[group]][["ttest_raw_pval"]], na.rm=T),
               "ttest_adj_pval" = max(plot_df[,ttest_adj_pval_log10], group_y_max[[group]][["ttest_adj_pval"]], na.rm=T),
               "wilcox_raw_pval" = max(plot_df[,wilcox_raw_pval_log10], group_y_max[[group]][["wilcox_raw_pval"]], na.rm=T),
               "wilcox_adj_pval" = max(plot_df[,wilcox_adj_pval_log10], group_y_max[[group]][["wilcox_adj_pval"]], na.rm=T),
               "cohend_abs" = max(plot_df[,cohend_abs], group_y_max[[group]][["cohend_abs"]], na.rm=T))
    group_x_max[[group]] <- list("logfc" = max(plot_df[,logfc], group_x_max[[group]][["logfc"]], na.rm=T))
    group_x_min[[group]] <- list("logfc" = min(plot_df[,logfc], group_x_min[[group]][["logfc"]], na.rm=T))
  }
  
  plot_df[ttest_raw_pval >= sig_lvl, ttest_raw_cat:=ns_string]
  plot_df[ttest_raw_pval < sig_lvl, ttest_raw_cat:=s_string]
  plot_df[logfc >= log_fc_thres & ttest_raw_pval < sig_lvl, ttest_raw_cat:=up_string]
  plot_df[logfc <= -log_fc_thres & ttest_raw_pval < sig_lvl, ttest_raw_cat:=down_string]
  
  plot_df[ttest_adj_pval >= sig_lvl, ttest_adj_cat:=ns_string]
  plot_df[ttest_adj_pval < sig_lvl, ttest_adj_cat:=s_string]
  plot_df[logfc >= log_fc_thres & ttest_adj_pval < sig_lvl, ttest_adj_cat:=up_string]
  plot_df[logfc <= -log_fc_thres & ttest_adj_pval < sig_lvl, ttest_adj_cat:=down_string]
  
  plot_df[wilcox_raw_pval >= sig_lvl, wilcox_raw_cat:=ns_string]
  plot_df[wilcox_raw_pval < sig_lvl, wilcox_raw_cat:=s_string]
  plot_df[logfc >= log_fc_thres & wilcox_raw_pval < sig_lvl, wilcox_raw_cat:=up_string]
  plot_df[logfc <= -log_fc_thres & wilcox_raw_pval < sig_lvl, wilcox_raw_cat:=down_string]
  
  plot_df[wilcox_adj_pval >= sig_lvl, wilcox_adj_cat:=ns_string]
  plot_df[wilcox_adj_pval < sig_lvl, wilcox_adj_cat:=s_string]
  plot_df[logfc >= log_fc_thres & wilcox_adj_pval < sig_lvl, wilcox_adj_cat:=up_string]
  plot_df[logfc <= -log_fc_thres & wilcox_adj_pval < sig_lvl, wilcox_adj_cat:=down_string]

  plot_df[cohend_abs >= effect_size_thres, cohend_abs_cat:=e_string]
  plot_df[cohend_abs < effect_size_thres, cohend_abs_cat:=ne_string]
  plot_df[logfc >= log_fc_thres & cohend_abs >= effect_size_thres, cohend_abs_cat:=e_up_string]
  plot_df[logfc <= -log_fc_thres & cohend_abs >= effect_size_thres, cohend_abs_cat:=e_down_string]
  
  colors_ttest_raw = colors[names(colors) %in% plot_df$ttest_raw_cat]
  colors_ttest_adj = colors[names(colors) %in% plot_df$ttest_adj_cat]
  colors_wilcox_raw = colors[names(colors) %in% plot_df$wilcox_raw_cat]
  colors_wilcox_adj = colors[names(colors) %in% plot_df$wilcox_adj_cat]
  colors_cohend_estimate = colors[names(colors) %in% plot_df$cohend_abs_cat]
  sub_plot_df = plot_df[!is.na(ttest_raw_pval_log10) & !is.na(logfc) & !is.na(cohend_abs)]

  sub_plot_df$ttest_raw_cat_show = sub_plot_df$RNA
  sub_plot_df$ttest_adj_cat_show = sub_plot_df$RNA
  sub_plot_df$wilcox_raw_cat_show = sub_plot_df$RNA
  sub_plot_df$wilcox_adj_cat_show = sub_plot_df$RNA
  sub_plot_df$cohend_abs_cat_show = sub_plot_df$RNA
  mask_ttest_raw = ((sub_plot_df$ttest_raw_cat != up_string) & (sub_plot_df$ttest_raw_cat != down_string)) 
  mask_ttest_adj = ((sub_plot_df$ttest_adj_cat != up_string) & (sub_plot_df$ttest_adj_cat != down_string)) 
  mask_wilcox_raw = ((sub_plot_df$wilcox_raw_cat != up_string) & (sub_plot_df$wilcox_raw_cat != down_string)) 
  mask_wilcox_adj = ((sub_plot_df$wilcox_adj_cat != up_string) & (sub_plot_df$wilcox_adj_cat != down_string)) 
  mask_cohend_abs = ((sub_plot_df$cohend_abs_cat != e_up_string) & (sub_plot_df$cohend_abs_cat != e_down_string)) 
  sub_plot_df[mask_ttest_raw,]$ttest_raw_cat_show = ""
  sub_plot_df[mask_ttest_adj,]$ttest_adj_cat_show = ""
  sub_plot_df[mask_wilcox_raw,]$wilcox_raw_cat_show = ""
  sub_plot_df[mask_wilcox_adj,]$wilcox_adj_cat_show = ""
  sub_plot_df[mask_cohend_abs,]$cohend_abs_cat_show = ""

  p1 = ggplot(sub_plot_df, aes(x=logfc, y=ttest_raw_pval_log10, color=ttest_raw_cat)) + 
    geom_point(size = 0.75) +
    geom_hline(yintercept = -log10(sig_lvl), color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    geom_text_repel(aes(label = ttest_raw_cat_show),
                    max.time = repel_max_time, max.iter = 1e6,
                    max.overlaps = Inf, 
                    box.padding = 0.1, 
                    point.patting = 0.1, 
                    force_pull = 2, 
                    force = 100, 
                    min.segment.length = 0, segment.size = 0.25,
                    seed = snakemake@params$parameters_porps$set_seed,
                    size = 6/11.04*3.88, family = snakemake@params$plots_props$font_family, color = "black") +
    xlab("Fold change (log2)") + ylab("Raw p-value\nfrom t-test (-log10)") +
    scale_color_manual(name="", values=colors_ttest_raw) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          plot.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size_header), 
          legend.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size),
          legend.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size, face = "bold")) + 
    guides(col=guide_legend(nrow=length(colors_ttest_raw)))

  p2 = ggplot(sub_plot_df, aes(x=logfc, y=ttest_adj_pval_log10, color=ttest_adj_cat)) + 
    geom_point(size = 0.75) +
    geom_hline(yintercept = -log10(sig_lvl), color="lightgrey", linetype="dashed", size = 0.25) +
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    geom_text_repel(aes(label = ttest_adj_cat_show),
                    max.time = repel_max_time, max.iter = 1e6,
                    max.overlaps = Inf,
                    box.padding = 0.1, 
                    point.patting = 0.1, 
                    force_pull = 2, 
                    force = 100, 
                    min.segment.length = 0, segment.size = 0.25,
                    seed = snakemake@params$parameters_porps$set_seed,
                    size = 6/11.04*3.88, family = snakemake@params$plots_props$font_family, color = "black") +
    xlab("Fold change (log2)") + ylab("Adj. p-value\nfrom t-test (-log10)") +
    scale_color_manual(name="", values=colors_ttest_adj) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          plot.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size_header), 
          legend.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size),
          legend.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size, face = "bold")) + 
    guides(col=guide_legend(nrow=length(colors_ttest_adj)))

  p3 = ggplot(sub_plot_df, aes(x=logfc, y=wilcox_raw_pval_log10, color=wilcox_raw_cat)) + 
    geom_point(size = 0.75) + 
    geom_hline(yintercept = -log10(sig_lvl), color="lightgrey", linetype="dashed", size = 0.25) +
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    geom_text_repel(aes(label = wilcox_raw_cat_show),
                    max.time = repel_max_time, max.iter = 1e6,
                    max.overlaps = Inf,
                    box.padding = 0.1, 
                    point.patting = 0.1, 
                    force_pull = 2, 
                    force = 100, 
                    min.segment.length = 0, segment.size = 0.25,
                    seed = snakemake@params$parameters_porps$set_seed,
                    size = 6/11.04*3.88, family = snakemake@params$plots_props$font_family, color = "black") +
    xlab("Fold change (log2)") + ylab("Raw p-value\nfrom Wilcoxon test(-log10)") +
    scale_color_manual(name="", values=colors_wilcox_raw) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          plot.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size_header), 
          legend.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size),
          legend.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size, face = "bold")) + 
    guides(col=guide_legend(nrow=length(colors_wilcox_raw)))

  p4 = ggplot(sub_plot_df, aes(x=logfc, y=wilcox_adj_pval_log10, color=wilcox_adj_cat)) + 
    geom_point(size = 0.75) + 
    geom_hline(yintercept = -log10(sig_lvl), color="lightgrey", linetype="dashed", size = 0.25) +
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    geom_text_repel(aes(label = wilcox_adj_cat_show),
                    max.time = repel_max_time, max.iter = 1e6,
                    max.overlaps = Inf,
                    box.padding = 0.1, 
                    point.patting = 0.1, 
                    force_pull = 2, 
                    force = 100, 
                    min.segment.length = 0, segment.size = 0.25,
                    seed = snakemake@params$parameters_porps$set_seed,
                    size = 6/11.04*3.88, family = snakemake@params$plots_props$font_family, color = "black") +
    xlab("Fold change (log2)") + ylab("Adj. p-value\nfrom Wilcoxon test (-log10)") + 
    scale_color_manual(name="", values=colors_wilcox_adj) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          plot.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size_header), 
          legend.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size),
          legend.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size, face = "bold")) + 
    guides(col=guide_legend(nrow=length(colors_wilcox_adj)))

  p5 = ggplot(sub_plot_df, aes(x=logfc, y=cohend_abs, color=cohend_abs_cat)) + 
    geom_point(size = 0.75) + 
    geom_hline(yintercept = abs(effect_size_thres), color="lightgrey", linetype="dashed", size = 0.25) +
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    geom_text_repel(aes(label = cohend_abs_cat_show),
                    max.time = repel_max_time, max.iter = 1e6,
                    max.overlaps = Inf, 
                    box.padding = 0.1, 
                    point.patting = 0.1, 
                    force_pull = 2, 
                    force = 100,
                    min.segment.length = 0, segment.size = 0.25,
                    seed = snakemake@params$parameters_porps$set_seed,
                    size = 6/11.04*3.88, family = snakemake@params$plots_props$font_family, color = "black") +
    xlab("Fold change (log2)") + ylab("Effect size (Cohen's d)") +
    scale_color_manual(name="", values=colors_cohend_estimate) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          plot.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size_header), 
          legend.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size),
          legend.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size, face = "bold")) + 
    guides(col=guide_legend(nrow=length(colors_cohend_estimate)))
  
  p6 = ggplot(sub_plot_df, aes(x=logfc, y=ttest_raw_pval_log10, color=ttest_raw_cat)) + 
    geom_point(size = 0.75) +
    geom_hline(yintercept = -log10(sig_lvl), color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    xlab("Fold change (log2)") + ylab("Raw p-value\nfrom t-test (-log10)") +
    scale_color_manual(name="", values=colors_ttest_raw) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          plot.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size_header), 
          legend.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size),
          legend.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size, face = "bold")) + 
    guides(col=guide_legend(nrow=length(colors_ttest_raw)))
  
  p7 = ggplot(sub_plot_df, aes(x=logfc, y=ttest_adj_pval_log10, color=ttest_adj_cat)) + 
    geom_point(size = 0.75) +
    geom_hline(yintercept = -log10(sig_lvl), color="lightgrey", linetype="dashed", size = 0.25) +
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    xlab("Fold change (log2)") + ylab("Adj. p-value\nfrom t-test (-log10)") +
    scale_color_manual(name="", values=colors_ttest_adj) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          plot.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size_header), 
          legend.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size),
          legend.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size, face = "bold")) + 
    guides(col=guide_legend(nrow=length(colors_ttest_adj)))
  
  p8 = ggplot(sub_plot_df, aes(x=logfc, y=wilcox_raw_pval_log10, color=wilcox_raw_cat)) + 
    geom_point(size = 0.75) + 
    geom_hline(yintercept = -log10(sig_lvl), color="lightgrey", linetype="dashed", size = 0.25) +
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    xlab("Fold change (log2)") + ylab("Raw p-value\nfrom Wilcoxon test(-log10)") +
    scale_color_manual(name="", values=colors_wilcox_raw) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          plot.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size_header), 
          legend.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size),
          legend.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size, face = "bold")) + 
    guides(col=guide_legend(nrow=length(colors_wilcox_raw)))
  
  p9 = ggplot(sub_plot_df, aes(x=logfc, y=wilcox_adj_pval_log10, color=wilcox_adj_cat)) + 
    geom_point(size = 0.75) + 
    geom_hline(yintercept = -log10(sig_lvl), color="lightgrey", linetype="dashed", size = 0.25) +
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    xlab("Fold change (log2)") + ylab("Adj. p-value\nfrom Wilcoxon test (-log10)") + 
    scale_color_manual(name="", values=colors_wilcox_adj) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          plot.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size_header), 
          legend.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size),
          legend.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size, face = "bold")) + 
    guides(col=guide_legend(nrow=length(colors_wilcox_adj)))
  
  p10 = ggplot(sub_plot_df, aes(x=logfc, y=cohend_abs, color=cohend_abs_cat)) + 
    geom_point(size = 0.75) + 
    geom_hline(yintercept = abs(effect_size_thres), color="lightgrey", linetype="dashed", size = 0.25) +
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    xlab("Fold change (log2)") + ylab("Effect size (Cohen's d)") +
    scale_color_manual(name="", values=colors_cohend_estimate) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          axis.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size), 
          plot.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size_header), 
          legend.text = element_text(family = snakemake@params$plots_props$font_family, size = font_size),
          legend.title = element_text(family = snakemake@params$plots_props$font_family, size = font_size, face = "bold")) + 
    guides(col=guide_legend(nrow=length(colors_cohend_estimate)))
  
  if(is.null(group)) {
    # 4 x 4
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.ttest.raw.labels.png", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p1, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi, bg = "transparent")
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.ttest.raw.labels.svg", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p1, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi, bg = "transparent")
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.ttest.adj.labels.png", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p2, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi, bg = "transparent")
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.ttest.adj.labels.svg", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p2, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi, bg = "transparent")
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.wilcox.raw.labels.png", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p3, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi, bg = "transparent")
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.wilcox.raw.labels.svg", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p3, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi, bg = "transparent")
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.wilcox.adj.labels.png", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p4, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi, bg = "transparent")
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.wilcox.adj.labels.svg", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p4, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi, bg = "transparent")
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.cohend_estimate.labels.png", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p5, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi, bg = "transparent")
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.cohend_estimate.labels.svg", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p5, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi, bg = "transparent")
  
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.ttest.raw.png", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p6, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi)
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.ttest.raw.svg", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p6, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi)
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.ttest.adj.png", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p7, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi)
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.ttest.adj.svg", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p7, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi)
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.wilcox.raw.png", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p8, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi)
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.wilcox.raw.svg", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p8, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi)
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.wilcox.adj.png", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p9, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi)
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.wilcox.adj.svg", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p9, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi)
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.cohend_estimate.png", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p10, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi)
    ggsave(file = sprintf("%s/%s__%s_vs_%s%s%s.cohend_estimate.svg", output_folder, prop, c$g1, c$g2, paired_string, sub_string), plot = p10, width = 4, height = 4, unit = plots_props$image_units, dpi = plots_props$dpi)
  } else {
      group_plots[[group]] <- append(group_plots[[group]], list(p1, p2, p3, p4, p5))
  }
}

for(i in 1:length(comps)){
  c = comps[[i]]
  if("volcano" %in% names(c)) {
    group <- c$volcano
    prop = names(comps)[i]
    paired_string = ifelse(!is.null(c$paired), sprintf("_paired_%s", c$paired), "")
    sub_string = ifelse(!is.null(c$subset), sprintf("_%s=%s", c$subset[1], c$subset[2]), "")
    p1 <- group_plots[[group]][[1]] + ylim(NA, group_y_max[[group]][["ttest_raw_pval"]]) + xlim(group_x_min[[group]][["logfc"]], group_x_max[[group]][["logfc"]])
    p2 <- group_plots[[group]][[2]] + ylim(NA, group_y_max[[group]][["ttest_adj_pval"]]) + xlim(group_x_min[[group]][["logfc"]], group_x_max[[group]][["logfc"]])
    p3 <- group_plots[[group]][[3]] + ylim(NA, group_y_max[[group]][["wilcox_raw_pval"]]) + xlim(group_x_min[[group]][["logfc"]], group_x_max[[group]][["logfc"]])
    p4 <- group_plots[[group]][[4]] + ylim(NA, group_y_max[[group]][["wilcox_adj_pval"]]) + xlim(group_x_min[[group]][["logfc"]], group_x_max[[group]][["logfc"]])
    p5 <- group_plots[[group]][[5]] + ylim(NA, group_y_max[[group]][["cohend_abs"]]) + xlim(group_x_min[[group]][["logfc"]], group_x_max[[group]][["logfc"]])
    save_plot(file.path(snakemake@output[[1]], sprintf("%s__%s_vs_%s%s%s.ttest.raw.pdf", prop, c$g1, c$g2, paired_string, sub_string)), p1, base_aspect_ratio = 0.9, base_height=4.8, device=cairo_pdf)
    save_plot(file.path(snakemake@output[[1]], sprintf("%s__%s_vs_%s%s%s.ttest.adj.pdf", prop, c$g1, c$g2, paired_string, sub_string)), p2, base_aspect_ratio = 0.9, base_height=4.8, device=cairo_pdf)
    save_plot(file.path(snakemake@output[[1]], sprintf("%s__%s_vs_%s%s%s.wilcox.raw.pdf", prop, c$g1, c$g2, paired_string, sub_string)), p3, base_aspect_ratio = 0.9, base_height=4.8, device=cairo_pdf)
    save_plot(file.path(snakemake@output[[1]], sprintf("%s__%s_vs_%s%s%s.wilcox.adj.pdf", prop, c$g1, c$g2, paired_string, sub_string)), p4, base_aspect_ratio = 0.9, base_height=4.8, device=cairo_pdf)
    save_plot(file.path(snakemake@output[[1]], sprintf("%s__%s_vs_%s%s%s.cohend_estimate.pdf", prop, c$g1, c$g2, paired_string, sub_string)), p5, base_aspect_ratio = 0.9, base_height=4.8, device=cairo_pdf)
    group_plots[[group]] <- group_plots[[group]][-(1:5)]
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

