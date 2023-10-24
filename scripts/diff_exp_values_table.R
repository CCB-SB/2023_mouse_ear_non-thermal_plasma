suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(openxlsx))

set.seed(snakemake@params$parameters_props$set_seed)

print("supplementary plot vulcano")


# ------------------------------- Input data -----------------------------------
diff_exp_logs = snakemake@params$diff_exp_logs
comps = snakemake@params$comparisons
results_folder = snakemake@params$results_folder

# load data
diff_exp = fread(sprintf("%s/matrices/diff_exp/diff_exp%s.csv", results_folder, diff_exp_logs), sep='\t', header=T)
annot = fread(snakemake@input$annot, sep='\t', colClasses=c(ID="character"))
interesting_miRNA = fread(snakemake@input$miRNA_list)

# save folder
output_folder = sprintf("%s/matrices/supp_volcano/", results_folder)
dir.create(output_folder, recursive=TRUE)


# ----------------------------------- 4 ears -----------------------------------
# ------------------------- plot fc versus expr median -------------------------
  
plot_df_all = c()
selected_colnames = c("RNA")
for(i in 1:length(comps)) {
  # which comp do we investigate and which groups are considered
  prop = names(comps)[i]
  c = comps[[i]]
  c$g1 = c$g[1]
  c$g2 = c$g[2]
  
  prop = names(comps)[i]
  
  paired_string = ifelse(!is.null(c$paired), sprintf("_paired_%s", c$paired), "")
  sub_string = ifelse(!is.null(c$subset), sprintf("_%s=%s", c$subset[1], c$subset[2]), "")
  
  ttest_adj_pval = sprintf("ttest_adjp_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)
  fc = sprintf("fc_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)
  auc = sprintf("auc_value_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)
  cohensd = sprintf("cohend_estimate_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)
  
  plot_df = data.table(RNA=diff_exp[["RNA"]],
                       ttest_adj_pval=diff_exp[[ttest_adj_pval]],
                       fc=diff_exp[[fc]],
                       auc=diff_exp[[auc]],
                       cohensd=diff_exp[[cohensd]]
                       )
  
  selected_colnames = append(selected_colnames, c(ttest_adj_pval, fc, auc, cohensd))
  
  fwrite(plot_df, sprintf("%s/supp_table_%s__%s_vs_%s%s%s.csv", output_folder, prop, c$g1, c$g2, paired_string, sub_string), sep = "\t", row.names = FALSE)
  write.xlsx(plot_df, sprintf("%s/supp_table__%s__%s_vs_%s%s%s.xlsx", output_folder, prop, c$g1, c$g2, paired_string, sub_string), colNames = TRUE, rowNames = FALSE, append = FALSE)
}

# save all
fwrite(diff_exp[, ..selected_colnames], sprintf("%s/supp_table_all.csv", output_folder), sep = "\t", row.names = FALSE)
write.xlsx(diff_exp[, ..selected_colnames], sprintf("%s/supp_table_all.xlsx", output_folder), colNames = TRUE, rowNames = FALSE, append = FALSE)


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

  
