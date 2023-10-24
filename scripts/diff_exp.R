suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(effsize))
suppressPackageStartupMessages(library(pROC))

setDTthreads(1)

log_list = snakemake@params$log_list

for (logs in log_list) {
  logs_value = gsub("_log", "", logs)
  # Read data
  expr = fread(snakemake@input[[sprintf("expr%s", logs)]], sep='\t', header=T)
  annot = fread(snakemake@input$annot, sep='\t', colClasses=c(ID="character"))
  # Get ids
  rna_ids = expr[[colnames(expr)[1]]]
  # Remove unmapped samples
  expr[, ( setdiff(names(expr), annot$ID) ) := NULL] 

  # drop unexpressed 
  to_keep = rowSums(expr != 0) > 0
  rna_ids = rna_ids[to_keep]
  expr = expr[to_keep]
  # Parse parameters
  comps_fc_versus_expr = snakemake@params$comparisons_fc_versus_expr
  comps_volcano = snakemake@params$comparisons_volcano
  comps_heatmap_fc = snakemake@params$comparisons_heatmap_fc
  comps_comparisons_supp_table = snakemake@params$comparisons_supp_table
  comps = c(comps_volcano, comps_fc_versus_expr, comps_heatmap_fc, comps_comparisons_supp_table)

  adjustment = snakemake@params$adjustment

  shapiro_rawp = apply(expr, 1, function(x) {
    if(length(unique(x)) == 1){
      return(NA)
    }
    return(shapiro.test(x)$p.value)
  })
  shapiro_adjp = p.adjust(shapiro_rawp, method=adjustment)

  result = data.frame(RNA=rna_ids, shapiro_rawp, shapiro_adjp)

  get_ids <- function(mat, annot, c, prop) {
    if(!is.null(c$subset)){
      sel1 = annot[[prop]] == c$g1 & annot[[c$subset[1]]] == c$subset[2]
      sel2 = annot[[prop]] == c$g2 & annot[[c$subset[1]]] == c$subset[2]
    } else {
      sel1 = annot[[prop]] == c$g1
      sel2 = annot[[prop]] == c$g2
    }
    
    if(!is.null(c$paired)){
      pair_prop_1 = annot[sel1][[c$paired]]
      pair_prop_2 = annot[sel2][[c$paired]]
      # keep only ids that are in both, and arrange them in the same order (i.e. order 2 same as 1)
      g1_ids = annot[sel1]$ID[!is.na(match(pair_prop_1, pair_prop_2))]
      g2_ids = annot[sel2]$ID[na.omit(match(pair_prop_1, pair_prop_2))]
    } else {
      g1_ids = annot[sel1]$ID
      g2_ids = annot[sel2]$ID
    }
    
    g1_ids = g1_ids[g1_ids %in% colnames(mat)]
    g2_ids = g2_ids[g2_ids %in% colnames(mat)]

    return(list(g1_ids, g2_ids))
}

  compute_measures <- function(x, g1_ids, g2_ids, c, prop, paired_string, sub_string) {
    sub1 = x[g1_ids]
    sub2 = x[g2_ids]

    result = list()
    result[[sprintf("mean_%s__%s%s%s", prop, c$g1, paired_string, sub_string)]] = mean(sub1)
    result[[sprintf("mean_%s__%s%s%s", prop, c$g2, paired_string, sub_string)]] = mean(sub2)
    result[[sprintf("median_%s__%s%s%s", prop, c$g1, paired_string, sub_string)]] = median(sub1)
    result[[sprintf("median_%s__%s%s%s", prop, c$g2, paired_string, sub_string)]] = median(sub2)
    result[[sprintf("num_%s__%s%s%s", prop, c$g1, paired_string, sub_string)]] = length(sub1)
    result[[sprintf("num_%s__%s%s%s", prop, c$g2, paired_string, sub_string)]] = length(sub2)
    return(result)
  }

  compute_other = function(x, g1_ids, g2_ids, test_func, getter, ...) {
    sub1 = x[g1_ids]
    sub2 = x[g2_ids]

    return(getter(test_func(sub1, sub2, ...)))
  }

  compute_auc <- function(sub1, sub2, ...){
    return(auc(c(rep.int(2, length(sub1)), rep.int(1, length(sub2))), c(sub1, sub2), quiet=T, direction="<"))
  }

  compute_p_value = function(x, g1_ids, g2_ids, test_func, paired, ...) {
    sub1 = x[g1_ids]
    sub2 = x[g2_ids]

    if(length(unique(sub1)) == 1 ||
        length(unique(sub2)) == 1){
      p_val = NA
    } else {
      p_val = test_func(sub1, sub2, paired=!is.null(paired), ...)$p.value
    }
    return(p_val)
  }

  for(i in 1:length(comps)){
    c = comps[[i]]
    # this is only needed if comparisons are loaded from the volcano part of config.yaml
    c$g1 = c$g[1]
    c$g2 = c$g[2]

    prop = names(comps)[i]

    paired_string_measures = ifelse(!is.null(c$paired), sprintf("_paired_%s_%s_vs_%s", c$paired, c$g1, c$g2), "")
    paired_string = ifelse(!is.null(c$paired), sprintf("_paired_%s", c$paired), "")
    sub_string = ifelse(!is.null(c$subset), sprintf("_%s=%s", c$subset[1], c$subset[2]), "")

    ids <- get_ids(expr, annot, c, prop)
    g1_ids <- ids[[1]]
    g2_ids <- ids[[2]]

    res = apply(expr, 1, compute_measures, g1_ids, g2_ids, c, prop, paired_string_measures, sub_string)
    res_df = data.table(do.call(rbind, res))
    res_df <- sapply(res_df, unlist)
    inter <- intersect(colnames(result), colnames(res_df))
    if(length(inter) > 0) {
      res_df <- res_df[,!(colnames(res_df) %in% inter)]
    }
    result <- cbind(result, res_df)

   if(as.numeric(logs_value) != 0){
      result[[sprintf("fc_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]] = as.numeric(logs_value)^(result[[sprintf("median_%s__%s%s%s", prop, c$g1, paired_string_measures, sub_string)]] - result[[sprintf("median_%s__%s%s%s", prop, c$g2, paired_string_measures, sub_string)]])
    } else {
      result[[sprintf("fc_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]] = result[[sprintf("median_%s__%s%s%s", prop, c$g1, paired_string_measures, sub_string)]] / result[[sprintf("median_%s__%s%s%s", prop, c$g2, paired_string_measures, sub_string)]]
    }

    ttest = apply(expr, 1, compute_p_value, g1_ids, g2_ids, t.test, c$paired)
    result[, (sprintf("%s_rawp_%s__%s_vs_%s%s%s", "ttest", prop, c$g1, c$g2, paired_string, sub_string))] = ttest
    result[, (sprintf("%s_adjp_%s__%s_vs_%s%s%s", "ttest", prop, c$g1, c$g2, paired_string, sub_string))] = p.adjust(ttest, method=adjustment)

    wilcox = apply(expr, 1, compute_p_value, g1_ids, g2_ids, wilcox.test, c$paired)
    result[, (sprintf("%s_rawp_%s__%s_vs_%s%s%s", "wilcox", prop, c$g1, c$g2, paired_string, sub_string))] = wilcox
    result[, (sprintf("%s_adjp_%s__%s_vs_%s%s%s", "wilcox", prop, c$g1, c$g2, paired_string, sub_string))] = p.adjust(wilcox, method=adjustment)

    cohend <- apply(expr, 1, compute_other, g1_ids, g2_ids, cohen.d, function(x) {return(x$estimate)}, paired=!is.null(c$paired))
    result[, (sprintf("%s_estimate_%s__%s_vs_%s%s%s", "cohend", prop, c$g1, c$g2, paired_string, sub_string))] = cohend

    area_u_c <- apply(expr, 1, compute_other, g1_ids, g2_ids, compute_auc, function(x) {return(as.numeric(x))})
    result[, (sprintf("%s_value_%s__%s_vs_%s%s%s", "auc", prop, c$g1, c$g2, paired_string, sub_string))] = area_u_c
  }

  rownames(result) = NULL

  result = result[order(result[,12]),]

  results_folder = snakemake@params$results_folder
  output_folder = sprintf("%s/matrices/diff_exp", results_folder)
  dir.create(output_folder, recursive=TRUE)
  write.table(result, sprintf("%s/diff_exp%s.csv", output_folder, logs), sep='\t', row.names = F)
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])