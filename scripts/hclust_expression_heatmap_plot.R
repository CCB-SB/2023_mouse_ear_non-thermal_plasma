suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_porps$set_seed)

print("clustering")

# load data
tbl = fread(snakemake@input$expr_log10, sep='\t', header=T)
rna_ids = tbl[[colnames(tbl)[1]]]

annot = fread(snakemake@input$annot, sep='\t', colClasses=c(ID="character"))

prop = snakemake@params$prop
additional_props = snakemake@params$additional_props
params = snakemake@params$params
colors = snakemake@params$colors
plot_props = snakemake@params$plot_props
top_list = snakemake@params$top_list


results_folder = snakemake@params$results_folder
prefix = sprintf("%s/heatmap_complex/expr_hclust_complete", results_folder)
dir.create(prefix, recursive=TRUE)

scale = function (x, rows, columns) {
  if(rows){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    x = (x - m)/s
  }
  if(columns){
    x = t(x)
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    x = (x - m)/s
    x = t(x)
  }
  return(x)
}


cv = function(mat){
  apply(mat, 1, function(x) sd(x, na.rm=T)) / rowMeans(mat, na.rm=T)
}


plot_heatmap = function(df, row_names, plot_props, width, height, filename){
  df = as.matrix(df)
  rownames(df) = row_names
  
  if(params$scale != "none"){
    scale_row = params$scale %in% c("row", "both")
    scale_col = params$scale %in% c("column", "both")
    df = scale(df, scale_row, scale_col)
    df[is.na(df)] = 0
  }
  
  annotation_col = data.frame(dummy=1:ncol(df))
  annotation_col[[prop]] = as.character(annot[[prop]][match(colnames(df), annot$ID)])
  rownames(annotation_col) = colnames(df)
  annotation_col$dummy = NULL
  
  annotation_colors = list()
  ucolors <- unlist(colors[[prop]])
  if(is.vector(ucolors) && all(annotation_col[[prop]] %in% names(ucolors))) {
    annotation_colors[[prop]] = ucolors
  }

  if(!is.null(additional_props)){
    if(is.vector(additional_props)){
      for(p in additional_props){
        annotation_col[[p]] = as.character(annot[[p]][match(colnames(df), annot$ID)])
        ucolors <- unlist(colors[[p]])
        if(is.vector(ucolors) && all(annotation_col[[p]] %in% names(ucolors))) {
          annotation_colors[[p]] = ucolors
        }
      }
    } else {
      annotation_col[[additional_props]] = as.character(annot[[additional_props]][match(colnames(df), annot$ID)])
      annotation_colors[[additional_props]] = unlist(colors[[additional_props]])
    }
  }
  
  show_rownames = nrow(df) < 15

  column_ha = HeatmapAnnotation(df = annotation_col, col = annotation_colors, gp = gpar(col = "white"), simple_anno_size = unit(0.3, "cm"), 
                                annotation_name_side = "left", annotation_name_gp = gpar(fontsize = plot_props$font_size, fontface = "bold", fontfamily=plot_props$font_family))
  
  p = ComplexHeatmap::Heatmap(df, col = viridis(100),
                              rect_gp = gpar(col = "white", lwd = .05),
                              height = unit(height - 3.75, "cm"), width = unit(width - 1.5, "cm"),
                              row_names_side = "left", row_names_gp = gpar(fontsize = 6, fontfamily=plot_props$font_family), column_names_side = "bottom", column_names_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family),
                              cluster_columns = function(m) hclust(dist(m), method = "complete"), clustering_distance_columns = "euclidean", #show_columns_dend = FALSE
                              row_dend_reorder = FALSE, cluster_rows = FALSE,
                              show_column_names = FALSE,
                              show_row_names = show_rownames,
                              show_heatmap_legend = FALSE,
                              top_annotation = column_ha
  )
  legend_ticks = seq(from = min(t(df)), to = max(t(df)), length.out = 5)
  legend_tick_labels = round(legend_ticks, digits = 2)
  legend = ComplexHeatmap::Legend(title = "Standardised expr. value (rpmmm, log10)", col_fun = circlize::colorRamp2(legend_ticks, viridis(5)), at = legend_tick_labels, direction="horizontal", 
                                  title_gp = gpar(fontsize = plot_props$font_size, fontface = "bold", fontfamily=plot_props$font_family), labels_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family))
  
  png(sprintf("%s.png",filename), width = width, height = height, units = plot_props$image_units, res = plot_props$dpi)
  draw(p, show_annotation_legend = FALSE, padding = unit(c(0.15, 0.25, -1.2, 0.25), "cm"))
  draw(legend, x = unit(1.3, "cm"), y = unit(.05, "cm"), just = c("left", "bottom"))
  dev.off()
  svglite(sprintf("%s.svg",filename), width = width / 2.54, height = height / 2.54)
  draw(p, show_annotation_legend = FALSE, padding = unit(c(0.15, 0.25, -1.2, 0.25), "cm"))
  draw(legend, x = unit(1.3, "cm"), y = unit(.05, "cm"), just = c("left", "bottom"))
  dev.off()
}

for (tp in unique(annot$Time)) {
  group_tp = annot[annot$Time == tp]$ID
  tmp = colnames(tbl) %in% group_tp
  tbl_2 = tbl[,..tmp]

  # keep only expression
  expr = tbl_2[, sapply(tbl_2, is.numeric), with=F]
  # complete 
  plot_heatmap(expr, NULL, plot_props, width = plot_props$image_width*2/3, height = plot_props$image_height, filename = file.path(prefix, sprintf("clustering_%s", tp)))

  # only expressed
  expr.filter = rowSums(expr != 0) > 0
  plot_heatmap(expr[expr.filter], NULL, plot_props, width = plot_props$image_width*2/3, height = plot_props$image_height, filename = file.path(prefix, sprintf("clustering.expressed_%s", tp)))
}

# keep only expression
expr = tbl[, sapply(tbl, is.numeric), with=F]

# complete 
plot_heatmap(expr, NULL, plot_props, width = plot_props$image_width, height = plot_props$image_height, filename = file.path(prefix, "clustering"))

# only expressed
expr.filter = rowSums(expr != 0) > 0
plot_heatmap(expr[expr.filter], NULL, plot_props, width = plot_props$image_width, height = plot_props$image_height, filename = file.path(prefix, "clustering.expressed"))

if("min_expr" %in% names(params)) {
  tbl.raw = fread(snakemake@input$raw, sep='\t')
  stopifnot(assertthat::are_equal(rna_ids, tbl.raw[[colnames(tbl.raw)[1]]]))
  expr.raw = tbl.raw[, sapply(tbl.raw, is.numeric), with=F]
  min_expr.filter = rowSums(expr.raw >= params$min_expr) != 0
  plot_heatmap(expr[min_expr.filter], NULL, plot_props, width = plot_props$image_width, height = plot_props$image_height, filename = file.path(prefix, sprintf("clustering.min_%d", params$min_expr)))
}

if(params$ranking == "cv"){
  expr.filter.cv = abs(cv(expr[expr.filter])) 
  ranked.filter.cv = frank(-expr.filter.cv, ties.method = "first")
  if("min_expr" %in% names(params)) {
    min_expr.filter.cv = abs(cv(expr[min_expr.filter]))
    ranked.min_expr.filter.cv = frank(-min_expr.filter.cv, ties.method = "first")
  }
} else if(params$ranking == "variance"){
  expr.filter.cv = abs(apply(expr, 1, function(x) var(x, na.rm = T)))
  ranked.filter.cv = frank(-expr.filter.cv, ties.method = "first")
  if("min_expr" %in% names(params)) {
    min_expr.filter.cv = abs(apply(expr[min_expr.filter], 1, function(x) var(x, na.rm = T)))
    ranked.min_expr.filter.cv = frank(-min_expr.filter.cv, ties.method = "first")
  }
}

for(i in top_list){
  if(i <= 50) {
    row_names = rna_ids[expr.filter][ranked.filter.cv <= i]
    if("min_expr" %in% names(params)) {
      row_names.min_expr = rna_ids[min_expr.filter][ranked.min_expr.filter.cv <= i]
    }
  } else {
    row_names = NULL
    row_names.min_expr = NULL
  }
  plot_heatmap(expr[expr.filter][ranked.filter.cv <= i], row_names, plot_props, width = plot_props$image_width, height = plot_props$image_height, filename = file.path(prefix, sprintf("clustering.top_%d", i)))
  if("min_expr" %in% names(params)) {
    plot_heatmap(expr[min_expr.filter][ranked.min_expr.filter.cv <= i], row_names.min_expr, plot_props, width = plot_props$image_width, height = plot_props$image_height, filename = file.path(prefix, sprintf("clustering.min_%d.top_%d", params$min_expr, i)))
  }

  for(tp in unique(annot$Time)){
    group_tp = annot[annot$Time == tp,]$ID
    tmp = colnames(expr) %in% group_tp
    expr_2 = expr[,..tmp]
    plot_heatmap(expr_2[expr.filter][ranked.filter.cv <= i], row_names, plot_props=plot_props, width = plot_props$image_width*2/3, height = plot_props$image_height, filename=file.path(prefix, sprintf("clustering.top_%d_Time_%s", i, tp)))
    if("min_expr" %in% names(params)) {
      plot_heatmap(expr_2[min_expr.filter][ranked.min_expr.filter.cv <= i], row_names.min_expr, plot_props, width=plot_props$image_width*2/3, height=plot_props$image_height, filename = file.path(prefix, sprintf("clustering.min_%d.top_%d_Time_%s", params$min_expr, i, tp)))
  }
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
