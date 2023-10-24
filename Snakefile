from os.path import join

from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.4.0")


##### load config #####

configfile: "config.yaml"

result_files = [join(config["results_folder"], "steps/{}.done".format(m)) for m in config['modules'] if m != "common"]

##### target rules #####
rule all:
    input: result_files


##### setup singularity #####
rule rpmmm_norm:
    input: "data/expression/raw/GSE236788_miRNA_expression_raw.tsv"
    output:
        result=join(config["results_folder"], "steps/rpmmm_norm.done"),
        expr_rpmm="data/expression/norm/miRNA_quantification_rpmmm_norm.csv"
    run:
        import pandas as pd
        from pathlib import Path
        tbl = pd.read_csv(input[0], sep='\t')
        tbl.iloc[:,1:] = tbl.iloc[:,1:] / tbl.groupby("miRNA").mean().sum() * 1e6
        folder = Path("data/expression/norm")
        folder.mkdir(parents=True, exist_ok=True)
        tbl.to_csv("data/expression/norm/miRNA_quantification_rpmmm_norm.csv", sep='\t', index=False)
        with open(output[0], "w") as file:
            file.write("")

rule log2_transform:
    input:
        expr="data/expression/norm/miRNA_quantification_rpmmm_norm.csv",
        result=join(config["results_folder"], "steps/rpmmm_norm.done")
    output:
        result=join(config["results_folder"], "steps/log2_transform.done"),
        expr_rpmmm_log2="data/expression/norm/miRNA_quantification_rpmmm_norm_log2.csv"
    params: 
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/log2_transform.R"

rule log10_transform:
    input:
        expr="data/expression/norm/miRNA_quantification_rpmmm_norm.csv",
        result=join(config["results_folder"], "steps/rpmmm_norm.done")
    output:
        result=join(config["results_folder"], "steps/log10_transform.done"),
        expr_rpmmm_log10="data/expression/norm/miRNA_quantification_rpmmm_norm_log10.csv"
    params: 
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/log10_transform.R"

rule create_raw_detection_matrix:
    input: "data/expression/raw/GSE236788_miRNA_expression_raw.tsv"
    output:
        result=join(config["results_folder"], "steps/create_raw_detection_matrix.done"),
        detect="data/expression/raw/raw_detection_miRNA.csv",
    params: min_detection=5
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/create_raw_detection_matrix.R"

rule filter_expression_files_per_group:
    input:
        result=join(config["results_folder"], "steps/rpmmm_norm.done"),
        result_log2=join(config["results_folder"], "steps/log2_transform.done"),
        result_log10=join(config["results_folder"], "steps/log10_transform.done"),
        expr="data/expression/norm/miRNA_quantification_rpmmm_norm.csv",
        expr_log2="data/expression/norm/miRNA_quantification_rpmmm_norm_log2.csv",
        expr_log10="data/expression/norm/miRNA_quantification_rpmmm_norm_log10.csv",
        detect="data/expression/raw/raw_detection_miRNA.csv",
        annot=config["annot"],
    output: 
        result=join(config["results_folder"], "steps/filter_expression_files_per_group.done"),
        expr_log2=config["expr_log2"],
        expr_log10=config["expr_log10"],
        expr=config["expr"]
    params:
        detection_rate=(100 / 100),
        filter_variable="Group",
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/filter_expression_files_per_group.R"

rule alignment_stats_bar_plot:
    input:
        annot=config["annot"],
        mapping_info=config["mapping_info"],
    output: result=join(config["results_folder"], "steps/alignment_stats_bar_plot.done")
    params: 
        plots_props=config["common"].get("plots", None),
        parameters_props=config["common"].get("parameters", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/alignment_stats_bar_plot.R"

rule mapping_stats_bar_plot:
    input:
        annot=config["annot"],
        cleaning_stats=config["cleaning_stats"],
        rna_composition_mapped=config["rna_composition_mapped"],
        rna_composition_total=config["rna_composition_total"],
    output: result=join(config["results_folder"], "steps/mapping_stats_bar_plot.done")
    params: 
        plots_props=config["common"].get("plots", None),
        parameters_props=config["common"].get("parameters", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/mapping_stats_bar_plot.R"

###### 
rule diff_exp:
    input:
        expr_log2=config["expr_log2"],
        expr_log10=config["expr_log10"],
        annot=config["annot"]
    output: result=join(config["results_folder"], "steps/diff_exp.done")
    params: 
        comparisons_fc_versus_expr=config["fc_versus_expr_scatter_plot"].get("comparisons", None),
        comparisons_volcano=config["volcano_scatter_plot"].get("comparisons", None),
        comparisons_heatmap_fc=config["comp_diff_exp_values_heatmap_plot"].get("comparisons", None),
        comparisons_supp_table=config["diff_exp_values_table"].get("comparisons", None),
        log_list=config["diff_exp"].get("log_list", None),
        adjustment=config["diff_exp"]["adjustment"],
        results_folder=config["results_folder"]
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/diff_exp.R"

rule pca_scatter_plot:
    input: 
        expr=config["expr"],
        annot=config["annot"],
    output: result=join(config["results_folder"], "steps/pca_scatter_plot.done")
    params: 
        properties=config["pca_scatter_plot"].get("properties", None),
        plots_props=config["common"].get("plots", None),
        parameters_props=config["common"].get("parameters", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/pca_scatter_plot.R"

rule hclust_expression_heatmap_plot:
    input: 
        expr_log10=config["expr_log10"],
        annot=config["annot"],
    output: result=join(config["results_folder"], "steps/hclust_expression_heatmap_plot.done")
    params: 
        prop=config["hclust_expression_heatmap_plot"].get("prop", None),
        additional_props=config["hclust_expression_heatmap_plot"].get("additional_props", None),
        params=config["hclust_expression_heatmap_plot"].get("hclust_params", None),
        top_list=config["hclust_expression_heatmap_plot"].get("top_list", None),
        plot_props=config["common"].get("plots", None),
        parameters_props=config["common"].get("parameters", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/hclust_expression_heatmap_plot.R"

rule fc_versus_expr_scatter_plot:
    input: 
        result=join(config["results_folder"], "steps/diff_exp.done"),
        annot=config["annot"]
    output: result=join(config["results_folder"], "steps/fc_versus_expr_scatter_plot.done")
    params: 
        diff_exp_logs=config["fc_versus_expr_scatter_plot"].get("diff_exp_logs", None),
        comparisons=config["fc_versus_expr_scatter_plot"].get("comparisons", None),
        test_values_th=config["fc_versus_expr_scatter_plot"].get("thresholds", None),
        plots_props=config["common"].get("plots", None),
        parameters_props=config["common"].get("parameters", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/fc_versus_expr_scatter_plot.R"

rule pvca_bar_plot:
    input: 
        expr_log2=config["expr_log2"],
        expr_log10=config["expr_log10"],
        annot=config["annot"],
    output: result=join(config["results_folder"], "steps/pvca_bar_plot.done")
    params: 
        props=config["pvca_bar_plot"].get("properties", None),
        parameters=config["pvca_bar_plot"].get("parameters", None),
        plots_props=config["common"].get("plots", None),
        parameters_props=config["common"].get("parameters", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/pvca_bar_plot.R"

rule volcano_scatter_plot:
    input: result=join(config["results_folder"], "steps/diff_exp.done"),
    output: result=join(config["results_folder"], "steps/volcano_scatter_plot.done")
    params: 
        diff_exp_logs=config["fc_versus_expr_scatter_plot"].get("diff_exp_logs", None),
        comparisons=config["volcano_scatter_plot"].get("comparisons", None),
        parameters=config["volcano_scatter_plot"].get("parameters", None),
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/volcano_scatter_plot.R"

rule mirna_expression_median_violin_plot:
    input: 
        expr=config["expr_log10"],
        annot=config["annot"],
        miRNA_list=config["miRNA_list"]
    output: result=join(config["results_folder"], "steps/mirna_expression_median_violin_plot.done")
    params: 
        comparisons=config["mirna_expression_median_violin_plot"].get("comparisons", None),
        plots_props=config["common"].get("plots", None),
        parameters_props=config["common"].get("parameters", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/mirna_expression_median_violin_plot.R"

rule comp_diff_exp_values_heatmap_plot:
    input: 
        result=join(config["results_folder"], "steps/diff_exp.done"),
        miRNA_list=config["miRNA_list"]
    output: result=join(config["results_folder"], "steps/comp_diff_exp_values_heatmap_plot.done")
    params: 
        diff_exp_logs=config["fc_versus_expr_scatter_plot"].get("diff_exp_logs", None),
        comparisons=config["comp_diff_exp_values_heatmap_plot"].get("comparisons", None),
        test_values=config["comp_diff_exp_values_heatmap_plot"].get("test_values", None),
        plots_props=config["common"].get("plots", None),
        parameters_props=config["common"].get("parameters", None),
        results_folder=config["results_folder"]
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/comp_diff_exp_values_heatmap_plot.R"

rule correlation_comp_time_heatmap_plot_table:
    input: 
        expr=config["expr"],
        annot=config["annot"],
        miRNA_list=config["miRNA_list"]
    output: result=join(config["results_folder"], "steps/correlation_comp_time_heatmap_plot_table.done")
    params:
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/correlation_comp_time_heatmap_plot_table.R"

rule mirna_expression_over_time_box_plot:
    input: 
        expr_log10=config["expr_log10"],
        annot=config["annot"],
        miRNA_list=config["miRNA_list"]
    output: result=join(config["results_folder"], "steps/mirna_expression_over_time_box_plot.done")
    params:
        polynom_degree=config["mirna_expression_over_time_box_plot"].get("polynom_degree", None),
        plots_props=config["common"].get("plots", None),
        parameters_props=config["common"].get("parameters", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/mirna_expression_over_time_box_plot.R"

rule pathway_analysis_results_heatmap_plots:
    input: miRNA_list=config["miRNA_list"]
    output: result=join(config["results_folder"], "steps/pathway_analysis_results_heatmap_plots.done")
    params:
        filter_params=config["pathway_analysis_results_heatmap_plots"].get("filter_params", None),
        plot_props=config["common"].get("plots", None),
        parameters_props=config["common"].get("parameters", None),
        results_folder=config["results_folder"]
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/pathway_analysis_results_heatmap_plots.R"

rule halushka_violin_plots:
    input: 
    output: result=join(config["results_folder"], "steps/halushka_violin_plots.done")
    params: 
        plots_props=config["common"].get("plots", None),
        parameters_props=config["common"].get("parameters", None),
        results_folder=config["results_folder"]
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/halushka_violin_plots.R"

rule diff_exp_values_table:
    input: 
        result=join(config["results_folder"], "steps/diff_exp.done"),
        annot=config["annot"],
        miRNA_list=config["miRNA_list"]
    output: result=join(config["results_folder"], "steps/diff_exp_values_table.done")
    params: 
        diff_exp_logs=config["diff_exp_values_table"].get("diff_exp_logs", None),
        comparisons=config["diff_exp_values_table"].get("comparisons", None),
        results_folder=config["results_folder"]
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/diff_exp_values_table.R"

rule make_publication_folder:
    input:
        join(config["results_folder"], "steps/alignment_stats_bar_plot.done"),
        join(config["results_folder"], "steps/mapping_stats_bar_plot.done"),
        join(config["results_folder"], "steps/pca_scatter_plot.done"),
        join(config["results_folder"], "steps/hclust_expression_heatmap_plot.done"),
        join(config["results_folder"], "steps/fc_versus_expr_scatter_plot.done"),
        join(config["results_folder"], "steps/pvca_bar_plot.done"),
        join(config["results_folder"], "steps/volcano_scatter_plot.done"),
        join(config["results_folder"], "steps/mirna_expression_median_violin_plot.done"),
        join(config["results_folder"], "steps/comp_diff_exp_values_heatmap_plot.done"),
        join(config["results_folder"], "steps/correlation_comp_time_heatmap_plot_table.done"),
        join(config["results_folder"], "steps/mirna_expression_over_time_box_plot.done"),
        join(config["results_folder"], "steps/pathway_analysis_results_heatmap_plots.done"),
        join(config["results_folder"], "steps/halushka_violin_plots.done"),
        join(config["results_folder"], "steps/diff_exp_values_table.done")
    output: result=join(config["results_folder"], "steps/make_publication_folder.done")
    params:
    conda: "envs/mouse_cold_p.yml"
    script: "scripts/extract_paper_figures.py"

