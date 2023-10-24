#!/usr/bin/env python3

"""
Extract the relevant figures for the publication.
In addition, rename them to the appropriate figure name.
Finally, remove the source figure folder.

Usage:
Run this script from the folder where the "results" folder is located!
"""

from pathlib import Path
import shutil

# import file_mapping.py
# from .file_mapping import file_mapping

file_mapping = {
    "fig_1b": "results/bar_plot/alignment/mapping_stats",
    "fig_1c": "results/pca/Time_Status/pca_cropped",
    "fig_1d": "results/heatmap_complex/expr_hclust_complete/clustering.top_100", 
    "fig_1e": "results/scatter_plot/logfc_vs_expr_coloured=Time", 
    "fig_1f": "results/pvca/pvca", 
    "fig_2b": ["results/volcano/Status__Treatment_vs_Control.cohend_estimate", "results/volcano/Status__Treatment_vs_Control.ttest.adj.labels", "results/volcano/Status__Treatment_vs_Control_Time=120m.cohend_estimate", "results/volcano/Status__Treatment_vs_Control_Time=120m.ttest.adj.labels"], 
    "fig_2c": ["results/volcano/Tissue__left ear_vs_right ear_Status=Treatment.cohend_estimate", "results/volcano/Tissue__left ear_vs_right ear_Status=Treatment.ttest.adj.labels", "results/volcano/Tissue__left ear_vs_right ear_Status_Time=Treatment_120m.cohend_estimate", "results/volcano/Tissue__left ear_vs_right ear_Status_Time=Treatment_120m.ttest.adj.labels"], 
    "fig_2d": ["results/volcano/Status__Treatment_vs_Control_Tissue=left ear.cohend_estimate", "results/volcano/Status__Treatment_vs_Control_Tissue=left ear.ttest.adj.labels", "results/volcano/Status__Treatment_vs_Control_Tissue_Time=left ear_120m.cohend_estimate", "results/volcano/Status__Treatment_vs_Control_Tissue_Time=left ear_120m.ttest.adj.labels"],
    "fig_2e": "results/violin_plots/expr/Group/mmu-miR-144-3p_log10", 
    "fig_2f": "results/violin_plots/expr/Group/mmu-miR-223-3p_log10", 
    "fig_3a": ["results/heatmap_complex/ttest_adjp", "results/heatmap_complex/fc", "results/heatmap_complex/auc_value"], 
    "fig_3b": "results/heatmap_complex/corr_plots/Time_Status_filtered", 
    "fig_3c": "results/heatmap_complex/corr_plots/Timepoint_Group_filtered", 
    "fig_3d": "results/box_plot_l=2/mmu-miR-144-3p", 
    "fig_3e": "results/box_plot_l=2/mmu-miR-223-3p", 
    "fig_3f": "results/box_plot_l=2/mmu-miR-451a", 
    "fig_3g": "results/box_plot_l=2/mmu-miR-142a-5p", 
    "fig_4a": "results/heatmap_complex/pathway_pvalues/kegg_prediction_union", 
    "fig_4d": "results/halushka_plots/hsa-miR-144-3p_sorted_by_median_num_cell_type=50", 
    "fig_4e": "results/halushka_plots/hsa-miR-223-3p_sorted_by_median_num_cell_type=50", 
    "fig_4f": "results/halushka_plots/hsa-miR-451a_sorted_by_median_num_cell_type=50", 
    "fig_4g": "results/halushka_plots/hsa-miR-142-5p_sorted_by_median_num_cell_type=50", 
    "supp_fig_1a": "results/bar_plot/mapping/rna_composition_mapped", 
    "supp_fig_2a": ["results/volcano/Tissue__left ear_vs_right ear_Status=Control.cohend_estimate", "results/volcano/Tissue__left ear_vs_right ear_Status=Control.ttest.adj.labels", "results/volcano/Tissue__left ear_vs_right ear_Status_Time=Control_120m.cohend_estimate", "results/volcano/Tissue__left ear_vs_right ear_Status_Time=Control_120m.ttest.adj.labels"], 
    "supp_fig_2b": ["results/volcano/Status__Treatment_vs_Control_Tissue=right ear.cohend_estimate", "results/volcano/Status__Treatment_vs_Control_Tissue=right ear.ttest.adj.labels", "results/volcano/Status__Treatment_vs_Control_Tissue_Time=right ear_120m.cohend_estimate", "results/volcano/Status__Treatment_vs_Control_Tissue_Time=right ear_120m.ttest.adj.labels"], 
    "supp_tab_1": "results/matrices/supp_volcano/supp_table_all", 
    "supp_tab_2": "results/matrices/supp_corr/supp_table_all",
    }

# relevant file types
#file_types = ["png", "svg", "csv"]
file_types = ["png", "csv", "xlsx"]

# snakemake output folder
source_figures_folder = Path("results/")
# check if it is really a folder
if source_figures_folder.is_dir():
    print(f"{source_figures_folder} exists and is a folder.")
else:
    print(f"{source_figures_folder} is not a folder ... aborting!")
    raise SystemExit("Bye!")

# create output folder if non-existant
pub_figures_folder = Path("publication/")
pub_figures_folder.mkdir(parents=True, exist_ok=True)

# loop over all figures
for pub_figure, source_paths in file_mapping.items():
    # if the dictionary item is not a list (multiple sub-plots), then make one with 1 item
    if not (isinstance(source_paths, list) or isinstance(source_paths, tuple)):
        source_paths = [source_paths]
    for i, source_path in enumerate(source_paths):
        source_path = Path(source_path)
        for file_type in file_types:
            # build full filename with ending and check if it exists
            source_full_path = Path(f"{source_path}.{file_type}")
            if not source_full_path.is_file():
                # print(f"{source_full_path} does not exist ... skipping!")
                continue
            # only add number if len(source_paths) > 1
            if len(source_paths) > 1:
                target_full_path = pub_figures_folder / Path(f"{pub_figure}_{i+1:.0f}.{file_type}")
            else:
                target_full_path = pub_figures_folder / Path(f"{pub_figure}.{file_type}")
            shutil.copyfile(source_full_path, target_full_path)
            print(f"creating {target_full_path} from {source_full_path}")

# remove figures folder
print(f"removing {source_figures_folder}!")
shutil.rmtree(source_figures_folder)


with open(snakemake.output[0], "w") as file:
    file.write("")