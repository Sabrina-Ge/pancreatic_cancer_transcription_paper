import pandas as pd
import scanpy as sc
import numpy as np
import re
import matplotlib.pyplot as plt
import scanpy.external as sce
import scrublet as scr

PERCENT_MITO_CUTOFF = 15
N_NEIGH = 50
MIN_GENES = 400
RES = [0.8, 1.0, 1.2]

mtx_paths_file = "source/all_mat_dir_full_path.txt"
gene_sets_files = "data/scRNAseq_genesets.csv"

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Read in sample paths
h5_paths = [line.rstrip() for line in open(mtx_paths_file)]

# Collect samples and sample ids, names
adatas = [sc.read_10x_mtx(path) for path in h5_paths]
sample_data = [re.findall("/[a-zA-Z0-9_]+/[a-zA-Z0-9_]+/cellranger/", path)[0] for path in h5_paths]
sample_ids = [data.split("/")[2] for data in sample_data]
sample_names = [data.split("/")[1] for data in sample_data]
sample_names = [sample_name.split("_")[1] if sample_name.startswith("PDA_") else sample_name for sample_name in sample_names]

# Special cases
sample_names = [re.sub("^0", "", sample_name) for sample_name in sample_names]
sample_names = [re.sub("_CMP$", "", sample_name) if sample_name != "85948_CMP" else sample_name for sample_name in sample_names]

adata_dict = dict(zip(sample_ids, adatas))
sample_id_to_name = dict(zip(sample_ids, sample_names))

# Merge samples
merged = sc.concat(adata_dict, label="orig_ident", index_unique="_")
merged.var_names_make_unique()
merged.obs["Sample"] = pd.Series([sample_id_to_name[id] for id in merged.obs["orig_ident"]]).values

# Doublet prediction 
doublet_data = {}
for sample_id in merged.obs["Sample"].unique(): # per Sample vs per orig.ident?
    adata = merged[merged.obs["Sample"] == sample_id,] 
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.05)  # 5% doublet rate for input 10000 cells
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    doublet_data[sample_id] = pd.DataFrame({"doublet_score": doublet_scores, "predicted_doublets": predicted_doublets}, index = adata.obs.index) 
    scrub.plot_histogram()
    plt.savefig(f"scrublet_hist_{sample_id}.png")

all_doublet_data = pd.concat(doublet_data.values())
merged.obs["doublet_score"] = all_doublet_data["doublet_score"] 
merged.obs["predicted_doublets"] = all_doublet_data["predicted_doublets"] 

# QC and Filtering
sc.pl.highest_expr_genes(merged, n_top=20, )
plt.savefig("highest_expr_genes.png")

merged.write("merged_pre_filtering.h5ad")

sc.pp.filter_cells(merged, min_genes=MIN_GENES)
merged.var['mt'] = merged.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(merged, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(merged, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
plt.savefig("qc_violin.pdf")
sc.pl.scatter(merged, x='total_counts', y='pct_counts_mt')
plt.savefig("qc_scatter1.pdf")
sc.pl.scatter(merged, x='total_counts', y='n_genes_by_counts')
plt.savefig("qc_scatter2.pdf")

merged = merged[merged.obs.pct_counts_mt < PERCENT_MITO_CUTOFF, :]

# Save pre-normalization (RAW) counts
merged.layers["raw_counts"] = merged.X.copy()

# Normalize
sc.pp.normalize_total(merged, target_sum=1e4)
sc.pp.log1p(merged)
sc.pp.highly_variable_genes(merged)
#sc.pp.highly_variable_genes(merged, flavor="seurat_v3")

# SAM
sam_obj = sce.tl.sam(merged, inplace=True, k=N_NEIGH)

# Clustering 
for res in RES:
    sam_obj.clustering(param=res, method="leiden")
    sam_obj.adata.obs.rename(columns={"leiden_clusters": "leiden_clusters_" + str(res)}, inplace=True)

# Scoring
modules = pd.read_csv(gene_sets_files)
[sc.tl.score_genes(merged, gene_list=modules.iloc[i]["genes"].split(" "), score_name=modules.iloc[i]["Geneset"]) for i in range(modules.shape[0])]
for filter_term in ["Moffitt", "Notta", "Collisson", "Bailey", "Raghavan", "Hwang"]:
    merged.obs[filter_term + "_subtype"] = merged.obs.filter(regex=filter_term + "_").idxmax(1)


merged.obs["Celltype_raw"] = merged.obs.filter(regex="Markers_").idxmax(1)
for res in RES:
    merged.obs["Celltype_by_cluster_" + str(res)] = merged.obs.groupby("leiden_clusters_" + str(res)).transform("mean").filter(regex="Markers_").idxmax(1)

# replace Markers_* with just *
start = len("Markers_")
merged.obs["Celltype_raw"] = merged.obs["Celltype_raw"].str[start:]
for res in RES:
    merged.obs["Celltype_by_cluster_" + str(res)] = merged.obs["Celltype_by_cluster_" + str(res)].str[start:]

for res in RES:
    sc.pl.umap(merged, color="Celltype_by_cluster_" + str(res))
    plt.savefig("qc_umap_celltype_" + str(res) + ".pdf")

# Differential expression
sc.tl.rank_genes_groups(merged, 'leiden_clusters_0.8', method='wilcoxon', pts=True)

# Save complete anndata object
merged.write("merged.h5ad")

# Pancreatic only subset
merged_panc = merged[merged.obs["Celltype_by_cluster_0.8"].isin(["Malignant", "Duct", "Acinar", "Endocrine"])]
sc.pp.highly_variable_genes(merged_panc)
sam_obj_panc = sce.tl.sam(merged_panc, inplace=True, k=N_NEIGH)
sam_obj_panc.clustering(param=0.8, method="leiden")
sc.tl.rank_genes_groups(merged_panc, 'leiden_clusters', method='wilcoxon', pts=True)

merged_panc.write("merged_panc.h5ad")

# Malignant only subset
merged_mal = merged[merged.obs["Celltype_by_cluster_0.8"].isin(["Malignant"])]
sc.pp.highly_variable_genes(merged_mal)
sam_obj_mal = sce.tl.sam(merged_mal, inplace=True, k=N_NEIGH)
sam_obj_mal.clustering(param=0.8, method="leiden")
sc.tl.rank_genes_groups(merged_mal, 'leiden_clusters', method='wilcoxon', pts=True)

merged_mal.write("merged_mal.h5ad")

# Single sample re-do UMAP
for sample in merged.obs["Sample"].unique():
    try:
        cells_to_subset = merged.obs[merged.obs["Sample"] == sample].index
        sample_obj = merged[cells_to_subset]
        sam_sample_obj = sce.tl.sam(sample_obj, inplace=True, k=N_NEIGH)
        sam_sample_obj.clustering(param=0.8, method="leiden")
        sc.tl.rank_genes_groups(sample_obj, 'leiden_clusters', method='wilcoxon', pts=True)
        sample_obj.write("sample_{}.h5ad".format(sample))
        cells_to_subset = merged_mal.obs[merged_mal.obs["Sample"] == sample].index
        sample_obj = merged_mal[cells_to_subset]
        sam_sample_obj = sce.tl.sam(sample_obj, inplace=True, k=N_NEIGH)
        sam_sample_obj.clustering(param=0.8, method="leiden")
        sc.tl.rank_genes_groups(sample_obj, 'leiden_clusters', method='wilcoxon', pts=True)
        sample_obj.write("sample_{}_mal.h5ad".format(sample))
    except Exception as e:
        print(f"Error with sample {sample}: {e}")
        continue

print("DONE")
