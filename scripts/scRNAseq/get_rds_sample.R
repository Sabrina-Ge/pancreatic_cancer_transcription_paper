library(Seurat)
library(Matrix)
library(hdf5r)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
file <- args[1]

sample <- gsub(".h5ad", "", file)
subtype_ss <- read.csv("merged_subtype_scores.csv", header=TRUE, row.names=1)
celltype_final <- read.csv("merged_celltype_final.csv", header=TRUE, row.names=1)

adata <- H5File$new(filename=file, mode = "r")

cell_names <- adata[["obs"]][["_index"]][]
gene_names <- gsub("_", "-", adata[["var"]][["_index"]][])

# data counts
x <- adata[["X"]]
count_mat <- sparseMatrix(i = x[["indices"]][], p = x[["indptr"]][], x = x[["data"]][], 
                          index1 = FALSE, dims=c(length(gene_names), length(cell_names)))
dimnames(count_mat) <- list(gene_names, cell_names)
obj <- CreateSeuratObject(CreateAssayObject(data = count_mat), project=sample)

# raw counts
y <- adata[["layers"]][["raw_counts"]]
raw_count_mat <- sparseMatrix(i = y[["indices"]][], p = y[["indptr"]][], x = y[["data"]][], 
                              index1 = FALSE, dims=c(length(gene_names), length(cell_names)))
dimnames(raw_count_mat) <- list(gene_names, cell_names)
obj@assays$RNA@counts <- raw_count_mat

# Add DimReduc data
pcs <- t(adata[["obsm"]][["X_pca"]][,])
dimnames(pcs) <- list(cell_names, paste0("PC_", 1:ncol(pcs)))
obj[["PCA"]] <- CreateDimReducObject(embeddings=pcs, assay="RNA", key="PC_")

umap <- t(adata[["obsm"]][["X_umap"]][,])
dimnames(umap) <- list(cell_names, paste0("UMAP_", 1:ncol(umap)))
obj[["UMAP"]] <- CreateDimReducObject(embeddings=umap, assay="RNA", key="UMAP_")

# Add SAM cluster metadata
metadata_names <- names(adata[["obs"]])[!grepl("^__", names(adata[["obs"]]))]
clust_meta_df <- do.call(cbind.data.frame, lapply(metadata_names, function(x) adata[["obs"]][[x]][]))
colnames(clust_meta_df) <- metadata_names
rownames(clust_meta_df) <- as.character(clust_meta_df[["_index"]])
clust_meta_df[["_index"]] <- NULL

# Convert categoricals to strings
for (categorical in names(adata[["obs"]][["__categories"]])) {
  dict <- adata[["obs"]][["__categories"]][[categorical]][]
  clust_meta_df[[categorical]] <- dict[clust_meta_df[[categorical]] + 1]
}

# Convert some integers to categoricals
clust_meta_df[, "leiden_clusters"] <- as.factor(clust_meta_df[, "leiden_clusters"])

obj <- AddMetaData(obj, metadata=clust_meta_df)

# Grab DE genes into misc
de_adata <- adata[["uns"]][["rank_genes_groups"]]

clusters <- colnames(de_adata[["pvals"]][])
n_genes_per_cluster <- nrow(de_adata[["pvals"]][])

seurat_format_df <- cbind.data.frame(p_val = unlist(de_adata[["pvals"]][]), 
                                     avg_log2FC = unlist(de_adata[["logfoldchanges"]][]), 
                                     pct.1 = 1, 
                                     pct.2 = 1, 
                                     p_val_adj = unlist(de_adata[["pvals_adj"]][]), 
                                     cluster = rep(clusters, each=n_genes_per_cluster), 
                                     gene = unlist(de_adata[["names"]][]),
                                     z_score = unlist(de_adata[["scores"]][]))

# Use same cutoffs as Seurat::FindAllMarkers:
# abs(logfc) > 0.25, pct > 0.1
# also, pval_adj < 0.05
seurat_format_df_filtered <- seurat_format_df[abs(seurat_format_df$avg_log2FC) > 0.25 & 
                                                seurat_format_df$pct.1 > 0.1 & 
                                                seurat_format_df$p_val_adj < 0.05, ]

obj@misc[["de_output"]] <- seurat_format_df_filtered

obj <- AddMetaData(obj, metadata=subtype_ss[rownames(obj@meta.data), ])
obj <- AddMetaData(obj, metadata=celltype_final[rownames(obj@meta.data), , drop=FALSE])

# Close H5 file
adata$close_all()
saveRDS(obj, sprintf("%s.rds", sample))

