library(Seurat)
library(Matrix)
library(hdf5r)
library(data.table)

source("singler_subtype.R")
singler_reference_file <- "references/bulk_reference.RData"
sample_metadata_file <- "data/scRNAseq_metadata.csv"
sample_metadata <- read.csv(sample_metadata_file, fileEncoding="UTF-8-BOM", stringsAsFactors = FALSE)

for (file in c("merged.h5ad", "merged_panc.h5ad", "merged_mal.h5ad")) {
  
  file <- "merged_panc.h5ad"
  sample <- gsub(".h5ad", "", file)
  
  adata <- H5File$new(filename=file, mode = "r")
  
  cell_names <- adata[["obs"]][["_index"]][]
  gene_names <- gsub("_", "-", adata[["var"]][["_index"]][])
  
  # data counts
  x <- adata[["X"]]
  count_mat <- sparseMatrix(i = x[["indices"]][], p = x[["indptr"]][], x = x[["data"]][], 
                            index1 = FALSE, dims=c(length(gene_names), length(cell_names)))
  dimnames(count_mat) <- list(gene_names, cell_names)
  obj <- CreateSeuratObject(CreateAssayObject(data = count_mat), project="MergedPDA")
  
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
  if (!("Celltype_by_cluster" %in% colnames(clust_meta_df))) {
    clust_meta_df[, "Celltype_by_cluster"] <- clust_meta_df[, "Celltype_by_cluster_0.8"]
  }
  
  
  # Convert some integers to categoricals
  if (!("leiden_clusters" %in% colnames(clust_meta_df))) {
    for (res in c("0.8", "1.0", "1.2")) {
      clust_meta_df[, paste0("leiden_clusters_", res)] <- as.factor(clust_meta_df[, paste0("leiden_clusters_", res)])
    }
    clust_meta_df[, "leiden_clusters"] <- clust_meta_df[, "leiden_clusters_0.8"]
  } else {
    clust_meta_df[, "leiden_clusters"] <- as.factor(clust_meta_df[, "leiden_clusters"])
  }
  
  obj <- AddMetaData(obj, metadata=clust_meta_df)
  
  # Sample in the obj metadata, and sample in the metadata file should be an exact match
  rownames(sample_metadata) <- sample_metadata$sample
  sample_metadata_to_add <- sample_metadata[obj@meta.data$Sample, ]
  rownames(sample_metadata_to_add) <- rownames(obj@meta.data)
  
  obj <- AddMetaData(obj, metadata=sample_metadata_to_add)
  
  # Grab DE genes into misc
  de_adata <- adata[["uns"]][["rank_genes_groups"]]
  
  clusters <- colnames(de_adata[["pvals"]][])
  n_genes_per_cluster <- nrow(de_adata[["pvals"]][])
  
  # ordered by different genes, not by the ones expected
  pct_1_geneorder <- de_adata[["pts"]][["_index"]][]
  pct_1_list <- lapply(clusters, function(i){ 
    unordered <- de_adata[["pts"]][[i]][] 
    names(unordered) <- pct_1_geneorder
    return(unordered[de_adata[["names"]][][[i]]])
  })
  pct_1 <- do.call(c, pct_1_list)
  
  pct_2_geneorder <- de_adata[["pts"]][["_index"]][]
  pct_2_list <- lapply(clusters, function(i){ 
    unordered <- de_adata[["pts"]][[i]][] 
    names(unordered) <- pct_2_geneorder
    return(unordered[de_adata[["names"]][][[i]]])
  })
  pct_2 <- do.call(c, pct_2_list)
  
  seurat_format_df <- cbind.data.frame(p_val = unlist(de_adata[["pvals"]][]), 
                                       avg_log2FC = unlist(de_adata[["logfoldchanges"]][]), 
                                       pct.1 = pct_1, 
                                       pct.2 = pct_2, 
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
  
  fwrite(seurat_format_df_filtered, sprintf("%s_de_by_leiden_clusters.csv", sample))
  
  # Variable features
  VariableFeatures(obj) <- adata[["var"]][["_index"]][][adata[["var"]][["highly_variable"]][]]
  
  # Close H5 file
  adata$close_all()
  
  obj <- SetIdent(obj, value="leiden_clusters")
  
  if (sample == "merged") {
    mat <- obj@assays$RNA@data
    out <- singler_subtype(mat, reference_file=singler_reference_file)[[1]]
    obj <- AddMetaData(obj, out)
    
    data_to_save <- out
    
    write.csv(data_to_save, sprintf("merged_subtype_scores.csv"), quote=FALSE)
  } else {
    subtype_scores <- read.csv("merged_subtype_scores.csv", row.names=1)
    obj <- AddMetaData(obj, subtype_scores)
  }
  
  saveRDS(obj, sprintf("%s.rds", sample))
  
}