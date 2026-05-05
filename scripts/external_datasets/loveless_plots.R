library(Seurat)
library(data.table)
library(ggplot2)
library(patchwork)
library(sagesalad)

genesets <- fread("data/scRNAseq_genesets.csv")
genesets_list <- strsplit(genesets$genes, " ")
names(genesets_list) <- genesets$Geneset

tumour <- readRDS("external/loveless_2025_pdac_atlas/TumorDuctal_01242024.rds")
ductal <- readRDS("external/loveless_2025_pdac_atlas/Updated_DuctalClusters.rds")

new_genesets <- genesets_list[startsWith(names(genesets_list), "DE_")]
ductal <- AddModuleScore(ductal, features=new_genesets, name=names(new_genesets))
colnames(ductal@meta.data)[match(paste0(names(new_genesets), 1:length(names(new_genesets))), colnames(ductal@meta.data))] <- names(new_genesets)

tumour <- AddMetaData(tumour, ductal[[names(new_genesets)]])

source("singler_subtype.R")
singler_reference_file <- "references/bulk_reference.RData"

mat <- ductal@assays$RNA@data
out <- singler_subtype(mat, reference_file=singler_reference_file)[[1]]
ductal <- AddMetaData(ductal, out)
tumour <- AddMetaData(tumour, out)

png.date("loveless_UMAP_basal.png", width = 14, height=7)
FeaturePlot(tumour, features = c("DE_Basal1", "DE_Basal2"), raster = F, col = c("yellow", "blue"),
            min.cutoff=0, max.cutoff=max(tumour$DE_Basal2)) &
  theme(aspect.ratio = 1,
        axis.ticks=element_blank(),
        axis.text=element_blank())
dev.off()  

png.date("loveless_UMAP_basal_blend_ordered.png", width = 15, height=5)
FeaturePlot(tumour, features = c("DE_Basal1", "DE_Basal2"), order = T, blend = T , raster = F, col=c("#7570B3", "#DB6002")) &
  theme(aspect.ratio = 1,
        axis.ticks=element_blank(),
        axis.text=element_blank())
dev.off()

png.date("loveless_UMAP_basal_blend_ordered_tiny.png", width = 9, height=3)
FeaturePlot(tumour, features = c("DE_Basal1", "DE_Basal2"), order = T, blend = T , raster = F, col=c("#7570B3", "#DB6002")) &
  theme(aspect.ratio = 1)
dev.off()

