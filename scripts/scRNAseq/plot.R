library(Seurat)
library(Matrix)
library(hdf5r)
library(ggplot2)
library(sagesalad)
library(data.table)
library(RColorBrewer)
library(patchwork)
library(Hmisc)
library(corrplot)
library(ggpubr)
library(ggrepel)
library(pals)
library(vegan)
library(scCustomize)
library(ggmosaic)

obj <- readRDS("merged_celltype_final.rds")
sample <- "merged"
dir.create(sample, showWarnings = FALSE)

sample_metadata_file <- "data/scRNAseq_metadata.csv"
metadata <- fread(sample_metadata_file)

genesets <- fread("data/scRNAseq_genesets.csv")
genesets_list <- strsplit(genesets$genes, " ")
names(genesets_list) <- genesets$Geneset

celltype_colour_file <- "data/celltype_colours.csv"
celltype_colours <- fread(celltype_colour_file)
cell_plot_order <- celltype_colours$celltype
cell_plot_colours <- celltype_colours$hexcode
names(cell_plot_colours) <- cell_plot_order

Moffitt_colours = c("Basal-like" = "#BEC100", "Classical" = "#B0ABFF")
Moffitt_subtype_colors <- c("Moffitt_Basal-like" = "#BEC100", "Moffitt_Classical" = "#B0ABFF")

singler_subtype_colors <- c("#1C6CAB", "#814C42","#A4C0E5", "#FF7311", "grey40")
names(singler_subtype_colors) <- c("Classical1", "Basal1", "Classical2", "Basal2", "Mixed")
singler_subtype_colors <- singler_subtype_colors[c("Basal1", "Basal2", "Mixed", "Classical1", "Classical2")]

singler_subtype_colors_malignant <- c("#1C6CAB", "#814C42","#A4C0E5", "#FF7311", "grey40", "grey90")
names(singler_subtype_colors_malignant) <- c("Classical1", "Basal1", "Classical2", "Basal2", "Mixed", "non-malignant")
singler_subtype_colors_malignant <- singler_subtype_colors_malignant[c("Basal1", "Basal2", "Mixed", "Classical1", "Classical2", "non-malignant")]

sample_colours <- unique(c(alphabet(), alphabet2()))[1:length(unique(obj$Sample))]
names(sample_colours) <- unique(obj$Sample)
organ_colours <- c(Pa="#ffed6f", Lu="#80b1d3", Lv="#fb8072")
type_colours <- c(primary="#ccebc5", metastatic="#fccde5")


#### Set final celltype ####

# Fix 102409 label (informed by Numbat)
celltype_final <- ifelse(obj$Celltype_by_cluster == "Hepatocyte" & obj$Sample == "102409", "Malignant", celltype_final)

# Set COMP_0158_P to Neuroendocrine by udpated patient info
celltype_final <- ifelse(obj$Celltype_by_cluster == "Malignant" & obj$Sample == "COMP_0158_P", "Neuroendocrine", obj$Celltype_by_cluster)

# Assign Duct1 and Duct2
obj$Duct1vDuct2 <- ifelse(obj$DCZhou_Ductlike1 > obj$DCZhou_Ductlike2, "Duct1", "Duct2")
celltype_final <- ifelse(celltype_final == "Duct", obj$Duct1vDuct2, celltype_final)

celltype_final <- factor(as.character(celltype_final), levels=cell_plot_order)

obj$Celltype_final <- celltype_final
write.csv(obj[["Celltype_final"]], "merged_celltype_final.csv", quote=FALSE)

obj$Celltype_final_no_neuroendocrine <- ifelse(obj$Celltype_final == "Neuroendocrine", "Malignant", as.character(obj$Celltype_final))
obj$Celltype_final_no_neuroendocrine <- factor(as.character(obj$Celltype_final_no_neuroendocrine), levels=cell_plot_order)


#### Stats ####

# avg reads per cell = total_counts
# avg genes per cell = n_genes
# avg reads per malignant cell
# avg genes per malignant cell
# avg mito percent = total_counts_mt - > pct_counts_mt
# avg mito percent per malignant cell
# cells of each celltype

dt <- as.data.table(obj@meta.data)
summary <- dt[ , .(mean_genes=mean(n_genes), mean_counts=mean(total_counts), mean_mt=mean(pct_counts_mt), 
                   median_genes=median(n_genes), median_counts=median(total_counts), median_mt=median(pct_counts_mt)), by="Sample"]
fwrite(summary, sprintf("%s/summary_all_stats.csv", sample))

dt_mal <- dt[Celltype_final=="Malignant", ]
summary_mal <- dt_mal[ , .(mean_genes=mean(n_genes), mean_counts=mean(total_counts), mean_mt=mean(pct_counts_mt), 
                           median_genes=median(n_genes), median_counts=median(total_counts), median_mt=median(pct_counts_mt)), by="Sample"]
fwrite(summary_mal, sprintf("%s/summary_malignant_stats.csv", sample))

pdf.date(sprintf("%s/%s_scatter_counts_vs_genes.pdf", sample, sample), width=7, height=7)
FeatureScatter(obj, "total_counts", "n_genes") + theme(aspect.ratio = 1)
FeatureScatter(obj, "total_counts", "n_genes", group.by="Celltype_final", cols=cell_plot_colors) + theme(aspect.ratio = 1)
dev.off()

subtype_data <- as.data.table(obj[[c("Sample", "subtype", "Celltype_final")]])
subtype_data <- subtype_data[Celltype_final == "Malignant", ]
summary_df <- as.data.frame.matrix(table(subtype_data[, .(Sample, subtype)]))
summary_dt <- as.data.table(summary_df, keep.rownames="Sample")

summary_dt[, total:= rowSums(.SD[, .(Basal1, Basal2, Classical1, Classical2, Mixed)])]

prop_cols <- paste0(c("Basal1", "Basal2", "Classical1", "Classical2", "Mixed"), "_pct")
summary_dt[, (prop_cols) := lapply(.SD, function(x) sprintf("%.2f", x/total * 100)), .SDcols=c("Basal1", "Basal2", "Classical1", "Classical2", "Mixed")]
fwrite(summary_dt, datestamp(sprintf("%s/%s_count_pct_subtype_summary.csv", sample, sample)))


#### Pseudobulk summary ####

source("singler_subtype.R")

subtype_data <- as.data.table(obj[[c("Sample", "subtype", "Celltype_final")]])
subtype_data <- subtype_data[Celltype_final == "Malignant", ]
summary_df <- as.data.frame.matrix(table(subtype_data[, .(Sample, subtype)]))
summary_dt <- as.data.table(summary_df, keep.rownames="Sample")

Idents(obj) <- "Sample"
sample_mat_all <- AggregateExpression(obj, assay="RNA", slot="counts")
norm_sample_mat <- log(apply(sample_mat_all$RNA, 2, function(x) x * 1e6 / sum(x)) + 1) # logCPM

singler_out <- singler_subtype(norm_sample_mat, reference_file="bulk_reference.RData")[[1]]
rownames(singler_out) <- gsub("-", "_", gsub("^g", "", rownames(singler_out)))
summary_dt$singler_by_sample <- singler_out[summary_dt$Sample, "subtype"]

fwrite(summary_dt, sprintf("%s/%s_count_subtype_summary.csv", sample, sample))


#### Plotting: Duct1 vs Duct2 ####

duct_de <- FindMarkers(obj, ident.1="Duct1", ident.2="Duct2", group.by="Celltype_final")
write.csv(duct_de, sprintf("%s/%s_duct1_duct2_de.csv", sample, sample))

duct1_filtered <- duct_de[duct_de$p_val_adj < 0.05 & duct_de$avg_log2FC > 2,]
duct2_filtered <- duct_de[duct_de$p_val_adj < 0.05 & duct_de$avg_log2FC < -2,]

duct1_filtered <- duct1_filtered[order(duct1_filtered$p_val_adj), ]
duct2_filtered <- duct2_filtered[order(duct2_filtered$p_val_adj), ]

n <- 25
subtype_genesets <- c("DE_Basal1", "DE_Basal2", "DE_Classical1", "DE_Classical2")
ducttypes <- list(duct1=duct1_filtered, duct2=duct2_filtered)
intersecting_genes <- matrix(nrow=4, ncol=2, dimnames=list(subtype_genesets, names(ducttypes)))
for (subtype_geneset in subtype_genesets) {
  for (ducttype_name in names(ducttypes)) {
    intersecting_genes[subtype_geneset, ducttype_name] <- paste0(intersect(genesets_list[[subtype_geneset]], rownames(ducttypes[[ducttype_name]])[1:n]), collapse=" ")
  }
}
write.csv(intersecting_genes, sprintf("%s/%s_interesecting_genes_duct1_duct2_vs_subtypes.csv", sample, sample))

top_genes <- c(rownames(duct1_filtered)[1:n], rownames(duct2_filtered)[1:n])

scaled_obj <- ScaleData(obj, features=top_genes)
duct_obj <- subset(scaled_obj, Celltype_final %in% c("Duct1", "Duct2"))
duct_obj$Celltype_final <- factor(as.character(duct_obj$Celltype_final))
png.date(sprintf("%s/%s_heatmap_duct1_duct2.png", sample, sample), width=7, height=14)
DoHeatmap(duct_obj, features=top_genes, group.by="Celltype_final", raster=FALSE) + 
  scale_fill_gradient2(low="plum4", mid="black", high="yellow", midpoint=0, na.value = "white")
dev.off()

tosti_cp <- readRDS("external_data/tosti/chronic_pancreatitis/tosti_cp.RDS")
tosti_cp <- ScaleData(tosti_cp, features=top_genes)
tosti_cp_sub <- subset(tosti_cp, Cluster %in% c("Ductal", "MUC5B+ Ductal"))
tosti_cp_sub$Cluster <- factor(as.character(tosti_cp_sub$Cluster))
png.date(sprintf("%s/%s_heatmap_tosti_cp_duct.png", sample, sample), width=7, height=14)
DoHeatmap(tosti_cp_sub, features=top_genes, group.by="Cluster", raster=FALSE)  + 
  scale_fill_gradient2(low="plum4", mid="black", high="yellow", midpoint=0, na.value = "white")
dev.off()



#### Plotting: Assessing celltype genes ####

genes_to_plot <- unlist(genesets_list[grep("Markers_", names(genesets_list))])
names(genes_to_plot) <- NULL
genes_to_plot <- unique(genes_to_plot)

pdf.date(sprintf("%s/%s_dotplot_marker_geneset_ungrouped.pdf", sample, sample), height=7, width=12)
print(DotPlot(obj, features=genes_to_plot, group.by="Celltype_final") +
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)))
dev.off()

genes_to_plot_grouped <- genesets_list[grep("Markers_", names(genesets_list))]
genes_to_plot_grouped$Markers_Myofibroblast <- genes_to_plot_grouped$Markers_Myofibroblast[genes_to_plot_grouped$Markers_Myofibroblast != "DCN"]
names(genes_to_plot_grouped) <- gsub("^Markers_", "", names(genes_to_plot_grouped))
genes_to_plot_grouped <- genes_to_plot_grouped[c('Malignant', 'Duct', 'Acinar', 'Endocrine', 'Fibroblast', 'Myofibroblast', 'Endothelial', 
                                                 'BPlasmaCell', 'Macrophage', 'Mast', 'TCell', 'NKCell', 'Hepatocyte', 'Lung')]

pdf.date(sprintf("%s/%s_dotplot_marker_geneset_grouped.pdf", sample, sample), height=7, width=14)
print(DotPlot(obj, features=genes_to_plot_grouped, group.by="Celltype_final") +
        ylab("Cell types") + xlab("") + 
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), 
              strip.text.x.top=element_text(angle=90, hjust=0, vjust=0.5)))
dev.off()


markers <- obj@misc$de_output
setDT(markers)
top_dt <- markers[order(cluster, -avg_log2FC)][, head(.SD, 20) , by=cluster]

obj <- ScaleData(obj, features=unique(c(VariableFeatures(obj), top_dt$gene)))

cluster_naming <- unique(obj[[c("Celltype_final", "leiden_clusters")]])
cluster_naming <- cluster_naming[cluster_naming$Celltype_final != "Malignant", ]
celltypes_to_investigate <- unique(cluster_naming$Celltype_final[duplicated(cluster_naming$Celltype_final)])

for (celltype in celltypes_to_investigate) {
  
  clusters <- cluster_naming[cluster_naming$Celltype_final == celltype, ]$leiden_cluster
  top_sub <- top_dt[cluster %in% clusters]$gene
  
  pdf.date(sprintf("%s/%s_heatmap_stroma_markers_%s.pdf", sample, celltype), height=10, width=10)
  print(DoHeatmap(subset(obj, subset = leiden_clusters %in% clusters), features=top_sub, slot="data") + ggtitle(celltype) + 
          scale_fill_gradientn(colours = rev(brewer.pal(8, "RdBu")), na.value="white"))
  print(DoHeatmap(subset(obj, subset = leiden_clusters %in% clusters), features=top_sub, slot="scale.data") + ggtitle(celltype))
  dev.off()
  
}


#### Plotting: Scatterplots ####

obj_mal <- subset(obj, subset = Celltype_final == "Malignant")
data <- obj_mal[[c("Basal1_score", "Basal2_score",  "Classical1_score", "Classical2_score", "subtype")]]
numeric_data <- as.matrix(data[, c("Basal1_score", "Basal2_score",  "Classical1_score", "Classical2_score")])
axis_ranges <- range(numeric_data) * 1.1

pdf.date(sprintf("%s/%s_densityplot_singler_scores_all_samples.pdf", sample, sample), width=10, height=10)
ggplot(data, aes(x=Basal2_score, y=Basal1_score)) + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  xlim(axis_ranges[1], axis_ranges[2]) + ylim(axis_ranges[1], axis_ranges[2]) +
  theme(aspect.ratio = 1) +
  theme_bw()
ggplot(data, aes(x=Classical1_score, y=Basal2_score)) + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  xlim(axis_ranges[1], axis_ranges[2]) + ylim(axis_ranges[1], axis_ranges[2]) +
  theme(aspect.ratio = 1) +
  theme_bw()
ggplot(data, aes(x=Classical1_score, y=Basal1_score)) + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  xlim(axis_ranges[1], axis_ranges[2]) + ylim(axis_ranges[1], axis_ranges[2]) +
  theme(aspect.ratio = 1) +
  theme_bw()
dev.off()

library(ggpp)
ch_plot <- ggplot(data, aes(x=Basal2_score, y=Basal1_score)) + 
  geom_point(shape=".") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  stat_quadrant_counts(aes(label = after_stat(pc.label)), digits=4) +
  xlim(axis_ranges[1], axis_ranges[2]) + ylim(axis_ranges[1], axis_ranges[2]) +
  theme(aspect.ratio = 1) +
  theme_bw()
ic_plot <- ggplot(data, aes(x=Classical1_score, y=Basal2_score)) + 
  geom_point(shape=".") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  stat_quadrant_counts(aes(label = after_stat(pc.label)), digits=4) +
  xlim(axis_ranges[1], axis_ranges[2]) + ylim(axis_ranges[1], axis_ranges[2]) +
  theme(aspect.ratio = 1) +
  theme_bw()
ih_plot <- ggplot(data, aes(x=Classical1_score, y=Basal1_score)) + 
  geom_point(shape=".") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  stat_quadrant_counts(aes(label = after_stat(pc.label)), digits=4) +
  xlim(axis_ranges[1], axis_ranges[2]) + ylim(axis_ranges[1], axis_ranges[2]) +
  theme(aspect.ratio = 1) +
  theme_bw()
ch_plot2 <- ggplot(data, aes(x=Classical1_score, y=Classical2_score)) + 
  geom_point(shape=".") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  stat_quadrant_counts(aes(label = after_stat(pc.label)), digits=4) +
  xlim(axis_ranges[1], axis_ranges[2]) + ylim(axis_ranges[1], axis_ranges[2]) +
  theme(aspect.ratio = 1) +
  theme_bw()
ic_plot2 <- ggplot(data, aes(x=Classical2_score, y=Basal2_score)) + 
  geom_point(shape=".") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  stat_quadrant_counts(aes(label = after_stat(pc.label)), digits=4) +
  xlim(axis_ranges[1], axis_ranges[2]) + ylim(axis_ranges[1], axis_ranges[2]) +
  theme(aspect.ratio = 1) +
  theme_bw()
ih_plot2 <- ggplot(data, aes(x=Classical2_score, y=Basal1_score)) + 
  geom_point(shape=".") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  stat_quadrant_counts(aes(label = after_stat(pc.label)), digits=4) +
  xlim(axis_ranges[1], axis_ranges[2]) + ylim(axis_ranges[1], axis_ranges[2]) +
  theme(aspect.ratio = 1) +
  theme_bw()
png.date(sprintf("%s/%s_scatterplot_singler_scores_all_samples.png", sample, sample), width=12, height=8)
print((ih_plot | ic_plot | ch_plot) / (ih_plot2 | ic_plot2 | ch_plot2))
dev.off()


#### Plotting: UMAPs ####

for (var in colnames(obj@meta.data)) {
  png.date(sprintf("%s/%s_UMAP_raw_%s.png", sample, sample, var))
  if (is.numeric(obj@meta.data[[var]])) {
    print(FeaturePlot(obj, features=var, raster=FALSE, cols=c("yellow", "blue")))
  } else {
    print(DimPlot(obj, group.by=var, raster=FALSE) + ggtitle(var))
  }
  dev.off()
}

obj$subtype_mal <- ifelse(obj$Celltype_final  == "Malignant", obj$subtype, "non-malignant")
png.date(sprintf("%s/%s_UMAP_subtype_mal.png", sample, sample), width=10, height=7)
print(DimPlot(obj, group.by="subtype_mal", raster=FALSE, cols=singler_subtype_colors_malignant, na.value="gray90") + 
        theme(aspect.ratio = 1,
              axis.ticks=element_blank(),
              axis.text=element_blank()))
dev.off()

png.date(sprintf("%s/%s_UMAP_Celltype_final.png", sample, sample), width=10, height=7)
print(DimPlot(obj, group.by="Celltype_final", raster=FALSE, cols=cell_plot_colours) + 
        theme(aspect.ratio = 1,
              axis.ticks=element_blank(),
              axis.text=element_blank()))
dev.off()

png.date(sprintf("%s/%s_UMAP_leiden_clusters_wider.png", sample, sample), width=11, height=7)
print(DimPlot(obj, group.by="leiden_clusters", raster=FALSE, label=TRUE) + 
        theme(aspect.ratio = 1,
              axis.ticks=element_blank(),
              axis.text=element_blank()))
dev.off()

png.date(sprintf("%s/%s_UMAP_Sample_wider.png", sample, sample), width=11, height=7)
print(DimPlot(obj, group.by="Sample", raster=FALSE) + ggtitle("Sample") + 
        theme(aspect.ratio = 1,
              axis.ticks=element_blank(),
              axis.text=element_blank()))
dev.off()

png.date(sprintf("%s/%s_UMAP_Sample_split.png", sample, sample), width=20, height=20)
print(DimPlot(obj, group.by="Sample", split.by="Sample", raster=FALSE, ncol=5))
dev.off()

select_genes <- c("AGR2", "EPCAM", "CFTR", "KRT19", "KRT17", "IL2RG", "SOX9",
                  "MUC5AC", "MUC5B", "S100P", "S100A10", "FOXQ1", "ONECUT2",
                  "IL18", "IL1RN", "IL1A", "TIMP1", "OIT1")
for (gene in select_genes) {
  png.date(sprintf("%s/%s_UMAP_%s.png", sample, sample, gene))
  print(FeaturePlot(obj, gene, raster=FALSE) + 
          theme(aspect.ratio = 1,
                axis.ticks=element_blank(),
                axis.text=element_blank()))
  dev.off()
}

obj_panc <- readRDS("merged_panc.rds")
obj_panc <- AddMetaData(obj_panc, obj[["Celltype_final"]])
obj_sub_exo_endo <- subset(obj_panc, Celltype_final %in% c("Duct1", "Duct2", "Acinar", "Endocrine"))
pdf.date(sprintf("%s/%s_UMAP_celltype_final_subexoendo.pdf", sample, sample), width=10, height=7)
print(DimPlot(obj_sub_exo_endo, group.by="Celltype_final", cols=cell_plot_colours, raster=FALSE) + 
        theme(aspect.ratio = 1,
              axis.ticks=element_blank(),
              axis.text=element_blank()))

dev.off()
pdf.date(sprintf("%s/%s_UMAP_celltype_final_panc.pdf", sample, sample), width=10, height=7)
print(DimPlot(obj_panc, group.by="Celltype_final", cols=cell_plot_colours, raster=FALSE) + 
        theme(aspect.ratio = 1,
              axis.ticks=element_blank(),
              axis.text=element_blank()))

dev.off()
obj_panc <- AddMetaData(obj_panc, obj[["Celltype_final_no_neuroendocrine"]])
png.date(sprintf("%s/%s_UMAP_celltype_final_no_neuroendocrine_panc.png", sample, sample), width=10, height=7)
print(DimPlot(obj_panc, group.by="Celltype_final_no_neuroendocrine", cols=cell_plot_colours, raster=FALSE) + 
        theme(aspect.ratio = 1,
              axis.ticks=element_blank(),
              axis.text=element_blank()))

dev.off()

pdf.date(sprintf("%s/%s_UMAP_subtype_score_subexoendo.pdf", sample, sample), width=10, height=7)
for (score in c("Basal1_score", "Basal2_score", "Classical1_score", "Classical2_score")) {
  print(FeaturePlot(obj_sub_exo_endo, score, raster=FALSE) + 
          scale_color_gradient(low="yellow", high="blue") +
          theme(aspect.ratio = 1,
                axis.ticks=element_blank(),
                axis.text=element_blank()))
}
dev.off()
pdf.date(sprintf("%s/%s_UMAP_subtype_score_subexoendo_samescale.pdf", sample, sample), width=10, height=7)
for (score in c("Basal1_score", "Basal2_score", "Classical1_score", "Classical2_score")) {
  print(FeaturePlot(obj_sub_exo_endo, score, raster=FALSE) + 
          scale_color_gradient(limits=c(-0.15, 0.15), low="yellow", high="blue") +
          theme(aspect.ratio = 1,
                axis.ticks=element_blank(),
                axis.text=element_blank()))
}
dev.off()


#### Plotting: Diversity grouped plots ####

subtype_data <- as.data.table(obj[[c("Sample", "subtype", "Celltype_final")]])
subtype_data <- subtype_data[Celltype_final == "Malignant", ]
melted_data <- subtype_data[, .(count=.N), by=.(Sample, subtype)]

diversity_data <- melted_data
diversity_data[, total:=sum(count), by=Sample]
diversity_data <- diversity_data[total > 50]
diversity_data_short <- dcast(diversity_data, Sample + total ~ subtype, value.var="count")

subtype_names <- c("Basal1", "Basal2", "Classical1", "Classical2", "Mixed")
diversity_data_short[, (subtype_names) := lapply(.SD, nafill, fill=0), .SDcols=subtype_names]
diversity_data_short$Injury <- diversity_data_short$Classical1 + diversity_data_short$Classical2
diversity_data_short$Invasive <- diversity_data_short$Basal1 + diversity_data_short$Basal2

diversity_data_short$Injury_prop <- diversity_data_short$Injury / diversity_data_short$total
diversity_data_short$Invasive_prop <- diversity_data_short$Invasive / diversity_data_short$total

prop_mat <- data.frame(diversity_data_short[, .(Injury_prop, Invasive_prop)], row.names=diversity_data_short$Sample)
set.seed(0)
km_out <- kmeans(prop_mat, centers = 4)

write.csv(as.data.frame(km_out$cluster), sprintf("%s/%s_injury_v_invasive_grouping.csv", sample, sample))

cluster_named <- c("high", "low_classic", "moderate", "low_basal")[km_out$cluster]
names(cluster_named) <- names(km_out$cluster)
cluster_named <- factor(cluster_named, levels=c("low_basal", "high", "moderate", "low_classic"))

diversity_data_short$cluster <- as.character(km_out$cluster[diversity_data_short$Sample])
diversity_data_short$cluster_named <- cluster_named[diversity_data_short$Sample]

pdf.date(sprintf("%s/%s_scatter_prop_injury_vs_invasive.pdf", sample, sample), width=5, height=5)
ggplot(diversity_data_short, aes(x=Injury_prop, y=Invasive_prop, colour=cluster_named)) + 
  geom_point() +
  theme_bw() + 
  scale_colour_manual(values = c("#C7652C", "grey50", "#768695", "#6096C8")) +
  theme(panel.grid = element_blank()) + 
  coord_fixed()
dev.off()

cluster_melted_data <- merge(melted_data, diversity_data_short)
cluster_melted_data$Sample <- factor(cluster_melted_data$Sample, levels=names(sort(km_out$cluster)))
cluster_melted_data$subtype <- factor(cluster_melted_data$subtype, levels=names(singler_subtype_colors))

barplot_count <- ggplot(cluster_melted_data, aes(x=Sample, y=count, fill=subtype)) + 
  geom_col() +
  scale_fill_manual(values=singler_subtype_colors) +
  scale_y_continuous(expand=expansion(mult=0.01)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_grid(cols=vars(cluster_named), scales="free", space="free")
barplot_prop <- ggplot(cluster_melted_data, aes(x=Sample, y=count, fill=subtype)) + 
  geom_col(position="fill") +
  scale_fill_manual(values=singler_subtype_colors) +
  scale_y_continuous(expand=expansion(mult=0.01)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_grid(cols=vars(cluster_named), scales="free", space="free")
pdf.date(sprintf("%s/%s_barplot_subtype_kmeans_cluster_ordered.pdf", sample, sample), width=11, height=4)
print(barplot_count)
print(barplot_prop)
dev.off()


# Tile plot
tile_data <- unique(obj@meta.data[, c("Sample", "type")])
setDT(tile_data)
tile_data$cluster <- km_out$cluster[tile_data$Sample]
tile_data$cluster_named <- cluster_named[tile_data$Sample]

tile_data$Sample <- factor(as.character(tile_data$Sample), levels=levels(cluster_melted_data$Sample))
tile_data <- tile_data[!is.na(Sample)]
colnames(tile_data) <- c("Sample",  "Type", "cluster", "cluster_named")

tile_data_melt <- melt(tile_data, id.vars=c("Sample", "cluster", "cluster_named"))

tileplot <- ggplot(tile_data_melt, aes(x=Sample, y=variable)) + 
  geom_tile(aes(fill=value, width=0.8, height=0.8), colour="black", size=0.25) + 
  scale_fill_manual(values=c(primary="#ccebc5", metastatic="#fccde5"),
                    na.value="lightgrey") +
  labs(x="", y="") +
  scale_x_discrete(position="bottom") +
  theme_grey() +
  theme(plot.background=element_blank(), 
        panel.background=element_blank(), 
        panel.border=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + 
  facet_grid(cols=vars(cluster_named), scales="free", space="free")
pdf.date(sprintf("%s/%s_tileplot_subtype_kmeans_cluster_ordered.pdf", sample, sample), width=13, height=1.8)
print(tileplot)
dev.off()

pdf.date(sprintf("%s/%s_barplot_subtype_kmeans_cluster_ordered_with_tileplot.pdf", sample, sample), width=20, height=11)
barplot_count / tileplot + plot_layout(heights = c(6, 1))
barplot_prop / tileplot  + plot_layout(heights = c(6, 1))
dev.off()

summary_dt <- fread(sprintf("%s/%s_count_subtype_summary.csv", sample, sample))
tile_data_pseudobulk <- merge(tile_data, summary_dt[, .(Sample, singler_by_sample)], by="Sample")

tile_data_pseudobulk_melt <- melt(tile_data_pseudobulk, id.vars=c("Sample", "cluster", "cluster_named"))

tileplot <- ggplot(tile_data_pseudobulk_melt, aes(x=Sample, y=variable)) + 
  geom_tile(aes(fill=value, width=0.8, height=0.8), colour="black", size=0.25) + 
  scale_fill_manual(values=c(primary="#ccebc5", metastatic="#fccde5",
                             Basal1="#814C42", Basal2="#FF7311", hybrid="grey40", Classical1="#1C6CAB", Classical2="#A4C0E5", Mixed="grey40"),
                    na.value="lightgrey") +
  labs(x="", y="") +
  scale_x_discrete(position="bottom") +
  theme_grey() +
  theme(plot.background=element_blank(), 
        panel.background=element_blank(), 
        panel.border=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + 
  facet_grid(cols=vars(cluster_named), scales="free", space="free")
pdf.date(sprintf("%s/%s_tileplot_subtype_kmeans_cluster_ordered.pdf", sample, sample), width=13, height=4)
print(tileplot)
dev.off()


# Barplot with all samples 

subtype_data_all <- as.data.table(obj[[c("Sample", subtype_metadata_col, "Celltype_final")]])
subtype_data_all <- subtype_data_all[Celltype_final %in% c("Malignant", "Neuroendocrine"), ]
melted_data_all <- subtype_data_all[, .(count=.N), by=.(Sample, subtype=get(subtype_metadata_col))]

diversity_data_all <- melted_data_all
diversity_data_all[, total:=sum(count), by=Sample]

diversity_data_all_short <- dcast(diversity_data_all, Sample + total ~ subtype, value.var="count")

subtype_names <- c("Basal1", "Basal2", "Classical1", "Classical2", "Mixed")
diversity_data_all_short[, (subtype_names) := lapply(.SD, nafill, fill=0), .SDcols=subtype_names]

diversity_data_all_short$Injury <- diversity_data_all_short$Classical1 + diversity_data_all_short$Classical2
diversity_data_all_short$Invasive <- diversity_data_all_short$Basal1 + diversity_data_all_short$Basal2

diversity_data_all_short$Injury_prop <- diversity_data_all_short$Injury / diversity_data_all_short$total
diversity_data_all_short$Invasive_prop <- diversity_data_all_short$Invasive / diversity_data_all_short$total

diversity_data_all_short$cluster <- as.character(km_out$cluster[diversity_data_all_short$Sample])
diversity_data_all_short$cluster_named <- cluster_named[diversity_data_all_short$Sample]

cluster_melted_data_all <- merge(melted_data_all, diversity_data_all_short)
cluster_melted_data_all$subtype <- factor(cluster_melted_data_all$subtype, levels=names(singler_subtype_colors))

barplot_count <- ggplot(cluster_melted_data_all, aes(x=Sample, y=count, fill=subtype)) + 
  geom_col() +
  scale_fill_manual(values=singler_subtype_colors) +
  scale_y_continuous(expand=expansion(mult=0.01)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid = element_blank(),
        strip.background = element_rect(fill="white")) + 
  facet_grid(cols=vars(cluster_named), scales="free", space="free")
barplot_prop <- ggplot(cluster_melted_data_all, aes(x=Sample, y=count, fill=subtype)) + 
  geom_col(position="fill") +
  scale_fill_manual(values=singler_subtype_colors) +
  scale_y_continuous(expand=expansion(mult=0.01)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid = element_blank(),
        strip.background = element_rect(fill="white")) + 
  facet_grid(cols=vars(cluster_named), scales="free", space="free")
pdf.date(sprintf("%s/%s_barplot_subtype_kmeans_cluster_ordered_all.pdf", sample, sample), width=11, height=4)
print(barplot_count)
print(barplot_prop)
dev.off()


# Incorporating in numbat (data from numbat_post_analysis.R) + met/ primary data 
clone_and_phenotype_diversity_data <- fread(sprintf("%s/%s_clones_per_sample_data.csv", sample, sample))
clone_and_phenotype_diversity_data$num_malignant_clones <- clone_and_phenotype_diversity_data$num_clone - 1
clone_and_phenotype_diversity_data$injury_v_invasive_grouping <- km_out$cluster[clone_and_phenotype_diversity_data$Sample]
clone_and_phenotype_diversity_data[Sample == "89221", "num_malignant_clones"] <- 1

fwrite(clone_and_phenotype_diversity_data, datestamp("merged_clones_and_injury_v_invasive_grouping_per_sample_data.csv"))

distinct_cnv_paths <- list.files(path = "numbat",
                                 pattern = "sample_.*_NB_distinct_cnv.csv",
                                 full.names = TRUE)
distinct_cnv <- lapply(distinct_cnv_paths, fread)
names(distinct_cnv) <- sapply(distinct_cnv_paths, function(path) {
  gsub(".*/sample_(.*)_NB_distinct_cnv.csv", "\\1", path)
})
distinct_cnv_counts <- sapply(distinct_cnv, nrow)
names(distinct_cnv_counts) <- names(distinct_cnv)

combined_data <- merge(tile_data, clone_and_phenotype_diversity_data, all = TRUE)
combined_data$distinct_cnv_counts <- distinct_cnv_counts[as.character(combined_data$Sample)]

combined_data_filtered <- combined_data[!is.na(cluster) & !is.na(num_clone)]
combined_data_filtered$type <- factor(as.character(combined_data_filtered$type), levels=c("primary", "metastatic"))

pdf.date(sprintf("%s/%s_boxplot_subtype_kmeans_cluster_distinct_cnv_num_malignant_clones.pdf", sample, sample), width=4, height=4)
set.seed(0)
ggplot(combined_data_filtered, aes(x=cluster_named, y=distinct_cnv_counts)) +
  stat_boxplot(geom = "errorbar", width=0.4) +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.2, alpha = 0.5) +
  stat_compare_means() +
  theme_bw() +
  theme(panel.grid = element_blank(), aspect.ratio = 1)
ggplot(combined_data_filtered, aes(x=type, y=distinct_cnv_counts)) +
  stat_boxplot(geom = "errorbar", width=0.4) +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.2, alpha = 0.5) +
  stat_compare_means() +
  theme_bw() +
  theme(panel.grid = element_blank(), aspect.ratio = 1)
ggplot(combined_data_filtered, aes(x=cluster_named, y=num_malignant_clones)) +
  stat_boxplot(geom = "errorbar", width=0.4) +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.2, alpha = 0.5) +
  scale_y_continuous(breaks = 0:9) +
  stat_compare_means() +
  theme_bw() +
  theme(panel.grid = element_blank(), aspect.ratio = 1)
ggplot(combined_data_filtered, aes(x=type, y=num_malignant_clones)) +
  stat_boxplot(geom = "errorbar", width=0.4) +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.2, alpha = 0.5) +
  scale_y_continuous(breaks = 0:9) +
  stat_compare_means() +
  theme_bw() +
  theme(panel.grid = element_blank(), aspect.ratio = 1)
dev.off()

combined_data_filtered$num_malignant_clones_factor <- factor(combined_data_filtered$num_malignant_clones)
combined_data_filtered$num_malignant_clones_factor_reduced <- ifelse(combined_data_filtered$num_malignant_clones >= 3, "3+", as.character(combined_data_filtered$num_malignant_clones))
combined_data_filtered$num_malignant_clones_factor_reduced <- factor(as.character(combined_data_filtered$num_malignant_clones_factor_reduced))

out_all <- fisher.test(table(combined_data_filtered[, c("Type", "num_malignant_clones_factor")]))
out_reduced <- fisher.test(table(combined_data_filtered[, c("Type", "num_malignant_clones_factor_reduced")]))
pdf.date(sprintf("%s/%s_barplot_type_num_malignant_clones.pdf", sample, sample), width=8, height=2.5)
ggplot(combined_data_filtered, aes(x=type, fill=num_malignant_clones_factor)) + 
  geom_bar(position="fill") + 
  theme_bw() + 
  labs(caption = sprintf("Fisher's exact test, p=%f", out_all$p.value)) + 
  theme(panel.grid = element_blank(), aspect.ratio = 1)
ggplot(combined_data_filtered, aes(x=type, fill=num_malignant_clones_factor_reduced)) + 
  geom_bar(position="fill") + 
  scale_fill_manual(values=c("1"="#e31a1c", "2" = "#33a02c", "3+" = "#1f78b4")) + 
  theme_bw() + 
  labs(caption = sprintf("Fisher's exact test, p=%f", out_reduced$p.value)) + 
  theme(panel.grid = element_blank(), aspect.ratio = 1)
dev.off()

pdf.date(sprintf("%s/%s_violinplot_subtype_kmeans_cluster_distinct_cnv_num_malignant_clones.pdf", sample, sample), width=4, height=4)
set.seed(0)
ggplot(combined_data_filtered, aes(x=cluster_named, y=distinct_cnv_counts)) + 
  geom_violin() + 
  geom_jitter(height = 0, width = 0.2, alpha = 0.5) + 
  stat_compare_means() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), aspect.ratio = 1)
ggplot(combined_data_filtered, aes(x=type, y=distinct_cnv_counts)) + 
  geom_violin() + 
  geom_jitter(height = 0, width = 0.2, alpha = 0.5) + 
  stat_compare_means() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), aspect.ratio = 1)
ggplot(combined_data_filtered, aes(x=cluster_named, y=num_malignant_clones)) + 
  geom_violin() + 
  geom_jitter(height = 0, width = 0.2, alpha = 0.5) + 
  scale_y_continuous(breaks = 0:9) + 
  stat_compare_means() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), aspect.ratio = 1)
ggplot(combined_data_filtered, aes(x=type, y=num_malignant_clones)) + 
  geom_violin() + 
  geom_jitter(height = 0, width = 0.2, alpha = 0.5) + 
  scale_y_continuous(breaks = 0:9) + 
  stat_compare_means() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), aspect.ratio = 1)
dev.off()

pdf.date(sprintf("%s/%s_barplot_subtype_kmeans_cluster_num_malignant_clones.pdf", sample, sample), width=4, height=4)
ggplot(combined_data_filtered, aes(x=num_malignant_clones)) + 
  geom_bar() + 
  theme_bw() + 
  scale_x_continuous(breaks = 0:9) + 
  facet_wrap(~cluster_named, ncol=1)
ggplot(combined_data_filtered, aes(x=num_malignant_clones)) + 
  geom_bar() + 
  theme_bw() + 
  scale_x_continuous(breaks = 0:9) + 
  facet_wrap(~type, ncol=1)
dev.off()

combined_data_filtered$clone_group <- ifelse(combined_data_filtered$num_malignant_clones > 1, "2+", "1")
pdf.date(sprintf("%s/%s_barplot_subtype_kmeans_cluster_clone_group.pdf", sample, sample), width=4, height=4)
ggplot(combined_data_filtered, aes(x=clone_group, fill = cluster_named )) + 
  geom_bar() + 
  theme_bw() 
ggplot(combined_data_filtered, aes(x=clone_group, fill = type)) + 
  geom_bar() + 
  theme_bw()
ggplot(combined_data_filtered, aes(fill=clone_group, x = cluster_named )) + 
  geom_bar() + 
  scale_fill_manual(values = c("1"="grey50", "2+" = "blue")) + 
  theme_bw() 
ggplot(combined_data_filtered, aes(fill=clone_group, x = type)) + 
  geom_bar() + 
  scale_fill_manual(values = c("1"="grey50", "2+" = "blue")) + 
  theme_bw()
dev.off()

pdf.date(sprintf("%s/%s_mosaic_subtype_kmeans_cluster_clone_group.pdf", sample, sample), width=4, height=4)
ggplot(combined_data_filtered) +
  geom_mosaic(aes(x = product(cluster_named), fill = clone_group)) +
  geom_mosaic_text(aes(x = product(cluster_named), fill = clone_group, label = after_stat(.wt))) + 
  scale_fill_manual(values = c("1"="grey50", "2+" = "blue")) +
  labs(caption=paste0("Fisher's Exact Test, p=", 
                      round(fisher.test(table(combined_data_filtered[, c("cluster_named", "clone_group")]))$p.value, 4))) + 
  theme_classic() + 
  theme(aspect.ratio = 1)
ggplot(combined_data_filtered) +
  geom_mosaic(aes(x = product(type), fill = clone_group)) +
  geom_mosaic_text(aes(x = product(type), fill = clone_group, label = after_stat(.wt))) + 
  scale_fill_manual(values = c("1"="grey50", "2+" = "blue")) + 
  labs(caption=paste0("Fisher's Exact Test, p=", 
                      round(fisher.test(table(combined_data_filtered[, c("type", "clone_group")]))$p.value, 4))) + 
  theme_classic() + 
  theme(aspect.ratio = 1)
dev.off()

tile_data$low_classic_other <- ifelse(tile_data$cluster_named == "low_classic", "low_classic", "other")
pdf.date(sprintf("%s/%s_mosaic_subgroup_type.pdf", sample, sample), width=4, height=4)
ggplot(tile_data) +
  geom_mosaic(aes(x = product(cluster_named), fill = Type)) +
  geom_mosaic_text(aes(x = product(cluster_named), fill = Type, label = after_stat(.wt))) + 
  scale_fill_manual(values = c(primary="#ccebc5", metastatic="#fccde5")) +
  labs(caption=paste0("Fisher's Exact Test, p=", 
                      round(fisher.test(table(tile_data[, c("cluster_named", "Type")]))$p.value, 4))) + 
  theme_classic() + 
  theme(aspect.ratio = 1)
ggplot(tile_data) +
  geom_mosaic(aes(x = product(low_classic_other), fill = Type)) +
  geom_mosaic_text(aes(x = product(low_classic_other), fill = Type, label = after_stat(.wt))) + 
  scale_fill_manual(values = c(primary="#ccebc5", metastatic="#fccde5")) +
  labs(caption=paste0("Fisher's Exact Test, p=", 
                      round(fisher.test(table(tile_data[, c("low_classic_other", "Type")]))$p.value, 4))) + 
  theme_classic() + 
  theme(aspect.ratio = 1)
dev.off()

pdf.date(sprintf("%s/%s_scatterplot_subtype_kmeans_cluster_distinct_cnv.pdf", sample, sample), width=6, height=4)
ggplot(combined_data_filtered, aes(x=num_malignant_clones, y=distinct_cnv_counts)) + 
  geom_jitter(height = 0, width=0.1) + 
  scale_x_continuous(breaks=1:10)
ggplot(combined_data_filtered, aes(x=num_malignant_clones, y=distinct_cnv_counts, colour=Tissue)) + 
  geom_jitter(height = 0, width=0.1) + 
  scale_x_continuous(breaks=1:10)
ggplot(combined_data_filtered, aes(x=num_malignant_clones, y=distinct_cnv_counts, colour=cluster_named)) + 
  geom_jitter(height = 0, width=0.1) + 
  scale_x_continuous(breaks=1:10)
dev.off()


#### Plotting: Prop cycling/emt by subtype ####

obj_mal <- subset(obj, Celltype_final == "Malignant")

# Barplots
ccs_name <- "Connor_CC"
proliferating_cutoff <- mean(obj@meta.data[, ccs_name]) + sd(obj@meta.data[, ccs_name])
obj$proliferating <- ifelse(obj@meta.data[, ccs_name] > proliferating_cutoff, "proliferating", "non-proliferating")
obj$EMT_status <- ifelse(obj$Hallmark_EMT > 0, "EMT", "non-EMT")

obj_mal <- subset(obj, Celltype_final == "Malignant")

bar_data <- obj_mal[[c("subtype", "proliferating")]]
bar_data$proliferating <- factor(bar_data$proliferating, levels=c("non-proliferating", "proliferating"))
bar_data$subtype <- factor(bar_data$subtype , levels=c("Basal1", "Basal2", "Mixed", "Classical1", "Classical2"))
pdf.date(sprintf("%s/%s_barplot_proliferating_by_subtype_horizontal.pdf", sample, sample), width=6, height=3)
ggplot(bar_data, aes(x=subtype, fill=proliferating)) + 
  geom_bar() + 
  scale_fill_manual(values=c("non-proliferating"="lightgrey", "proliferating"="#0BC182")) + 
  coord_flip() + 
  theme_bw()
ggplot(bar_data, aes(x=subtype, fill=proliferating)) + 
  geom_bar(position="fill") + 
  scale_fill_manual(values=c("non-proliferating"="lightgrey", "proliferating"="#0BC182")) + 
  ylab("proportion") + 
  coord_flip() + 
  theme_bw(base_size=14)
dev.off()
pdf.date(sprintf("%s/%s_barplot_proliferating_by_subtype.pdf", sample, sample), width=4, height=5)
ggplot(bar_data, aes(x=subtype, fill=proliferating)) + 
  geom_bar() + 
  scale_fill_manual(values=c("non-proliferating"="lightgrey", "proliferating"="#0BC182")) + 
  theme_bw(base_size=14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggplot(bar_data, aes(x=subtype, fill=proliferating)) + 
  geom_bar(position="fill") + 
  scale_fill_manual(values=c("non-proliferating"="lightgrey", "proliferating"="#0BC182")) + 
  ylab("proportion") + 
  scale_y_continuous(expand = expansion(add = c(0.01, 0.01))) + 
  theme_bw(base_size=14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

bar_data_prolif <- bar_data[bar_data$proliferating == "proliferating", ]
pdf.date(sprintf("%s/%s_barplot_proliferating_by_subtype_no_fill.pdf", sample, sample), width=4, height=5)
ggplot(bar_data_prolif, aes(x=subtype, fill=proliferating)) + 
  geom_bar(fill="#0BC182") + 
  scale_fill_manual(values=c("non-proliferating"="lightgrey", "proliferating"="#0BC182")) + 
  ylab("proportion") + 
  scale_y_continuous(expand = expansion(add = c(0.01, 0.01))) + 
  theme_bw(base_size=14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

bar_data_dt <- as.data.table(bar_data)
bar_data_counts <- merge(bar_data_dt[, .(count = .N), by = .(subtype, proliferating)], bar_data_dt[, .(total_count = .N), by = subtype])
bar_data_counts[, prop := count/ total_count]
bar_data_counts_prolif <- bar_data_counts[proliferating == "proliferating"]
pdf.date(sprintf("%s/%s_barplot_proliferating_by_subtype_prop.pdf", sample, sample), width=4, height=5)
ggplot(bar_data_counts_prolif, aes(x=subtype, y=prop, fill=subtype)) + 
  geom_col() + 
  scale_fill_manual(values=singler_subtype_colors) + 
  ylab("proportion") + 
  scale_y_continuous(expand = expansion(add = c(0.001, 0.01))) + 
  theme_bw(base_size=14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
  )
dev.off()
bar_data_dt$Basal2vrest <- ifelse(bar_data_dt$subtype == "Basal2", "Basal2", "other")
bar_data_mat <- table(bar_data_dt[, c("Basal2vrest", "proliferating")])
print(fisher.test(bar_data_mat))
print(fisher.test(bar_data_mat)$p.value)

bar_data <- obj_mal[[c("subtype", "EMT_status")]]
bar_data$EMT_status <- factor(bar_data$EMT_status, levels=c("non-EMT", "EMT"))
bar_data$subtype <- factor(bar_data$subtype , levels=c("Basal1", "Basal2", "Mixed", "Classical1", "Classical2"))
pdf.date(sprintf("%s/%s_barplot_emtstatus_by_subtype.pdf", sample, sample), width=6, height=3)
ggplot(bar_data, aes(x=subtype, fill=EMT_status)) + 
  geom_bar() + 
  scale_fill_manual(values=c("non-EMT"="lightgrey", "EMT"="#BBA620")) + 
  coord_flip() + 
  theme_bw()
ggplot(bar_data, aes(x=subtype, fill=EMT_status)) +   
  geom_bar(position="fill") + 
  scale_fill_manual(values=c("non-EMT"="lightgrey", "EMT"="#BBA620")) + 
  ylab("proportion") + 
  coord_flip() + 
  theme_bw()
dev.off()


#### Plotting: Program correlation #####

stat_sets <- list(
  emt=c(colnames(obj@meta.data)[grep("EMT", colnames(obj@meta.data))],
        "Basal1_score", "Basal2_score",  "Classical1_score", "Classical2_score"),
  cellcycle=c("Connor_CC", "Seurat_CC_G2M", "Seurat_CC_S", "Gavish_MP01_CellCycleG2M", "Gavish_MP02_CellCycleG1S",
              "Basal1_score", "Basal2_score",  "Classical1_score", "Classical2_score"),
  hwang=c(colnames(obj@meta.data[startsWith(colnames(obj@meta.data), "Hwang_")]), 
          "Basal1_score", "Basal2_score",  "Classical1_score", "Classical2_score")
)

for (stat_set_name in names(stat_sets)) {
  
  stat_set <- stat_sets[[stat_set_name]]
  
  mat_df <- obj[[c("Celltype_final", stat_set)]]
  mat <- as.matrix(mat_df[mat_df$Celltype_final=="Malignant", -1])
  corr_data <- rcorr(mat)
  
  corr_data_P_cleaned <- corr_data$P
  corr_data_P_cleaned[is.na(corr_data_P_cleaned)] <- 0
  
  pdf.date(sprintf("%s/%s_corrplot_mal_%s_vs_subtype.pdf", sample, sample, stat_set_name), width=10, height=10)
  corrplot(corr_data$r, order="hclust", col=rev(colorRampPalette(colors=brewer.pal(11, "PiYG"))(200)),
           p.mat=corr_data_P_cleaned, insig="blank", tl.col = "black")
  corrplot(corr_data$r, method="number", order="hclust", col=rev(colorRampPalette(colors=brewer.pal(11, "PiYG"))(200)),
           p.mat=corr_data_P_cleaned, insig="blank", tl.col = "black")
  dev.off()
  
  pdf.date(sprintf("%s/%s_corrplot_mal_%s_vs_subtype_small.pdf", sample, sample, stat_set_name), width=5, height=5)
  corrplot(corr_data$r, order="hclust", col=rev(colorRampPalette(colors=brewer.pal(11, "PiYG"))(200)),
           p.mat=corr_data_P_cleaned, insig="blank", tl.col = "black")
  corrplot(corr_data$r, method="number", order="hclust", col=rev(colorRampPalette(colors=brewer.pal(11, "PiYG"))(200)),
           p.mat=corr_data_P_cleaned, insig="blank", tl.col = "black")
  dev.off()
  
  pdf.date(sprintf("%s/%s_corrplot_mal_%s_vs_subtype_large.pdf", sample, sample, stat_set_name), width=20, height=20)
  corrplot(corr_data$r, order="hclust", col=rev(colorRampPalette(colors=brewer.pal(11, "PiYG"))(200)),
           p.mat=corr_data_P_cleaned, insig="blank", tl.col = "black")
  corrplot(corr_data$r, method="number", order="hclust", col=rev(colorRampPalette(colors=brewer.pal(11, "PiYG"))(200)),
           p.mat=corr_data_P_cleaned, insig="blank", tl.col = "black")
  dev.off()
  
}

# Everything, sorted by correlation with the important subtypes
pdf.date(sprintf("%s/%s_corrplot_mal_everything_vs_subtype_ordered.pdf", sample, sample), width=20, height=20)
for (subtype_score in c("Basal1_score", "Basal2_score",  "Classical1_score", "Classical2_score")) {
  corr_order <- names(sort(corr_data$r[, subtype_score], decreasing=TRUE))
  print(corrplot(corr_data$r[corr_order, corr_order], order="original", col=rev(colorRampPalette(colors=brewer.pal(11, "PiYG"))(200)),
                 p.mat=corr_data_P_cleaned[corr_order, corr_order], insig="blank", tl.col = "black"))
}
dev.off()


#### Plotting: Subtype score continuum ####

library(ComplexHeatmap)

subtype_names <- c("Basal1_score", "Basal2_score",  "Classical1_score", "Classical2_score")

top_anno_names_cont <- c("TP63", "KRT5", "KRT6A", "S100A2", "OVOL1", 
                         "CDX2", "LAMC2", "VIM", "FOXJ1", "SNAI2", "CDH2", 
                         "CDH1", "EPCAM", "GATA6", 
                         "Connor_CC",
                         "Hallmark_EMT",
                         "Dilly_pEMT",
                         "Gavish_MP12_EMTI",
                         "Gavish_MP13_EMTII",
                         "Gavish_MP14_EMTIII",
                         "Gavish_MP15_EMTIV",
                         "Puram_pEMT",
                         "Klomp_PDAC_siKRAS_KRASi_iKras_UP",
                         "Singh_KRASDependency_addiction",
                         "Raghavan_scBasal",
                         "Raghavan_scClassical",
                         "Raghavan_IC")

top_anno_names <- c("subtype",  "type", "Sample", top_anno_names_cont)
tag <- sample
subtype_dataset <- FetchData(obj, c("Celltype_final", subtype_names, top_anno_names))
setDT(subtype_dataset, keep.rownames="cell")
subtype_dataset$subtype <- factor(as.character(subtype_dataset$subtype), levels=rev(names(singler_subtype_colors)))

scale_cols <- c(subtype_names, top_anno_names_cont)
subtype_dataset_scaled <- apply(as.matrix(subtype_dataset[, ..scale_cols]), 2, function(x) {(x - min(x))/(max(x) - min(x))})
colnames(subtype_dataset_scaled) <- paste0(colnames(subtype_dataset_scaled), "_scaled")

subtype_dataset <- cbind(subtype_dataset, subtype_dataset_scaled)

col_fun_scaled <- circlize::colorRamp2(c(0, 1), c("#f7f7f7", "#b2182b"))
anno_fun_scaled <- circlize::colorRamp2(c(0, 1), c("#f7f7f7", "purple"))

col_list_scaled <- c(list(singler_subtype_colors), rep(list(anno_fun_scaled), length(top_anno_names_cont)))
names(col_list_scaled) <- c("subtype", paste0(top_anno_names_cont, "_scaled"))

max_lim <- max(apply(subtype_dataset[, ..subtype_names], 2, max))
min_lim <- min(apply(subtype_dataset[, ..subtype_names], 2, min))
col_fun <- circlize::colorRamp2(c(min_lim, 0, max_lim), c("royalblue",  "#f7f7f7", "#b2182b"))

max_lim <- max(apply(subtype_dataset[, ..top_anno_names_cont], 2, max))
min_lim <- min(apply(subtype_dataset[, ..top_anno_names_cont], 2, min))
anno_fun <- circlize::colorRamp2(c(min_lim, 0, 1, max_lim), c("#d1e5f0",  "#f7f7f7", "#fddbc7", "#b2182b"))

col_list <- c(list(singler_subtype_colors), list(type_colours), list(sample_colours), rep(list(anno_fun), length(top_anno_names_cont)))
names(col_list) <- top_anno_names

# PCA for ordering
pca_out <- princomp(subtype_dataset[, ..subtype_names])
subtype_dataset$pc1 <- pca_out$scores[, 1]

# Quick density plot
pdf.date(sprintf("%s/%s_density_subtype_pc1.pdf", sample, tag), height=5, width=12)
ggplot(subtype_dataset, aes(x=pc1, fill=subtype)) +
  scale_fill_manual(values=singler_subtype_colors) +
  geom_density()
ggplot(subtype_dataset, aes(x=pc1, fill=Sample)) +
  scale_fill_manual(values=sample_colours) +
  geom_density()
dev.off()

subtype_dataset_sub <- subtype_dataset[Celltype_final == "Malignant", ..subtype_names]
subtype_order_score <- order(subtype_dataset[Celltype_final == "Malignant", pc1])

combined_plot_anno_names <- c(top_anno_names, "Sample")
anno_data_sub <- subtype_dataset[Celltype_final == "Malignant", ..combined_plot_anno_names][subtype_order_score, ]
subtype_dataset_sub_mat <- t(as.matrix(subtype_dataset_sub))[, subtype_order_score]

png.date(sprintf("%s/%s_heatmap_subtype_continuous_combined.png", sample, tag), height=12, width=20)
print(Heatmap(subtype_dataset_sub_mat, 
              col = col_fun,
              name = 'Merged',
              top_annotation = HeatmapAnnotation(df=anno_data_sub, col=col_list, show_legend=c(rep(TRUE, length(col_list)), FALSE)),
              cluster_columns = FALSE, 
              cluster_rows = FALSE,
              heatmap_legend_param = list(direction="horizontal")))
dev.off()

combined_plot_anno_names <- c("subtype", paste0(top_anno_names_cont, "_scaled"), "Sample")
anno_data_sub <- subtype_dataset[Celltype_final == "Malignant", ..combined_plot_anno_names][subtype_order_score, ]
subtype_dataset_sub_mat <- t(as.matrix(subtype_dataset_sub))[, subtype_order_score]

png.date(sprintf("%s/%s_heatmap_subtype_continuous_combined_scaled.png", sample, tag), height=12, width=20)
print(Heatmap(subtype_dataset_sub_mat, 
              col = col_fun,
              name = 'Merged',
              top_annotation = HeatmapAnnotation(df=anno_data_sub, col=col_list_scaled, show_legend=c(rep(TRUE, length(col_list)), FALSE)),
              cluster_columns = FALSE, 
              cluster_rows = FALSE))
dev.off()

# Barplot summary
nbins <- 100
subtype_dataset_ordered <- subtype_dataset[Celltype_final == "Malignant", ][order(pc1)]
subtype_dataset_ordered$even_bin <- sort(rep(1:nbins, length.out=nrow(subtype_dataset_ordered)))

even_bin_summary <- subtype_dataset_ordered[, lapply(.SD, mean, na.rm=TRUE), by=even_bin, .SDcols=c(subtype_names, top_anno_names_cont, paste0(top_anno_names_cont, "_scaled"))]
even_bin_counts <- as.data.frame.matrix(table(subtype_dataset_ordered[, .(even_bin, subtype)]))
even_bin_summary_counts <- cbind(even_bin_summary, even_bin_counts)

subtype_dataset_ordered$post_neoadjuvant <- ifelse(subtype_dataset_ordered$post_neoadjuvant == 1, "yes", "no")

pdf.date(sprintf("%s/%s_barplot_subtype_continuous_%sbins.pdf", sample, tag, nbins), height=5, width=12)
ggplot(subtype_dataset_ordered, aes(x=even_bin, fill=subtype)) + geom_bar() + scale_fill_manual(values=singler_subtype_colors) + 
  theme_bw() + scale_x_continuous(expand=c(0, 0), limits=c(1, nbins), breaks=c(1, 97, nbins), minor_breaks=c(1, seq(from=10, to=100, by=10)))
ggplot(subtype_dataset_ordered, aes(x=even_bin, fill=type)) + geom_bar() + scale_fill_manual(values=type_colours) + 
  theme_bw() + scale_x_continuous(expand=c(0, 0), limits=c(1, nbins), breaks=c(1, 97, nbins), minor_breaks=c(1, seq(from=10, to=100, by=10)))
ggplot(subtype_dataset_ordered, aes(x=even_bin, fill=post_neoadjuvant)) + geom_bar() +  
  theme_bw() + scale_x_continuous(expand=c(0, 0), limits=c(1, nbins), breaks=c(1, 97, nbins), minor_breaks=c(1, seq(from=10, to=100, by=10)))
ggplot(subtype_dataset_ordered, aes(x=even_bin, fill=Sample)) + geom_bar() + scale_fill_manual(values=sample_colours) + 
  theme_bw() + scale_x_continuous(expand=c(0, 0), limits=c(1, nbins), breaks=c(1, 97, nbins), minor_breaks=c(1, seq(from=10, to=100, by=10)))

for (feature in c(subtype_names, top_anno_names_cont, paste0(top_anno_names_cont, "_scaled"), names(singler_subtype_colors))) {
  print(ggplot(even_bin_summary_counts, aes(x=even_bin, y=.data[[feature]])) + geom_col() + 
          theme_bw() + scale_x_continuous(expand=c(0, 0), limits=c(1, nbins), breaks=c(1, 97, nbins), minor_breaks=c(1, seq(from=10, to=100, by=10))))
}
dev.off()

# Genes scaled line plot
line_style_groupings <- c(CDH1="early", EPCAM="early", GATA6="early", 
                          CDH2="middle", VIM="middle", LAMC2="middle", SNAI2="middle", FOXJ1="middle",
                          KRT5="late", KRT6A="late", S100A2="late", TP63="late", OVOL1="late", TP73="late")
line_style_groupings <- factor(line_style_groupings, levels=unique(line_style_groupings))
custom_line_paras <- c("even_bin", names(line_style_groupings))
even_bin_summary_melt <- melt(even_bin_summary_counts[, ..custom_line_paras], id.vars="even_bin", value.name="mean_score", variable.name="program")
even_bin_summary_melt$type <- line_style_groupings[as.character(even_bin_summary_melt$program)]

seed(1)
grouped_plot_gene_colors <- sample(scales::hue_pal()(length(names(line_style_groupings))))
names(grouped_plot_gene_colors) <- names(line_style_groupings)

grouped_plots <- lapply(levels(even_bin_summary_melt$type), function(x) {
  even_bin_summary_melt_sub <- even_bin_summary_melt[type==x]
  ggplot(even_bin_summary_melt_sub, aes(x=even_bin, y=mean_score, color=program)) + 
    geom_line() + 
    theme_bw() + 
    scale_color_manual(values=grouped_plot_gene_colors) + 
    scale_x_continuous(expand=c(0, 0), limits=c(1, nbins), breaks=c(1, 97, nbins), minor_breaks=c(1, seq(from=10, to=100, by=10)))
})

bin_scaled <- even_bin_summary_counts[, ..custom_line_paras]
bin_scaled <- cbind(bin_scaled[, 1], apply(as.matrix(bin_scaled[, 2:ncol(bin_scaled)]), 2, function(x) {(x - min(x))/(max(x) - min(x))}))
bin_scaled_melt <- melt(bin_scaled, id.vars="even_bin", value.name="scaled_mean_score", variable.name="program")
bin_scaled_melt$type <- line_style_groupings[as.character(bin_scaled_melt$program)]

grouped_plots_scaled <- lapply(levels(even_bin_summary_melt$type), function(x) {
  bin_scaled_melt_sub <- bin_scaled_melt[type==x]
  ggplot(bin_scaled_melt_sub, aes(x=even_bin, y=scaled_mean_score, color=program)) + 
    geom_line() + 
    theme_bw() + 
    scale_color_manual(values=grouped_plot_gene_colors) + 
    scale_y_continuous(breaks=c(0, 1)) + 
    scale_x_continuous(expand=c(0.005, 0.005), limits=c(1, nbins), breaks=c(1, 97, nbins), minor_breaks=c(1, seq(from=10, to=100, by=10)))
})

pdf.date(sprintf("%s/%s_lineplot_subtype_continuous_grouped_%sbins.pdf", sample, tag, nbins), height=5, width=12)
wrap_plots(grouped_plots, ncol=1) + plot_layout(guides = 'collect', axes="collect_x")
wrap_plots(grouped_plots_scaled, ncol=1) + plot_layout(axes="collect")
ggplot(even_bin_summary_melt, aes(x=even_bin, y=mean_score, color=program)) +
  geom_line() +
  facet_wrap(~type, scales="free_y", ncol=1) +
  theme_bw() 
dev.off()


# Other programs scaled line plot
line_style_groupings <- c(rep("EMT", 5), rep("pEMT", 2), rep("KRAS", 2), rep("Raghavan", 3), "Cell Cycle")
names(line_style_groupings) <- c("Gavish_MP12_EMTI", "Gavish_MP13_EMTII", "Gavish_MP14_EMTIII", "Gavish_MP15_EMTIV", "Hallmark_EMT", 
                                 "Puram_pEMT", "Dilly_pEMT", 
                                 "Singh_KRASDependency_addiction", "Klomp_PDAC_siKRAS_KRASi_iKras_UP", 
                                 "Raghavan_scClassical", "Raghavan_IC", "Raghavan_scBasal",
                                 "Connor_CC")
line_style_groupings <- factor(line_style_groupings, levels=unique(line_style_groupings))
custom_line_paras <- c("even_bin", names(line_style_groupings))
even_bin_summary_melt <- melt(even_bin_summary_counts[, ..custom_line_paras], id.vars="even_bin", value.name="mean_score", variable.name="program")
even_bin_summary_melt$type <- line_style_groupings[as.character(even_bin_summary_melt$program)]

seed(1)
grouped_plot_gene_colors <- sample(scales::hue_pal()(length(names(line_style_groupings))))
names(grouped_plot_gene_colors) <- names(line_style_groupings)

grouped_plots <- lapply(levels(even_bin_summary_melt$type), function(x) {
  even_bin_summary_melt_sub <- even_bin_summary_melt[type==x]
  ggplot(even_bin_summary_melt_sub, aes(x=even_bin, y=mean_score, color=program)) + 
    geom_line() + 
    theme_bw() + 
    scale_color_manual(values=grouped_plot_gene_colors) + 
    scale_x_continuous(expand=c(0, 0), limits=c(1, nbins), breaks=c(1, 97, nbins), minor_breaks=c(1, seq(from=10, to=100, by=10)))
})

bin_scaled <- even_bin_summary_counts[, ..custom_line_paras]
bin_scaled <- cbind(bin_scaled[, 1], apply(as.matrix(bin_scaled[, 2:ncol(bin_scaled)]), 2, function(x) {(x - min(x))/(max(x) - min(x))}))
bin_scaled_melt <- melt(bin_scaled, id.vars="even_bin", value.name="scaled_mean_score", variable.name="program")
bin_scaled_melt$type <- line_style_groupings[as.character(bin_scaled_melt$program)]

grouped_plots_scaled <- lapply(levels(even_bin_summary_melt$type), function(x) {
  bin_scaled_melt_sub <- bin_scaled_melt[type==x]
  ggplot(bin_scaled_melt_sub, aes(x=even_bin, y=scaled_mean_score, color=program)) + 
    geom_line() + 
    theme_bw() + 
    scale_color_manual(values=grouped_plot_gene_colors) + 
    scale_y_continuous(breaks=c(0, 1)) + 
    scale_x_continuous(expand=c(0.005, 0.005), limits=c(1, nbins), breaks=c(1, 97, nbins), minor_breaks=c(1, seq(from=10, to=100, by=10)))
})

pdf.date(sprintf("%s/%s_lineplot_subtype_continuous_grouped_%sbins.pdf", sample, tag, nbins), height=7, width=12)
wrap_plots(grouped_plots, ncol=1) + plot_layout(guides = 'collect', axes="collect_x")
wrap_plots(grouped_plots_scaled, ncol=1) + plot_layout(axes="collect")
ggplot(even_bin_summary_melt, aes(x=even_bin, y=mean_score, color=program)) +
  geom_line() +
  facet_wrap(~type, scales="free_y", ncol=1) +
  theme_bw() 
dev.off()



pdf.date(sprintf("%s/%s_heatmap_subtype_split_all.pdf", sample, sample), height=4, width=9)
for (plot_sample in unique(subtype_dataset$Sample)) {
  tryCatch({
    subtype_dataset_sub <- subtype_dataset[Sample == plot_sample & Celltype_final == "Malignant", ..subtype_names]
    subtype_order_score <- order(subtype_dataset[Sample == plot_sample & Celltype_final == "Malignant", pc1])
    anno_data_sub <- subtype_dataset[Sample == plot_sample & Celltype_final == "Malignant", ..top_anno_names][subtype_order_score, ]
    subtype_dataset_sub_mat <- t(as.matrix(subtype_dataset_sub))[, subtype_order_score]
    
    print(Heatmap(subtype_dataset_sub_mat, 
                  col = col_fun,
                  name = plot_sample,
                  column_split = anno_data_sub$subtype,
                  top_annotation = HeatmapAnnotation(df=anno_data_sub, col=col_list), 
                  cluster_columns = FALSE, 
                  cluster_rows = FALSE))
  }, error = function(e) {message("Problems with plotting ", plot_sample)})
}
dev.off()


#### Plotting: Subtypes across pancreas lineage ####

obj$pancreas_grouping <- ifelse(as.character(obj$Celltype_final) %in% c("Acinar", "Duct1", "Duct2", "Endocrine", "Malignant"), as.character(obj$Celltype_final), "other")
obj$pancreas_grouping <- ifelse(obj$pancreas_grouping == "Malignant",  obj$subtype, obj$pancreas_grouping)
obj$pancreas_grouping <- factor(obj$pancreas_grouping, levels= c("other", "Endocrine", "Acinar", "Duct1", "Duct2",
                                                                 "Classical2", "Classical1", "Mixed", 
                                                                 "Basal2", "Basal1"))
pancreas_grouping_colours <- c(singler_subtype_colors, other="gray90", cell_plot_colours)

vln_plot_data <- FetchData(obj, c("pancreas_grouping", "Classical1_score", "Classical2_score", "Basal1_score", "Basal2_score"))

vln_plot_data_melt <- melt(vln_plot_data, variable.name="subtype", value.name="score")
setDT(vln_plot_data_melt)
vln_plot_data_melt <- vln_plot_data_melt[, .(mean_score=mean(score), median_score=median(score)), by=.(subtype, pancreas_grouping)]

scaled_vln_plot_data <- as.matrix(vln_plot_data[, c("Classical1_score", "Classical2_score", "Basal1_score", "Basal2_score")])
rownames(scaled_vln_plot_data) <- vln_plot_data[, "pancreas_grouping"]
scaled_vln_plot_data <- scale(scaled_vln_plot_data)
scaled_vln_plot_data_melt <- melt(scaled_vln_plot_data, varnames=c("pancreas_grouping", "subtype"), value.name="score")
setDT(scaled_vln_plot_data_melt)
scaled_vln_plot_data_melt <- scaled_vln_plot_data_melt[, .(mean_scaled_score=mean(score), median_scaled_score=median(score)), by=.(subtype, pancreas_grouping)]
scaled_vln_plot_data_melt$pancreas_grouping <- factor(as.character(scaled_vln_plot_data_melt$pancreas_grouping), levels=levels(vln_plot_data_melt$pancreas_grouping))

vln_plot_data_melt[, scaled_mean_score := (mean_score - min(mean_score)) / (max(mean_score) - min(mean_score)), by=subtype]

pdf.date(sprintf("%s/%s_heatmap_subtype_scores_over_panc_lineage.pdf", sample, sample), height=4.5, width=10)
ggplot(vln_plot_data_melt, aes(x=pancreas_grouping, y=subtype, fill=scaled_mean_score)) + 
  geom_tile()  + 
  scale_fill_gradient(low="white", high="red")
ggplot(vln_plot_data_melt, aes(x=pancreas_grouping, y=subtype, fill=scaled_mean_score)) + 
  geom_tile()  + 
  scale_fill_gradient2(low="#f7fcfd", mid="#66c2a4", high="#00441b", midpoint=0.5)
ggplot(vln_plot_data_melt, aes(x=pancreas_grouping, y=subtype, fill=scaled_mean_score)) + 
  geom_tile()  + 
  scale_fill_gradient2(low="#f7fcfd", mid="#ccece6", high="#00441b", midpoint=0.5)
ggplot(vln_plot_data_melt, aes(x=pancreas_grouping, y=subtype, fill=scaled_mean_score)) + 
  geom_tile()  + 
  scale_fill_gradient2(low="#2166ac", mid="#f7f7f7", high="#b2182b", midpoint=0.5)
ggplot(vln_plot_data_melt, aes(x=pancreas_grouping, y=subtype, fill=scaled_mean_score)) + 
  geom_tile()  + 
  scale_fill_gradientn(colours=c("#2166ac", "#d1e5f0", "#f7f7f7", "#fddbc7", "#b2182b"), breaks=seq(0, 1, length.out=5))
ggplot(vln_plot_data_melt, aes(x=pancreas_grouping, y=subtype, fill=scaled_mean_score)) + 
  geom_tile()  + 
  scale_fill_gradient2(low="plum4", mid="black", high="yellow", midpoint=0.5)

vln_plot_data_melt_final <- vln_plot_data_melt
vln_plot_data_melt_final$subtype <- gsub("_score$", "", vln_plot_data_melt_final$subtype)
pdf.date(sprintf("%s/%s_heatmap_subtype_scores_over_panc_lineage_final.pdf", sample, sample), height=2, width=8)
ggplot(vln_plot_data_melt_final, aes(x=pancreas_grouping, y=subtype, fill=scaled_mean_score)) + 
  geom_tile() + 
  coord_fixed() +
  scale_fill_gradient2(low="plum4", mid="black", high="yellow", midpoint=0.5) + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,colour="black"), 
        axis.text.y = element_text(colour="black"), 
        axis.ticks = element_blank(),
        panel.background = element_blank())
dev.off()

pdf.date(sprintf("%s/%s_vlnplot_subtype_scores_over_panc_lineage_refmean.pdf", sample, sample), height=4.5, width=7.5)
for (y_var in c("Classical1_score", "Classical2_score", "Basal1_score", "Basal2_score")) {
  print(ggplot(vln_plot_data, aes(x=pancreas_grouping, y=.data[[y_var]])) + 
          geom_violin(aes(fill=pancreas_grouping)) +
          geom_boxplot(width=0.25, outlier.shape=NA) + 
          geom_hline(yintercept=0) +
          scale_fill_manual(values=pancreas_grouping_colours) + 
          stat_compare_means(ref.group = ".all.", label="p.signif", method="wilcox", paired=FALSE, method.args=list(alternative="greater"),
                             symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))) + 
          theme_bw() + labs(caption="alternative: greater; ns: p > 0.05; *: p <= 0.05; **: p <= 0.01; ***: p <= 0.001; ****: p <= 0.0001") + 
          theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
          NoLegend())
}
dev.off() 

pdf.date(sprintf("%s/%s_vlnplot_subtype_scores_over_panc_lineage.pdf", sample, sample), height=4.5, width=7.5)
for (y_var in c("Classical1_score", "Classical2_score", "Basal1_score", "Basal2_score")) {
  print(ggplot(vln_plot_data, aes(x=pancreas_grouping, y=.data[[y_var]])) + 
          geom_violin(aes(fill=pancreas_grouping)) +
          geom_boxplot(width=0.25, outlier.shape=NA) + 
          geom_hline(yintercept=0) +
          scale_fill_manual(values=pancreas_grouping_colours) + 
          stat_compare_means(ref.group = "other", label="p.signif", method="wilcox", paired=FALSE, method.args=list(alternative="greater"),
                             symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))) + 
          theme_bw() + labs(caption="alternative: greater; ns: p > 0.05; *: p <= 0.05; **: p <= 0.01; ***: p <= 0.001; ****: p <= 0.0001") + 
          theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
          NoLegend())
}
dev.off()

pdf.date(sprintf("%s/%s_boxplot_subtype_scores_over_panc_lineage_stats.pdf", sample, sample), height=6, width=9)
for (y_var in c("Classical1_score", "Classical2_score", "Basal1_score", "Basal2_score")) {
  print(ggsummarystats(
    vln_plot_data, x = "pancreas_grouping", y = y_var, fill="pancreas_grouping",
    ggfunc = ggboxplot, digits=6, palette=pancreas_grouping_colours))
}
dev.off()

vln_plot_data_melt <- melt(vln_plot_data, variable.name="subtype", value.name="score")
pdf.date(sprintf("%s/%s_vlnplot_subtype_scores_over_panc_lineage_wrapped.pdf", sample, sample), height=4.5, width=7.5)
  print(ggplot(vln_plot_data_melt, aes(x=pancreas_grouping, y=score)) +
          geom_hline(yintercept=0) +
          geom_violin(aes(fill=pancreas_grouping)) +
          scale_fill_manual(values=pancreas_grouping_colours) + 
          geom_boxplot(width=0.25, outlier.shape=NA) +
          theme_bw() + 
          theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
          NoLegend() +
          facet_wrap(vars(subtype), ncol=2, scales="free_y"))
  print(ggplot(vln_plot_data_melt, aes(x=pancreas_grouping, y=score)) +
          geom_hline(yintercept=0) +
          geom_violin(aes(fill=pancreas_grouping)) +
          scale_fill_manual(values=pancreas_grouping_colours) + 
          geom_boxplot(width=0.25, outlier.shape=NA) +
          theme_bw() + 
          theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
          NoLegend() +
          facet_wrap(vars(subtype), ncol=2))
dev.off() 

vln_plot_data <- FetchData(obj, c("SST", "GCG", "INS", "pancreas_grouping"))
pdf.date(sprintf("%s/%s_vlnplot_endocrine_genes_over_panc_lineage.pdf", sample, sample), height=7, width=7)
for (y_var in c("SST", "GCG", "INS")) {
  print(ggplot(vln_plot_data, aes(x=pancreas_grouping, y=.data[[y_var]])) + 
          geom_violin(aes(fill=pancreas_grouping)) +
          geom_boxplot(width=0.25, outlier.shape=NA) +  
          theme_bw())
}
dev.off() 

Idents(obj) <- "pancreas_grouping"
pdf.date(sprintf("%s/%s_dotplot_endocrine_genes_over_panc_lineage.pdf", sample, sample), height=3, width=7)
DotPlot(obj, group.by="pancreas_grouping", features=c("SST", "GCG", "INS", "CHGA", "MAFA", "ARX", "PDX1"))
dev.off() 


#### Plotting: KRAS status ####

# Data from Vartrix runs
root_path <- "source/KRAS_count_mats/"
cell_end_to_file <- c("85948" = "85948_KRAS_amp_mat.rds",
                      "91412" = "91412_KRAS_amp_mat.rds",
                      "94930_3p" = "94930_KRAS_amp_mat.rds",
                      "96460" = "96460_KRAS_amp_mat.rds",
                      "87235_CD45m" = "87235_KRAS_amp_mat.rds",
                      "91610" = "91610_KRAS_amp_mat.rds",
                      "95092" = "95092_KRAS_amp_mat.rds",
                      "G9903_CD45m" = "G9903_KRAS_amp_mat.rds",
                      "87784_CD45p" = "87784_KRAS_amp_mat.rds",
                      "91706" = "91706_KRAS_amp_mat.rds",
                      "95373" = "95373_KRAS_amp_mat.rds",
                      "PDA_115053_Pa_P_rna10x" = "115053_kras_amp_mat.rds",
                      "PDA_053403_Pa_P_rna10x" = "53403_kras_amp_mat.rds",
                      "97189" = "97189_kras_amp_mat.rds",
                      "PDA_139313_Pa_P_rna10x" = "139313_kras_amp_mat.rds",
                      "PDA_064630_Pa_P_rna10x" = "64630_kras_amp_mat.rds")

kras_list <- list()
for (cell_end in names(cell_end_to_file)) {
  file <- paste0(root_path, cell_end_to_file[[cell_end]])
  data <- readRDS(file)
  
  rownames(data) <- paste0(rownames(data), "_", cell_end)
  kras_list[[cell_end]] <- data
}

new_kras_data <- as.data.frame(do.call(rbind, kras_list))

obj <- AddMetaData(obj, new_kras_data)
obj$num_mut_cat <- ifelse(obj$num_mut >= 2, "2+", as.character(obj$num_mut))
obj$num_wt_cat <- ifelse(obj$num_wt >= 2, "2+", as.character(obj$num_wt))

obj$KRAS_mut_data <- !is.na(obj$num_mut)

obj$num_mut_bool <- obj$num_mut > 0
obj$num_wt_bool <- obj$num_wt > 0
obj$KRAS_status <- ifelse(obj$num_wt_bool & !obj$num_mut_bool, "wt only", ifelse(obj$num_mut_bool, "mut", NA))
obj$KRAS_status <- factor(as.character(obj$KRAS_status), levels=c("wt only", "mut"))

png.date(sprintf("%s/%s_UMAP_KRAS_mut_num.png", sample, sample), width=16, height=7)
FeaturePlot(obj, features=c("num_mut", "num_wt"), raster=FALSE, order=TRUE) * 
  theme(aspect.ratio = 1,
        axis.ticks=element_blank(),
        axis.text=element_blank())
dev.off()

png.date(sprintf("%s/%s_UMAP_KRAS_mut_cat.png", sample, sample), width=16, height=7)
DimPlot(obj, group.by=c("num_mut_cat", "num_wt_cat"), cols=c("0"="grey70", "1"="blue", "2+"="red"), 
        order=TRUE, raster=FALSE, na.value="grey90") * 
  theme(aspect.ratio = 1,
        axis.ticks=element_blank(),
        axis.text=element_blank())
dev.off()

png.date(sprintf("%s/%s_UMAP_KRAS_status.png", sample, sample), width=7, height=7)
DimPlot(obj, group.by="KRAS_status", cols=c("wt only"="blue", "mut"="red"), 
        order=TRUE, raster=FALSE, na.value="grey90") * 
  theme(aspect.ratio = 1,
        axis.ticks=element_blank(),
        axis.text=element_blank())
dev.off()


barplot_data <- reshape2::melt(table(obj[[c("num_mut_cat", "Celltype_final")]]))
barplot_data_sub <- barplot_data[barplot_data$Celltype_final %in% c("Acinar", "Duct1", "Duct2", "Malignant"),]

pdf.date(sprintf("%s/%s_barplot_KRAS_mut_cat.pdf", sample, sample), width=5, height=3)
ggplot(barplot_data, aes(x=Celltype_final, y=value, fill=num_mut_cat)) + 
  geom_col(position="fill") + 
  scale_y_continuous(expand=expansion(mult=0.01)) + 
  scale_fill_manual(values=c("0"="grey90", "1"="#FFD3C2", "2+"="red")) +
  labs(fill="# KRAS mut") +
  theme_bw() + 
  xlab("Celltype") + ylab("proportion") +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5), 
        panel.grid=element_blank())
ggplot(barplot_data_sub, aes(x=Celltype_final, y=value, fill=num_mut_cat)) + 
  geom_col(position="fill") + 
  labs(fill="# KRAS mut") +
  xlab("Celltype") + ylab("proportion") +
  scale_y_continuous(expand=expansion(mult=0.01)) + 
  scale_fill_manual(values=c("0"="grey90", "1"="#FFD3C2", "2+"="red")) + 
  theme_bw() + theme(panel.grid=element_blank())
dev.off()


mut_wt_data <- as.data.table(obj[[c("num_mut", "num_wt", "Celltype_final")]])
total_num_wt_data <- mut_wt_data[, .(total_mut=sum(num_mut, na.rm=TRUE), total_wt=sum(num_wt, na.rm=TRUE), 
                                     cells_with_data=sum(!is.na(num_mut)), cells_without_data=sum(is.na(num_mut)), total_cells=.N), 
                                 by=Celltype_final]
fwrite(total_num_wt_data, sprintf("%s/%s_KRAS_mut_counts.csv", sample, sample))

total_num_wt_data_melt <- melt(total_num_wt_data[, .(Celltype_final, total_mut, total_wt)], id.vars="Celltype_final", value.name="count")

pdf.date(sprintf("%s/%s_barplot_KRAS_mut_vs_wt.pdf", sample, sample), width=5, height=3)
ggplot(total_num_wt_data_melt, aes(x=Celltype_final, y=count, fill=variable)) + 
  geom_col() + 
  scale_fill_manual(values=c(total_mut="red", total_wt="blue")) + 
  theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5))
ggplot(total_num_wt_data_melt, aes(x=Celltype_final, y=count, fill=variable)) + 
  geom_col(position="fill") + 
  scale_fill_manual(values=c(total_mut="red", total_wt="blue")) + 
  ylab("proportion") +
  theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5)) 
ggplot(total_num_wt_data_melt, aes(x=Celltype_final, y=count, fill=variable)) + 
  geom_col(position="dodge") + 
  scale_fill_manual(values=c(total_mut="red", total_wt="blue")) + 
  ylab("proportion") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5),
                     panel.grid = element_blank()) 
dev.off()

celltypes_cutoff <- as.character(total_num_wt_data[total_mut + total_wt > 250, Celltype_final])

total_num_wt_data_melt_sub <- total_num_wt_data_melt[Celltype_final %in% celltypes_cutoff]
pdf.date(sprintf("%s/%s_barplot_KRAS_mut_vs_wt_cutoff.pdf", sample, sample), width=5, height=3)
ggplot(total_num_wt_data_melt_sub, aes(x=Celltype_final, y=count, fill=variable)) + 
  geom_col() + 
  scale_fill_manual(values=c(total_mut="red", total_wt="blue")) + 
  theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5))
ggplot(total_num_wt_data_melt_sub, aes(x=Celltype_final, y=count, fill=variable)) + 
  geom_col(position="fill") + ylab("proportion") +
  scale_fill_manual(values=c(total_mut="red", total_wt="blue")) + 
  theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5)) 
dev.off()

# For single sample KRAS plotting
sample_by_kras <- table(obj[[c("Sample", "KRAS_status")]])
sample_by_kras <- sample_by_kras[sample_by_kras[, "wt only"] > 0, ]
saveRDS(sample_by_kras, "sample_by_kras.rds")
saveRDS(obj[[c("num_mut", "num_wt", "KRAS_status", "Celltype_final")]], "kras_status_by_cell.rds")





