library(ComplexHeatmap)
library(circlize)
library(readxl)
library(data.table)
library(sagesalad)

stamp <- "carpenter_nanostring"

sample_metadata <- read_xlsx("external/carpenter_2023_nanostring/GSE226829_PostQC_segments_annotation.xlsx")
sample_metadata$sample_name <- gsub("-", ".", sample_metadata$`library name`)
sample_metadata$cell_type <- factor(sample_metadata$cell_type, levels = c("Acinar", "Duct", "ADM", "PanIN", "Glandular_Tumor", "PoorlyDiffTumor"))
sample_metadata <- sample_metadata[order(sample_metadata$DiseaseStatus, sample_metadata$patient_id), ]

norm_expression <- read_xlsx("external/carpenter_2023_nanostring/GSE226829_PostQC_logQ3_normalized_counts.xlsx")
norm_expr_mat <- as.matrix(norm_expression[, 2:ncol(norm_expression)])
rownames(norm_expr_mat) <- norm_expression[[1]]

all(sample_metadata$sample_name %in% colnames(norm_expr_mat))
all(colnames(norm_expr_mat) %in% sample_metadata$sample_name)

genesets_df <- fread("data/scRNAseq_genesets.csv")
geneset_list <- strsplit(genesets_df$genes, " ")
names(geneset_list) <- genesets_df$Geneset

disease_status_colours <- c(Healthy="grey90", Tumor="#5542FA")
patient_id_colours <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f")
names(patient_id_colours) <- unique(sample_metadata$patient_id)


#### Plot heatmap ####

sample_metadata_sub <- as.data.frame(sample_metadata[, c("cell_type", "patient_id", "DiseaseStatus")])
rownames(sample_metadata_sub) <- sample_metadata[["sample_name"]]

duct_de <- fread("data/merged_celltype_final_duct1_duct2_de.csv")

duct1_filtered <- duct_de[duct_de$p_val_adj < 0.05 & duct_de$avg_log2FC > 2,]
duct2_filtered <- duct_de[duct_de$p_val_adj < 0.05 & duct_de$avg_log2FC < -2,]

duct1_filtered <- duct1_filtered[order(duct1_filtered$p_val_adj), ]
duct2_filtered <- duct2_filtered[order(duct2_filtered$p_val_adj), ]

query_genes_raw <- list(acinar = c("PRSS1", "PRSS2", "PNLIP"), duct1 = duct1_filtered[[1]][1:25], duct2 = duct2_filtered[[1]][1:25])
query_genes_valid <- lapply(query_genes_raw, function(x) {x[x %in% rownames(norm_expr_mat)]})

unscaled_mat <- log2(norm_expr_mat[unlist(query_genes_valid), rownames(sample_metadata_sub)])
col_fun1 <- colorRamp2(quantile(unscaled_mat, c(0.01,  0.99)), c("cornsilk",  "darkred"))
col_fun2 <- colorRamp2(quantile(unscaled_mat, c(0.01, 0.5, 0.99)), c("#FF00FF", "black", "yellow"))

pdf.date(sprintf("%s_heatmap_de_acinarductgenes.pdf", stamp), width = 10, height = 10)
Heatmap(unscaled_mat,
        name = "log Q3 norm counts",
        col = col_fun1,
        top_annotation = HeatmapAnnotation(df = sample_metadata_sub[, c("patient_id", "DiseaseStatus")],
                                           col = list(DiseaseStatus=disease_status_colours, patient_id=patient_id_colours)),
        column_split = sample_metadata_sub$cell_type,
        cluster_columns = TRUE, 
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 9),
        row_title_gp = gpar(fontsize = 9),
        row_split = rep(names(query_genes_valid), sapply(query_genes_valid, length)),
        cluster_rows = FALSE)
Heatmap(unscaled_mat,
        name = "log Q3 norm counts",
        col = col_fun2,
        top_annotation = HeatmapAnnotation(df = sample_metadata_sub[, c("patient_id", "DiseaseStatus")],
                                           col = list(DiseaseStatus=disease_status_colours, patient_id=patient_id_colours)),
        column_split = sample_metadata_sub$cell_type,
        cluster_columns = TRUE, 
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 9),
        row_title_gp = gpar(fontsize = 9),
        row_split = rep(names(query_genes_valid), sapply(query_genes_valid, length)),
        cluster_rows = FALSE)
Heatmap(t(scale(t(unscaled_mat))),
        name = "z-score",
        top_annotation = HeatmapAnnotation(df = sample_metadata_sub[, c("patient_id", "DiseaseStatus")],
                                           col = list(DiseaseStatus=disease_status_colours, patient_id=patient_id_colours)),
        column_split = sample_metadata_sub$cell_type,
        cluster_columns = TRUE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 9),
        row_title_gp = gpar(fontsize = 9),
        row_split = rep(names(query_genes_valid), sapply(query_genes_valid, length)),
        cluster_rows = FALSE)
dev.off()

for (patient in unique(sample_metadata_sub$patient_id)) {
        sample_metadata_sub_patient <- sample_metadata_sub[sample_metadata_sub$patient_id == patient, ]
        unscaled_mat_patient <- unscaled_mat[, rownames(sample_metadata_sub_patient)]
        
        pdf.date(sprintf("%s_heatmap_de_acinarductgenes_%s.pdf", stamp, patient), width = 8, height = 8)
        print(Heatmap(unscaled_mat_patient,
                name = "log Q3 norm counts",
                col = col_fun1,
                top_annotation = HeatmapAnnotation(df = sample_metadata_sub_patient[, c("patient_id", "DiseaseStatus")],
                                                   col = list(DiseaseStatus=disease_status_colours, patient_id=patient_id_colours)),
                column_split = sample_metadata_sub_patient$cell_type,
                cluster_columns = TRUE, 
                cluster_column_slices = FALSE,
                show_column_names = FALSE,
                row_names_gp = gpar(fontsize = 9),
                row_title_gp = gpar(fontsize = 9),
                row_split = rep(names(query_genes_valid), sapply(query_genes_valid, length)),
                cluster_rows = FALSE))
        print(Heatmap(unscaled_mat_patient,
                name = "log Q3 norm counts",
                col = col_fun2,
                top_annotation = HeatmapAnnotation(df = sample_metadata_sub_patient[, c("patient_id", "DiseaseStatus")],
                                                   col = list(DiseaseStatus=disease_status_colours, patient_id=patient_id_colours)),
                column_split = sample_metadata_sub_patient$cell_type,
                cluster_columns = TRUE, 
                cluster_column_slices = FALSE,
                show_column_names = FALSE,
                row_names_gp = gpar(fontsize = 9),
                row_title_gp = gpar(fontsize = 9),
                row_split = rep(names(query_genes_valid), sapply(query_genes_valid, length)),
                cluster_rows = FALSE))
        print(Heatmap(t(scale(t(unscaled_mat_patient))),
                name = "z-score",
                top_annotation = HeatmapAnnotation(df = sample_metadata_sub_patient[, c("patient_id", "DiseaseStatus")],
                                                   col = list(DiseaseStatus=disease_status_colours, patient_id=patient_id_colours)),
                column_split = sample_metadata_sub_patient$cell_type,
                cluster_columns = TRUE,
                cluster_column_slices = FALSE,
                show_column_names = FALSE,
                row_names_gp = gpar(fontsize = 9),
                row_title_gp = gpar(fontsize = 9),
                row_split = rep(names(query_genes_valid), sapply(query_genes_valid, length)),
                cluster_rows = FALSE))
        dev.off()
}
