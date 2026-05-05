library(Seurat)
library(Matrix)
library(hdf5r)
library(ggplot2)
library(sagesalad)
library(data.table)
library(Hmisc)
library(corrplot)
library(RColorBrewer)
library(patchwork)
library(numbat)
library(princomp)
library(igraph)
library(ggraph)
library(tidygraph)
library(scatterpie)z

Moffitt_colours_malignant <- c("Basal-like" = "#BEC100",
                               "Classical" = "#B0ABFF",
                               "non-malignant" = "lightgrey")

singler_subtype_colors <- c("#1C6CAB", "#814C42","#A4C0E5", "#FF7311", "grey40")
names(singler_subtype_colors) <- c("Classical1", "Basal1", "Classical2", "Basal2", "Mixed")
singler_subtype_colors <- singler_subtype_colors[c("Basal1", "Basal2", "Mixed", "Classical1", "Classical2")]

singler_subtype_colors_malignant <- c("#1C6CAB", "#814C42","#A4C0E5", "#FF7311", "grey40", "grey90")
names(singler_subtype_colors_malignant) <- c("Classical1", "Basal1", "Classical2", "Basal2", "Mixed", "non-malignant")
singler_subtype_colors_malignant <- singler_subtype_colors_malignant[c("Basal1", "Basal2", "Mixed", "Classical1", "Classical2", "non-malignant")]

genesets <- fread("data/scRNAseq_genesets.csv")
genesets_list <- strsplit(genesets$genes, " ")
names(genesets_list) <- genesets$Geneset

marker_genes_to_plot <- unlist(genesets_list[grep("^Markers_", names(genesets_list))], use.names = FALSE)
names(marker_genes_to_plot) <- rep(names(genesets_list[grep("^Markers_", names(genesets_list))]), lengths(genesets_list[grep("^Markers_", names(genesets_list))]))
marker_genes_to_plot <- marker_genes_to_plot[!duplicated(marker_genes_to_plot)]

celltype_order <- gsub("Markers_", "", names(genesets_list[grep("^Markers_", names(genesets_list))]))

celltype_colour_file <- "data/celltype_colours.csv"
celltype_colours <- fread(celltype_colour_file)
celltype_order <- celltype_colours$celltype
cell_plot_colours <- celltype_colours$hexcode
names(cell_plot_colours) <- celltype_order

cell_plot_colours_updated <- cell_plot_colours
cell_plot_colours_updated["Basal_97727"] <- "brown" 

dir.create("numbat", showWarnings = FALSE)

all_files <- list.files(pattern="^sample_.*\\.rds$")

for (file in all_files[1:length(all_files)]) {
  tryCatch({
    sample <- gsub(".rds", "", file)
    obj <- readRDS(file)
    
    dir.create(sample, showWarnings = FALSE)
    
    obj$Celltype_by_cluster <- obj$Celltype_by_cluster_0.8
    
    celltype_final <- read.csv("merged_celltype_final.csv", row.names = 1)
    obj <- AddMetaData(obj, celltype_final, "Celltype_final")
    obj$Celltype_final <- factor(as.character(obj$Celltype_final), levels=celltype_order)
    
    obj$subtype_mal <- ifelse(obj$Celltype_final  == "Malignant", obj$subtype, "non-malignant")
    obj$subtype_mal <- factor(as.character(obj$subtype_mal), levels=names(singler_subtype_colors_malignant))     
    
    # Manual update based on Numbat review
    if (sample == "sample_97727") { 
      obj$Celltype_final_updated <- ifelse(obj$subtype_mal == "Basal1", "Basal_97727", as.character(obj$Celltype_final))
      obj$Celltype_final_updated <- factor(as.character(obj$Celltype_final_updated), levels=names(cell_plot_colours_updated))
      
      basal_97727_cells <- colnames(obj)[obj$subtype_mal == "Basal1"]
      celltype_final[basal_97727_cells, ] <- "Basal_97727"
      write.csv(celltype_final,  "merged_celltype_final_updated_97727.csv")
    }
    
    #### Plotting: UMAP ####

    for (var in colnames(obj@meta.data)) {
      png.date(sprintf("%s/%s_UMAP_raw_%s.png", sample, sample, var))
      if (is.numeric(obj@meta.data[[var]])) {
        print(FeaturePlot(obj, features=var, raster=FALSE) +
                theme(aspect.ratio = 1,
                      axis.ticks=element_blank(),
                      axis.text=element_blank()))
      } else {
        print(DimPlot(obj, group.by=var, raster=FALSE) + ggtitle(var) +
                theme(aspect.ratio = 1,
                      axis.ticks=element_blank(),
                      axis.text=element_blank()))
      }
      dev.off()
    }

    
    #### Plotting: Select genes UMAP ####
    
    select_genes <- c("EPCAM", "CFTR", "KRT19", "KRT17", "IL2RG",
                      "MUC5AC", "S100P", "S100A10", "FOXQ1", "ONECUT2",
                      "IL18", "IL1RN", "IL1A", "TIMP1", "OIT1")
    select_features <- select_genes[select_genes %in% rownames(obj)]
    for (gene in select_features) {
      png.date(sprintf("%s/%s_UMAP_raw_%s.png", sample, sample, gene))
      print(FeaturePlot(obj, gene, raster=FALSE, col=c("yellow", "blue")) +
              theme(aspect.ratio = 1,
                    axis.ticks=element_blank(),
                    axis.text=element_blank()))
      dev.off()
    }
    for (gene in select_features) {
      png.date(sprintf("%s/%s_UMAP_raw_ordered_%s.png", sample, sample, gene), width=5, height=5)
      print(FeaturePlot(obj, gene, raster=FALSE, order=TRUE, col=c("yellow", "blue")) +
              theme(aspect.ratio = 1,
                    axis.ticks=element_blank(),
                    axis.text=element_blank()))
      dev.off()
    }
  
    
    #### Plotting: Dotplot ####
    
    select_genes <- c("EPCAM", "CFTR", "KRT19", "KRT17", "IL2RG",
                      "MUC5AC", "S100P", "S100A10", "FOXQ1", "ONECUT2",
                      "IL18", "IL1RN", "IL1A", "TIMP1", "OIT1")
    select_features <- select_genes[select_genes %in% rownames(obj)]
    pdf.date(sprintf("%s/%s_dotplot_selected_genes.pdf", sample, sample), height=7, width=20)
    print(DotPlot(obj, features=select_features, group.by="Celltype_final") + RotatedAxis() +
            theme(strip.text.x=element_text(angle = 90, hjust=0, vjust = 0.5)))
    dev.off()
    
    if (sample == "sample_97727") {
      obj$Celltype_final_updated_simplified <- ifelse(obj$Celltype_final_updated %in% c("Malignant", "Basal_97727"), as.character(obj$Celltype_final_updated), "other")
      
      pdf.date(sprintf("%s/%s_dotplot_selected_genes_updated_celltype.pdf", sample, sample), height=7, width=7)
      print(DotPlot(obj, features=select_features, group.by="Celltype_final_updated") + RotatedAxis() +
              theme(strip.text.x=element_text(angle = 90, hjust=0, vjust = 0.5)))
      print(DotPlot(obj, features=select_features, group.by="Celltype_final_updated_simplified") + RotatedAxis() +
              theme(strip.text.x=element_text(angle = 90, hjust=0, vjust = 0.5)))
      dev.off()
    }

    
    #### Plotting: specific UMAPs ####

    png.date(sprintf("%s/%s_UMAP_subtype_mal.png", sample, sample), width=10, height=7)
    print(DimPlot(obj, group.by="subtype_mal", raster=FALSE, cols=singler_subtype_colors_malignant) +
            theme(aspect.ratio = 1,
                  axis.ticks=element_blank(),
                  axis.text=element_blank()))
    dev.off()
    
    scores <- obj[[c("Basal1_score", "Basal2_score",  "Classical1_score", "Classical2_score")]]
    max_lim <- max(scores)
    min_lim <- min(scores)
    for (score in c("Basal1_score", "Basal2_score",  "Classical1_score", "Classical2_score")) {
      png.date(sprintf("%s/%s_UMAP_subtype_cont_%s.png", sample, sample, score), width=6, height=4)
      print(FeaturePlot(obj, features = score, raster=FALSE) +
              scale_colour_gradient2(low="yellow", mid="grey", high="purple", limits=c(min_lim, max_lim), midpoint=0) + 
              theme(aspect.ratio = 1,
                    axis.ticks=element_blank(),
                    axis.text=element_blank()))
      dev.off()
    }
    
    png.date(sprintf("%s/%s_UMAP_subtype_mal.png", sample, sample), width=10, height=7)
    print(DimPlot(obj, group.by="subtype_mal", raster=FALSE, cols=singler_subtype_colors_malignant) +
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

    if (sample == "sample_97727") {
      png.date(sprintf("%s/%s_UMAP_Celltype_final_updated.png", sample, sample), width=10, height=7)
      print(DimPlot(obj, group.by="Celltype_final_updated", raster=FALSE, cols=cell_plot_colours_updated) +
              theme(aspect.ratio = 1,
                    axis.ticks=element_blank(),
                    axis.text=element_blank()))
      dev.off()
    }

    png.date(sprintf("%s/%s_UMAP_leiden_clusters_wider.png", sample, sample), width=11, height=7)
    print(DimPlot(obj, group.by="leiden_clusters", raster=FALSE, label=TRUE) +
            theme(aspect.ratio = 1,
                  axis.ticks=element_blank(),
                  axis.text=element_blank()))
    dev.off()
    
    obj$Moffitt_subtype_malignant <- ifelse(obj$Celltype_by_cluster == "Malignant", gsub("^Moffitt_", "", obj$Moffitt_subtype), "non-malignant")
    obj$Moffitt_Classical_malignant <- ifelse(obj$Celltype_by_cluster == "Malignant", obj$Moffitt_Classical, NA)
    obj$`Moffitt_Basal.like_malignant` <- ifelse(obj$Celltype_by_cluster == "Malignant", obj$`Moffitt_Basal-like`, NA)

    pdf.date(sprintf("%s/%s_UMAP_Moffitt_subtype_malignant.pdf", sample, sample))
    print(DimPlot(obj, group.by="Moffitt_subtype_malignant", raster=FALSE, cols=Moffitt_colours_malignant))
    print(FeaturePlot(obj, features="Moffitt_Classical_malignant") + scale_colour_gradient(low="yellow", high="blue", na.value="lightgrey"))
    print(FeaturePlot(obj, features="Moffitt_Basal.like_malignant") + scale_colour_gradient(low="yellow", high="blue", na.value="lightgrey"))
    dev.off()
    
    if (!endsWith(sample, "_mal")) {

      numbat_base_dir <- "numbat/max_iter_5_no_converge"

      i <- 5
      while (!file.exists(sprintf("%s/%s/joint_post_%s.tsv", numbat_base_dir, sample, i)) & i > 0) {
        i <- i - 1
      }

      cat("Using iteration", i, "\n")

      nb <- Numbat$new(out_dir = sprintf("%s/%s", numbat_base_dir, sample), i=i)

      annot_metadata <- obj[[c("Celltype_final", "subtype_mal", "leiden_clusters")]]
      annot_metadata$cell <- gsub("_.*", "", rownames(annot_metadata))

      fwrite(nb$joint_post %>% filter(p_cnv > 0.99) %>% distinct(seg, cnv_state, CHROM) %>% arrange(CHROM, seg),
             file = sprintf("numbat/%s_NB_distinct_cnv.csv", sample))

      joint_plot_data <- nb$joint_post[, c("cell", "cnv_state", "CHROM", "seg_start", "seg_end", "p_cnv")]
      joint_plot_data <- joint_plot_data[joint_plot_data$p_cnv > 0.99]
      anno_data <- obj[[c("subtype_mal", "Celltype_final", "total_counts", "Basal1_score", "Basal2_score", "Classical1_score", "Classical2_score")]]
      if (sample == "sample_97727") {
        anno_data <- obj[[c("subtype_mal", "Celltype_final", "Celltype_final_updated", "total_counts", "Basal1_score", "Basal2_score", "Classical1_score", "Classical2_score")]]
      }
      anno_data$subtype_mal <- factor(anno_data$subtype_mal, levels=names(singler_subtype_colors_malignant))

      anno_data$cell <- gsub("_.*$", "", rownames(anno_data))

      pca_out <- princomp(anno_data[, c("Basal1_score", "Basal2_score", "Classical1_score", "Classical2_score")])
      anno_data$pc1 <- pca_out$scores[, 1]
      anno_data <- anno_data[order(anno_data$Celltype_final, anno_data$pc1), ]
      anno_data$cell_number <- 1:nrow(anno_data)
      anno_data$cell_factor <- factor(anno_data$cell, levels=anno_data$cell)

      heatmap_plot_data <- merge(joint_plot_data, anno_data, by="cell", all=TRUE)
      heatmap_plot_data <- merge(heatmap_plot_data, nb$clone_post[, c("cell", "clone_opt")], by="cell", all.x=TRUE)

      # populate cells with no cnv
      ref_neutral <- nb$segs_consensus[cnv_state == "neu"][1]
      heatmap_plot_data[is.na(cnv_state),  CHROM := ref_neutral$CHROM ]
      heatmap_plot_data[is.na(cnv_state),  seg_start := ref_neutral$seg_start ]
      heatmap_plot_data[is.na(cnv_state),  seg_end := ref_neutral$seg_end ]

      annot_data_melt <- melt(heatmap_plot_data[, c("cell_factor", "clone_opt", "Basal1_score", "Basal2_score", "Classical1_score", "Classical2_score")],
                              id.vars=c("cell_factor", "clone_opt"))
      clone_heatmap <- ggplot(heatmap_plot_data, aes(x=seg_start, xend=seg_end, y=cell_factor, yend=cell_factor, color=cnv_state)) +
        geom_segment(linewidth=0.1) +
        geom_segment(inherit.aes = FALSE,
                     aes(x = seg_start, xend = seg_end, y = 1, yend = 1),
                     data = nb$segs_consensus, linewidth = 0, color = 'white', alpha = 0) +
        facet_grid(cols=vars(CHROM), rows=vars(clone_opt), scales="free", space="free", switch="y") +
        theme(panel.spacing = unit(0, 'mm'),
              panel.border = element_rect(linewidth = 0.5, color = 'gray', fill = NA),
              strip.background = element_blank(),
              axis.text = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks = element_blank(),
              plot.margin = margin(0,0,0,0, unit = 'mm'),
              axis.line = element_blank(),
              legend.box.background = element_blank(),
              legend.background = element_blank(),
              legend.margin = margin(0,0,0,0),
              strip.text.y.left=element_blank()) +
        scale_x_continuous(expand = expansion(0)) +
        scale_color_manual(values = c('amp' = 'darkred', 'del' = 'darkblue', 'bamp' = 'pink', 'loh' = 'darkgreen', 'bdel' = 'blue'),
                           labels = c('amp' = 'AMP', 'del' = 'DEL', 'bamp' = 'BAMP', 'loh' = 'CNLoH', 'bdel' = 'BDEL'),
                           na.translate=FALSE)
      score_bar <- ggplot(annot_data_melt, aes(x=variable, y=cell_factor, fill=value)) +
        geom_tile() +
        facet_grid(rows=vars(clone_opt), scales="free", space="free", switch="y") +
        scale_fill_gradient2(high="#b2182b", low="#d1e5f0", mid="#f7f7f7") +
        theme(panel.spacing = unit(0, 'mm'),
              panel.border = element_rect(linewidth = 0.5, color = 'gray', fill = NA),
              strip.background = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=rel(0.3)),
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks = element_blank(),
              plot.margin = margin(0,0,0,0, unit = 'mm'),
              axis.line = element_blank(),
              legend.box.background = element_blank(),
              legend.background = element_blank(),
              legend.margin = margin(0,0,0,0),
              strip.text.y.left=element_blank())
      annotation_bar <- ggplot(heatmap_plot_data, aes(x=1, y=cell_factor, fill=subtype_mal)) +
        geom_tile() +
        scale_fill_manual(values=singler_subtype_colors_malignant) +
        facet_grid(rows=vars(clone_opt), scales="free", space="free", switch="y") +
        scale_x_continuous(expand = expansion(0)) +
        theme(panel.spacing = unit(0, 'mm'),
              panel.border = element_rect(linewidth = 0.5, color = 'gray', fill = NA),
              strip.background = element_blank(),
              axis.text = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks = element_blank(),
              plot.margin = margin(0,0,0,0, unit = 'mm'),
              axis.line = element_blank(),
              legend.box.background = element_blank(),
              legend.background = element_blank(),
              legend.margin = margin(0,0,0,0))
      celltype_bar <- ggplot(heatmap_plot_data, aes(x=1, y=cell_factor, fill=Celltype_final)) +
        geom_tile() +
        scale_fill_manual(values=cell_plot_colours) +
        facet_grid(rows=vars(clone_opt), scales="free", space="free", switch="y") +
        scale_x_continuous(expand = expansion(0)) +
        theme(panel.spacing = unit(0, 'mm'),
              panel.border = element_rect(linewidth = 0.5, color = 'gray', fill = NA),
              strip.background = element_blank(),
              axis.text = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks = element_blank(),
              plot.margin = margin(0,0,0,0, unit = 'mm'),
              axis.line = element_blank(),
              legend.box.background = element_blank(),
              legend.background = element_blank(),
              legend.margin = margin(0,0,0,0))
      if (sample == "sample_97727") {
        celltype_bar <- ggplot(heatmap_plot_data, aes(x=1, y=cell_factor, fill=Celltype_final_updated)) +
          geom_tile() +
          scale_fill_manual(values=cell_plot_colours_updated) +
          facet_grid(rows=vars(clone_opt), scales="free", space="free", switch="y") +
          scale_x_continuous(expand = expansion(0)) +
          theme(panel.spacing = unit(0, 'mm'),
                panel.border = element_rect(linewidth = 0.5, color = 'gray', fill = NA),
                strip.background = element_blank(),
                axis.text = element_blank(),
                axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                axis.ticks = element_blank(),
                plot.margin = margin(0,0,0,0, unit = 'mm'),
                axis.line = element_blank(),
                legend.box.background = element_blank(),
                legend.background = element_blank(),
                legend.margin = margin(0,0,0,0))
      }

      pdf.date(sprintf("numbat/%s_NB_clone_heatmap.pdf", sample), width=8.5, height=4)
      print(celltype_bar + annotation_bar + score_bar + clone_heatmap + plot_layout(widths=c(0.25, 0.25, 1, 20), guides="collect"))
      dev.off()
      
      pdf.date(sprintf("numbat/%s_NB_clone_heatmap_taller.pdf", sample), width=8.5, height=10)
      print(celltype_bar + annotation_bar + score_bar + clone_heatmap + plot_layout(widths=c(0.25, 0.25, 1, 20), guides="collect"))
      dev.off()

      heatmap_plot_data_mal <- heatmap_plot_data[heatmap_plot_data$Celltype_final == "Malignant" & heatmap_plot_data$clone_opt != 1, ]
      
      annot_data_melt_mal <- melt(heatmap_plot_data_mal[, c("cell_factor", "clone_opt", "Basal1_score", "Basal2_score", "Classical1_score", "Classical2_score")],
                              id.vars=c("cell_factor", "clone_opt"))
      clone_heatmap <- ggplot(heatmap_plot_data_mal, aes(x=seg_start, xend=seg_end, y=cell_factor, yend=cell_factor, color=cnv_state)) +
        geom_segment(linewidth=0.1) +
        geom_segment(inherit.aes = FALSE,
                     aes(x = seg_start, xend = seg_end, y = 1, yend = 1),
                     data = nb$segs_consensus, linewidth = 0, color = 'white', alpha = 0) +
        facet_grid(cols=vars(CHROM), rows=vars(clone_opt), scales="free", space="free", switch="y") +
        theme(panel.spacing = unit(0, 'mm'),
              panel.border = element_rect(linewidth = 0.5, color = 'gray', fill = NA),
              strip.background = element_blank(),
              axis.text = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks = element_blank(),
              plot.margin = margin(0,0,0,0, unit = 'mm'),
              axis.line = element_blank(),
              legend.box.background = element_blank(),
              legend.background = element_blank(),
              legend.margin = margin(0,0,0,0),
              strip.text.y.left=element_blank()) +
        scale_x_continuous(expand = expansion(0)) +
        scale_color_manual(values = c('amp' = 'darkred', 'del' = 'darkblue', 'bamp' = 'pink', 'loh' = 'darkgreen', 'bdel' = 'blue'),
                           labels = c('amp' = 'AMP', 'del' = 'DEL', 'bamp' = 'BAMP', 'loh' = 'CNLoH', 'bdel' = 'BDEL'),
                           na.translate=FALSE)
      score_bar <- ggplot(annot_data_melt_mal, aes(x=variable, y=cell_factor, fill=value)) +
        geom_tile() +
        facet_grid(rows=vars(clone_opt), scales="free", space="free", switch="y") +
        scale_fill_gradient2(high="#b2182b", low="#d1e5f0", mid="#f7f7f7") +
        theme(panel.spacing = unit(0, 'mm'),
              panel.border = element_rect(linewidth = 0.5, color = 'gray', fill = NA),
              strip.background = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=rel(0.3)),
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks = element_blank(),
              plot.margin = margin(0,0,0,0, unit = 'mm'),
              axis.line = element_blank(),
              legend.box.background = element_blank(),
              legend.background = element_blank(),
              legend.margin = margin(0,0,0,0),
              strip.text.y.left=element_blank())
      annotation_bar <- ggplot(heatmap_plot_data_mal, aes(x=1, y=cell_factor, fill=subtype_mal)) +
        geom_tile() +
        scale_fill_manual(values=singler_subtype_colors_malignant) +
        facet_grid(rows=vars(clone_opt), scales="free", space="free", switch="y") +
        scale_x_continuous(expand = expansion(0)) +
        theme(panel.spacing = unit(0, 'mm'),
              panel.border = element_rect(linewidth = 0.5, color = 'gray', fill = NA),
              strip.background = element_blank(),
              axis.text = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks = element_blank(),
              plot.margin = margin(0,0,0,0, unit = 'mm'),
              axis.line = element_blank(),
              legend.box.background = element_blank(),
              legend.background = element_blank(),
              legend.margin = margin(0,0,0,0))

      pdf.date(sprintf("numbat/%s_NB_clone_heatmap_mal_only_strict.pdf", sample), width=8, height=4)
      print(annotation_bar + score_bar + clone_heatmap + plot_layout(widths=c(0.25, 1, 20), guides="collect"))
      dev.off()

      # Based on numbat code for plot_mut_history()
      mut_g <- nb$mut_graph
      mut_df <- as_tbl_graph(mut_g)
      mut_df <- mut_df %>% activate(edges) %>%
        mutate(n_mut = unlist(purrr::map(stringr::str_split(to_label, ','), length))) %>%
        mutate(length = n_mut)

      pdf.date(sprintf("numbat/%s_NB_blank_mut_dendrogram.pdf", sample), height=4, width=8)
      print(ggraph(mut_df,
                   layout = 'dendrogram',
                   length = length) +
              geom_edge_elbow(aes(label = stringr::str_trunc(to_label, 20, side = 'center')),
                              vjust = -1,
                              hjust = -0.2,
                              arrow = arrow(length = unit(2, "mm")),
                              end_cap = circle(6, 'mm'),
                              start_cap = circle(6, 'mm'),
                              label_size = 4) +
              theme_void() +
              scale_x_continuous(expand = expansion(0.2)) +
              guides(color = 'none') +
              # geom_node_text(aes(label = clone), size = 6) +
              coord_flip() + scale_y_reverse())
      dev.off()

      pdf.date(sprintf("numbat/%s_NB_blank_mut_dendrogram_full_label.pdf", sample), height=3, width=7)
      print(ggraph(mut_df,
                   layout = 'dendrogram',
                   length = length) +
              geom_edge_elbow(aes(label = to_label),
                              vjust = -1,
                              hjust = -0.2,
                              arrow = arrow(length = unit(2, "mm")),
                              end_cap = circle(6, 'mm'),
                              start_cap = circle(6, 'mm'),
                              label_size = 4) +
              theme_void() +
              scale_x_continuous(expand = expansion(0.2)) +
              guides(color = 'none') +
              coord_flip() + scale_y_reverse())
      dev.off()

      pie_data_discrete <- as.data.frame(table(heatmap_plot_data[, c("subtype_mal", "clone_opt")]))
      pie_data_discrete_no_mal_strict <- pie_data_discrete[pie_data_discrete$clone_opt!= 1 & pie_data_discrete$subtype_mal != "non-malignant", ]
      pie_data_cont_raw <- heatmap_plot_data[, c("Basal1_score", "Basal2_score", "Classical1_score", "Classical2_score", "clone_opt")]
      pie_data_cont_raw <- melt(pie_data_cont_raw, id.vars="clone_opt", variable.name="subtype_score", value="score")
      pie_data_cont <- pie_data_cont_raw[, .(mean_score = mean(score)), by=.(clone_opt, subtype_score)][order(clone_opt, subtype_score)]
      pie_data_cont$subtype_score <- gsub("_score$", "", pie_data_cont$subtype_score)
      pie_data_cont$clone_opt <- as.factor(pie_data_cont$clone_opt)

      pdf.date(sprintf("numbat/%s_NB_clone_charts_for_dendrogram.pdf", sample), height=8, width=8)
      print(ggplot(pie_data_discrete, aes(x="", y=Freq, fill=subtype_mal)) +
              geom_col(position="fill") +
              scale_fill_manual(values=singler_subtype_colors_malignant) +
              coord_polar("y", start=0) +
              facet_wrap(~clone_opt) +
              theme_void() +
              guides(color = 'none'))
      print(ggplot(pie_data_discrete_no_mal_strict, aes(x="", y=Freq, fill=subtype_mal)) +
              geom_col(position="fill") +
              scale_fill_manual(values=singler_subtype_colors_malignant) +
              coord_polar("y", start=0) +
              facet_wrap(~clone_opt) +
              theme_void() +
              guides(color = 'none'))
      print(ggplot(pie_data_discrete, aes(x=clone_opt, y=Freq, fill=subtype_mal)) +
              geom_col(position="fill") +
              scale_fill_manual(values=singler_subtype_colors_malignant) +
              theme_void() +
              guides(color = 'none'))
      print(ggplot(pie_data_cont, aes(x=clone_opt, y=mean_score, fill=subtype_score)) +
              geom_col(position="dodge") +
              scale_fill_manual(values=singler_subtype_colors_malignant) +
              theme_void() +
              guides(color = 'none'))
      print(ggplot(pie_data_cont, aes(x=clone_opt, y=subtype_score, fill=subtype_score, alpha=mean_score)) +
              geom_tile(color="white") +
              scale_fill_manual(values=singler_subtype_colors_malignant) +
              theme_void() +
              guides(color = 'none') +
              coord_fixed())
      dev.off()

      pdf.date(sprintf("numbat/%s_NB_clone_chart.pdf", sample), height=7, width=10)
      print(ggplot(pie_data_cont, aes(x=clone_opt, y=mean_score, fill=subtype_score)) +
              geom_col(position="dodge") +
              geom_hline(yintercept=0, linetype='dotted') +
              scale_fill_manual(values=singler_subtype_colors_malignant) +
              theme_classic())
      print(ggplot(pie_data_cont, aes(x=clone_opt, y=mean_score, fill=subtype_score)) +
              geom_col(position="dodge") +
              geom_hline(yintercept=0, linetype='dotted') +
              facet_wrap(~subtype_score, ncol=1, scales="free") +
              scale_fill_manual(values=singler_subtype_colors_malignant) +
              theme_bw())
      print(ggplot(pie_data_cont, aes(x=clone_opt, y=mean_score, fill=subtype_score)) +
              geom_col(position="dodge") +
              geom_hline(yintercept=0, linetype='dotted') +
              facet_wrap(~subtype_score, ncol=2, scales="free") +
              scale_fill_manual(values=singler_subtype_colors_malignant) +
              theme_bw())
      dev.off()

      ##
      n_clones <- length(unique(nb$clone_post$clone_opt))

      pdf.date(sprintf("numbat/%s_NB_bulk_clones.pdf", sample), height=1.5 * n_clones, width=10)
      print(plot_bulks(nb$bulk_clones,
                       min_LLR = 10,
                       legend = TRUE))
      dev.off()

      pdf.date(sprintf("numbat/%s_NB_bulk_clones_wide.pdf", sample), height=1.5 * n_clones, width=30)
      print(plot_bulks(nb$bulk_clones,
                       min_LLR = 10,
                       legend = TRUE) &
              theme(legend.direction="horizontal",
                    legend.position="bottom"))
      dev.off()



      pdf.date(sprintf("numbat/%s_NB_builtin_plots.pdf", sample), height=6, width=8)
      print(nb$plot_phylo_heatmap(annot=annot_metadata, annot_bar_width=1,
                                  annot_scale=scale_fill_manual(values=c(cell_plot_colours, singler_subtype_colors_malignant))))
      print(plot_bulks(nb$bulk_clones,
                       min_LLR = 10, # filtering CNVs by evidence
                       legend = TRUE))
      print(nb$plot_clone_profile())
      print(nb$plot_sc_tree(
        label_size = 3,
        branch_width = 0.5,
        tip_length = 0.5,
        tip = TRUE))
      print(nb$plot_mut_history())
      dev.off()

      numbat_metadata <- nb$clone_post[, c("cell", "clone_opt",  "p_opt", "p_cnv", "p_cnv_x", "p_cnv_y", "compartment_opt")]

      cell_to_rowname_cell <- rownames(annot_metadata)
      names(cell_to_rowname_cell) <- annot_metadata$cell
      rownames(numbat_metadata) <- cell_to_rowname_cell[numbat_metadata$cell]
      numbat_metadata$clone_opt <- factor(numbat_metadata$clone_opt)

      obj <- AddMetaData(obj, numbat_metadata)

      write.csv(obj[[c("clone_opt", "cell")]], sprintf("numbat/%s_clone_data.csv", sample))

      pdf.date(sprintf("numbat/%s_NB_UMAP_posterior_prob.pdf", sample), width=28, height=7)
      print(FeaturePlot(obj, features=c("p_cnv", "p_cnv_x", "p_cnv_y", "p_opt"), ncol=4) * theme(aspect.ratio = 1,
                                                                                                 axis.ticks=element_blank(),
                                                                                                 axis.text=element_blank()))
      dev.off()

      clone_colours <- c("grey90", sort(brewer.pal(n_clones - 1, "Set2")))
      names(clone_colours) <- 1:n_clones

      pdf.date(sprintf("numbat/%s_NB_UMAP_clones.pdf", sample), width=7, height=7)
      print(DimPlot(obj, group.by="clone_opt") +
              scale_colour_manual(values=clone_colours) +
              theme(aspect.ratio = 1,
                    axis.ticks=element_blank(),
                    axis.text=element_blank()))
      print(DimPlot(obj, group.by="Celltype_final") +
              scale_color_manual(values=cell_plot_colours) +
              theme(aspect.ratio = 1,
                    axis.ticks=element_blank(),
                    axis.text=element_blank()))
      dev.off()

    } else {
    numbat_clone_file <- sprintf("numbat/%s_clone_data.csv", gsub("_mal$", "", sample))
    if (file.exists(numbat_clone_file)) {
      clone_data <- read.csv(numbat_clone_file, row.names = 1)
      clone_data$clone_opt <- factor(clone_data$clone_opt)

      obj <- AddMetaData(obj, clone_data)


      n_clones <- length(unique(clone_data$clone_opt))
      clone_colours1 <- c("grey90", "#6A0FFF", "#39FF1C", "#F5BC16", "#FF0DF6", "#FF1500")[1:n_clones]
      names(clone_colours1) <- 1:n_clones

      clone_colours2 <- c("grey90", "#FF1500", "#F5BC16", "#39FF1C", "#6A0FFF",  "#FF0DF6")[1:n_clones]
      names(clone_colours2) <- 1:n_clones

      clone_colours3 <- c("grey90", "#ffc107",  "#ff5722", "#009e73", "#c2185b", "#f959b6")[1:n_clones]
      names(clone_colours2) <- 1:n_clones

      clone_colours <- c("grey90", sort(brewer.pal(n_clones - 1, "Set2")))
      names(clone_colours) <- 1:n_clones


      pdf.date(sprintf("numbat/%s_NB_UMAP_clones.pdf", sample), width=7, height=7)
      print(DimPlot(obj, group.by="clone_opt") +
              scale_colour_manual(values=clone_colours1) +
              theme(aspect.ratio = 1,
                    axis.ticks=element_blank(),
                    axis.text=element_blank()))
      print(DimPlot(obj, group.by="clone_opt") +
              scale_colour_manual(values=clone_colours2) +
              theme(aspect.ratio = 1,
                    axis.ticks=element_blank(),
                    axis.text=element_blank()))
      print(DimPlot(obj, group.by="clone_opt") +
              scale_colour_manual(values=clone_colours3) +
              theme(aspect.ratio = 1,
                    axis.ticks=element_blank(),
                    axis.text=element_blank()))
      print(DimPlot(obj, group.by="clone_opt") +
              scale_colour_manual(values=clone_colours) +
              theme(aspect.ratio = 1,
                    axis.ticks=element_blank(),
                    axis.text=element_blank()))
      print(DimPlot(obj, group.by="subtype_mal") +
              scale_colour_manual(values=singler_subtype_colors_malignant) +
              theme(aspect.ratio = 1,
                    axis.ticks=element_blank(),
                    axis.text=element_blank()))
      dev.off()

      pie_data_cont_raw <- obj[[c("Basal1_score", "Basal2_score", "Classical1_score", "Classical2_score", "clone_opt")]]
      setDT(pie_data_cont_raw)
      pie_data_cont_raw <- melt(pie_data_cont_raw, id.vars="clone_opt", variable.name="subtype_score", value="score")
      pie_data_cont <- pie_data_cont_raw[, .(mean_score = mean(score)), by=.(clone_opt, subtype_score)][order(clone_opt, subtype_score)]
      pie_data_cont$subtype_score <- gsub("_score$", "", pie_data_cont$subtype_score)
      pie_data_cont$clone_opt <- as.factor(pie_data_cont$clone_opt)

      pie_data_cont <- pie_data_cont[clone_opt != 1]

      pdf.date(sprintf("numbat/%s_NB_barplot_clone_by_subtype.pdf", sample), width=8, height=3)
      print(ggplot(pie_data_cont, aes(x=clone_opt, y=mean_score, fill=subtype_score)) +
              geom_col(position="dodge") +
              geom_hline(yintercept=0, linetype='dotted') +
              scale_fill_manual(values=singler_subtype_colors_malignant) +
              ylab("mean subtype score") + xlab("Clone") +
              theme_classic())
      dev.off()


    } else {
      cat("No numbat file ", numbat_clone_file, "\n")
    }
    }
    
    }, error = function(e) {message("Problems plotting ", file)})
}


sample_by_kras <- readRDS("sample_by_kras.rds")
kras_status_by_cell <- readRDS("kras_status_by_cell.rds")

for (single_sample in rownames(sample_by_kras)) {
  obj_sample <- readRDS(sprintf("sample_%s.rds", single_sample))
  obj_sample <- AddMetaData(obj_sample, kras_status_by_cell)
  
  pdf.date(sprintf("%s/%s_UMAP_KRAS_status.pdf", sample, single_sample), width=5, height=5)
  print(DimPlot(obj_sample, group.by="KRAS_status", cols=c("wt only"="blue", "mut"="red"), 
                order=TRUE, raster=FALSE, na.value="grey90") * 
          theme(aspect.ratio = 1,
                axis.ticks=element_blank(),
                axis.text=element_blank()))
  print(DimPlot(obj_sample, group.by="Celltype_final", cols=cell_plot_colours, 
                order=FALSE, raster=FALSE, na.value="grey90") * 
          theme(aspect.ratio = 1,
                axis.ticks=element_blank(),
                axis.text=element_blank()))
  dev.off()
  
  mal_file <- sprintf("sample_%s_mal.rds", single_sample)
  if (file.exists(mal_file)) {
    obj_sample_mal <- readRDS(mal_file)
    obj_sample_mal <- AddMetaData(obj_sample_mal, obj[[c("num_mut", "num_wt", "KRAS_status")]])
    
    pdf.date(sprintf("%s/%s_mal_UMAP_KRAS_status.pdf", sample, single_sample), width=5, height=5)
    print(DimPlot(obj_sample_mal, group.by="KRAS_status", cols=c("wt only"="blue", "mut"="red"), 
                  order=TRUE, raster=FALSE, na.value="grey90") * 
            theme(aspect.ratio = 1,
                  axis.ticks=element_blank(),
                  axis.text=element_blank()))
    print(DimPlot(obj_sample_mal, group.by="subtype", cols=singler_subtype_colors, 
                  order=FALSE, raster=FALSE, na.value="grey90") * 
            theme(aspect.ratio = 1,
                  axis.ticks=element_blank(),
                  axis.text=element_blank()))
    dev.off()
  }
  
}
