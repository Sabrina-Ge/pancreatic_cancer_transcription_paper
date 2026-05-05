library(Seurat)
library(ggplot2)
library(sagesalad)
library(data.table)
library(RColorBrewer)

rds_file <- "merged.rds"
base_dir <-  "numbat"

obj <- readRDS(rds_file)

singler_subtype_colors <- c("#1C6CAB", "#814C42","#A4C0E5", "#FF7311", "grey40")
names(singler_subtype_colors) <- c("Classical1", "Basal1", "Classical2", "Basal2", "Mixed")
singler_subtype_colors <- singler_subtype_colors[c("Basal1", "Basal2", "Mixed", "Classical1", "Classical2")]

singler_subtype_colors_malignant <- c("#1C6CAB", "#814C42","#A4C0E5", "#FF7311", "grey40", "grey90")
names(singler_subtype_colors_malignant) <- c("Classical1", "Basal1", "Classical2", "Basal2", "Mixed", "non-malignant")
singler_subtype_colors_malignant <- singler_subtype_colors_malignant[c("Basal1", "Basal2", "Mixed", "Classical1", "Classical2", "non-malignant")]

celltype_colour_file <- "data/celltype_colours.csv"
celltype_colours <- fread(celltype_colour_file)
cell_plot_order <- celltype_colours$celltype
cell_plot_colours <- celltype_colours$hexcode
names(cell_plot_colours) <- cell_plot_order

#### Extract numbat clone data ####

# Files output by plot_sample.R
all_numbat_clone_files <- list.files(path="numbat", 
                                     pattern="sample_.*_clone_data.csv", 
                                     full.names=TRUE)
clone_data_list <- lapply(all_numbat_clone_files, read.csv, row.names=1)
clone_data <- do.call(rbind.data.frame, clone_data_list)

write.csv(clone_data, datestamp("merged_cell_to_clone.csv"))

obj <- AddMetaData(obj, clone_data)

#### Tables ####

# clones per sample w/ met vs primary info, num cells total
clone_per_sample_data <- obj[[c("clone_opt", "Sample", "type")]]
setDT(clone_per_sample_data)

clone_per_sample <- clone_per_sample_data[, .(num_clone=max(clone_opt)), by=.(Sample, type)]

fwrite(clone_per_sample, datestamp("merged_clones_per_sample_data.csv"))

# num cells per clone per sample
num_per_clone <- clone_per_sample_data[, .N, by=.(Sample, clone_opt)]
num_per_clone <- num_per_clone[order(Sample, clone_opt)]

fwrite(num_per_clone, datestamp("merged_clone_counts.csv"))

# intersect between clone 1 and non-malignant for samples
normal_overlap_data <- obj[[c("clone_opt", "Sample", "Celltype_final")]]
normal_overlap_data$normal_by_clone <- ifelse(normal_overlap_data$clone_opt == 1, "non_aneuploid", "aneuploid")
normal_overlap_data$normal_by_celltype <- ifelse(normal_overlap_data$Celltype_final == "Malignant", "malignant", "non_malignant")

setDT(normal_overlap_data)
table(normal_overlap_data[, .(normal_by_clone, normal_by_celltype)])

dated_file <- datestamp("merged_clone_celltype_normal_overlap.csv")
fwrite(table(normal_overlap_data[, .(normal_by_clone, normal_by_celltype)]), dated_file)

for (smpl in unique(normal_overlap_data$Sample)) {
  sub_tbl <- table(normal_overlap_data[Sample == smpl, .(normal_by_clone, normal_by_celltype)])
  if (length(sub_tbl) > 0) {
    write(c("", smpl), dated_file, append=TRUE)
    fwrite(sub_tbl, dated_file, append=TRUE)
  }
}

#### Plots ####

obj_mal <- subset(obj, Celltype_final == "Malignant")

samples_with_many_clones <- clone_per_sample[num_clone > 2, Sample]

all_clone_data_raw <- obj_mal[[c("clone_opt", "Sample", "subtype")]]
setDT(all_clone_data_raw)
all_clone_data <- all_clone_data_raw[clone_opt != 1 & Sample %in% samples_with_many_clones]
all_clone_data$early_vs_late <- ifelse(all_clone_data$clone_opt == 2, "early", "late")

counts_for_prop_limits <- as.data.frame(table(all_clone_data[, .(Sample, early_vs_late)]))
samples_to_remove <- counts_for_prop_limits[counts_for_prop_limits$Freq < 10, "Sample"]

all_clone_data <-  all_clone_data[!Sample %in% samples_to_remove]

pdf.date(sprintf("%s/merged_subtype_vs_clone.pdf", base_dir), width=7.5, height=5)
ggplot(all_clone_data, aes(x=early_vs_late, fill=subtype)) + 
  geom_bar() +
  scale_fill_manual(values=singler_subtype_colors_name) + 
  theme_bw() + theme(panel.grid=element_blank(), strip.background = element_rect(fill="white", colour="black")) + 
  facet_wrap(~Sample, nrow=3)
ggplot(all_clone_data, aes(x=early_vs_late, fill=subtype)) + 
  geom_bar(position="fill") +
  scale_fill_manual(values=singler_subtype_colors_name) + 
  theme_bw() + theme(panel.grid=element_blank(), strip.background = element_rect(fill="white", colour="black")) + 
  facet_wrap(~Sample, nrow=3)
dev.off()

all_clone_data$grouping <- ifelse(all_clone_data$subtype %in% c("Basal1", "Basal2"), "Basal", all_clone_data$subtype)
all_clone_data$grouping <- ifelse(all_clone_data$subtype %in% c("Classical1", "Classical2"), "Classical", all_clone_data$grouping)

all_clone_data[, clone_count := .N, by=.(Sample, early_vs_late)]

all_clone_prop <- unique(all_clone_data[, .(subtype_count=.N, subtype_prop=.N/clone_count), by=.(Sample, early_vs_late, subtype)])

for (subtype in c("Classical1", "Classical2", "Basal1", "Basal2")){
  
  subtype_clone_prop <- all_clone_prop[subtype == subtype]
  
  subtype_clone_prop_short <- dcast(subtype_clone_prop[, .(Sample, early_vs_late, subtype_prop)], Sample ~ early_vs_late)
  subtype_clone_prop_short$early[is.na(subtype_clone_prop_short$early)] <- 0  
  subtype_clone_prop_short$late[is.na(subtype_clone_prop_short$late)] <- 0 
  
  subtype_clone_prop_short$diff <- subtype_clone_prop_short$late - subtype_clone_prop_short$early
  subtype_clone_prop_short$direction <- "similar"
  subtype_clone_prop_short$direction <- ifelse(subtype_clone_prop_short$diff > 0.15, "increase", subtype_clone_prop_short$direction)
  subtype_clone_prop_short$direction <- ifelse(subtype_clone_prop_short$diff < -0.15, "decrease", subtype_clone_prop_short$direction)
  
  subtype_clone_prop_plot <- merge(subtype_clone_prop, subtype_clone_prop_short)
  subtype_clone_prop_plot <- merge(melt(subtype_clone_prop_short, measure.vars=c("early", "late"), variable.name="early_vs_late", value.name="subtype_prop"),
                                   subtype_clone_prop, all.x=TRUE)
  subtype_clone_prop_plot$Sample <- reorder(subtype_clone_prop_plot$Sample, subtype_clone_prop_plot$diff)
  
  
  pdf.date(sprintf("%s/merged_scatter_subtype_vs_clone_%s.pdf", base_dir, subtype), width=4, height=4)
  print(ggplot(subtype_clone_prop_plot, aes(x=early_vs_late, y=subtype_prop)) + 
          geom_point() + geom_line(aes(group=Sample)) +
          ylab(sprintf("%s proportion", subtype)) +
          theme_bw() + theme(panel.grid=element_blank(), strip.background = element_rect(fill="white", colour="black")))
  print(ggplot(subtype_clone_prop_plot, aes(x=early_vs_late, y=subtype_prop, colour=diff)) + 
          geom_point() + geom_line(aes(group=Sample)) +
          scale_colour_gradient2(low="blue", mid="lightgrey", high="red") +
          ylab(sprintf("%s proportion", subtype))  +
          theme_bw() + theme(panel.grid=element_blank(), strip.background = element_rect(fill="white", colour="black")))
  print(ggplot(subtype_clone_prop_plot, aes(x=early_vs_late, y=subtype_prop, colour=direction)) + 
          geom_point() + geom_line(aes(group=Sample)) +
          scale_colour_manual(values = c(similar="lightgrey", increase="salmon2", decrease="cornflowerblue")) +
          ylab(sprintf("%s proportion", subtype))  + xlab("clones") +
          theme_bw() + theme(panel.grid=element_blank(), strip.background = element_rect(fill="white", colour="black"))) 
  dev.off()
  
}

