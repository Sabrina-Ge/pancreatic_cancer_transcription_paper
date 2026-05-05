library(numbat)
library(Seurat)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

sample <- args[1]
outdir <- sprintf("./sample_%s", sample)
dir.create(outdir, showWarnings=FALSE)

obj <- readRDS("merged.rds")
sample_cells <- rownames(obj@meta.data)[obj@meta.data$sample == sample]
count_mat <- obj@assays$RNA@counts[, sample_cells]

df_allele <- fread(sprintf('numbat/allele_counts_v20240109/%s_allele_counts.tsv.gz', sample))

colnames(count_mat) <- gsub("_.*", "", colnames(count_mat))

run_numbat(count_mat, ref_hca, df_allele, genome="hg38", t=1e-5, ncores=1, max_iter=5, plot=TRUE, out_dir=outdir)

cat(sample, " - numbat.R run complete")
