#!/usr/bin/env Rscript

# Usage: kras_parse_MM_to_rds.R sample_name directory

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 2) {
   stop("Incorrect number of arguments")
}

cat(paste0("Running kras_parse_MM_to_rds.R ", args[1], args[2], "\n"))

sample <- args[1]
setwd(args[2])

barcodes <- scan("barcodes.tsv", what="character")
vars <- scan("features_vartrix.tsv", what="character")

mut_data <- t(as.matrix(Matrix::readMM(file="mut_vartrix.tsv")))
rownames(mut_data) <- barcodes
colnames(mut_data) <- vars

mutation_data <- cbind(vars, colSums(mut_data))
colnames(mutation_data) <- c("mutation", "count")
write.csv(mutation_data, sprintf("kras_%s_mutation_total_counts.csv", sample), quote=FALSE, row.names=FALSE)

mut_col_to_keep <- which.max(colSums(mut_data))

wt_data <- t(as.matrix(Matrix::readMM(file="wt_vartrix.tsv")))

data <- cbind(mut_data[, mut_col_to_keep], wt_data[, mut_col_to_keep])
colnames(data) <- c("num_mut", "num_wt")
rownames(data) <- barcodes

saveRDS(data, paste0(sample, "_kras_amp_mat.rds"))

cat("Complete")
