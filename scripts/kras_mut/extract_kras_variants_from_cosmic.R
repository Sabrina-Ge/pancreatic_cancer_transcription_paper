library(data.table)
library(stringr)

cosmic <- fread("data/hg38_cosmic91.txt.gz")
colnames(cosmic) <- c("chr", "start", "stop", "ref", "alt", "info")

kras_chr <- "12"
kras_coords <- c("25204789", "25250936")

cosmic_kras <- cosmic[chr == kras_chr & start >= kras_coords[1] & stop <= kras_coords[2]]
cosmic_kras_panc <- cosmic_kras[grep("pancreas", cosmic_kras$info)]

cosmic_kras_panc$panc_count <- as.numeric(str_extract(cosmic_kras_panc$info, "[=,](\\d+)\\(pancreas\\)", group=1))
cosmic_kras_panc <- cosmic_kras_panc[order(-panc_count)]
cosmic_kras_panc_top

cosmic_kras_panc_top <- cosmic_kras_panc[panc_count >= 10]
cosmic_kras_panc_top

sink("cosmic_kras_panc_atleast10count.vcf")
cat("##fileformat=VCFv4.2\n")
cat("##reference=data/cellranger-GRCh38-3.0.0-KRAS/fasta/genome.fa\n")
cat("##contig=<ID=12,length=133275309>\n")
cat('##ALT=<ID=*,Description="Represents allele(s) other than observed.">\n')
cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")
for (i in 1:nrow(cosmic_kras_panc_top)) {
  cat(cosmic_kras_panc_top[i, chr], "\t", cosmic_kras_panc_top[i, start], "\t.\t", 
      cosmic_kras_panc_top[i, ref], "\t", cosmic_kras_panc_top[i, alt], "\t0\t.\t.\t.\n", sep="")
}
sink()

