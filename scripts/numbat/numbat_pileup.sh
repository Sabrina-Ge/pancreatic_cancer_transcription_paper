module load R/4.3.0
module load miniconda
startconda
conda activate tools/numbat/cellsnp_lite_conda
module load samtools

ncores=1
sample=$1
bam_path="$2"
barcode_path="$3"

Rscript ~/R/4.3/numbat/bin/pileup_and_phase.R \
    --label $sample \
    --samples $sample \
    --bams "$bam_path" \
    --barcodes "$barcode_path" \
    --outdir ./sample_$sample \
    --gmap tools/numbat/eagle/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
    --snpvcf tools/numbat/references/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
    --paneldir tools/numbat/references/1000G_hg38 \
    --ncores $ncores
