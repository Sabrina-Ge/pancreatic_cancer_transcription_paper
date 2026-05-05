#!/bin/bash

# Run script in folder with directories of the samples

module load R/4.3.0 

sample_file=$1
kras_file="cosmic_kras_panc_atleast10count.vcf"

for sample in `cat $sample_file`; do
  cd $sample
  mkdir kras_vartrix/
  gunzip -c $sample/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > kras_vartrix/barcodes.tsv
  echo "Running Vartrix on $sample with vcf file $kras_file"
  cd kras_vartrix/
  tools/vartrix/vartrix-v1.1.3-x86_64-linux/vartrix -v "$kras_file" \
    --bam ../$sample/outs/possorted_genome_bam.bam \
    --fasta data/cellranger-GRCh38-3.0.0-KRAS/fasta/genome.fa \
    --cell-barcodes barcodes.tsv \
    --umi \
    --scoring-method coverage \
    --out-matrix mut_vartrix.tsv \
    --ref-matrix wt_vartrix.tsv \
    --out-variants features_vartrix.tsv \
    --log-level info

  kras_parse_MM_to_rds.R $sample .

  cd ../../ 
done
