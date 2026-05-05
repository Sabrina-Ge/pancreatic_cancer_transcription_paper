#!/bin/bash

module load cellranger/3.1.0

transcriptome=tools/cellranger_ref/GRCh38-3.0.0

for id in $(cat KRAS_samples.txt); do
  cd $id
  cellranger count --id=$id \
                 --transcriptome=$transcriptome \
                 --fastqs=fastq \
                 --nosecondary \
                 --localcores=1 \
                 --mempercore=64
  cd ..              
done 
