bam_paths_file="data/numbat_sample_to_bam_path.txt"
barcode_paths_file="data/numbat_sample_to_barcodes_path.txt"

echo Using samples from $bam_paths_file and $barcode_paths_file

while read sample_bam bam_path <&3 && read sample_barcode barcode_path <&4 ; do
  if [[ $sample_bam != $sample_barcode ]] ; then
    echo Mismatch in samples in reference files: $sample_bam $sample_barcode
    continue
  fi
  echo Submitting job with $sample_bam $bam_path $barcode_path  
  qsub -P pcsi -cwd -b y -j y -o ./ -l h_vmem=64G -l h_rt=2:0:0:0 -N nb_$sample_bam numbat_pileup.sh $sample_bam $bam_path $barcode_path 
done 3<"$bam_paths_file" 4<"$barcode_paths_file"
