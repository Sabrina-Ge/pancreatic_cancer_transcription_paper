bam_paths_file="data/numbat_sample_to_bam_path.txt"

echo Using samples from $bam_paths_file

while read sample_bam bam_path; do
  echo Submitting numbat.R job with $sample_bam
  qsub -P pcsi -cwd -b y -j y -o ./logs/ -l h_vmem=64G -l h_rt=2:0:0:0 -N nbR_$sample_bam "module load R/4.3.0;Rscript numbat.R $sample_bam" 
done <"$bam_paths_file"
