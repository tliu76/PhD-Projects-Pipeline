

cd /scratch2/tanxinli/wgs.smoking.out/TM_871719_untrimmed

fastq1=(`find . -name "TM_871719_HJHGYCCXY_*_1.fq" | sed -n ${SLURM_ARRAY_TASK_ID}p`)
DIR=$(dirname ${fastq1})
base=$(basename ${fastq1} _1.fq)
fastq2="${DIR}/$(basename ${fastq1} _1.fq)_2.fq"

output_forward_paired="$base.paired_1.fq"
output_forward_unpaired="$base.unpaired_1.fq"
output_reverse_paired="$base.paired_2.fq"
output_reverse_unpaired="$base.unpaired_2.fq"

java -jar trimmomatic-0.39.jar PE -threads 20 -phred33 $fastq1 $fastq2 /scratch2/tanxinli/wgs.smoking.out/TM_871719_trimmed/$output_forward_paired /scratch2/tanxinli/wgs.smoking.out/TM_871719_trimmed/$output_forward_unpaired /scratch2/tanxinli/wgs.smoking.out/TM_871719_trimmed/$output_reverse_paired /scratch2/tanxinli/wgs.smoking.out/TM_871719_trimmed/$output_reverse_unpaired SLIDINGWINDOW:4:15 MINLEN:36