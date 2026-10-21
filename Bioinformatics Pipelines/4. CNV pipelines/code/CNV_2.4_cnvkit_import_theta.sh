

srun hostname

export PATH="/project/desmith_488/tanxinli/SINGULARITY_CACHEDIR/:$PATH"

for i in {1..23} # remove sample90707 with low purity(low acuracy), remove sample90581 no results
do
   echo "Welcome $i times"

cd /scratch2/tanxinli/theta2.out/

## input:
input_theta=(`ls *.BEST.results | sed -n ${i}p`)
base=$(basename ${input_theta} .BEST.results)
prefix="sample"
sample_name=${base#"$prefix"}
tumor_cns="T$sample_name.aligned.duplicates_marked.recalibrated.cns"

# From an existing CNVkit reference
singularity exec /project/desmith_488/tanxinli/SINGULARITY_CACHEDIR/cnvkit_latest.sif cnvkit.py import-theta /scratch2/tanxinli/cnvkit.out/$tumor_cns $input_theta -d /scratch2/tanxinli/theta2.out/

done