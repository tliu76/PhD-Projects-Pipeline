#!/bin/bash
#SBATCH --job-name=assoc_chr
#SBATCH --mail-type=END
#SBATCH --account=desmith_128
#SBATCH --partition=epyc-64
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --array=1-22
#SBATCH --output=/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/revised/log/assoc_chr_%A_%a.out
#SBATCH --error=/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/revised/log/assoc_chr_%A_%a.err


cd /project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/revised/result

module load ver/2506 gcc/14.3.0
module load r/4.5.0

# 让 BLAS/LAPACK 用多线程
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

CHR=${SLURM_ARRAY_TASK_ID}

echo "Running single-variant assoc for chr ${CHR} ..."

PHENOS=(  "n.all_INT"                    "n.nonig_INT"                  "n.ig_INT"                     "n.total.least.rag_INT"       
"n.total.both.rag_INT"         "n.nonig.least.rag_INT"        "n.ig.least.rag_INT"           "n.nonig.both.rag_INT"        
 "n.ig.both.rag_INT"            "prop_total_least_rag.sv_INT"  "prop_ig_least_rag.sv_INT"     "prop_nonig_least_rag.sv_INT" 
 "n_raghit_only_one_bkp_INT"    "n.ig.rag_only_INT"            "n.nonig.rag_only_INT"         "prop_raghit_only_one_bkp_INT"
 "prop_ig_rag_only_INT"         "prop_nonig_rag_only_INT"      "n.all.del_INT"                "n.ig.del_INT"                
 "n.nonig.del_INT"              "n.total.least.rag.del_INT"    "n.ig.least.rag.del_INT"       "n.nonig.least.rag.del_INT"   
 "prop_total_least_rag.del_INT" "prop_ig_least_rag.del_INT"    "prop_nonig_least_rag.del_INT" "n_raghit_only_del_INT"       
 "n.ig.rag_only.del_INT"        "n.nonig.rag_only.del_INT"     "prop_raghit_only_del_INT"     "prop_ig_rag_only_del_INT"    
 "prop_nonig_rag_only_del_INT"  "n.all.dup_INT"                "n.ig.dup_INT"                 "n.nonig.dup_INT"             
 "n.total.least.rag.dup_INT"    "n.ig.least.rag.dup_INT"       "n.nonig.least.rag.dup_INT"    "prop_total_least_rag.dup_INT"
 "prop_ig_least_rag.dup_INT"    "prop_nonig_least_rag.dup_INT" "n_raghit_only_dup_INT"        "n.ig.rag_only.dup_INT"       
 "n.nonig.rag_only.dup_INT"     "prop_raghit_only_dup_INT"     "prop_ig_rag_only_dup_INT"     "prop_nonig_rag_only_dup_INT" 
 "n.all.inv_INT"                "n.ig.inv_INT"                 "n.nonig.inv_INT"              "n.total.least.rag.inv_INT"   
 "n.ig.least.rag.inv_INT"       "n.nonig.least.rag.inv_INT"    "prop_total_least_rag.inv_INT" "prop_ig_least_rag.inv_INT"   
 "prop_nonig_least_rag.inv_INT" "n_raghit_only_inv_INT"        "n.ig.rag_only.inv_INT"        "n.nonig.rag_only.inv_INT"    
 "prop_raghit_only_inv_INT"     "prop_ig_rag_only_inv_INT"     "prop_nonig_rag_only_inv_INT"  "n.all.tra_INT"               
 "n.ig.tra_INT"                 "n.nonig.tra_INT"              "n.total.least.rag.tra_INT"    "n.ig.least.rag.tra_INT"      
"n.nonig.least.rag.tra_INT"    "prop_total_least_rag.tra_INT" "prop_ig_least_rag.tra_INT"    "prop_nonig_least_rag.tra_INT"
    
)
R_SCRIPT="/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/revised/6.single_assoc.R" 


# 循环遍历每一个 phenotype
for pheno in "${PHENOS[@]}"; do
    
    echo "------------------------------------------------------"
    echo "Running GWAS for phenotype: ${pheno} on Chr ${CHR}..."
    echo "------------------------------------------------------"
    
    # 运行 R 脚本，传入 染色体号 和 phenotype名
    Rscript ${R_SCRIPT} ${CHR} ${pheno}
    
    # 检查上一步是否成功
    if [ $? -eq 0 ]; then
        echo "Successfully finished ${pheno} on Chr ${CHR}"
    else
        echo "ERROR: Failed running ${pheno} on Chr ${CHR}"
    fi

done
