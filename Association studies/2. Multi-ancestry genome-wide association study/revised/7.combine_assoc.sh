#!/bin/bash
#SBATCH --job-name=combine_phenos
#SBATCH --mail-type=END
#SBATCH --account=desmith_128
#SBATCH --partition=epyc-64
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --array=0-71
#SBATCH --output=/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/revised/log/combine_%A_%a.out
#SBATCH --error=/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/revised/log/combine_%A_%a.err

# 注意：合并文件通常不需要 20个CPU 和 60G 内存，我稍微调低了资源以更容易排队，
# 如果你的文件非常巨大，请改回原来的数值。

# 1. 定义表型数组
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
 
  
# 2. 获取当前 Job Array ID 对应的表型
# SLURM_ARRAY_TASK_ID 会从 0 变到 71
CURRENT_PHENO=${PHENOS[$SLURM_ARRAY_TASK_ID]}

echo "Running combination for phenotype: $CURRENT_PHENO"

# 工作目录
cd /project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/result/null_model

module load ver/2506 gcc/14.3.0
module load r/4.5.0

# 3. 运行 R 脚本，并将表型名称作为参数传递
Rscript /project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/revised/7.combine_assoc.R "$CURRENT_PHENO"


