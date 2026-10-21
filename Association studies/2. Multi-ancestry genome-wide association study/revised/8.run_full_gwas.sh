#!/bin/bash
#SBATCH --job-name=plot_assoc
#SBATCH --mail-type=END
#SBATCH --account=desmith_128
#SBATCH --partition=epyc-64
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --array=1-72
#SBATCH --output=/scratch1/tanxinli/output/plot_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/plot_%A_%a.err



module load ver/2506 gcc/14.3.0
module load r/4.5.0
# glpk is needed by igraph
module load glpk/5.0 || true

# ensure glpk/igraph shared objects are findable
export LD_LIBRARY_PATH="/home1/tanxinli/R_lib:/home1/tanxinli/R_lib/lib64:${LD_LIBRARY_PATH}"
if [[ -n "${EBROOTGLPK}" ]]; then
  export LD_LIBRARY_PATH="${EBROOTGLPK}/lib:${LD_LIBRARY_PATH}"
fi

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
 
 

# 2. 获取当前表型
# SLURM_ARRAY_TASK_ID 是 1 到 72
# 数组索引是 0 到 71，所以要减 1
IDX=${SLURM_ARRAY_TASK_ID}
TRAIT=${PHENOS[$IDX-1]}


Rscript /project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/revised/8.assoc_multi.R "${TRAIT}"
