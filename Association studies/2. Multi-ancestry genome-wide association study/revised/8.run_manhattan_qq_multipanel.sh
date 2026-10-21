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
set -euo pipefail

BASE_DIR="/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/result/null_model"
RSCRIPT="/project2/desmith_128/tanxin/scratch1_new/mp2prt_gwas/filtered/code/genesis/0323_SV/revised/8.run_all_phenotypes_multipanel.R"

# -------- phenotype list --------
PHENOS=(
 "n.all_INT"                    "n.nonig_INT"                  "n.ig_INT"                     "n.total.least.rag_INT"       
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


PHENO="${PHENOS[$((SLURM_ARRAY_TASK_ID - 1))]}"

Rscript "${RSCRIPT}" "${PHENO}"




