#!/bin/bash
#SBATCH --job-name=mergeplot
#SBATCH --mail-type=END
#SBATCH --account=desmith_128
#SBATCH --partition=epyc-64
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --array=1
#SBATCH --output=/scratch1/tanxinli/output/mergeplot_%A_%a.out
#SBATCH --error=/scratch1/tanxinli/error/mergeplot_%A_%a.err



module load imagemagick 2>/dev/null || echo "Trying to use system ImageMagick..."


BASE_DIR="/scratch1/tanxinli/mp2prt_gwas/filtered/result/pcrelate_test10/result0116"
OUTPUT_PDF_DIR="${BASE_DIR}/Summary_Plots"

mkdir -p ${OUTPUT_PDF_DIR}

PHENOS=(
  "n.all.del_INT"
  "n.ig.del_INT"
  "n.nonig.del_INT"
  "n.total.least.rag.del_INT"
  "n.total.both.rag.del_INT"
  "n.nonig.least.rag.del_INT"
  "n.nonig.both.rag.del_INT"
  "n.ig.least.rag.del_INT"
  "n.ig.both.rag.del_INT"
  "prop.total.least.rag.del_INT"
  "prop.total.both.rag.del_INT"
  "prop.nonig.least.rag.del_INT"
  "prop.nonig.both.rag.del_INT"
  "prop.ig.least.rag.del_INT"
  "prop.ig.both.rag.del_INT"
)




LIST_QQ=""
LIST_MAN=""

for TRAIT in "${PHENOS[@]}"; do
    IMG_QQ="${BASE_DIR}/${TRAIT}/plots/qqplot_${TRAIT}.png"
    IMG_MAN="${BASE_DIR}/${TRAIT}/plots/manhattan_${TRAIT}.png"

    if [[ -f "$IMG_QQ" ]]; then
        LIST_QQ="$LIST_QQ $IMG_QQ"
    else
        echo "Warning: Missing QQ plot for $TRAIT"
    fi

    if [[ -f "$IMG_MAN" ]]; then
        LIST_MAN="$LIST_MAN $IMG_MAN"
    else
        echo "Warning: Missing Manhattan plot for $TRAIT"
    fi
done

convert $LIST_QQ -quality 100 "${OUTPUT_PDF_DIR}/Summary_QQ_Plots.pdf"

convert $LIST_MAN -quality 100 "${OUTPUT_PDF_DIR}/Summary_Manhattan_Plots.pdf"

