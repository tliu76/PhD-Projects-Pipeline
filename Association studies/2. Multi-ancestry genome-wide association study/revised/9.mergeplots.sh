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

# ====================================================
# 这是一个后期处理脚本，请在所有画图任务完成后运行一次
# 运行方式: bash 8.merge_plots.sh
# ====================================================

# 1. 加载 ImageMagick (大多数集群都有，如果没有请联系管理员或使用 R 合并)
# 尝试加载常见的名称，如果不确定可以用 'module spider imagemagick' 查一下
module load imagemagick 2>/dev/null || echo "Trying to use system ImageMagick..."

# 2. 定义目录和表型 (直接复制之前的列表)
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

echo "开始合并图片..."

# ====================================================
# 模式 A: 所有 QQ 图合成一个 PDF，所有 Manhattan 合成一个 PDF
# ====================================================

# 准备文件列表
LIST_QQ=""
LIST_MAN=""

for TRAIT in "${PHENOS[@]}"; do
    # 定位之前的图片路径
    IMG_QQ="${BASE_DIR}/${TRAIT}/plots/qqplot_${TRAIT}.png"
    IMG_MAN="${BASE_DIR}/${TRAIT}/plots/manhattan_${TRAIT}.png"

    # 检查文件是否存在
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

# 使用 convert 命令合并 (Quality 100 保证清晰度)
echo "正在生成 Summary_QQ_Plots.pdf ..."
convert $LIST_QQ -quality 100 "${OUTPUT_PDF_DIR}/Summary_QQ_Plots.pdf"

echo "正在生成 Summary_Manhattan_Plots.pdf ..."
convert $LIST_MAN -quality 100 "${OUTPUT_PDF_DIR}/Summary_Manhattan_Plots.pdf"


# ====================================================
# 模式 B: (可选) 生成“左边QQ-右边Manhattan”的对照报告
# ====================================================
echo "正在生成 Combined_Report.pdf (Side-by-Side) ..."

TEMP_PAGES=""

for TRAIT in "${PHENOS[@]}"; do
    IMG_QQ="${BASE_DIR}/${TRAIT}/plots/qqplot_${TRAIT}.png"
    IMG_MAN="${BASE_DIR}/${TRAIT}/plots/manhattan_${TRAIT}.png"
    
    PAGE_OUT="${OUTPUT_PDF_DIR}/temp_${TRAIT}.png"
    
    if [[ -f "$IMG_QQ" && -f "$IMG_MAN" ]]; then
        # +append 表示水平拼接 (QQ在左，Manhattan在右)
        # -resize 统一高度，防止图片大小不一
        convert "$IMG_QQ" "$IMG_MAN" +append "${PAGE_OUT}"
        TEMP_PAGES="$TEMP_PAGES ${PAGE_OUT}"
    fi
done

# 将所有拼接好的单页合并成一个 PDF
convert $TEMP_PAGES -quality 100 "${OUTPUT_PDF_DIR}/All_Phenotypes_Report.pdf"

# 清理临时文件
rm ${OUTPUT_PDF_DIR}/temp_*.png

echo "------------------------------------------------"
echo "全部完成！结果保存在: ${OUTPUT_PDF_DIR}"
echo "1. Summary_QQ_Plots.pdf       (所有QQ图)"
echo "2. Summary_Manhattan_Plots.pdf (所有Manhattan图)"
echo "3. All_Phenotypes_Report.pdf   (左右对照合并版)"
echo "------------------------------------------------"