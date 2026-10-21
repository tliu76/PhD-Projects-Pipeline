#!/usr/bin/env Rscript
.libPaths("/home1/tanxinli/R_lib")

# Install required packages if needed
if (!require("png", quietly = TRUE)) install.packages("png")
if (!require("grid", quietly = TRUE)) install.packages("grid")
if (!require("gridExtra", quietly = TRUE)) install.packages("gridExtra")

library(png)
library(grid)
library(gridExtra)
library(ggplot2)

base_dir <- "/scratch1/tanxinli/mp2prt_gwas/filtered/result/pcrelate_test10/result0116"
outdir <- file.path(base_dir, "multipanel")
dir.create(outdir, showWarnings = FALSE)

phenos <- c(
  "n.total.least.rag.del_INT",
  "n.nonig.least.rag.del_INT",
  "n.ig.least.rag.del_INT",
  "prop.total.least.rag.del_INT",
  "prop.nonig.least.rag.del_INT",
  "prop.ig.least.rag.del_INT"
)

plot_subdir <- "plots2"  # 如果你用的是 plots，就改成 "plots"

# Function to create multi-panel figure
create_multipanel <- function(file_paths, output_path) {
  # Read all PNG files
  imgs <- lapply(file_paths, readPNG)
  
  # Convert to grobs
  grobs <- lapply(imgs, function(img) {
    rasterGrob(img, interpolate = TRUE)
  })
  
  # Arrange in 2x3 grid
  combined <- grid.arrange(
    grobs[[1]], grobs[[2]], grobs[[3]],
    grobs[[4]], grobs[[5]], grobs[[6]],
    ncol = 3, nrow = 2
  )
  
  # Save to file
  ggsave(
    filename = output_path,
    plot = combined,
    width = 18,
    height = 12,
    dpi = 300,
    bg = "white"
  )
}

# =========================
# 1. Manhattan multi-panel
# =========================
man_files <- file.path(
  base_dir,
  phenos,
  plot_subdir,
  paste0("manhattan_", phenos, ".png")
)

cat("Reading Manhattan plots:\n")
print(man_files)

create_multipanel(
  man_files,
  file.path(outdir, "Figure_Manhattan_6phenotypes.png")
)

# =========================
# 2. QQ multi-panel
# =========================
qq_files <- file.path(
  base_dir,
  phenos,
  plot_subdir,
  paste0("qqplot_", phenos, ".png")
)

cat("Reading QQ plots:\n")
print(qq_files)

create_multipanel(
  qq_files,
  file.path(outdir, "Figure_QQ_6phenotypes.png")
)

cat(">>> Multi-panel figures saved to:", outdir, "\n")