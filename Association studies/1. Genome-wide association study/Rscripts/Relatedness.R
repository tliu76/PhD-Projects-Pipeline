.libPaths("/home1/tanxinli/R_lib")
library(ggplot2)
pdf("hist_relatedness.pdf")
#relatedness = read.table("zoom_pihat.genome", header=T)

relatedness = read.table("relatedness.genome", header=T)

ggplot(relatedness, aes(x = PI_HAT)) +
  geom_histogram(
    bins = 10,                # equivalent to breaks = 10
    fill = "skyblue",
    color = "white"
  ) +
  labs(
    title = "Histogram of PI_HAT",
    x = "PI_HAT",
    y = "Frequency"
  ) +
  theme_minimal()
dev.off() 

pdf("high_relatedness.pdf")
highrelate = read.table("high_relatedness_file", header=T)
ggplot(highrelate, aes(x = PI_HAT)) +
  geom_histogram(
    bins = 10,                # equivalent to breaks = 10
    fill = "skyblue",
    color = "white"
  ) +
  labs(
    title = "Histogram of PI_HAT > 0.2",
    x = "PI_HAT",
    y = "Frequency"
  ) +
  theme_minimal()
  
dev.off() 
