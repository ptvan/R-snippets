library(viridis)
library(ComplexHeatmap)
library(cowplot)
library(ggridges)
library(survminer)
library(ggpubr)
library(venn)

# viridis provides colorblindness-safe palettes 
ggplot(mtcars, aes(wt, mpg)) + 
  geom_point(size=4, aes(colour = factor(cyl))) +
  scale_color_viridis(discrete=TRUE) +
  theme_bw()

# ComplexHeatmap allows partitioning and more robust annotation


# ggridges (formerly ggjoy) plots ridgelines


# ggsurvplot from survminer


# ggviolin from ggpubr


# venn supports up to 7-set Venn diagrams
