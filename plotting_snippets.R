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

# ggpubr provides p-values for comparing groups with ggboxplot
ggboxplot(ToothGrowth, 
          x="dose", 
          y="len", 
          color="dose",
          shape="dose",
          add="jitter"
          ) +
    stat_compare_means(comparisons= list(c("0.5","1"), c("1","2"), c("0.5","2"))) 

# ggviolin from ggpubr, can show significant pvals with */**/***
ggviolin(ToothGrowth, 
         x = "dose", 
         y = "len", 
         fill = "dose",
         add = "boxplot"
         , add.params = list(fill = "white")
         ) +
stat_compare_means(comparisons=list(c("0.5","1"), c("1","2"), c("0.5","2"))
                   ,label="p.signif"
                   ) 

# ComplexHeatmap allows partitioning and more robust annotation
nr1 <- 4; nr2 <- 8; nr3 <- 6; nr <- nr1 + nr2 + nr3
nc1 <- 6; nc2 <- 8; nc3 <- 10; nc <- nc1 + nc2 + nc3

mat <- cbind(rbind(matrix(rnorm(nr1*nc1, mean = 1,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc1, mean = 0,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc1, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc2, mean = 0,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc2, mean = 1,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc2, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc3, mean = 0.5, sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc3, mean = 0.5, sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc3, mean = 1,   sd = 0.5), nr = nr3))
)
mat <- mat[sample(nr, nr), sample(nc, nc)]
rownames(mat) <- paste0("row", seq_len(nr))
colnames(mat) <- paste0("column", seq_len(nc))

Heatmap(mat, name = "mat", row_km = 2)

# ggridges (formerly ggjoy) plots ridgelines
ggplot(iris, aes(x = Sepal.Length, y = Species)) + 
  geom_density_ridges(scale = 0.9)

# ggsurvplot from survminer, support risk tables, confidence intervals and p-vals
fit <- survfit(Surv(time, status) ~ sex, data = lung)

ggsurvplot(fit,  
           size = 1,  
           linetype = "strata", 
           break.time.by = 250, 
           palette = c("#E7B800", "#2E9FDF"), 
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE
)


# venn supports up to 7-set Venn diagrams
x <- as.data.frame(matrix(sample(0:1, 150, replace = TRUE), ncol = 5))
venn(x, snames = "A, B, C, D, E", zcolor = "red, blue, green, yellow, purple")
