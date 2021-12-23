#####################################################
#       Prj: GCP and MCP                            #
#       Assignments: PCA and heatmap                #
#       Author: Shawn Wang @ shawnwang@henu.edu.cn  #
#       Date: Nov 28,2020                           #
#       Location: HENU, Kaifeng, China              #
#####################################################
# PCA ---------------------------------------------------------------------
x <- read.csv("ExpMat.csv")
head(x)

library(tidyverse)

x %>% 
  column_to_rownames("Name") %>% 
  as.matrix()-> mat
group = data.frame(
  row.names  = colnames(x)[-1],
  group = rep(c("GCP","MCP"),each = 11)
)

library(PCAtools)

pca_result <- pca(mat = mat,metadata = group,scale = T,removeVar = 0.1)
pdf("Final_File 2/PCA_merge_pos.neg.pdf",height = 8,width = 7)
biplot(pca_result,colby  = "group",
       colkey = c('GCP' = 'forestgreen', 'MCP' = 'purple'),
       # ellipse config
       ellipse = TRUE,
       ellipseFill = TRUE,
       ylim = c(-100,100),
      # ellipseLineSize = 1.0,
       legendPosition = 'bottom',
       legendLabSize = 12, legendIconSize = 8.0,
       title = "PCA biplot",hline = 0,vline = 0,
       subtitle = "Normlized peak result")
dev.off()


# Heatmap ---------------------------------------------------------------------


library(circlize)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
ht = ComplexHeatmap::Heatmap(
  matrix = t(scale(t(mat))),
  col =  colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),## set boundary of legend
  cluster_columns = F,
  show_row_names = F,
  name = " ",
  column_split = rep(c("GCP","MCP"),each = 11)
)

pdf("Final_File 2/Heatmap.pdf",height = 8,width = 7)
draw(ht)
dev.off()