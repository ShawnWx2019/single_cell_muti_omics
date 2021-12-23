######################################################
#      Prj:   single cell multi omics                #
#      Assignment:   GO enrichment MCP vs GCP        #
#      Author: Shawn Wang   @ shawnwang@henu.edu.cn  #
#      Date:    Aug 02, 2021                         #
######################################################
a <- read.delim("~/15.PostDoc/02.Project/11.Others/foldchange.txt",header = T,sep = "\t")
geneid = a$GeneID

log2fc = log2(a$Fold.change)
names(log2fc) = a$GeneID
#BiocManager::install("ggnewscale")
#BiocManager::install("enrichplot")
#BiocManager::install("ggupset")
BiocManager::install("pathview")
library(clusterProfiler)
library(org.At.tair.db)
library(enrichplot)

## enrich by clusterprofiler
ekegg <- enrichKEGG(geneid,organism = "ath",,pvalueCutoff = 0.05)
ego <- enrichGO(gene = geneid,OrgDb = "org.At.tair.db",keyType = "TAIR",ont = "BP",pvalueCutoff = 0.05)
egosimp <- simplify(ego,cutoff=0.7,by="p.adjust",select_fun=min,measure="Wang")
## output table
tbl.kegg <- as.data.frame(ekegg)
write.table(ekegg,file = "~/15.PostDoc/02.Project/11.Others/kegg_result.xls",row.names = F,sep = "\t",quote = F)

## barplot and dotplot
dot_plot = dotplot(ekegg,showCategory = 24)
bar_plot = barplot(ekegg,showCategory = 24)
net_plot = cnetplot(ekegg,showCategory = 24,foldChange = log2fc)
heat_plot = heatplot(ekegg,foldChange = log2fc)
upset_plot = upsetplot(ekegg)
browseKEGG(ekegg,"ath00020")
ggsave(dot_plot,filename = "~/15.PostDoc/02.Project/11.Others/dot_plot.pdf",width = 10,height = 12)
ggsave(bar_plot,filename = "~/15.PostDoc/02.Project/11.Others/bar_plot.pdf",width = 10,height = 12)
ggsave(net_plot,filename = "~/15.PostDoc/02.Project/11.Others/net_plot.pdf",width = 48,height = 48)
ggsave(net_plot,filename = "~/15.PostDoc/02.Project/11.Others/net_plot_overview.pdf",width = 10,height = 10)
ggsave(heat_plot,filename = "~/15.PostDoc/02.Project/11.Others/heat_plot.pdf",width = 80,height = 10,limitsize = FALSE)
ggsave(upset_plot,filename = "~/15.PostDoc/02.Project/11.Others/up_set.pdf",width = 20,height = 12)
data("geneList")
## GO
dot_plot = dotplot(egosimp,showCategosimpry = 30)
bar_plot = barplot(egosimp,showCategosimpry = 30)
net_plot = cnetplot(egosimp,showCategosimpry = 30,foldChange = log2fc)
heat_plot = heatplot(egosimp,foldChange = log2fc,showCategory = 10)
upset_plot = upsetplot(egosimp)
ggsave(dot_plot,filename = "~/15.PostDoc/02.Project/11.Others/GO_dot_plot.pdf",width = 10,height = 12)
ggsave(bar_plot,filename = "~/15.PostDoc/02.Project/11.Others/GO_bar_plot.pdf",width = 10,height = 12)
ggsave(net_plot,filename = "~/15.PostDoc/02.Project/11.Others/GO_net_plot.pdf",width = 48,height = 48)
ggsave(net_plot,filename = "~/15.PostDoc/02.Project/11.Others/GO_net_plot_overview.pdf",width = 10,height = 10)
ggsave(heat_plot,filename = "~/15.PostDoc/02.Project/11.Others/GO_heat_plot.pdf",width = 80,height = 10,limitsize = FALSE)
ggsave(upset_plot,filename = "~/15.PostDoc/02.Project/11.Others/GO_up_set.pdf",width = 20,height = 12)

gosim.tbl <- data.frame(egosimp)
write.table(gosim.tbl,"~/15.PostDoc/02.Project/11.Others/GO_result.xls",row.names = F,sep = "\t")
