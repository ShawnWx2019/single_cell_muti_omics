#####################################################
#       Prj: GCP and MCP                            #
#       Assignments: Prot link meta                 #
#       Author: Shawn Wang @ shawnwang@henu.edu.cn  #
#       Date: Nov 28,2020                           #
#       Location: HENU, Kaifeng, China              #
#####################################################

# Flavonoid biosynthesis pathway ------------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(rvest))
suppressMessages(library(ComplexHeatmap))

## api of ath 2 ath00941
url = "http://rest.kegg.jp/link/ath/ath00941"
## rvest as web spider to get details and convert to dataframe.
ath00941 = read_html(url) %>% 
  html_text() %>% 
  read.table(text = ., sep = "\t") %>% 
  rename("Ath_PathwayID" = "V1","GeneID" = "V2") %>% 
  mutate(
    GeneID = gsub(pattern = "ath:",replacement = "",x = GeneID)
  )

## fetch the details from uniprot
uniprot_spider = function(IDlist){
  final_list = list()
  i = 1
  for (i in 1:length(IDlist)) {
    tryCatch({
      Sys.sleep(1)
      GeneID = IDlist[i]
      url = paste0("https://www.uniprot.org/uniprot/?query=",GeneID,"&fil=organism%3A%22Arabidopsis+thaliana+%28Mouse-ear+cress%29+%5B3702%5D%22&sort=score")
      a = read_html(url) %>% 
        html_nodes('table.grid') %>% 
        html_table() %>% 
        as.data.frame() %>% 
        select("Entry","Entry.name","Protein.names","Gene.names","Organism") %>% 
        mutate(Gene_ID = GeneID)
      final_list[[i]] = a[2,]
      print(paste0(GeneID," Succeed!"))
    },
    error = function(e) {
      print(paste0(GeneID," Failed!"))
    }
    )
  }
  final = bind_rows(final_list)
  return(final)
}


ID_list = ath00941$GeneID

uniprot_result = uniprot_spider(IDlist = ID_list)

## import protein profile matrix
options(scipen = 9)
all_mat <- read.csv("../Protein_mat_all.csv")

## get proteins abundance which in Flavonoid biosynthesis pathway
Flo_mat <- inner_join(
  uniprot_result,
  all_mat,by = c("Gene_ID" = "Accession")
)

## make heatmap matrix
Flo_mat %>% select(Fold.change,P.value) %>% 
  mutate(Fold.change = gsub(pattern = "#DIV/0!",replacement = Inf,Fold.change)) %>% 
  mutate_if(is.character,as.numeric) %>% 
  mutate(
    `log2(ratio)` = case_when(
      Fold.change < 2 & Fold.change > 1 ~ "2 fold up-regular",
      Fold.change < 4 & Fold.change > 2 ~ "2 to 4 fold up-regular",
      Fold.change < 10 & Fold.change > 4 ~ "4 to 10 fold up-regular",
      Fold.change > 10 | Fold.change == Inf ~ "More than 10 fold up-regular",
      Fold.change < 1  ~ "Down regular"
    ),
    P.value = case_when(
      P.value < 0.001 ~ "***",
      P.value < 0.01 & P.value > 0.001 ~ "**",
      P.value < 0.05 & P.value > 0.01 ~ "*",
      TRUE ~ " "
    )
  ) %>% 
  select(`log2(ratio)`,P.value) %>% 
  mutate(gid = rownames(mat2)) %>% 
  column_to_rownames("gid")-> anno_mat

## log2 value heatmap
log2_mat = data.frame(
  row.names = rownames(anno_mat),
  `log2(ratio)` =anno_mat$`log2(ratio)`
) %>% 
  as.matrix()

## protein heatmap main body
p_mat = data.frame(
  row.names = rownames(anno_mat),
  P.value = anno_mat$P.value
) %>% 
  as.matrix()
## log2fc heatmap
log2fc_anno = Heatmap(
  as.matrix(log2_mat),
  name = "log2(GCP/MCP) ",
  show_column_names = T,
  col = structure(c("red","salmon","pink","blue"),names=c("More than 10 fold up-regular","4 to 10 fold up-regular","2 fold up-regular","Down regular")),
  cell_fun = function(j,i,x,y,w,h,col){
    grid.text(p_mat[i,j],x,y)
  }
)
## merge 2 heatmap
ht_list1 = ht_list+log2fc_anno
draw(ht_list1)
