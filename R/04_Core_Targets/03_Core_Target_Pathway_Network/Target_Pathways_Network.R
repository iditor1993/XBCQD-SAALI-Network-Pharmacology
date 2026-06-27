# setwd() removed - use RStudio Project or set working directory to project root                                #设置工作目录

# 清除当前环境中的所有变量
rm(list = ls())


library(dplyr)
library(openxlsx)
library(tidyverse)  




#------------------------创建文件夹、定义输入\输出文件路径----------------------
Filter_fold="results/04_Core_Targets/04_Core_Target_Pathway_Network"


#创建文件夹路径Filter_fold
if (!dir.exists(Filter_fold)) {
  dir.create(Filter_fold, recursive = TRUE)
  cat("文件夹已创建:", Filter_fold, "\n")
} else {
  cat("文件夹已存在:", Filter_fold, "\n")
}



#输入文件路径
Final_core_target_path <-"results/04_Core_Targets/01_Core_Targets.csv"
KEGG_result_path <- "results/03_Intersection_Targets/02_GO_KEGG_Enrichment/02_KEGG/KEGG_enrich.csv"
Tatget_PPi_path <- "results/04_Core_Targets/03_Core_PPI/Target_PPi.csv"


#输出文件路径
type_path <- "results/04_Core_Targets/04_Core_Target_Pathway_Network/type.xlsx"
veen_cytoscape_path <- "results/04_Core_Targets/04_Core_Target_Pathway_Network/Target_Pathways_network.xlsx"
#-------------------------------------------------------------------------------






#---------------------------整理出kegg_target_network---------------------------
kegg_enrich <- read.csv(KEGG_result_path, stringsAsFactors = FALSE)    #读取文件

#读取核心靶点
Final_core_target <- read.csv(Final_core_target_path,header=TRUE, check.names = FALSE)

#筛选出pvalue<=0.05显著性前30的kegg通路
kegg_enrich_filter <- kegg_enrich[kegg_enrich$pvalue <= 0.05, ]
kegg_top30 <- kegg_enrich_filter %>%
  dplyr::arrange(pvalue) %>%       #升序排序
  dplyr::slice(1:30)           #选出前30


# 分割Gene_symbol列，并重塑数据框  
kegg_target_network <- kegg_top30%>%  
  rowwise() %>%  
  mutate(Gene_symbol = strsplit(as.character(Gene_symbol), "/")) %>%  
  unnest(Gene_symbol) %>%  
  ungroup() %>%  
  dplyr::select(KEGGID, Gene_symbol)

kegg_target_network <- merge(kegg_target_network,Final_core_target, by = "Gene_symbol")%>%  
  dplyr::select(KEGGID, Gene_symbol)

names(kegg_target_network)[names(kegg_target_network) == "KEGGID"] <- "from_node"  #修改列名
names(kegg_target_network)[names(kegg_target_network) == "Gene_symbol"] <- "to_node"  #修改列名

#写入文件
write.xlsx(kegg_target_network,veen_cytoscape_path,rowNames = FALSE)
#-------------------------------------------------------------------------------






#--------------------------------整理出type文件---------------------------------
type_kegg <- kegg_target_network %>%
             dplyr::select(node=from_node)%>%
             distinct()
type_kegg$type <- "kegg"


type_gene <- kegg_target_network %>%
             dplyr::select(node=to_node)%>%
             distinct()
type_gene$type <- "target"


#合并3个type
data_frames <- list()    # 初始化一个空的列表来存储数据框
#合并数据框
data_frames[[length(data_frames) + 1]] <- type_kegg 
data_frames[[length(data_frames) + 1]] <- type_gene 
type_combine <- do.call(rbind, data_frames)          #合并数据框

write.xlsx(type_combine, type_path, rowNames = FALSE)
#-------------------------------------------------------------------------------


