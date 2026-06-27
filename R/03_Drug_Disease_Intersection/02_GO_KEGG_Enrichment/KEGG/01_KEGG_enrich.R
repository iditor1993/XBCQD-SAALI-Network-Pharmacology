# setwd() removed - use RStudio Project or set working directory to project root                                #设置工作目录

# 清除当前环境中的所有变量
rm(list = ls())


library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readxl)
library(AnnotationDbi)





#------------------------创建文件夹、定义输入\输出文件路径----------------------
Filter_fold="results/03_Intersection_Targets/02_GO_KEGG_Enrichment/02_KEGG/kegg_plot"

#创建文件夹路径Filter_fold
if (!dir.exists(Filter_fold)) {
  dir.create(Filter_fold, recursive = TRUE)
  cat("文件夹已创建:", Filter_fold, "\n")
} else {
  cat("文件夹已存在:", Filter_fold, "\n")
}


#输入文件为交集target文件
Intersection_target_path <- "results/03_Intersection_Targets/Intersection_Targets.csv"

#GO富集结果输出文件路径
KEGG_result_path <- "results/03_Intersection_Targets/02_GO_KEGG_Enrichment/02_KEGG/KEGG_enrich.csv"  
#-------------------------------------------------------------------------------





#----------------------------读取交集靶点、进行ID转换---------------------------
#读取输入文件
Intersection_target <- read.csv(Intersection_target_path, stringsAsFactors = FALSE,)  


#进行ID转换
data_ENTREZID<-bitr(
  Intersection_target$Gene_symbol ,
  fromType = 'SYMBOL',
  toType = 'ENTREZID',
  OrgDb = 'org.Hs.eg.db'
)
#-------------------------------------------------------------------------------






#-----------------------------------KEGG富集-------------------------------------- 
options(timeout = 600)
  # 进行KEGG富集分析
  KEGG_result <- clusterProfiler::enrichKEGG(gene = data_ENTREZID$ENTREZID,
                                             organism = 'hsa',
                                             pAdjustMethod = 'BH', # 使用Benjamini-Hochberg方法进行p值校正
                                             universe     = NULL,
                                             pvalueCutoff = 1, # 设定p值阈值
                                             qvalueCutoff = 1, # 设定q值阈值
  )
  
  
  
  KEGG_result<- KEGG_result@result
  names(KEGG_result)[names(KEGG_result) == "ID"] <- "KEGGID"
  names(KEGG_result)[names(KEGG_result) == "category"] <- "Category"
  names(KEGG_result)[names(KEGG_result) == "subcategory"] <- "Subcategory"
  names(KEGG_result)[names(KEGG_result) == "p.adjust"] <- "padj"
  # 将富集分析的结果保存为新的CSV文件
  write.csv(KEGG_result, KEGG_result_path, row.names = FALSE)  
#-------------------------------------------------------------------------------
  
  
  
  
  
  
#----------------------------添加Gene_symbol列----------------------------------  
  # 将kegg_enrich中的geneID列根据"/"分割成多个geneID，并创建一个新的数据框来存储匹配后的symbol
  KEGG_result$geneID_split <- strsplit(KEGG_result$geneID, "/")
  KEGG_result_final <- data.frame(
    Category = KEGG_result$Category,
    Subcategory = KEGG_result$Subcategory,
    KEGGID = KEGG_result$KEGGID,
    Description = KEGG_result$Description,
    GeneRatio = KEGG_result$GeneRatio,
    BgRatio = KEGG_result$BgRatio,
    pvalue= KEGG_result$pvalue,
    padj = KEGG_result$padj,
    qvalue = KEGG_result$qvalue,
    geneID = KEGG_result$geneID,
    Gene_symbol = sapply(KEGG_result$geneID_split, function(x) paste(data_ENTREZID$SYMBOL[match(x, data_ENTREZID$ENTREZID)], collapse = "/")),
    Count = KEGG_result$Count
  )
  
  # 将进行ID转换后的数据表保存到原来的文件路径中
  write.csv(KEGG_result_final, KEGG_result_path, row.names = FALSE)
#-------------------------------------------------------------------------------  
  
  
  
  
  
  