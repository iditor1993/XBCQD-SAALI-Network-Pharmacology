# setwd() removed - use RStudio Project or set working directory to project root                                #设置工作目录


# 清除当前环境中的所有变量
rm(list = ls())



library(clusterProfiler)
library(org.Hs.eg.db)





#------------------------创建文件夹、定义输入\输出文件路径----------------------
Filter_fold="results/03_Intersection_Targets/02_GO_KEGG_Enrichment/01_GO"

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
GO_result_path <- "results/03_Intersection_Targets/02_GO_KEGG_Enrichment/01_GO/GO_enrich.csv"  
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
 





#-----------------------------------GO富集-------------------------------------- 
  go_result <- enrichGO(gene         = data_ENTREZID$ENTREZID,
                        OrgDb        = org.Hs.eg.db,
                        keyType      = 'ENTREZID',
                        pAdjustMethod = 'BH', # 使用Benjamini-Hochberg方法进行p值校正
                        ont = "all",
                        universe     = NULL,
                        pvalueCutoff = 1, # 设定p值阈值
                       qvalueCutoff = 1, # 设定q值阈值
  )
  go_result_CSV<- go_result@result  #将GO富集的结果转换为表格
  names(go_result_CSV)[names(go_result_CSV) == "ID"] <- "GOID"     #修改表头
  names(go_result_CSV)[names(go_result_CSV) == "p.adjust"] <- "padj"     #修改表头
  
  # 将富集分析的结果保存为新的CSV文件
  write.csv(go_result_CSV, GO_result_path, row.names = FALSE)  
#-------------------------------------------------------------------------------  
  
  
  
  
  
  
#----------------------------添加Gene_symbol列----------------------------------
  # 将go_enrich中的geneID列根据"/"分割成多个geneID，并创建一个新的数据框来存储匹配后的ENSEMBL ID
  go_result_CSV$geneID_split <- strsplit(go_result_CSV$geneID, "/")
  go_result_CSV_final <- data.frame(
    ONTOLOGY = go_result_CSV$ONTOLOGY,
    GOID = go_result_CSV$GOID,
    Description = go_result_CSV$Description,
    GeneRatio = go_result_CSV$GeneRatio,
    BgRatio = go_result_CSV$BgRatio,
    pvalue= go_result_CSV$pvalue,
    padj = go_result_CSV$padj,
    qvalue = go_result_CSV$qvalue,
    geneID = go_result_CSV$geneID,
    Gene_symbol = sapply(go_result_CSV$geneID_split, function(x) paste(data_ENTREZID$SYMBOL[match(x, data_ENTREZID$ENTREZID)], collapse = "/")),
    Count = go_result_CSV$Count
  )
  
  # 将进行ID转换后的数据表保存到原来的文件路径中
  write.csv(go_result_CSV_final, GO_result_path, row.names = FALSE)  
#-------------------------------------------------------------------------------  
  
  