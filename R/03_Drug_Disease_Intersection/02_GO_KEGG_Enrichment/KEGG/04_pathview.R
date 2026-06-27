# setwd() removed - use RStudio Project or set working directory to project root                                #设置工作目录


# 加载必要的库
library(clusterProfiler)
library(pathview)



#输入文件路径
deg_KEGG_all_path <- "results/03_Intersection_Targets/02_GO_KEGG_Enrichment/02_KEGG/KEGG_enrich.csv"

# 读取输入文件
data_KEGG_all <- read.csv(deg_KEGG_all_path, stringsAsFactors = FALSE,head=TRUE)






#------------------------创建文件夹、定义输入\输出文件路径----------------------
Filter_fold="results/03_Intersection_Targets/02.GO_KEGG富集分析/02.KEGG/kegg_plot/KEGG_TEMP"

#创建文件夹路径Filter_fold
if (!dir.exists(Filter_fold)) {
  dir.create(Filter_fold, recursive = TRUE)
  cat("文件夹已创建:", Filter_fold, "\n")
} else {
  cat("文件夹已存在:", Filter_fold, "\n")
}

#-------------------------------------------------------------------------------



#定义临时工作路径和kegg图缓存路径
work_path=paste0("results/03_Intersection_Targets/02_GO_KEGG_Enrichment/02_KEGG/kegg_plot")
kegg.dir_path="KEGG_TEMP"
setwd(work_path) 

# 初始化一个空向量来存储出现错误的KEGG ID  
kegg_all_error <- c()  




# 遍历每个通路
for (i in 1:nrow(data_KEGG_all)) {
  # 提取通路 ID 和基因列表
  kegg_id <- data_KEGG_all[i, "KEGGID"]
  kegg_gene_id <- data.frame(ENTREZID= unlist(strsplit(data_KEGG_all[i, "geneID"], "/")))
  
  kegg_gene_id$count <- 1
  
  geneList <- kegg_gene_id$count
  names(geneList)=kegg_gene_id$ENTREZID
  head(geneList)
  
  # 打开一个新的图形设备来保存图像
  png()
  
  # 使用 tryCatch 来尝试绘制通路图，并在发生错误时跳过
  tryCatch({
    # 使用 pathview 生成通路图并保存为 PDF 文件（注意：这里是保存为 PNG）
    pathview(gene.data = geneList,
             species = "hsa",
             pathway.id = kegg_id,
             kegg.dir =  kegg.dir_path, # 用于缓存下载的通路数据
             limit = list(gene=max(abs(geneList))),
             bins = list(gene = 10),  # 设置颜色分段
             both.dirs = list(gene = FALSE),
             same.layer = T)
    
    
  }, error = function(e) {
    # 打印错误信息并跳过当前循环
    cat(sprintf("Error plotting pathway %s: %s\n", kegg_id, e$message))
    # 将出错的kegg_id添加到kegg_error向量中
    kegg_all_error <<- c(kegg_all_error, kegg_id) 
    
  })
  dev.off()
}
#回到原始工作目录
# setwd() removed - use RStudio Project or set working directory to project root       




#如果存在通路不能正常添加颜色的则执行一下部分
temp_kegg_all_error <- "results/03_Intersection_Targets/02.GO_KEGG富集分析/02.KEGG/kegg_plot/KEGG_TEMP/kegg_all_error.csv"

if (length(kegg_all_error)!=0) {
  # 将向量转换为数据框，其中向量的每个元素都是一行  
  kegg_all_error_df <- data.frame(KEGGID = kegg_all_error) 
  
  # 写入临时文件  
  write.csv(kegg_all_error_df, temp_kegg_all_error, row.names = FALSE, quote = FALSE)
  
  # 读取临时CSV文件到数据框  
  temp_kegg_all_error <- read.csv(temp_kegg_all_error,header = TRUE) 
  
  # 遍历数据框的每一行  
  for (i in 1:nrow(temp_kegg_all_error)) {  
    # 提取当前行的某个列的值，这里以'KEGGID'列为例  
    current_kegg_id <- temp_kegg_all_error[i, "KEGGID"]  
    print(current_kegg_id)
    
    #拼接需要移动的文件路径
    kegg_path_in=paste0("results/03_Intersection_Targets/02.GO_KEGG富集分析/02.KEGG/kegg_plot/KEGG_TEMP/", current_kegg_id, ".png")
    kegg_path_out=paste0("results/03_Intersection_Targets/02.GO_KEGG富集分析/02.KEGG/kegg_plot/", current_kegg_id, ".pathview.png")
    
    #移动文件
    file.copy(from = kegg_path_in, to = kegg_path_out, overwrite = TRUE) 
    
    # 判断是否有未关闭的图形设备，# 如果有，则循环关闭所有图形设备
    open_devices <- dev.list()
    if (length(open_devices) > 0) {
      while (length(dev.list()) > 0) {
        dev.off()
      }} 
    
    
  }
  
  
} 


#删除KEGG_TEMP文件夹
folder_path <- "results/03_Intersection_Targets/02.GO_KEGG富集分析/02.KEGG/kegg_plot/KEGG_TEMP"   # 指定文件夹路径
unlink(folder_path, recursive = TRUE, force = TRUE)     # 删除文件夹及其所有内容
