# setwd() removed - use RStudio Project or set working directory to project root                                #设置工作目录

# 清除当前环境中的所有变量
rm(list = ls())


library(dplyr)




#------------------------创建文件夹、定义输入\输出文件路径----------------------
Filter_fold="results/03_Intersection_Targets/03_PPI_Network"

#创建文件夹路径Filter_fold
if (!dir.exists(Filter_fold)) {
  dir.create(Filter_fold, recursive = TRUE)
  cat("文件夹已创建:", Filter_fold, "\n")
} else {
  cat("文件夹已存在:", Filter_fold, "\n")
}



#输入文件路径
Intersection_target_path <- "results/03_Intersection_Targets/Intersection_Targets.csv"
protein_to_symbol_path <- "data/STRING/9606.protein.info.v12.0.txt"
protein_links_path <- "data/STRING/9606.protein.links.v12.0.txt"


#输出文件路径
Tatget_PPi_path <- "results/03_Intersection_Targets/03_PPI_Network/Target_PPi.csv"
Edge_Node_count_path <- "results/03_Intersection_Targets/03_PPI_Network/Edge_Node_count.csv"
#-------------------------------------------------------------------------------






#--------------------------------读取输入文件-----------------------------------
Intersection_target <- read.csv(Intersection_target_path, stringsAsFactors = FALSE)  
protein_to_symbol <- read.csv(protein_to_symbol_path, stringsAsFactors = FALSE, sep ="\t")
protein_links <- read.csv(protein_links_path, stringsAsFactors = FALSE, sep =" ")
#-------------------------------------------------------------------------------




#--------------------------------进行PPI分析------------------------------------
#将string_protein_添加到Intersection_target
Intersection_target <- merge(Intersection_target,
                             protein_to_symbol[, c("string_protein_id", "preferred_name")],
                    by.x = "Gene_symbol", by.y ="preferred_name")


#将from_symbol添加到data_link
data_links <- merge(protein_links,
                    Intersection_target,
                    by.x = "protein1", by.y ="string_protein_id")
names(data_links)[names(data_links) == "Gene_symbol"] <- "from_symbol"

#将to_symbol添加到data_link
data_links <- merge(data_links,
                    Intersection_target,
                    by.x = "protein2", by.y ="string_protein_id")
names(data_links)[names(data_links) == "Gene_symbol"] <- "to_symbol"

#修改列名
names(data_links)[names(data_links) == "protein1"] <- "from_protein"
names(data_links)[names(data_links) == "protein2"] <- "to_protein"
data_links <- data_links %>% dplyr::select (from_protein, to_protein, from_symbol, to_symbol, combined_score) %>%
  rowwise() %>%    #从这一行到最后，都是去除重复行的过程
  mutate(
    sorted_from = pmin(from_symbol, to_symbol),
    sorted_to = pmax(from_symbol, to_symbol)
  ) %>%
  ungroup() %>%
  distinct(sorted_from, sorted_to, .keep_all = TRUE)%>%
  dplyr::select(-sorted_from, -sorted_to)



#------------------------------设置置信度、保存---------------------------------
#如果edge太多，可设置combined_score的阈值进行筛选
data_links <- data_links[data_links$combined_score >= 400, ]    #保留combined_score >= 400的行

#写入文件
write.csv(data_links, Tatget_PPi_path, row.names = FALSE)  
#-------------------------------------------------------------------------------







#统计node数量
all_nodes <- c(data_links$from_symbol, data_links$to_symbol)  # 合并两列的所有元素
unique_nodes <- unique(all_nodes)   # 找出唯一元素
num_unique_nodes<- length(unique_nodes)    # 统计唯一元素的数量

#统计edge数量
num_edges <- nrow(data_links)


#创建数据框保存node和edge的统计数
count_df <- data.frame(Node_count = c(num_unique_nodes),
                       Edge_count= c(num_edges)
)

#写入文件
write.csv(count_df, Edge_Node_count_path, row.names = FALSE)  











