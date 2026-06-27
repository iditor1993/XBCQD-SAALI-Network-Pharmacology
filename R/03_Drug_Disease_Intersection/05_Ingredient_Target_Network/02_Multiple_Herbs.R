# setwd() removed - use RStudio Project or set working directory to project root                                #设置工作目录

# 清除当前环境中的所有变量
rm(list = ls())



library(dplyr)
library(openxlsx)



#------------------------创建文件夹、定义输入\输出文件路径----------------------
Filter_fold="results/03_Intersection_Targets/05_Ingredient_Target_Network"

#创建文件夹路径Filter_fold
if (!dir.exists(Filter_fold)) {
  dir.create(Filter_fold, recursive = TRUE)
  cat("文件夹已创建:", Filter_fold, "\n")
} else {
  cat("文件夹已存在:", Filter_fold, "\n")
}




#输入文件路径
ingredient_target_path <-"results/01_Drug_Targets/02_Database_Targets_Summary.csv"    #输出文件路径
Intersection_target_path <- "results/03_Intersection_Targets/Intersection_Targets.csv"
activate_compound_path <- "results/Drug_Active_Ingredients.csv"


#输出文件路径
type_path <- "results/03_Intersection_Targets/05_Ingredient_Target_Network/type.xlsx"
network_path <- "results/03_Intersection_Targets/05_Ingredient_Target_Network/Ingredint_target_network.xlsx"
#-------------------------------------------------------------------------------






#读取输入文件
ingredient_target <- read.csv(ingredient_target_path,header=TRUE, check.names = FALSE)
Intersection_target <- read.csv(Intersection_target_path,header=TRUE, check.names = FALSE)
activate_compound <- read.csv(activate_compound_path,header=TRUE, check.names = FALSE, fileEncoding = "GBK",colClasses = "character")


#ingredient_target筛选出交集靶点部分
ingredient_target <-merge(ingredient_target,Intersection_target, by = "Gene_symbol")


#生成ingredient_target_network文件
ingredient_target_net <- ingredient_target %>% dplyr::select(Mol_ID, Gene_symbol)%>%distinct()
names(ingredient_target_net)[names(ingredient_target_net) == "Mol_ID"] <- "from_node"
names(ingredient_target_net)[names(ingredient_target_net) == "Gene_symbol"] <- "to_node"



#生成  中药_ingredient_network文件
medicine_ingredient_net <- merge(activate_compound, ingredient_target_net,by.x = "Mol_ID",by.y = "from_node")%>% 
                          dplyr::select("Herb_name", "Mol_ID")%>%distinct()
names(medicine_ingredient_net)[names(medicine_ingredient_net) == "Herb_name"] <- "from_node"
names(medicine_ingredient_net)[names(medicine_ingredient_net) == "Mol_ID"] <- "to_node"


#生成最终的network文件
data_frames <- list()    # 初始化一个空的列表来存储数据框
data_frames[[length(data_frames) + 1]] <- ingredient_target_net
data_frames[[length(data_frames) + 1]] <- medicine_ingredient_net
combined_net <- do.call(rbind, data_frames)%>% 
  distinct() 



#保存为文件
write.xlsx(combined_net,network_path,rowNames = FALSE)


#生成typek文件
#ingredient的type
drug_list <- unique(medicine_ingredient_net$from_node)
head(drug_list)

data_frames<- list()

for (drug_name in drug_list) {
  
  print(drug_name)
  
  
  drug_unique_commpound <- medicine_ingredient_net %>%
    group_by(to_node) %>%
    filter(n() == 1 & from_node == drug_name) %>%
    ungroup()
  
  names(drug_unique_commpound)[names(drug_unique_commpound) == "from_node"] <- "type"
  names(drug_unique_commpound)[names(drug_unique_commpound) == "to_node"] <- "node"
  
  drug_unique_commpound <- drug_unique_commpound %>%dplyr::select(node,type)
  
  data_frames[[length(data_frames) + 1]] <- drug_unique_commpound
  
  
}

ingredient_type <- do.call(rbind, data_frames)

#target的type
target_type <- list()
target_type$node <- ingredient_target$Gene_symbol
target_type$type <- "target"
target_type <- as.data.frame(target_type)%>%
  distinct()


#Herb_name的type
medicine_type <- list()
medicine_type$node <- medicine_ingredient_net$from_node
medicine_type$type <- "medicine"
medicine_type <- as.data.frame(medicine_type)%>%
  distinct()



#合并数据框
data_frames <- list()
#将多个个数据库来源的target合并为一个列表里
data_frames[[length(data_frames) + 1]] <- ingredient_type
data_frames[[length(data_frames) + 1]] <- target_type
data_frames[[length(data_frames) + 1]] <- medicine_type
combined_type <- do.call(rbind, data_frames)%>%
  distinct()

#保存为type文件
write.xlsx(combined_type,type_path,rowNames = FALSE)


