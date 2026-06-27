# setwd() removed - use RStudio Project or set working directory to project root                                #设置工作目录

# 清除当前环境中的所有变量
rm(list = ls())


library(dplyr)
library(openxlsx)






#------------------------创建文件夹、定义输入\输出文件路径----------------------
Filter_fold="results/03_Intersection_Targets/04_Drug_Disease_Target_Network"

#创建文件夹路径Filter_fold
if (!dir.exists(Filter_fold)) {
  dir.create(Filter_fold, recursive = TRUE)
  cat("文件夹已创建:", Filter_fold, "\n")
} else {
  cat("文件夹已存在:", Filter_fold, "\n")
}


#输入文件路径
Drug_Disease_name_path <- "results/Drug_Disease_Names.csv"
Intersection_target_path <- "results/03_Intersection_Targets/Intersection_Targets.csv"


#输出文件路径
type_path <- "results/03_Intersection_Targets/04_Drug_Disease_Target_Network/type.xlsx"
veen_cytoscape_path <- "results/03_Intersection_Targets/04_Drug_Disease_Target_Network/Disease&Drug_Target_network.xlsx"
#-------------------------------------------------------------------------------





#--------------------------生成network文件--------------------------------------
Drug_Disease_name <- read.csv(Drug_Disease_name_path)  
Intersection_target <- read.csv(Intersection_target_path)


drug_name <- Drug_Disease_name$Drug_name
disease_name <- Drug_Disease_name$Disease_name


Intersection_target$to_node <-drug_name
durg_node <- Intersection_target


Intersection_target$to_node <-disease_name
disease_node <- Intersection_target



combined_node <- rbind(durg_node, disease_node) %>%
                               dplyr::select(from_node = Gene_symbol, to_node)


#写入文件
write.xlsx(combined_node,veen_cytoscape_path,rowNames = FALSE)
#-------------------------------------------------------------------------------











#-----------------------------生成type文件--------------------------------------
Intersection_target$type <- "target"

target_type <- Intersection_target %>%
               dplyr::select(node=Gene_symbol, type)


disease_drut_type <- data.frame(
  node = c(drug_name,disease_name),
  type = c("drug_name","disease_name")
)


combined_type <- rbind(target_type, disease_drut_type) %>%
  dplyr::select(node,type)


#写入文件
write.xlsx(combined_type, type_path, rowNames = FALSE)
#-------------------------------------------------------------------------------



