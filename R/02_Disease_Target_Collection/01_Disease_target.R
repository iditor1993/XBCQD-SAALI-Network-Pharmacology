# setwd() removed - use RStudio Project or set working directory to project root                                #设置工作目录

# 清除当前环境中的所有变量
rm(list = ls())

library(dplyr)
library(ggplot2)



#--------------------------------创建文件夹-------------------------------------
#定义01.Drug_target文件夹路径
Filter_fold="results/02_Disease_Targets"

#创建文件夹路径Filter_fold
if (!dir.exists(Filter_fold)) {
  dir.create(Filter_fold, recursive = TRUE)
  cat("文件夹已创建:", Filter_fold, "\n")
} else {
  cat("文件夹已存在:", Filter_fold, "\n")
}
#-------------------------------------------------------------------------------





#-------------------------------检查文件内容------------------------------------
# 获取文件夹中的所有文件名（包含完整路径）
folder_path <- "data/disease_targets"
file_list <- list.files(path = folder_path, full.names = TRUE)

#输出使用的数据库名称
database_list <- tools::file_path_sans_ext(basename(file_list))
cat("使用的数据库：\n",paste(database_list, collapse = "\n"))



#检查文件
for (database_file in file_list) {

  current_database_target <- read.csv(database_file)
  
  col_list <- c("Gene_symbol")
  
  # 检查col_list是否都存在于current_database_target中
  if(!all(col_list %in% colnames(current_database_target))) {
    cat("文件 ", database_file ," 读取有误，请检查文件格式，文件名，文件内容等\n")
    cat("务必确保文件格式是逗号分隔符的csv文件\n")
    cat("务必确保表格中包含Mol_ID、Gene_symbol这两列\n")
    stop("请安要求检查、修改文件")
  }
  
}
#-------------------------------------------------------------------------------







#-----------------------------读取、汇总文件------------------------------------
# 获取文件夹中的所有文件名（包含完整路径）
folder_path <- "data/disease_targets"
file_list <- list.files(path = folder_path, full.names = TRUE)
database_list <- tools::file_path_sans_ext(basename(file_list))


data_frames <- list()    # 初始化一个空的列表来存储数据框

#读取、汇总文件
for (database_name in database_list) {
  
  database_file <- paste0("疾病靶点预测/", database_name, ".csv")
  
  current_database_target <- read.csv(database_file)
  
  current_database_target <- current_database_target %>% 
                            dplyr::select(Gene_symbol)%>% 
                            distinct()
  
  current_database_target$Data_base <- database_name
  
  #将多个个数据库来源的target合并为一个列表里
  data_frames[[length(data_frames) + 1]] <- current_database_target
  
}


#汇总所有数据库
combined_df <- do.call(rbind, data_frames)
write.csv(combined_df,"results/02_Disease_Targets/02_Database_Targets_Summary.csv",row.names = FALSE)  #保留文件



#清除重复项
unique_target <-combined_df%>% 
                dplyr::select(Gene_symbol)%>% 
                distinct(Gene_symbol, .keep_all = TRUE)
write.csv(unique_target,"results/02_Disease_Targets/01_Unique_Disease_Targets.csv",row.names = FALSE)  #保留文件
#-------------------------------------------------------------------------------









#----------------------------数据库-靶点数量统计--------------------------------
Database_target <-combined_df%>% 
  dplyr::select(Data_base, Gene_symbol)%>% 
  distinct()

count_df <- Database_target  %>%
  count(Data_base, name = "Count")


ggplot(count_df, aes(x = Data_base, y = Count)) +
  geom_col(width = 0.5, fill = "#4682b4") + # 柱子填充颜色
  geom_text(aes(label = Count), # 添加文本图层
            vjust = -0.5, # 调整文本位置，使其位于柱子顶部
            color = "black",size = 2) + # 设置文本颜色，确保在深色背景上可见
  labs(x = "",y = "Count of targets") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA), # 添加外框线
        axis.text.x = element_text(hjust = 1, angle = 45, size = 6, colour = "black"),
        axis.text.y = element_text(hjust = 1, size = 6, colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(color = "black", size = 14),  #设置坐标轴标题格式
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 2, unit = "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("results/02_Disease_Targets/03_Database_Target_Stats.png", width = 8, height = 8, dpi = 1000, bg = "white")
ggsave("results/02_Disease_Targets/03_Database_Target_Stats.pdf", width = 8, height = 8, dpi = 1000, bg = "white")
#-------------------------------------------------------------------------------








