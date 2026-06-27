# setwd() removed - use RStudio Project or set working directory to project root                                #设置工作目录

# 清除当前环境中的所有变量
rm(list = ls())


library(VennDiagram)


#--------------------------------创建文件夹-------------------------------------
#定义01.Drug_target文件夹路径
Filter_fold="results/04_Core_Targets/02_Venn"

#创建文件夹路径Filter_fold
if (!dir.exists(Filter_fold)) {
  dir.create(Filter_fold, recursive = TRUE)
  cat("文件夹已创建:", Filter_fold, "\n")
} else {
  cat("文件夹已存在:", Filter_fold, "\n")
}
#-------------------------------------------------------------------------------






#-----------------------------读取、汇总文件------------------------------------
# 获取文件夹中的所有文件名（包含完整路径）
folder_path <- "results/04_Core_Targets/01_Core_Targets_Screening"
file_list <- list.files(path = folder_path, full.names = TRUE)
cytohubba_list <- tools::file_path_sans_ext(basename(file_list))


#定义空的向量,用于存储绘制韦恩图的Gene_symbol
gene_list <- list()


# 遍历core_methods_file的每一行
for (Method in cytohubba_list) {
  print(Method )

  Core_target_path <-paste0("results/04_Core_Targets/01.核心靶点筛选/",Method, ".csv")    #输出文件路径
  
  Core_target_list <- read.csv(Core_target_path,header=TRUE, check.names = FALSE, sep = "\t")
  
  gene_list[[Method]] <- Core_target_list$Gene_symbol 
}
#-------------------------------------------------------------------------------






#---------------------------------绘制韦恩--------------------------------------
#以下是绘制韦恩图
num_core_methods<-length(cytohubba_list)   #统计使用了多少种筛选核心靶点的方法
#韦恩图的输出路径
png_path="results/04_Core_Targets/02_Venn/Core_target_veen.png"

#定义需要的填充颜色
fill_color <- c('yellow','#c20df6','red','green','orange', 'blue')
fill_current <- fill_color[1:num_core_methods]  #截取当前差异组合数的颜色

#定义需要的标签字体颜色
label_color <- c('black', 'black','black','black','black')
label_color_current <- label_color[1:num_core_methods]  #截取当前差异组合数的颜色

#定义标签与圆圈之间的距离
dis <- c(0.05,0.05,0.05,0.05,0.05,0.05,0.05)
dis_current <- dis[1:num_core_methods]  #截取当前差异组合数的标签与圆圈之间的距离



data.list<-gene_list
venn.diagram(
  data.list,
  filename = png_path,#输出图片名
  cat.default.pos = "outer",#'text',标签在圆圈内,'outer'标签在圆圈外
  fill=fill_current,#圆圈填充颜色
  col = NA,#圆圈边框颜色
  cat.col = label_color_current,# "orange", "darkorchid4"),#标签字体颜色
  cat.cex = 1,#标签字体大小
  cat.dist=dis_current,#标签与圆圈之间的距离
  #cat.pos = c(180,180,180),#标签名称与圆圈之间的角度，0点到12点
  cat.fontface=1,#标签字体设置为斜体,1常规，2加粗，3斜体，4加粗+斜体
  label.col='black',#数字颜色
  cex=2,#数字大小
  fontface=0.5,#数字字体设置为常规
  Resolution=300,#输出图片分辨率
  imagetype="png",#输出图片类型
  scaled=F)#不根据比例显示大小
#-------------------------------------------------------------------------------






#---------------------------------提取表格--------------------------------------
#以下是提取各交集的gene列表
target_veen_path <-"results/04_Core_Targets/02_Venn/Core_target_veen.csv"

gene_vector <- unlist(gene_list)         #将gene_list列表合并成一个向量
unique_gene_vector<- unique(gene_vector)    #去除重复项
gene_list_frame <- data.frame(Gene_symbol = unique_gene_vector)   #转为数据框



# 遍历ingredient列的每一行
for (Method in cytohubba_list) {
  
  Core_target_path <-paste0("results/04_Core_Targets/01.核心靶点筛选/",Method, ".csv")    #输出文件路径
  
  Core_target_list <- read.csv(Core_target_path,header=TRUE, check.names = FALSE, sep = "\t")
  
  #找出靶点
  gene_list_frame[[as.character(Method)]] <- ""    #新建一个空白列  
  gene_list_frame[[as.character(Method)]]<- ifelse(gene_list_frame$Gene_symbol %in%Core_target_list$Gene_symbol , "TRUE", "FALSE")   #按照gene_id列进行匹配，判断是否存在于该集合中
  
}


#写入文件
write.csv(gene_list_frame,target_veen_path,row.names = FALSE)
#-------------------------------------------------------------------------------






#------------------------------提取交集核心靶点---------------------------------
#输出文件路径
final_target_path <- "results/04_Core_Targets/01_Core_Targets.csv"


i=1
# 遍历ingredient列的每一行
for (Method in cytohubba_list) {

  Core_target_path <-paste0("results/04_Core_Targets/01.核心靶点筛选/",Method, ".csv")    #输出文件路径
  
  Core_target_list <- read.csv(Core_target_path,header=TRUE, check.names = FALSE, sep = "\t")
  
  #通过表格合并相同值的方式筛选出交集靶点
  if (i==1) {
    Core_target_A <- Core_target_list
  }
  
  if (i>1) {
    
    Core_target_A <- merge(Core_target_A, Core_target_list, by ="Gene_symbol")
  }
  
  i=i+1
  
}


write.table(Core_target_A, final_target_path, row.names = FALSE, sep = "\t")
#-------------------------------------------------------------------------------



