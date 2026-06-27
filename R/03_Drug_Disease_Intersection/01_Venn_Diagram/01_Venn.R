# setwd() removed - use RStudio Project or set working directory to project root                                #设置工作目录

# 清除当前环境中的所有变量
rm(list = ls())

library(VennDiagram)




#--------------------------------创建文件夹-------------------------------------
Filter_fold="results/03_Intersection_Targets/01_Venn_Diagram"

#创建文件夹路径Filter_fold
if (!dir.exists(Filter_fold)) {
  dir.create(Filter_fold, recursive = TRUE)
  cat("文件夹已创建:", Filter_fold, "\n")
} else {
  cat("文件夹已存在:", Filter_fold, "\n")
}
#-------------------------------------------------------------------------------





#----------------------------------读取文件-------------------------------------
#读取靶点文件
Drug_target <- read.csv("results/01_Drug_Targets/01_Unique_Drug_Targets.csv", check.names = FALSE)  
Disease_target <- read.csv("results/02_Disease_Targets/01_Unique_Disease_Targets.csv", check.names = FALSE)  
Drug_Disease_name <- read.csv("results/Drug_Disease_Names.csv", check.names = FALSE)   
#-------------------------------------------------------------------------------




#------------------------------整理作图数据-------------------------------------
#定义空的向量
gene_list <- list()
group_list <- c()


#创建输入数据的列表
gene_list <-list(Drug_target$Gene_symbol,
                 Disease_target$Gene_symbol
)

#每个集合的名称
group_list <- c(Drug_Disease_name[1, "Drug_name"], Drug_Disease_name[1, "Disease_name"])
#-------------------------------------------------------------------------------





#-------------------------------绘制韦恩图--------------------------------------
#韦恩图的输出路径
png_path="results/03_Intersection_Targets/01_Venn_Diagram/target_veen.png"
pdf_path="results/03_Intersection_Targets/01_Venn_Diagram/target_veen.pdf"



data.list<-gene_list
venn.diagram(
  data.list,
  category.names = group_list,
  filename = png_path,          #输出图片名
  cat.default.pos = "outer",    #'text',标签在圆圈内,'outer'标签在圆圈外
  fill = c('yellow','#c20df6'), #圆圈填充颜色
  col = NA,                     #圆圈边框颜色
  cat.col = c('black', 'black'),#标签字体颜色
  cat.cex = 1.5,                #标签字体大小
  cat.dist=c(0.05,0.05),        #标签与圆圈之间的距离
  cat.pos = c(180,180),         #标签名称与圆圈之间的角度，0点到12点
  cat.fontface=1,               #标签字体设置为斜体,1常规，2加粗，3斜体，4加粗+斜体
  label.col='black',            #数字颜色
  cex=1.5,                      #数字大小
  fontface=1,                   #数字字体设置为常规
  Resolution=1000,               #输出图片分辨率
  imagetype="png",              #输出图片类型
  scaled=F)                     #不根据比例显示大小
#-------------------------------------------------------------------------------







#------------------------------提取交集表格-------------------------------------
#以下是提取各交集的gene列表
target_veen_path <-"results/03_Intersection_Targets/01_Venn_Diagram/target_veen.csv"

gene_vector <- unlist(gene_list)         #将gene_list列表合并成一个向量
unique_gene_vector<- unique(gene_vector)    #去除重复项
gene_list_frame <- data.frame(Gene_symbol = unique_gene_vector)   #转为数据框



#找出Drug_target的靶点
gene_list_frame[[group_list[1]]] <- ""      #新建一个空白列，列名为药物名称    
gene_list_frame[[group_list[1]]] <- ifelse(gene_list_frame$Gene_symbol %in% Drug_target$Gene_symbol , "TRUE", "FALSE")   #按照gene_id列进行匹配，判断是否存在于该集合中


#找出Disease_target的靶点
gene_list_frame[[group_list[2]]] <- ""      #新建一个空白列，列名为药物名称    
gene_list_frame[[group_list[2]]] <- ifelse(gene_list_frame$Gene_symbol %in% Disease_target$Gene_symbol , "TRUE", "FALSE")   #按照gene_id列进行匹配，判断是否存在于该集合中


#写入文件
write.csv(gene_list_frame,target_veen_path,row.names = FALSE)  
#-------------------------------------------------------------------------------





#——-----------------------------最终交集基因-----------------------------------
#输出文件路径
Intersection_target_path <- "results/03_Intersection_Targets/Intersection_Targets.csv"


#通过表格合并的方式获取交集基因列表
Intersection_target <- merge(Drug_target, Disease_target, by ="Gene_symbol")%>% 
  dplyr::select(Gene_symbol)


write.table(Intersection_target, Intersection_target_path, row.names = FALSE,sep = "\t")  
#-------------------------------------------------------------------------------


