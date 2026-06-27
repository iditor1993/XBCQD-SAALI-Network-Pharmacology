# setwd() removed - use RStudio Project or set working directory to project root                                #设置工作目录

# 清除当前环境中的所有变量
rm(list = ls())



library(tidyverse)
library(openxlsx)
library(dplyr)
library(ggsankeyfier) 
library(MetBrewer)
library(colorspace)
library(ggh4x)
library(ggfun)
library(ggnewscale)
library(grid)
library(gridExtra)
library(cowplot)



#------------------------创建文件夹、定义输入\输出文件路径----------------------
Filter_fold="results/04_Core_Targets/05_Sankey_Diagram"

#创建文件夹路径Filter_fold
if (!dir.exists(Filter_fold)) {
  dir.create(Filter_fold, recursive = TRUE)
  cat("文件夹已创建:", Filter_fold, "\n")
} else {
  cat("文件夹已存在:", Filter_fold, "\n")
}



#输入文件为KEGG富集分析的表格
KEGG_result_path <- "results/03_Intersection_Targets/02_GO_KEGG_Enrichment/02_KEGG/KEGG_enrich.csv"
Tatget_PPi_path <- "results/04_Core_Targets/03_Core_PPI/Target_PPi.csv"
core_target_path <- "results/04_Core_Targets/01_Core_Targets.csv"
drug_target_path <- "results/01_Drug_Targets/02_Database_Targets_Summary.csv"

#输出文件路径
sankey_png <- "results/04_Core_Targets/05_Sankey_Diagram/sankey.png"
sankey_pdf <- "results/04_Core_Targets/05_Sankey_Diagram/sankey.pdf"
kegg_png <- "results/04_Core_Targets/05_Sankey_Diagram/kegg_dot.png"
kegg_pdf <- "results/04_Core_Targets/05_Sankey_Diagram/kegg_dot.pdf"
plot_out_png <- "results/04_Core_Targets/05_Sankey_Diagram/sankey_kegg_merge.png"
plot_out_pdf <- "results/04_Core_Targets/05_Sankey_Diagram/sankey_kegg_merge.pdf"
#-------------------------------------------------------------------------------




#---------------------------------读取输入文件----------------------------------
Tatget_PPi <- read.csv(Tatget_PPi_path,header=TRUE, check.names = FALSE)
drug_target <- read.csv(drug_target_path, stringsAsFactors = FALSE)


ingredient_target <- drug_target %>%
                     dplyr::select(Mol_ID,Gene_symbol) %>%
                     distinct()

combined_df <- ingredient_target %>%
  dplyr::select(Mol_ID,Gene_symbol)
#-------------------------------------------------------------------------------




#-------------------------------整理作图数据------------------------------------
# 读取文件
kegg_enrich <- read.csv(KEGG_result_path, stringsAsFactors = FALSE)


# 定义一个函数来解析分数并转换成小数
parse_fraction <- function(fraction_str) {
  # 检查输入是否为分数形式（格式为numerator/denominator）
  if (grepl("/", fraction_str)) {
    # 使用strsplit分割字符串，得到分子和分母
    parts <- unlist(strsplit(fraction_str, "/"))
    # 将字符串转换为数值，并计算小数
    numerator <- as.numeric(parts[1])
    denominator <- as.numeric(parts[2])
    return(numerator / denominator)
  } else {
    # 如果不是分数格式，尝试直接转换为数值（或根据需要处理错误）
    return(as.numeric(fraction_str))
  }
}


#kegg的KEGG富集气泡图
# 计算-log10(pvalue)值并添加到数据框中
kegg_enrich$neg_log10_pvalue <- -log10(kegg_enrich$pvalue)

# 应用函数转换GeneRatio列
kegg_enrich$GeneRatio <- sapply(kegg_enrich$GeneRatio, parse_fraction)

#截取Descripton的前50个字符
kegg_enrich$Description <- substr(kegg_enrich$Description, 1, 50)       

#筛选出pvalue<0.05的通路
#kegg_enrich <- kegg_enrich[kegg_enrich$pvalue<0.05, ]

kegg_top30 <- kegg_enrich %>%
  dplyr::arrange(pvalue) %>%       #升序排序
  dplyr::slice(1:30)               #前30的通路


#通路对应的靶点
output_data_RNA <- kegg_top30 %>%  
  rowwise() %>%  
  mutate(Gene_symbol = strsplit(as.character(Gene_symbol), "/")) %>%  
  unnest(Gene_symbol) %>%  
  ungroup() %>%  
  dplyr::select(KEGGID, Gene_symbol)


#找出前30的核心靶点
target_node_count <- c(Tatget_PPi$from_symbol, Tatget_PPi$to_symbol)
node_counts <- table(target_node_count)   # 使用table()函数统计每个元素的数量
target_type <- as.data.frame(node_counts)   # 将表格转换为数据框
colnames(target_type) <- c("Gene_symbol", "count")  # 重命名列名
target_type <- target_type %>% arrange(desc(count))  %>% #按照count列进行降序排序
  dplyr::slice(1:30)           #找出Degree值前30的核心靶点


output_data_RNA  <- merge(output_data_RNA, target_type, by= "Gene_symbol") 
#-------------------------------------------------------------------------------






#----------------------------------桑基图---------------------------------------
#用于桑基图的数据表
sankey_df <- merge(combined_df, output_data_RNA, by = "Gene_symbol") %>%
  dplyr::select(Mol_ID, Gene_symbol, KEGGID)%>%
  dplyr::distinct()



Mol_ID_counts <- as.data.frame(table(sankey_df$Mol_ID))   
colnames(Mol_ID_counts ) <- c("Mol_ID", "count")  # 重命名列名
Mol_ID_counts <- Mol_ID_counts %>% arrange(desc(count))  %>% #按照count列进行降序排序
  dplyr::slice(1:30)           #每个分组选出前10

sankey_df <- merge(sankey_df, 
                   Mol_ID_counts %>% dplyr::select(Mol_ID),
                   by="Mol_ID")


df1 <- sankey_df %>% dplyr::select(1,2) %>% group_by(Mol_ID,Gene_symbol) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Mol_ID", "Gene_symbol"),
                      values_from = "n")
df2 <- sankey_df %>% dplyr::select(2,3) %>% group_by(Gene_symbol,KEGGID) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Gene_symbol", "KEGGID"),
                      values_from = "n")


# 确定需要的颜色数量
num_colors_needed <- length(unique(c(df1$node, df2$node)))

colors <- c("#FFFF66", "#FFCCFF", "#FFCC99", "#FFCC33", "#FF99FF", "#FF9999", "#FF9933", "#FF66FF", "#FF6666", "#FF33CC", "#FF3333", "#CCFFFF", "#CCFF99", "#CCFF33", "#CCCCFF", "#CCCC99", "#CCCC33", "#CC99FF", "#CC9999", "#CC9933", "#CC66FF", "#CC6699", "#99FFFF", "#99FF99", "#99FF33", "#99CCFF", "#99CC99", "#9999FF", "#999999", "#999966", "#999900", "#9966CC", "#996633", "#9933CC", "#9900FF", "#99FFFF", "#99FFCC", "#99FF66", "#99CC33", "#999999", "#9966CC", "#996633", "#00FF99", "#00CCCC", "#00CC66", "#FFFF66", "#FFCCFF", "#FFCC99", "#FFCC33", "#FF99FF", "#FF9999", "#FF9933", "#FF66FF", "#FF6666", "#FF33CC", "#FF3333", "#CCFFFF", "#CCFF99", "#CCFF33", "#CCCCFF", "#CCCC99","#FF9999", "#FF9933", "#FF66FF", "#FF6666", "#FF33CC", "#FF3333", "#CCFFFF", "#CCFF99", "#CCFF33", "#CCCCFF", "#CCCC99", "#CCCC33", "#CC99FF", "#CC9999", "#CC9933", "#CC66FF", "#CC6699", "#99FFFF", "#99FF99", "#99FF33", "#99CCFF", "#99CC99", "#9999FF", "#999999", "#999966", "#999900", "#9966CC", "#996633", "#9933CC", "#9900FF", "#99FFFF", "#99FFCC", "#99FF66", "#99CC33", "#999999", "#9966CC")

p1 <- ggplot(data=df1,aes(x = stage,y =n,group = node,
                          edge_id = edge_id,connector = connector))+
  # 绘制第 1，2 层级
  geom_sankeyedge(aes(fill = node),
                  position = position_sankey(order ="ascending",v_space="auto",
                                             width = 0.2))+
  geom_sankeynode(aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.2))+
  geom_text(data=df1 %>% filter(connector=="from"),
            aes(label = node),stat = "sankeynode",
            position = position_sankey(v_space ="auto",order="ascending",nudge_x=0),
            hjust=0.5,size=2.5,fontface="bold",color="black")+
  # 绘制第 2，3 层级
  geom_sankeyedge(data=df2,aes(fill = node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.2))+
  geom_sankeynode(data=df2,aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.2))+
  
  geom_text(data=df1 %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",#angle=-90,
            position = position_sankey(v_space ="auto",order="ascending", nudge_x=0),
            hjust=0.5,size=3,vjust=0.5,color="black",fontface="bold")+
  geom_text(data=df2 %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",#angle=-90,
            position = position_sankey(v_space ="auto",order="ascending",nudge_x=0),
            hjust=0.5,size=3,color="black",fontface="bold")+
  
  coord_cartesian(clip="off")+
  
  scale_fill_manual(values = colors[1:num_colors_needed]) +
  scale_color_manual(values = colors[1:num_colors_needed]) +
  theme_void()+
  theme(plot.margin = margin(1,0,0,0,unit = "cm"),
        axis.text.x=element_text(color="red",face="bold",size=10),
        legend.position="none")

ggsave(sankey_png, width = 7, height = 8, dpi = 300,bg = "white")
ggsave(sankey_pdf, width = 7, height = 8, dpi = 300,bg = "white")
#-------------------------------------------------------------------------------






#----------------------------------气泡图---------------------------------------
#以下是绘制气泡图
g <- ggplot_build(p1)  #获取桑基图里kegg标签的y轴刻度值
node_data <- g$data[[which(sapply(g$data, function(x) "y" %in% names(x)))[7]]]
node_data <- node_data %>% dplyr::select(label,y)

#将刻度值添加到富集分析表格中
enrich_df <- merge(kegg_top30,node_data,by.x = "KEGGID",by.y = "label")


# 定义一个函数来解析分数并转换成小数
parse_fraction <- function(fraction_str) {
  # 检查输入是否为分数形式（格式为numerator/denominator）
  if (grepl("/", fraction_str)) {
    # 使用strsplit分割字符串，得到分子和分母
    parts <- unlist(strsplit(fraction_str, "/"))
    # 将字符串转换为数值，并计算小数
    numerator <- as.numeric(parts[1])
    denominator <- as.numeric(parts[2])
    return(numerator / denominator)
  } else {
    # 如果不是分数格式，尝试直接转换为数值（或根据需要处理错误）
    return(as.numeric(fraction_str))
  }
}

# 计算-log10(pvalue)值并添加到数据框中
enrich_df$neg_log10_pvalue <- -log10(enrich_df$pvalue)

# 应用函数转换GeneRatio列
enrich_df$GeneRatio <- sapply(enrich_df$GeneRatio, parse_fraction)


plot_df <- enrich_df %>%
  dplyr::arrange(desc(y)) 



#以下是提取出Pvalue的最大值和最小值，用于图例的范围设置
MAX <- plot_df %>%         
  pull(neg_log10_pvalue) %>%            # 提取neg_log10_pvalue列
  max()                                 # 计算最大值

MIN <- plot_df %>%
  pull(neg_log10_pvalue) %>%            # 提取neg_log10_pvalue列
  min()



breaks <- seq(MIN, MAX, length.out = 4) # 这将生成一个包含5个元素的向量，用于4个分割点



p2<-plot_df %>%
  ggplot() + 
  geom_point(data = plot_df,
             aes(x = GeneRatio, y = y, fill = neg_log10_pvalue, size = Count), shape = 21, colour = "grey90") + 
  geom_segment(data = plot_df,
               aes(x = GeneRatio, y = y, yend = y, colour = neg_log10_pvalue),
               xend = -0.01) +
  scale_fill_gradient(low = "yellow", high = "red", name = "-log10(pvalue)",
                      labels = function(x) sprintf("%.1f", x), breaks = breaks) +  # 确保breaks已被定义
  scale_colour_gradient(low = "yellow", high = "red", name = "-log10(pvalue)",labels = function(x) sprintf("%.1f", x), breaks = breaks,guide = "none") +
  
  ggtitle(label = " ") +
  labs(x = "GeneRatio", y=element_blank()) +
  scale_size(range = c(3, 7), guide = guide_legend(override.aes = list(colour = "black"))) +
  theme_bw() + 
  theme(
    legend.background = element_roundrect(color = "#969696"),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.8, "lines"),
    legend.text = element_text(size = 8),
    panel.border = element_rect(linewidth = 0.5, color = "#000000"),
    axis.text = element_text(color = "#000000", size = 8),
    axis.ticks.length = unit(3, "pt"),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.title = element_text(color = "red", face="bold",size=12),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) 


ggsave(kegg_png,width = 4, height = 8, dpi = 300,bg = "white")
ggsave(kegg_pdf,width = 4, height = 8, dpi = 300,bg = "white")
#-------------------------------------------------------------------------------



