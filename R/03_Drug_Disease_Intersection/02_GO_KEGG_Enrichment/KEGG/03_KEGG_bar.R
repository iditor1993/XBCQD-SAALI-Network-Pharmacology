# setwd() removed - use RStudio Project or set working directory to project root                                #设置工作目录

# 清除当前环境中的所有变量
rm(list = ls())


library(gridExtra)
library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)
library(ggh4x)
library(ggfun)
library(ggnewscale)
library(grid)



#------------------------创建文件夹、定义输入\输出文件路径----------------------
#输入文件为kegg富集分析的表格
KEGG_result_path <- "results/03_Intersection_Targets/02_GO_KEGG_Enrichment/02_KEGG/KEGG_enrich.csv"

#保存KEGG富集柱形图的路径
KEGGbar_png <- "results/03_Intersection_Targets/02_GO_KEGG_Enrichment/02_KEGG/KEGG_bar.png"
KEGGbar_pdf <- "results/03_Intersection_Targets/02_GO_KEGG_Enrichment/02_KEGG/KEGG_bar.pdf"

#-------------------------------------------------------------------------------






#---------------------------读取文件、整理数据----------------------------------  
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
  # 读取文件
  kegg_enrich <- read.csv(KEGG_result_path, stringsAsFactors = FALSE)  
  
  # 计算-log10(pvalue)值并添加到数据框中
  kegg_enrich$neg_log10_pvalue <- -log10(kegg_enrich$pvalue)
  
  # 应用函数转换GeneRatio列
  kegg_enrich$GeneRatio <- sapply(kegg_enrich$GeneRatio, parse_fraction)
  
  #截取Descripton的前50个字符
  kegg_enrich$Description <- substr(kegg_enrich$Description, 1, 50)       
  
  #筛选出pvalue<0.05的通路
  kegg_enrich <- kegg_enrich[kegg_enrich$pvalue<0.05, ]
  
  kegg_top20 <- kegg_enrich %>%
    dplyr::arrange(pvalue) %>%       #升序排序
    dplyr::slice(1:30)           #每个分组选出前10
  
  plot_df <- kegg_top20 %>%
    dplyr::arrange(desc(GeneRatio)) %>%              #每一个Category都按照GeneRatio进行降序排序
    dplyr::mutate(Description = str_remove(Description, pattern = "\\(.*")) %>%     #删除Ddescription的（以及以后得内容
    dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered = T))
#-------------------------------------------------------------------------------
  
  
  
  
  
#-----------------------------------绘图----------------------------------------  
  
  #以下是提取出GeneRatio的最大值和最小值，用于图例的范围设置
  MAX <- plot_df %>%         
    pull(GeneRatio) %>%                      # 提取GeneRatio列
    max()                                 # 计算最大值
  
  MIN <- plot_df %>%
    pull(GeneRatio) %>%                      # 提取GeneRatio列
    min()
  
  
  
  breaks <- seq(MIN, MAX, length.out = 4) # 这将生成一个包含5个元素的向量，用于4个分割点
  
  
  
  p<-plot_df %>%
    ggplot() + 
    geom_col(data = plot_df,
             aes(x = neg_log10_pvalue, y = Description, fill = GeneRatio)) + 
    scale_fill_gradient(low = "yellow", high = "red", name = "GeneRatio",labels = function(x) sprintf("%.2e", x),breaks = breaks) + 
    ggtitle(label = "  ") + 
    labs(x = "-log10(pvalue)", y = "Description") + 
    scale_size(range = c(3, 7),guide = guide_legend(override.aes = list(colour = "#000000",fill= "#000000"))) +  #设置count图例的格式
    theme_bw() + 
    
    theme(
      legend.background = element_roundrect(color = "#969696"),    #设置图例外框线
      legend.title = element_text(size = 8),
      legend.key.size = unit(0.8, "lines"),   #设置图例图形大小
      legend.text = element_text(size = 8),   #设置图例文本大小
      panel.border = element_rect(size = 0.5, color = "#000000"),
      axis.text = element_text(color = "#000000", size = 12),
      axis.ticks.length = unit(3, "pt"),         #设置刻度线长度
      axis.ticks = element_line(color = "black", linewidth = 0.5),   #设置刻度线
      axis.title = element_text(color = "#000000", size = 15),
      plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
    )
  
  
  ggsave(KEGGbar_png,width = 10, height = 7, dpi = 300,bg = "white")   #在这里调整图片宽度、高度
  ggsave(KEGGbar_pdf,width = 10, height = 7, dpi = 300,bg = "white")
  #-------------------------------------------------------------------------------  

    
  
  
  