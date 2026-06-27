# setwd() removed - use RStudio Project or set working directory to project root                                #设置工作目录

# 清除当前环境中的所有变量
rm(list = ls())


# 加载必要的库
library(gridExtra)
library(dplyr)
library(readxl)
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
#输入文件为go富集分析的表格
GO_result_path <- "results/03_Intersection_Targets/02_GO_KEGG_Enrichment/01_GO/GO_enrich.csv"  
  
#输出文件路径
gobar_png <- "results/03_Intersection_Targets/02_GO_KEGG_Enrichment/01_GO/GO_bar.png"
gobar_pdf <- "results/03_Intersection_Targets/02_GO_KEGG_Enrichment/01_GO/GO_bar.pdf"
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
  
  
  #go的GO富集气泡图
  # 读取文件
  go_enrich <- read.csv(GO_result_path, stringsAsFactors = FALSE)  
  
  # 计算-log10(pvalue)值并添加到数据框中
  go_enrich$neg_log10_value <- -log10(go_enrich$pvalue)
  
  # 应用函数转换GeneRatio列
  go_enrich$GeneRatio <- sapply(go_enrich$GeneRatio, parse_fraction)
  
  #截取Descripton的前50个字符
  go_enrich$Description <- substr(go_enrich$Description, 1, 120)       
  
  #筛选出pvalue<0.05的通路
  #go_enrich <- go_enrich[go_enrich$pvalue<0.05, ]
  
  GO_top10 <- go_enrich %>%
    dplyr::group_by(ONTOLOGY) %>%    #设置分组
    dplyr::arrange(pvalue) %>%       #升序排序
    dplyr::slice(1:10) %>%           #每个分组选出前10
    dplyr::ungroup()                 #取消分组
  
  
  plot_df <- GO_top10 %>%
    dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF")), ordered = T)) %>%    #将ONTOLOGY设置为因子
    dplyr::arrange(ONTOLOGY, desc(GeneRatio)) %>%              #每一个ONTOLOGY都按照GeneRatio进行降序排序
    dplyr::mutate(Description = str_remove(Description, pattern = "\\(.*")) %>%     #删除DDescription的（以及以后得内容
    dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered = T))
#-------------------------------------------------------------------------------
  
  
  
  
  
#-----------------------------------绘图----------------------------------------  
  #以下是提取出每个ONTOLOGY的GeneRatio的最大值和最小值，用于图例的范围设置
  BP_MAX <- plot_df %>%
    filter(ONTOLOGY == "BP") %>%          # 筛选出ONTOLOGY为"BP"的行
    pull(GeneRatio) %>%                      # 提取GeneRatio列
    max()                                 # 计算最大值
  
  BP_MIN <- plot_df %>%
    filter(ONTOLOGY == "BP") %>%          # 筛选出ONTOLOGY为"BP"的行
    pull(GeneRatio) %>%                      # 提取GeneRatio列
    min()
  
  CC_MAX <- plot_df %>%
    filter(ONTOLOGY == "CC") %>%          # 筛选出ONTOLOGY为"CC"的行
    pull(GeneRatio) %>%                      # 提取GeneRatio列
    max()                                 # 计算最大值
  
  CC_MIN <- plot_df %>%
    filter(ONTOLOGY == "CC") %>%          # 筛选出ONTOLOGY为"CC"的行
    pull(GeneRatio) %>%                      # 提取GeneRatio列
    min()
  
  MF_MAX <- plot_df %>%
    filter(ONTOLOGY == "MF") %>%          # 筛选出ONTOLOGY为"MF"的行
    pull(GeneRatio) %>%                      # 提取GeneRatio列
    max()                                 # 计算最大值
  
  MF_MIN <- plot_df %>%
    filter(ONTOLOGY == "MF") %>%          # 筛选出ONTOLOGY为"MF"的行
    pull(GeneRatio) %>%                      # 提取GeneRatio列
    min()
  
  
  
  BP_breaks <- seq(BP_MIN, BP_MAX, length.out = 4) # 这将生成一个包含4个元素的向量，用于3个分割点
  CC_breaks <- seq(CC_MIN, CC_MAX, length.out = 4) # 这将生成一个包含4个元素的向量
  MF_breaks <- seq(MF_MIN, MF_MAX, length.out = 4) # 这将生成一个包含4个元素的向量
  
  
  
  plot_df %>%
    ggplot() + 
    geom_col(data = plot_df %>% dplyr::filter(ONTOLOGY == "MF"),
             aes(x = neg_log10_value, y = interaction(Description,ONTOLOGY), fill = GeneRatio)) + 
    scale_fill_gradient(low = "#b2eecf", high = "#41ae76", name = "MF GeneRatio",labels = function(x) sprintf("%.2e", x),breaks = MF_breaks) +       #添加一个颜色渐变图例 
    ggnewscale::new_scale_fill() +       #新增一个图层
    geom_col(data = plot_df %>% dplyr::filter(ONTOLOGY == "CC"),
             aes(x = neg_log10_value,, y = interaction(Description,ONTOLOGY), fill = GeneRatio)) + 
    scale_fill_gradient(low = "#90dee7", high = "#0a8a99", name = "CC GeneRatio",labels = function(x) sprintf("%.2e", x),breaks = CC_breaks) + 
    ggnewscale::new_scale_fill() +
    geom_col(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = neg_log10_value, y = interaction(Description,ONTOLOGY), fill = GeneRatio)) + 
    scale_fill_gradient(low = "#f2c7b9", high = "#e74716", name = "BP GeneRatio",labels = function(x) sprintf("%.2e", x),breaks = BP_breaks) + 
    guides(y = "axis_nested") + 
    ggtitle(label = "  ") + 
    labs(x = "-log10(pvalue)", y = "Description") + 
    scale_x_continuous(limits = c(0, max(plot_df$neg_log10_value, na.rm = TRUE)*1.1), expand = c(0, 0))+  # 设置x轴刻度范围为0~最大值*1.1，x=0端不扩展
    scale_size(range = c(2, 5),guide = guide_legend(override.aes = list(colour = "#000000",fill= "#000000"))) +  #设置count图例的格式
    theme_bw() + 
    
    theme(
      ggh4x.axis.nestline.y = element_line(size = 3, color = c("#41ae76", "#0a8a99", "#e74716")),   #设置ONTOLOGY分类标签的颜色
      ggh4x.axis.nesttext.y = element_text(colour = c("#41ae76", "#0a8a99", "#e74716")),
      legend.background = element_roundrect(color = "#969696"),    #设置图例外框线
      legend.title = element_text(size = 8),
      legend.key.size = unit(0.8, "lines"),   #设置图例图形大小
      legend.text = element_text(size = 8),   #设置图例文本大小
      panel.border = element_rect(size = 0.5, color = "#000000"),
      axis.text = element_text(color = "#000000", size = 12),
      axis.text.y = element_text(color = rep(c("#41ae76", "#0a8a99", "#e74716"), each = 10)),      #设置y轴标签的颜色
      axis.ticks.length = unit(3, "pt"),         #设置刻度线长度
      axis.ticks = element_line(color = "black", linewidth = 0.5),   #设置刻度线
      axis.title = element_text(color = "#000000", size = 15),
      plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
    )+
    scale_y_discrete(labels = function(y) str_wrap(y, width = 70))  # 设置y轴标签换行
  

  
  ggsave(gobar_png,width = 10, height = 7, dpi = 300,bg = "white")    #在这里调整图片宽度、高度
  ggsave(gobar_pdf,width = 10, height = 7, dpi = 300,bg = "white")
#-------------------------------------------------------------------------------  
  
  
  
  
  