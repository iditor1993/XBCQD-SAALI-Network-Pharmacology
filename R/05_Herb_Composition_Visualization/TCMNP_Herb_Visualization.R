library(devtools)
#安装TCMNP包
# 安装TCMNP
devtools::install_github("tcmlab/TCMNP",upgrade =FALSE,dependencies =TRUE)

#加载需要的R包
library(TCMNP)
library(tidyverse)
library(ggraph)
library(clusterProfiler,quietly=TRUE)
library(org.Hs.eg.db,quietly=TRUE)
library(DOSE,quietly=TRUE)
library(readxl)

#基础处方组成可视化
xbcqd_compostion <- data.frame(
  herb = c("shengshigao","shengdahuang","kuxingren","gualou"),
  weight = c(25,20,10,10)
)

tcm_comp(xbcqd_compostion)

#处方+药性+归经综合可视化
herb = c("shi gao","da huang","ku xing ren","gua lou")
xbcqd_compostion <- data.frame(
  herb = herb,
  weight = c(25,20,10,10),
  property=herb_pm[match(herb,herb_pm$Herb_name_pinyin),]$Property,
  flavor=herb_pm[match(herb,herb_pm$Herb_name_pinyin),]$Flavor,
  meridian=herb_pm[match(herb,herb_pm$Herb_name_pinyin),]$Meridian
)
tcm_compound(xbcqd_compostion)

#中药→成分→靶点 预测
#拼音名称检索
herbs2 <- c("shi gao","da huang","ku xing ren","gua lou")
xbcqd <-  herb_target(herbs2,type = "Herb_name_pin_yin")
head(xbcqd)

write_excel_csv(xbcqd, "data/herb_target_prediction/xbcqd_herb_target.csv")

