#在绘制气泡图和柱形图的脚本肿，ggh4x这个包最新版的去除掉了某些功能，
#导致绘图出bug，请安装0.3.0版本，方法如下：


#查看R包的版本,运行后出现版本号则代表已安装ggh4x包，未出现版本号则代表未安装ggh4x包
packageVersion("ggh4x") 



#如果版本号不是0.3.0
remove.packages("ggh4x")  #删除原有的R包



#指定安装0.3.0版本的ggh4x包
library(devtools)
devtools::install_version("ggh4x", version = "0.3.0")  



#再次检查R包的版本
packageVersion("ggh4x")  #查看R包的版本