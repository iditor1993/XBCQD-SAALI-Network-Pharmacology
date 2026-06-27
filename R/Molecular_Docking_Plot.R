library(ggplot2)
library(dplyr)
library(tidyr)

# 1. 准备数据（数据顺序与之前一致，包含您新增的 Licochalcone B）
data <- data.frame(
  Molecule = c(
    rep("25R-Timosaponin-AIII", 3), rep("Aloe Emodin", 3),
    rep("Catechin", 3), rep("Cucurbitacin C", 3),
    rep("Glabridin", 3), rep("Glycyrol", 3),
    rep("Phaseol", 3), rep("SMI-3", 3),
    rep("Licochalcone B", 3)
  ),
  Target = rep(c("FASN", "IL1B", "MAPK14"), 9),
  BindingEnergy = c(
    -9.8, -7.7, -8.5,  -8.3, -6.7, -8.6,  
    -8.0, -6.9, -9.5,  -8.3, -6.8, -8.4,  
    -10.3, -8.0, -10.8, -9.3, -7.1, -8.6,  
    -11.5, -7.3, -8.6, -9.3, -6.2, -8.9,  
    -9.0, -7.1, -8.2
  )
)

# 2. 设置因子水平（保证行和列顺序与图片一致）
data$Molecule <- factor(data$Molecule, 
                         levels = c("25R-Timosaponin-AIII", "Aloe Emodin", 
                                    "Catechin", "Cucurbitacin C", 
                                    "Glabridin", "Glycyrol", 
                                    "Phaseol", "SMI-3", 
                                    "Licochalcone B"))
data$Target <- factor(data$Target, levels = c("FASN", "IL1B", "MAPK14"))

# 3. 绘图
p <- ggplot(data, aes(x = Target, y = Molecule, fill = BindingEnergy)) +
  geom_tile(color = "white", size = 0.5) +
  
  # 添加数值标签
  geom_text(aes(label = sprintf("%.1f", BindingEnergy)), 
            color = "black", family = "serif", size = 5, fontface = "bold") +

  # 【关键修改】颜色映射：负值越大（如-11.5）越红，负值越小（如-6.2）越蓝，中间为米色
  scale_fill_gradientn(
    # 颜色顺序：红 -> 米色 -> 蓝（对应数值从低到高）
    colors = c("#d73027", "#fddbc7", "#f7f7f7", "#7facd6", "#2166ac"),
    # 锚点映射：-11.5 为红色，-8.5 为米色，-6.2 为蓝色
    values = scales::rescale(c(-11.5, -8.5, -6.2)),
    name = "Binding Energy (kcal/mol)"
  ) +
  
  labs(x = NULL, y = NULL) +
  
  guides(fill = guide_colorbar(title.position = "right", 
                               label.position = "right",
                               frame.colour = "black",
                               ticks.colour = "black")) +
  
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                               size = 14, color = "black", 
                               family = "serif", face = "bold"),
    axis.text.y = element_text(size = 14, color = "black", 
                               family = "serif", face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(family = "serif", size = 16, angle = 90, hjust = 0.5),
    legend.text = element_text(family = "serif", size = 14),
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10)
  )

# 4. 预览
print(p)

# 5. 导出图片（确保图例完整）
ggsave(filename = "heatmap_gradient_reversed.png", plot = p, width = 10, height = 7, dpi = 300)