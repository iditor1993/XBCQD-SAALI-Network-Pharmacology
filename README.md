[README.md](https://github.com/user-attachments/files/29419270/README.md)
# XBCQD-SAALI-Network-Pharmacology# XBCQD-SAALI Network Pharmacology

**Network Pharmacology Analysis of Xuanbai Chengqi Decoction (XBCQD) for Sepsis-Associated Acute Lung Injury (SAALI)**

---

## 项目简介 | Project Overview

本项目包含宣白承气汤（XBCQD）治疗脓毒症相关急性肺损伤（SAALI）的网络药理学分析完整 R 代码与运行结果。

This repository contains the complete R scripts and analysis results for the network pharmacology study of Xuanbai Chengqi Decoction (宣白承气汤) in the treatment of Sepsis-Associated Acute Lung Injury (SAALI).

---

## 分析流程 | Analysis Pipeline

| 步骤 | 内容 | 目录 |
|:---|:---|:---|
| 1 | 药物靶点收集 | `R/01_Drug_Target_Collection/` |
| 2 | 疾病靶点收集 | `R/02_Disease_Target_Collection/` |
| 3 | 药物-疾病交集靶点分析 | `R/03_Drug_Disease_Intersection/` |
| 4 | 核心靶点筛选与分析 | `R/04_Core_Targets/` |
| 5 | 中药复方可视化 | `R/05_Herb_Composition_Visualization/` |
| 6 | 分子对接图 | `R/Molecular_Docking_Plot.R` |

### 详细分析内容

- **01. 药物靶点汇总**: 整合多个数据库（TCMSP、SwissTargetPrediction 等）的药物活性成分靶点
- **02. 疾病靶点汇总**: 整合 GeneCards、OMIM、DisGeNET 等数据库的 SAALI 相关靶点
- **03. 交集靶点分析**:
  - Venn 图（药物 vs 疾病靶点交集）
  - GO 富集分析（柱状图 + 气泡图）
  - KEGG 通路富集分析（柱状图 + 气泡图 + Pathview）
  - PPI 蛋白互作网络（基于 STRING 数据库）
  - 药物-疾病-靶点网络图
  - 成分-靶点网络图
- **04. 核心靶点分析**:
  - 拓扑网络分析（Degree, Betweenness, Closeness）筛选核心靶点
  - 核心靶点 Venn 验证
  - 核心靶点 PPI 网络
  - 核心靶点-通路网络图
  - 成分-靶点-通路桑基图（Sankey Diagram）
- **05. 中药复方可视化**: 基于 TCMNP 包的中药成分、药性、归经可视化
- **分子对接**: 核心成分与核心靶点的分子对接结果可视化

---

## 项目结构 | Repository Structure

```
XBCQD-SAALI-Network-Pharmacology/
├── R/                          # R analysis scripts
│   ├── 01_Drug_Target_Collection/
│   ├── 02_Disease_Target_Collection/
│   ├── 03_Drug_Disease_Intersection/
│   │   ├── 01_Venn_Diagram/
│   │   ├── 02_GO_KEGG_Enrichment/
│   │   ├── 03_PPI_Network/
│   │   ├── 04_Drug_Disease_Target_Network/
│   │   └── 05_Ingredient_Target_Network/
│   ├── 04_Core_Targets/
│   │   ├── 01_Venn/
│   │   ├── 02_Core_PPI/
│   │   ├── 03_Core_Target_Pathway_Network/
│   │   └── 04_Sankey_Diagram/
│   ├── 05_Herb_Composition_Visualization/
│   └── Molecular_Docking_Plot.R
├── data/                       # Raw data and input files
│   ├── drug_targets/           # Drug target prediction data (TCMSP, etc.)
│   ├── disease_targets/      # Disease target data (GeneCards, etc.)
│   ├── STRING/                 # STRING database files
│   │   ├── 9606.protein.info.v12.0.txt
│   │   └── README.md         # Download link for protein.links file
│   ├── TCMSP/                  # TCMSP database files
│   ├── raw/                    # Other raw data (UniProt, GEO, etc.)
│   └── herb_target_prediction/ # TCMNP prediction results
├── results/                    # Analysis outputs
│   ├── 01_Drug_Targets/
│   ├── 02_Disease_Targets/
│   ├── 03_Intersection_Targets/
│   └── 04_Core_Targets/
├── docs/                       # Documentation and references
│   ├── notes/                  # Analysis notes
│   ├── tutorials/              # R troubleshooting tutorials
│   └── references/             # Reference papers
├── examples/                   # Example data for practice
├── .gitignore
├── LICENSE
└── README.md
```

---

## 依赖环境 | Dependencies

### R 包

```r
# Core packages
dplyr
ggplot2
VennDiagram
clusterProfiler
org.Hs.eg.db
DOSE
igraph
ggraph

# Network & visualization
pathview
STRINGdb

# TCM analysis
devtools::install_github("tcmlab/TCMNP")
```

### 外部数据库

| 数据库 | 版本 | 说明 |
|:---|:---|:---|
| STRING | v12.0 | 蛋白互作网络 |
| TCMSP | - | 中药系统药理学数据库 |
| KEGG | - | 通路富集分析 |
| GeneCards | - | 疾病靶点 |

---

## 使用说明 | Usage

### 1. 克隆仓库

```bash
git clone https://github.com/YOUR_USERNAME/XBCQD-SAALI-Network-Pharmacology.git
cd XBCQD-SAALI-Network-Pharmacology
```

### 2. 下载 STRING 数据库大文件

由于 `9606.protein.links.v12.0.txt` (602MB) 超过 GitHub 单文件限制，请从 STRING 官网下载：

- 下载地址: https://string-db.org/cgi/download?species_text=Homo+sapiens
- 文件: `9606.protein.links.v12.0.txt`
- 放置位置: `data/STRING/9606.protein.links.v12.0.txt`

### 3. 在 RStudio 中打开项目

双击 `XBCQD-SAALI-Network-Pharmacology.Rproj` 或在 RStudio 中通过 File → Open Project 打开。

### 4. 按顺序运行分析脚本

脚本已按分析流程编号，建议按顺序运行：

1. `R/01_Drug_Target_Collection/01_Drug_target.R`
2. `R/02_Disease_Target_Collection/01_Disease_target.R`
3. `R/03_Drug_Disease_Intersection/01_Venn_Diagram/01_Venn.R`
4. `R/03_Drug_Disease_Intersection/02_GO_KEGG_Enrichment/GO/01_GO_enrich.R`
5. ...（依此类推）

> **注意**: 所有 R 脚本已修改为使用**相对路径**。请确保在运行前将 R 工作目录设置为项目根目录。

---

## 中药复方组成 | Formula Composition

| 中药 | 剂量 (g) | 功效 |
|:---|:---|:---|
| 石膏 (Shi Gao) | 25 | 清热泻火 |
| 大黄 (Da Huang) | 20 | 泻下攻积 |
| 苦杏仁 (Ku Xing Ren) | 10 | 降气止咳 |
| 瓜蒌 (Gua Lou) | 10 | 清热化痰 |

---

## 引用 | Citation

如果您在研究中使用本项目的代码或数据，请引用：

> [Your Paper Citation Here]

---

## 许可 | License

本项目采用 [MIT License](LICENSE) 开源许可。

---

## 联系方式 | Contact

如有任何问题或建议，欢迎通过 GitHub Issues 提交。

