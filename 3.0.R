# 安装shiny包
install.packages("shiny")
# 安装和加载所需的包

install.packages("readxl")
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC")

library(ANCOMBC)
library(dplyr)
library(readxl)
library(nloptr)
library(dplyr)
source("F:/zby/tsaas/3 BAs/统计检验/ancom_bc/bc1.R")

# 从Excel文件中导入绝对丰度数据和元数据
# 假设您的文件分别命名为 "abundance_data.xlsx" 和 "metadata.xlsx"
# 并且数据在各自文件的第一个工作表上
abundance_table <- read_excel("multi-ab.xlsx")
meta_data <- read_excel("information.xlsx")

# 转换绝对丰度数据格式
species_names <- abundance_table$species
abundance_table <- as.matrix(abundance_table[,-1])
rownames(abundance_table) <- species_names

# 转换元数据格式
# 假设第一列是样本ID，第二列是分组变量
colnames(meta_data) <- c("sample_id", "group")

#two
# 数据预处理
pre.process <- feature_table_pre_process(feature.table = abundance_table,
                                         meta.data = meta_data,
                                         sample.var = "sample_id",
                                         group.var = "group",
                                         zero.cut = 0.90,
                                         lib.cut = 1000,
                                         neg.lb = TRUE)
feature.table <- pre.process$feature.table
group.name <- pre.process$group.name
group.ind <- pre.process$group.ind
struc.zero <- pre.process$structure.zeros

# 设置ANCOM-BC分析参数
grp.name <- group.name
grp.ind <- group.ind
adj.method <- "bonferroni"
tol.EM <- 1e-5
max.iterNum <- 100
perNum <- 1000
alpha <- 0.05

# 执行ANCOM-BC分析
out <- ANCOM_BC(feature.table = feature.table,
                grp.name = grp.name,
                grp.ind = grp.ind,
                struc.zero = struc.zero,
                adj.method = adj.method,
                tol.EM = tol.EM,
                max.iterNum = max.iterNum,
                perNum = perNum,
                alpha = alpha)

# 导出结果
res <- cbind(taxon = rownames(out$feature.table), out$res)
write.csv(res, "ANCOMBC_results.csv")

#multi
