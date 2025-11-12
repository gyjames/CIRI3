# 樃柦椷峴拞???庢??????揑嶲悢
args <- commandArgs(trailingOnly = TRUE)
# ???庢??????揑嶲悢

input_path <- args[1]
out_path <- args[2]
infor_path<- args[3]
BSJ_path<- args[4]
FSJ_path<- args[5]

infor<-read.table(infor_path,header = T)
BSJ_count<-read.table(BSJ_path,header = T,row.names = 1)
FSJ_count<-read.table(FSJ_path,header = T,row.names = 1)
ratio_count<-BSJ_count*2/(BSJ_count*2+FSJ_count)
ratio_count[is.nan(as.matrix(ratio_count))] <- 0

rMats_re <- read.delim(input_path)
rMats_re$FDR<-p.adjust(rMats_re$PValue)
rMats_re<-rMats_re[,c(1,8,9)]
ratio_count<-ratio_count[rMats_re$circRNA_ID,]



calc_FC <- function(expr, meta) {
  # expr: 基因表达矩阵 (行=基因, 列=样本)
  # meta: 样本信息表，第一列是样本名，第二列是分组信息
  stopifnot(all(meta[,1] %in% colnames(expr)))  # 检查样本是否匹配
  
  # 按样本表顺序重排表达矩阵
  expr <- expr[, meta[,1]]
  
  # 提取分组信息
  group <- as.factor(meta[,2])
  group_levels <- levels(group)
  
  if (length(group_levels) != 2) {
    stop("The number of groups is not 2. Please make sure that the second column of 'meta' contains exactly two groups.")
  }
  
  group1 <- group_levels[1]
  group2 <- group_levels[2]
  
  message(paste0("Comparing ", group2, " vs ", group1, " ..."))
  
  # 获取对应样本
  samples1 <- meta[group == group1, 1]
  samples2 <- meta[group == group2, 1]
  
  # 计算每组平均表达
  mean1 <- rowMeans(expr[, samples1, drop=FALSE], na.rm = TRUE)
  mean2 <- rowMeans(expr[, samples2, drop=FALSE], na.rm = TRUE)
  
  # 计算 Fold Change 和 log2FC（加伪计数避免除零）
  FC <- (mean1 + 1e-6) / (mean2 + 1e-6)
  log2FC <- log2(FC)
  
  # 计算Delta PSI
  DeltaPSI<-mean1 - mean2
  
  # 计算在多少样本中检测到
  PositiveRate<-rowMeans(expr > 0)
  
  # 输出结果
  result <- setNames(
    data.frame(
      rownames(expr),
      mean1,
      mean2,
      FC,
      log2FC,
      DeltaPSI,
      PositiveRate
    ),
    c("circRNA", paste0("Mean_", group1), paste0("Mean_", group2),
      "FC", "log2FC", "DeltaPSI", "PositiveRate")
  )
  
  return(result)
}


result<-calc_FC(ratio_count,infor)
result<-data.frame(result,rMats_re[,2:3])

write.table(result, file=out_path, quote = FALSE,sep = "\t")
