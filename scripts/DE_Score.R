# 从命令行中获取传递的参数
args <- commandArgs(trailingOnly = TRUE)
# 获取传递的参数
BSJ_count_Path <- args[1]
infor_Path <- args[2]
outPath <- args[3]
factor <- as.numeric(args[4])
pval <- as.numeric(args[5])
# 设置随机种子以确保结果可重复
set.seed(123)
de_score <- function(a, b, factor = 1, size = 100000, pval = 0.05) {
  x1 <- rgamma(size, shape = a + 1)
  x2 <- rgamma(size, shape = b + 1)
  
  z <- x1 / x2
  z <- sort(z)
  
  if (mean(z) * factor >= 1) {
    score <- max(log2(z[round(size * pval)] * factor), 0)
  } else {
    score <- min(log2(z[round(size * (1 - pval))] * factor), 0)
  }
  return(score)
}
#导入信息
BSJ_count<-read.delim(BSJ_count_Path)
infor<-read.delim(infor_Path)
BSJ_count<-BSJ_count[,c(0,which(infor$Class=="Case"),which(infor$Class=="Control"))+1]
for (i in 1:nrow(BSJ_count)) {
  tmp_de <- de_score(max(BSJ_count[i,2], 1), max(BSJ_count[i,3], 1), 1.0 / factor, size=100000, pval=pval)
  BSJ_count$DE_Score[i] <-tmp_de
}
##导出信息
write.table(BSJ_count,outPath,sep = "\t",row.names = F)




