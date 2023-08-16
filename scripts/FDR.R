# 樃柦椷峴拞???庢??????揑嶲悢
args <- commandArgs(trailingOnly = TRUE)
# ???庢??????揑嶲悢

input_path <- args[1]
out_path <- args[2]
rMats_re <- read.delim(input_path)
rMats_re$FDR<-p.adjust(rMats_re$PValue)
rMats_re<-rMats_re[,c(1,8,9)]
write.table(rMats_re, file=out_path, quote = FALSE,sep = "\t")
