
# ธฝ฿s???ๆ??????IQ
args <- commandArgs(trailingOnly = TRUE)
# ???ๆ??????IQ

lib_path <- args[1]
BSJ_path <- args[2]
gene_path <- args[3]
out_Path <- args[4]

library(edgeR)
#BSJ
lib_mtx <- read.delim(lib_path)
bsj_mtx <- read.delim(BSJ_path, row.names = 1)
gene_mtx <- read.delim(gene_path, row.names = 1)
colnames(lib_mtx)<-c("Sample","Path","Class","Num")
bsj_mtx<-bsj_mtx[,lib_mtx$Sample]
gene_mtx<-gene_mtx[,lib_mtx$Sample]

gene_DGE <- DGEList(counts = gene_mtx, group = lib_mtx$Class)

#??????countแI๎๖
gene_idx <- filterByExpr(gene_DGE)
gene_DGE <- gene_DGE[gene_idx, keep.lib.sizes=FALSE]
#???yป
gene_DGE <- calcNormFactors(gene_DGE)

#ช???
if ("Num" %in% colnames(lib_mtx)) {
  Num <- factor(lib_mtx$Num)
  treat <- factor(lib_mtx$Class, levels=c("Control", "Case"))
  design <- model.matrix(~Num + treat)
} else {
  treat <- factor(lib_mtx$Class, levels=c("Control", "Case"))
  design <- model.matrix(~treat)
}

circ_DGE <- DGEList(counts = bsj_mtx,
                    group = lib_mtx$Class,
                    lib.size = gene_DGE$samples[, "lib.size"],
                    norm.factors = gene_DGE$samples[, "norm.factors"])

#ฤZ๎๖\??????I???Ux
circ_DGE <- estimateDisp(circ_DGE, design, robust = TRUE)

#???๑???๖?????????ซอ^
circ_fit <- glmFit(circ_DGE, design)
circ_lrt <- glmLRT(circ_fit)

circ_df <- circ_lrt$table
circ_order <- order(circ_lrt$table$PValue)
circ_df$DE <- as.vector(decideTestsDGE(circ_lrt))
circ_df <- circ_df[circ_order, ]
circ_df$FDR <- p.adjust(circ_df$PValue, method="fdr")
write.table(circ_df, file=out_Path, quote = FALSE,sep = "\t")
