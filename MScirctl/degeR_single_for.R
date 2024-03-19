#!/usr/bin/Rscript

args <- commandArgs(T)
print(args[1])
print(args[2])

library("edgeR")
#library('ggplot2')


###读取数据
#setwd("G:/My_exercise/edgeR")
rawdata <- read.table(args[1],header=TRUE,row.names=1,check.names = FALSE)
head(rawdata)
#重命名列名
names(rawdata) <- c("F.1yr.OC.count","M.1yr.OC.count")
#进行分组
group <- factor(c("F.1yr.OC.count","M.1yr.OC.count"))

###过滤与标准化
y <- DGEList(counts=rawdata,genes=rownames(rawdata),group = group)
###TMM标准化
y<-calcNormFactors(y)
y$samples
###推测离散度,根据经验设置，若样本是人，设置bcv = 0.4，模式生物设置0.1.（这里没有经验，我就多试几个）
#bcv <- 0.1
bcv <- 0.2
#bcv <- 0.4
et <- exactTest(y, dispersion=bcv^2)
topTags(et)
summary(de <- decideTestsDGE(et))
###图形展示检验结果
png_draw <- tryCatch({png(args[2])},warning=function(w){print("Differential expression draw warning")},error=function(e){print("Differential expression draw error")})

detags <- rownames(y)[as.logical(de)];
plotSmear(et, de.tags=detags)
abline(h=c(-4, 4), col="blue");
dev.off()

###导出数据
DE <- et$table
operation_DE <- tryCatch({DE$significant <- as.factor(DE$PValue<0.05 & abs(DE$logFC) >1)},warning=function(w){print("Differential expression operation DE warning")},error=function(e){print("Differential expression operation DE error")})
write.table(DE,file=args[3],sep="\t",na="NA",quote=FALSE)
#write.csv(DE, "edgeR.F-vs-M.OC2.csv")

#DE2 <- topTags(et,20000)$table
#DE2$significant <- as.factor(DE2$PValue<0.05 & DE2$FDR<0.05 & abs(DE2$logFC) >1)
#write.csv(DE2, "F_1yr.OC-vs-M_1yr.OC3.csv")
