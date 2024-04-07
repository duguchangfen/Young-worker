#!/usr/bin/Rscript
library(ggplot2)
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
#png_draw <- tryCatch({png(args[2])},warning=function(w){print("Differential expression draw warning")},error=function(e){print("Differential expression draw error")})

#detags <- rownames(y)[as.logical(de)];
#plotSmear(et, de.tags=detags)
#abline(h=c(-4, 4), col="blue");
#dev.off()

###导出数据
DE <- et$table
operation_DE <- tryCatch({DE$significant <- as.factor(DE$PValue<0.05 & abs(DE$logFC) >1)},warning=function(w){print("Differential expression operation DE warning")},error=function(e){print("Differential expression operation DE error")})
write.table(DE,file=args[3],sep="\t",na="NA",quote=FALSE)

#volcano图
log2FC_cutoff = log2(1)
padj_cutoff = 0.05
#DEG_DEseq2 <- read.table("/home/chasing/Work/RNA-seq_pipeline/06_deseq/diffexpr_ColCon_vs_ColABA_0.01.txt", header = T, row.names = 1)
DEG_DEseq2 <- read.table(args[3], header = T, row.names = 1)
head(DEG_DEseq2)
##选取差异分析结果
need_DEG <- DEG_DEseq2[,c(1,3)] # 选取log2FoldChange, padj信息
head(need_DEG)
colnames(need_DEG) <- c('log2FoldChange','padj') 

need_DEG$significance  <- as.factor(ifelse(need_DEG$padj < padj_cutoff & abs(need_DEG$log2FoldChange) > log2FC_cutoff,
                                           ifelse(need_DEG$log2FoldChange > log2FC_cutoff ,'UP','DOWN'),'NOT'))
head(need_DEG)

title <- paste0(' Up :  ',nrow(need_DEG[need_DEG$significance =='UP',]) ,
                     '\n Down : ',nrow(need_DEG[need_DEG$significance =='DOWN',]),
                     '\n FoldChange >= ',round(2^log2FC_cutoff,3))

if(nrow(need_DEG[need_DEG$significance =='DOWN',])!=0) {colourss <- c('blue','grey','red')
} else {colourss <- c('grey','red')}

nrow(need_DEG[need_DEG$significance =='DOWN',])

g <- ggplot(data=need_DEG, 
            aes(x=log2FoldChange, y=-log10(padj), 
                color=significance)) +
  #点和背景
  geom_point(alpha=0.4, size=1) +
  theme_classic()+ #无网格线
  #坐标轴
  xlab("log2 ( FoldChange )") + 
  ylab("-log10 ( PValue )") +
  #标题文本
  ggtitle( title ) +
  #分区颜色                  
  scale_colour_manual(values = colourss)+ 
  #辅助线
  geom_vline(xintercept = c(-log2FC_cutoff,log2FC_cutoff),lty=4,col="grey",lwd=0.8) +
  geom_hline(yintercept = -log10(padj_cutoff),lty=4,col="grey",lwd=0.8) +
  #图例标题间距等设置
  theme(plot.title = element_text(hjust = 0.5), 
        plot.margin=unit(c(2,2,2,2),'lines'), #上右下左
        legend.title = element_blank(), #不显示图例标题
        legend.position="right")  #图例位置

ggsave(g,filename = args[2],width =8,height =7.5)
