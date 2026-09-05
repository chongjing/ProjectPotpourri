.libPaths("/rds/project/rds-FTKWLWDeHys/programs/R/R-3.6.0/lib/")
setwd("/rds/user/cx264/hpc-work/project/2.Jeongmin/01.Saskia//4.DEanalysis.FDR_0.05/04.pho2HP-myc_vs_mock/")
library(EnhancedVolcano)
library(gplots)
library(RColorBrewer)
library("ggplot2")
library(edgeR)
library(stringr)
library(statmod)
my_count_matrix <- read.csv("../../4.DEanalysis.FDR_0.01/04.pho2HP-myc_vs_mock/01.raw.count.pho2HP-myc_vs_mock.tab", sep = "\t", header=T)
head(my_count_matrix)
group <- c(1,1,1,1,2,2,2,2) ##first columns (1) are control, (2 is treatment)
Time <- factor(c("1","2","3","4","1","2","3","4"))
y <- DGEList(counts=my_count_matrix[,2:9], group=group,genes=my_count_matrix[,1])
keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y)
y <- calcNormFactors(y)
head(y)
design <- model.matrix(~Time+group)
rownames(design) <- colnames(y)
design
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
# [1] 0.1822498
fit <- glmQLFit(y, design, robust=TRUE)
qlf <- glmQLFTest(fit, coef=2:3)
topTags(qlf)
FDR <- p.adjust(qlf$table$PValue, method="BH")
sum(FDR < 0.05)
qlf <- glmQLFTest(fit)
topTags(qlf)
cpm(y)[rownames(topTags(qlf)),]
summary(decideTests(qlf,p.value=0.05))
#       group
#Down    72
#NotSig 23966
#Up      479

allTags2 <- topTags(qlf, n = nrow(qlf$genes), adjust.method = "BH", sort.by = "none", p.value = 1)
DEG_glmQLF <- as.data.frame(allTags2)
k1 <- (DEG_glmQLF$PValue < 0.05) & (DEG_glmQLF$logFC < -1)
k2 <- (DEG_glmQLF$PValue < 0.05) & (DEG_glmQLF$logFC > 1)
DEG_glmQLF$change <- ifelse(k1, "DOWN", ifelse(k2, "UP", "Not_significant"))
table(DEG_glmQLF$change)
DEG_glmQLF$CPM <- cpm(y)[rownames(qlf),]
write.csv(DEG_glmQLF, "1.4.QLF.DE.results.pho2HP-myc_vs_mock.FINAL.csv", quote = F)

#heatmap plot
deseq_results <- read.csv("1.4.QLF.DE.results.pho2HP-myc_vs_mock.FINAL.csv",sep=",", header = T)
diffGenes <- deseq_results[deseq_results$PValue < 0.05,][,c(2,9,10,11,12,13,14,15,16)]
clustRows <- hclust(as.dist(1-cor(t(diffGenes[,c(2,3,4,5,6,7,8,9)]),method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(diffGenes[,c(2,3,4,5,6,7,8,9)],method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9)
module.color <- module.color[as.vector(module.assign)]
myheatcolors1 <- bluered(75)
pdf("1.5.QLF.DiffGene.pho2HP-myc_vs_mock.heatmap.pdf")
heatmap.2(data.matrix(diffGenes[,c(2,3,4,5,6,7,8,9)]),
    Rowv=as.dendrogram(clustRows),
    Colv=as.dendrogram(clustColumns),
    RowSideColors=module.color,
    col=myheatcolors1, scale='row', labRow=NA,
    density.info="none", trace="none",
    cexRow=1, cexCol=1, margins=c(8,20), keysize=1, srtCol=45)
dev.off()
png("1.5.QLF.DiffGene.pho2HP-myc_vs_mock.heatmap.png",width = 636, height = 980)
heatmap.2(data.matrix(diffGenes[,c(2,3,4,5,6,7,8,9)]),
    Rowv=as.dendrogram(clustRows),
    Colv=as.dendrogram(clustColumns),
    RowSideColors=module.color,
    col=myheatcolors1, scale='row', labRow=NA,
    density.info="none", trace="none",
    cexRow=1, cexCol=1, margins=c(8,20), keysize=0.5, srtCol=45)
dev.off()

#volcano plot
pdf("1.6.QLF.DiffGene.pho2HP-myc_vs_mock.volcano.pdf", 7, 7)
EnhancedVolcano(deseq_results, lab = rownames(deseq_results), x = 'logFC', y = 'FDR', pCutoff = 0.05, 
    FCcutoff=1.0, ylim = c(0, 5), xlim = c(-10, 15), pointSize = 1.0, labSize = 0, 
    colAlpha = 1, legendIconSize=2.0, legendLabSize = 10) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
dev.off()
png("1.6.QLF.DiffGene.pho2HP-myc_vs_mock.volcano.png", 800, 800)
EnhancedVolcano(deseq_results, lab = rownames(deseq_results), x = 'logFC', y = 'FDR', pCutoff = 0.05, 
    FCcutoff=1.0, ylim = c(0, 5), xlim = c(-10, 15), pointSize = 1.0, labSize = 0, 
    colAlpha = 1, legendIconSize=2.0, legendLabSize = 10) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
dev.off()
q()

