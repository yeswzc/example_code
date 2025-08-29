setwd("/Users/wuz6/Documents/Project/08.spatialTranscriptome/01.RNA/01.edgeR.QCfirst/")
library(TCGAbiolinks)
library("survival")
library("survminer")
library(edgeR)
library(clusterProfiler);
library(dplyr)
library(pheatmap)
library(factoextra)

rm(list=ls())

red_blue_20 = c("#124984","#1f63a8","#2f79b5","#3f8ec0","#5fa5cd","#87beda","#a7d0e4",
                "#c5dfec","#dbeaf2","#edf2f5","#f9f0eb","#fbe3d4","#fbd0b9","#f7b799",
                "#f09c7b","#e17860","#d25849","#c13639","#ae172a","#8a0b25");

#c5 = msigdbr::msigdbr(species = "human", category = "C5")
db = msigdbr::msigdbr(species = "human", category = "C5") %>% dplyr::distinct(gene_symbol, gs_name) %>% as.data.frame()

###
lgg_clin <- TCGAquery_subtype(tumor = "lgg")
lgg_clin <- lgg_clin[lgg_clin$IDH.status == "Mutant" & lgg_clin$X1p.19q.codeletion == "non-codel",]
nrow(lgg_clin)

rna <- read.table("../../00.data/TCGALGG/TCGA.LGG.sampleMap_HiSeqV2.gz", head =T, row.names = 1, check.names = F)
rna[1:4,1:4]

lgg_clin[1:4,1:4]

colnames(rna) <- sapply(colnames(rna), function(x) substr(x, 1,12))
colnames(rna)
overlap <- intersect(colnames(rna), lgg_clin$patient)

lgg_clin <- lgg_clin[lgg_clin$patient %in% overlap,]
rna <- rna[,colnames(rna) %in% overlap]

rna <- DGEList(counts=rna)
rna <- cpm(rna, log = T)

plot(density(rna))
dim(rna)
dim(lgg_clin)

survival.data <- subset(lgg_clin, select = c("Survival..months.", "Vital.status..1.dead."))
head(survival.data)

###
gene.sd <- apply(rna, 1, sd)
gene.sd[1:4]

rna.top1k <- rna[order(gene.sd, decreasing = T)[1:1000],]

pc <- prcomp(t(rna.top1k))

pc.res <- data.frame(pc$x)
dim(pc.res)
pc.res[1:4,1:4]
ggplot(pc.res, aes(x = PC1, y = PC2)) + geom_point() + theme_bw()


pc2.top.genes = rownames(rna.top1k)[order(pc$rotation[,2], decreasing = T)[1:30]]
pc2.low.genes <- rownames(rna.top1k)[order(pc$rotation[,2], decreasing = F)[1:30]]


####
dsp.pc2.high <- read.table("01.CTA.PC2.high.genes.txt", head =F)[,1]
dsp.pc2.low <- read.table("01.CTA.PC2.low.genes.txt", head =F)[,1]
dsp.pc2.high <- dsp.pc2.high[dsp.pc2.high %in% rownames(rna)];length(dsp.pc2.high)
dsp.pc2.low <- dsp.pc2.low[dsp.pc2.low %in% rownames(rna)]; length(dsp.pc2.low)

gene.ann <- data.frame(class = c(rep("PC2 high", length(dsp.pc2.high)), 
                                 rep("PC2 low", length(dsp.pc2.low))))
row.names(gene.ann) <- c(dsp.pc2.high, dsp.pc2.low)

pheatmap(rna[c(dsp.pc2.high, dsp.pc2.low),], scale = "none",
         cluster_rows = T,
         annotation_row = gene.ann,cutree_cols = 3, show_colnames = F)


pc <- prcomp(t(rna[c(dsp.pc2.high, dsp.pc2.low),]))
pc.res <- data.frame(pc$x)
ggplot(pc.res, aes(x = PC1, y = PC2)) + geom_point() + theme_bw()



###CCGA
rna <- read.table("../../00.data/TCGALGG/CGGA.mRNAseq_693.RSEM-genes.20200506.txt", head =T, row.names = 1)
rna[1:4,1:4]
clin <- read.delim("../../00.data/TCGALGG/CGGA.mRNAseq_693_clinical.20200506.txt", head =T, stringsAsFactors = F, sep = "\t")
head(clin)

clin <- clin[clin$IDH_mutation_status =="Mutant" & clin$X1p19q_codeletion_status == "Non-codel",]

table(clin$PRS_type)
table(clin$Histology)
#table(clin$)

rna <- rna[, colnames(rna) %in% clin$CGGA_ID]
dim(rna)

rna <- DGEList(counts=rna)
rna <- cpm(rna, log = T)
#rna <- rna[grep("^MT-", rownames(rna), invert = T), ]

gene.sd <- apply(rna, 1, sd)
gene.sd[1:4]

rna.top1k <- rna[order(gene.sd, decreasing = T)[1:1000],]

pc <- prcomp(t(rna.top1k))

pc.res <- data.frame(pc$x[,1:4])
pc.res <- cbind(pc.res, clin[match(rownames(pc.res), clin$CGGA_ID),])

dim(pc.res)
pc.res[1:4,1:4]
ggplot(pc.res, aes(x = PC1, y = PC2, col = Grade)) + geom_point() + theme_bw()
ggplot(pc.res, aes(x = PC1, y = PC2, col = Histology)) + geom_point() + theme_bw()

ggplot(pc.res, aes(x = PC3, y = PC4, col = Histology)) + geom_point() + theme_bw()

fviz_pca_var(pc, col.var = "black", repel = TRUE)




###
dsp.pc1.high <- read.table("01.CTA.PC1.high.genes.txt", head =F)[,1]; length(dsp.pc1.high)
dsp.pc1.low <- read.table("01.CTA.PC1.low.genes.txt", head =F)[,1]; length(dsp.pc1.low)
dsp.pc1.high <- dsp.pc1.high[dsp.pc1.high %in% rownames(rna)];length(dsp.pc1.high) #0
dsp.pc1.low <- dsp.pc1.low[dsp.pc1.low %in% rownames(rna)]; length(dsp.pc1.low) #8

dsp.pc2.high <- read.table("01.CTA.PC2.high.genes.txt", head =F)[,1]
dsp.pc2.low <- read.table("01.CTA.PC2.low.genes.txt", head =F)[,1]
dsp.pc2.high <- dsp.pc2.high[dsp.pc2.high %in% rownames(rna)];length(dsp.pc2.high) #28
dsp.pc2.low <- dsp.pc2.low[dsp.pc2.low %in% rownames(rna)]; length(dsp.pc2.low) # 27

rna <- rna[c(dsp.pc2.high, dsp.pc2.low),]

gene.ann <- data.frame(PC2 = c(rep("High", length(dsp.pc2.high)), 
                                 rep("Low", length(dsp.pc2.low))))
row.names(gene.ann) <- c(dsp.pc2.high, dsp.pc2.low)


#gene.cluster <- hclust(dist(rna, method = "euclidean"), method = "ward.D")
gene.cluster <- hclust(as.dist(1 - cor(t(rna), method = "spearman")), method = "average")
sample.cluster <- hclust(dist(t(rna),method = "euclidean"), method = "ward.D")
#sample.cluster <- hclust(as.dist(1 - cor(rna, method = "spearman")), method = "average")

sample.ann <- data.frame(cluster = factor(cutree(sample.cluster, k = 3)),
                         Grade = clin$Grade[match(colnames(rna), clin$CGGA_ID)])

cluster.col <- RColorBrewer::brewer.pal(12, "Paired")[c(2,4,8)]
names(cluster.col) <- c('1', '2', '3')
pc.col <- c(red_blue_20[20], red_blue_20[1])
names(pc.col) <- c("High", "Low")
tumor.grade.col <- RColorBrewer::brewer.pal(8, "Greys")[c(4,6,8,1)]
names(tumor.grade.col) <- unique(sample.ann$Grade)
dev.off()

pdf("05.CGGA.dspPC2LowHighGene.heatmap.pdf", width = 8, height = 5,useDingbats = F)
pheatmap(rna, scale = "row",
         cluster_rows = gene.cluster,
         cluster_cols = sample.cluster,
         color = red_blue_20, show_colnames = F,
         annotation_col = sample.ann,
         annotation_row = gene.ann,
         annotation_colors = list(cluster = cluster.col, PC2 =pc.col, Grade = tumor.grade.col),
          cutree_cols = 3,
         fontsize_row = 5)


dev.off()

clin$cluster <- sample.ann$cluster[match(clin$CGGA_ID, rownames(sample.ann))]
fit <- survfit(Surv(OS/31, Censor..alive.0..dead.1.) ~ cluster, data = clin)
print(fit)

pdf("05.CGGA.dspPC2LowHighGeneCluster.survival.pdf", width = 10, height = 10, useDingbats = F)
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           palette = as.character(cluster.col),
           ggtheme = theme_bw() + theme(panel.grid = element_blank()) # Change ggplot2 theme
           ) 

dev.off()

clin1 <- clin[clin$cluster == 1 | clin$cluster==3,]

fit1 <- survfit(Surv(OS/31, Censor..alive.0..dead.1.) ~ cluster, data = clin1)
print(fit1)

pdf("05.CGGA.dspPC2LowHighGeneCluster13.survival.pdf", width = 10, height = 10, useDingbats = F)
ggsurvplot(fit1,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           palette = as.character(cluster.col[c(1,3)]),
           ggtheme = theme_bw() + theme(panel.grid = element_blank()), # Change ggplot2 theme
) 
dev.off()

identical(rownames(pc.res), names(sample.cluster$labels))
pc.res$cluster <- factor(cutree(sample.cluster, 3)[match(row.names(pc.res), sample.cluster$labels)])
ggplot(pc.res, aes(x = PC1, y = PC2, col = cluster)) + geom_point() + theme_bw()
