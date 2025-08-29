setwd("/Users/wuz6/Documents/Project/08.spatialTranscriptome/01.RNA/01.edgeR.QCfirst/")
library(TCGAbiolinks)
library("survival")
library("survminer")
library(edgeR)


rm(list=ls())

red_blue_20 = c("#124984","#1f63a8","#2f79b5","#3f8ec0","#5fa5cd","#87beda","#a7d0e4",
                "#c5dfec","#dbeaf2","#edf2f5","#f9f0eb","#fbe3d4","#fbd0b9","#f7b799",
                "#f09c7b","#e17860","#d25849","#c13639","#ae172a","#8a0b25");


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


####
programs <- read.csv("08.meta.programs.csv", head =T )
colnames(programs) <- paste0("MP", 1:4)
programs <- cbind(1, programs)
head(programs)
programs<- reshape2::melt(programs, id.vars = "1")
head(programs)
sum(is.na(programs[,3]))
sum(programs[,3] == "")
programs = programs[programs[,3] != "", -1]
head(programs)

programes1 <- programs[programs[,2] %in% rownames(rna), ]

gene.ann <- data.frame(program = programes1[,1])
row.names(gene.ann) <- programes1[,2]

pheatmap(rna[rownames(gene.ann),], scale = "row",
         cluster_rows = T,
         annotation_row = gene.ann,cutree_cols = 3, show_colnames = F)


###CCGA
rna <- read.table("../../00.data/TCGALGG/CGGA.mRNAseq_693.RSEM-genes.20200506.txt", head =T, row.names = 1)
rna[1:4,1:4]
clin <- read.delim("../../00.data/TCGALGG/CGGA.mRNAseq_693_clinical.20200506.txt", head =T, stringsAsFactors = F, sep = "\t")
head(clin)

clin <- clin[clin$IDH_mutation_status =="Mutant" & clin$X1p19q_codeletion_status == "Non-codel",]

rna <- rna[, colnames(rna) %in% clin$CGGA_ID]
dim(rna)
rna <- DGEList(counts=rna)
rna <- cpm(rna, log = T)

###
programes1 <- programs[programs[,2] %in% rownames(rna), ]
gene.ann <- data.frame(program = programes1[,1])
row.names(gene.ann) <- programes1[,2]

dev.off();dev.off()
... <- lapply(1:4, function(k){
rna1 <- rna[programes1[which(programes1[,1] == paste0("MP", k)),2],]
#gene.cluster <- hclust(dist(rna, method = "euclidean"), method = "ward.D")
gene.cluster <- hclust(as.dist(1 - cor(t(rna1), method = "spearman")), method = "average")
sample.cluster <- hclust(dist(t(rna1),method = "euclidean"), method = "ward.D")
#sample.cluster <- hclust(as.dist(1 - cor(rna, method = "spearman")), method = "average")

sample.ann <- data.frame(cluster = factor(cutree(sample.cluster, k = 2)))

cluster.col <- RColorBrewer::brewer.pal(12, "Paired")[c(4,8)]
names(cluster.col) <- c('1', '2')



clin$cluster <- sample.ann$cluster[match(clin$CGGA_ID, rownames(sample.ann))]
fit <- survfit(Surv(OS/31, Censor..alive.0..dead.1.) ~ cluster, data = clin)
print(fit)



pdf(paste0("08.metaPrograms.MP", k, ".CGGA.survival.pdf"), width = 10, height = 10, useDingbats = F)
print(pheatmap(rna1, scale = "row",
         cluster_rows = gene.cluster,
         cluster_cols = sample.cluster,
         color = red_blue_20, show_colnames = F,
         annotation_col = sample.ann,
         #annotation_row = gene.ann,
         annotation_colors = list(cluster = cluster.col, PC2 =pc.col),
         cutree_cols = 2,
         fontsize_row = 5))

print(ggsurvplot(fit,
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           palette = as.character(cluster.col),
           ggtheme = theme_bw() + theme(panel.grid = element_blank()) # Change ggplot2 theme
           ))

dev.off()

})

