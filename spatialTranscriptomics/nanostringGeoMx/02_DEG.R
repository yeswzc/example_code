setwd("/Users/wuz6/Documents/Project/08.spatialTranscriptome/01.RNA/02.DEG.samplewise/")
rm(list=ls())
library(edgeR)
library(ggrepel)

library(fgsea)

library(clusterProfiler);
library(dplyr)
c5 = msigdbr::msigdbr(species = "human", category = "C5")
h = msigdbr::msigdbr(species = "human", category = "H")
c5.fgsea.db = split(x = c5$gene_symbol, f = c5$gs_name)
names(c5.fgsea.db) = gsub("^GO\\w\\w\\_", "", names(c5.fgsea.db))
h.fgsea.db = split(x = h$gene_symbol, f = h$gs_name)
names(h.fgsea.db) <- gsub("^HALLMARK_", "", names(h.fgsea.db))

c5 = msigdbr::msigdbr(species = "human", category = "C5") %>% dplyr::distinct(gene_symbol, gs_name) %>% as.data.frame()
c2 = msigdbr::msigdbr(species = "human", category = "C2") %>% dplyr::distinct(gene_symbol, gs_name) %>% as.data.frame()
c6 = msigdbr::msigdbr(species = "human", category = "C6") %>% dplyr::distinct(gene_symbol, gs_name) %>% as.data.frame()
h = msigdbr::msigdbr(species = "human", category = "H") %>% dplyr::distinct(gene_symbol, gs_name) %>% as.data.frame()
#c5.fgsea <- split(x = c5.raw$gene_symbol, f = c5.raw$gs_name)
pair.col <- RColorBrewer::brewer.pal(12, "Paired")
loc1.col <- c("#808080", pair.col[6], "#D2691E", "#00FFFF")
names(loc1.col) <- c("GM", "BV", "TC", "IE")

red_blue_20 = c("#124984","#1f63a8","#2f79b5","#3f8ec0","#5fa5cd","#87beda","#a7d0e4",
                "#c5dfec","#dbeaf2","#edf2f5","#f9f0eb","#fbe3d4","#fbd0b9","#f7b799",
                "#f09c7b","#e17860","#d25849","#c13639","#ae172a","#8a0b25");
red_blue_50 = c("#0a3b70","#10457e","#15508d","#1b5a9c","#2065ab","#276eb0","#2e77b5",
                "#3480b9","#3b88be","#4291c2","#4f9bc7","#5fa5cd","#6eaed2","#7eb8d7",
                "#8dc2dc","#9bc9e0","#a7d0e4","#b3d6e8","#c0dceb","#cce2ef","#d5e7f1",
                "#ddebf2","#e4eef4","#ecf2f5","#f3f5f6","#f8f4f2","#f9efe9","#fae9df",
                "#fbe4d6","#fcdecd","#fcd7c2","#fbccb4","#f9c2a7","#f7b799","#f5ac8b",
                "#f2a17f","#ec9374","#e6866a","#e17860","#db6b55","#d55d4c","#ce4f45",
                "#c6413e","#bf3338","#b82531","#b1182b","#a21328","#930e26","#840924",
                "#760521");

sample.names = c(1399, 1905, 2023, 2751, 3984, 4461, 604, 6614, "8147M", "8147O")
sample.colors = RColorBrewer::brewer.pal(12, "Paired")[c(1:10)]
names(sample.colors) = sample.names

### Tumor grade
low.grade <- rep("low", 3)
names(low.grade) <- c(6614, 2751, 3984)
high.grade <- rep("high", 7)
names(high.grade) <- c(1905, 604, 1399, 2023,4461, "8147O", "8147M")
tumor.grade <- c(low.grade, high.grade)

grade.col <- c("#888888", "#000000")
names(grade.col) <- c("low", "high")

p.cutoff = 0.01

load("../01.edgeR.QCfirst/00.IDHastrocytoma.CTA.filtered.edgeR.obj.rda")

y <- y[ rownames(y) != "Negative Probe", ]

#cpm <- read.csv("../01.edgeR/00.IDHastrocytoma.CTA.cpm.filtered.csv", head =T, row.names = 1, check.names = F)
#head(cpm[,1:4])
#rownames(y)[1:4]
#y <- y[rownames(y) %in% rownames(cpm),colnames(y) %in% colnames(cpm)]
dim(y)

table(meta$annotation1)
colnames(y)[1:4]
y.ann.0 <- meta$annotation1[match(colnames(y), meta$ID)]
y.sample.ID <- sapply(colnames(y), function(x) unlist(strsplit(x, "\\-"))[1])
y.ann <- paste0(y.ann.0, ".", y.sample.ID)
y.ann[grep("GM", y.ann)] = "GM"
table(y.ann)

y.ann <- as.factor(y.ann)
design <- model.matrix(~0+y.ann)
colnames(design) <- gsub("^y.ann", "", colnames(design))

y <- estimateDisp(y, design)
fit <- glmQLFit(y,design)




#######################################################################################################################################
###1.GM - TC, sample-wise
table(y.ann)
TC.list <- grep("TC", unique(y.ann), value = T)
TC.list

qlf1 <- glmQLFTest(fit, contrast = makeContrasts(TC.1399-GM, levels=design))
rm(qlf1)

y$samples$group <- y.ann

TC.GM.DEG <- lapply(TC.list, function(id){
  sample.ID <- unlist(strsplit(id, "\\."))[2]
  message(id, ":", sample.ID)
  #qlf1 <- glmQLFTest(fit, contrast = makeContrasts(id-GM, levels=design))
  qlf1 <- exactTest(y, pair=c("GM", id))
  topTags(qlf1)
  deg <- qlf1$table
  deg$gene <- rownames(deg)
  deg$p.adjust = p.adjust(deg$PValue, method = "fdr")
  deg$sample = sample.ID
  
  write.csv(deg, file = paste0("01.", sample.ID,".TC-GM.DEG.csv"), row.names = F, quote = F)
  over.expressed <- rownames(deg)[deg$p.adjust<p.cutoff & deg$logFC > 0]; length(over.expressed)
  under.expressed <- rownames(deg)[deg$p.adjust<p.cutoff & deg$logFC < 0]; length(under.expressed)
  p.x <- deg$p.adjust[deg$p.adjust != 0]
  min.p <- min(p.x)
  deg$p.adjust[deg$p.adjust ==0] = min.p/10
  
  p <- ggplot(deg, aes(x = logFC, y = -log10(p.adjust), label = gene)) + 
    geom_point(size = 0.1, color = ifelse(deg$p.adjust <p.cutoff, "black", "gray")) + 
    geom_text_repel(data = deg[deg$gene %in% c(over.expressed, under.expressed),], 
                    min.segment.length = Inf, size = 1, seed = 123)+
    geom_hline(yintercept = -log10(p.cutoff), linetype = "dashed")+
    theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank())+
    scale_y_continuous(expand = expansion(mult = c(0, .1)))
    labs(title = paste0(sample.ID, ": TC v.s. GM"))
  
  pdf(paste0("01.", sample.ID,".TC-GM.DEG.pdf"), width = 4, height = 4, useDingbats = F)
  print(p)
  dev.off()
  deg$sample = sample.ID
  return(deg)
  
})
length(TC.GM.DEG)
TC.GM.DEG <- do.call(rbind, TC.GM.DEG)

up.gene <- table(TC.GM.DEG$gene, TC.GM.DEG$p.adjust <p.cutoff & TC.GM.DEG$logFC> 0)
down.gene <- table(TC.GM.DEG$gene, TC.GM.DEG$p.adjust <p.cutoff & TC.GM.DEG$logFC< 0)

up.gene <- up.gene[up.gene[,2] >= 10,]
down.gene <- down.gene[down.gene[,2] >= 10,]
nrow(up.gene)
nrow(down.gene)

###Plot heatmap of DEGs
cpm <- cpm(y, log = T)

gene.ann <- data.frame(DEG = c(rep("up", nrow(up.gene)), rep("down", nrow(down.gene))))
rownames(gene.ann) <- c(rownames(up.gene), rownames(down.gene))

roi.ann <- data.frame(location = meta$annotation1[match(colnames(cpm), meta$ID)] )
rownames(roi.ann) <- colnames(cpm)
roi.ann$sample <- sapply(rownames(roi.ann), function(x) unlist(strsplit(x, "\\-"))[1])
roi.ann$grade <- tumor.grade[roi.ann$sample]
table(roi.ann)

deg.col <- c(pair.col[c(6,2)]) #, "#808080")
names(deg.col) <- c("up", "down") #, "non-DEG")

annotation.col <- list(location = loc1.col, 
                       grade = grade.col,
                       sample = sample.colors,
                       DEG = deg.col)

###oder ROIs
gm.cluster <- hclust(as.dist(1- cor(cpm[rownames(gene.ann), which(roi.ann$location == "GM")], method = "spearman")), method = "average")
tc.cluster <- hclust(as.dist(1- cor(cpm[rownames(gene.ann), which(roi.ann$location == "TC")], method = "spearman")), method = "average")
ie.cluster <- hclust(as.dist(1- cor(cpm[rownames(gene.ann), which(roi.ann$location == "IE")], method = "spearman")), method = "average")
bv.cluster <- hclust(as.dist(1- cor(cpm[rownames(gene.ann), which(roi.ann$location == "BV")], method = "spearman")), method = "average")

roi.order <- c(gm.cluster$labels[gm.cluster$order],
               ie.cluster$labels[ie.cluster$order],
               tc.cluster$labels[tc.cluster$order],
               bv.cluster$labels[bv.cluster$order])

#cpm.z <- t(scale(t(cpm[rownames(gene.ann),roi.order]), scale = T, center = T))
roi.cluster <- hclust(as.dist(1 - cor(cpm[rownames(gene.ann),], method = "spearman")), method = "average")
#roi.cluster <- hclust(dist(t(cpm[rownames(gene.ann),]), method = "euclidean"), method = "ward.D2")
#roi.cluster <- hclust(dist(t(cpm.z), method = "euclidean"), method = "ward.D2")
#gene.cluster <- hclust(dist(cpm[all.deg,roi.order], method = "euclidean"), method = "ward.D2")
#gene.cluster <- hclust(dist(cpm.z, method = "euclidean"), method = "ward.D2")
gene.cluster <- hclust(as.dist(1 - cor(t(cpm[rownames(gene.ann),]), method = "spearman")), method = "average")

head(TC.GM.DEG)

mean.logFC <- aggregate(TC.GM.DEG$logFC, list(rownames(TC.GM.DEG)), mean)
mean.logFC <- mean.logFC[order(mean.logFC[,2]),]

mean.logFC <- mean.logFC[mean.logFC[,1] %in% rownames(gene.ann),]

pdf("01.DEG.heatmap.pdf", width = 12, height = 12, useDingbats = F)
pheatmap(
  cpm[rownames(gene.ann),],
  cluster_cols = roi.cluster,
  #cluster_cols = F,
  cluster_rows = gene.cluster,
  #cluster_rows = F,
  scale = "row",
  annotation_row = gene.ann, 
  annotation_col = roi.ann,
  annotation_colors = annotation.col,
  color = red_blue_50,
  show_rownames = T, show_colnames = F)
dev.off()


##############################################3
write.table(gene.ann, file = "01.TC-GM.allSampleConsistentDEG.tsv", quote = F, col.names = F)

over.expressed.GO = data.frame(enricher(gene = rownames(up.gene), TERM2GENE = c5, universe = row.names(cpm)))
over.expressed.H = data.frame(enricher(gene = rownames(up.gene), TERM2GENE = h, universe = row.names(cpm)))

###
under.expressed.GO = data.frame(enricher(gene = rownames(down.gene), TERM2GENE = c5, universe = row.names(cpm)))
under.expressed.H = data.frame(enricher(gene = rownames(down.gene), TERM2GENE = h, universe = row.names(cpm)))


write.csv(over.expressed.GO, file = "01.TC-GM.overexpressed.GO.csv", row.names = F)
write.csv(over.expressed.H, file = "01.TC-GM.overexpressed.H.csv", row.names = F)
write.csv(under.expressed.GO, file = "01.TC-GM.underexpressed.GO.csv", row.names = F)
write.csv(over.expressed.H, file = "01.TC-GM.underexpressed.H.csv", row.names = F)

over.expressed.H$ID = gsub("^HALLMARK_", "", over.expressed.H$ID)
under.expressed.H$ID = gsub("^HALLMARK_", "", under.expressed.H$ID)
over.expressed.GO$ID = gsub("^GO\\w\\w\\_", "", over.expressed.GO$ID)
under.expressed.GO$ID = gsub("^GO\\w\\w\\_", "", under.expressed.GO$ID)
#RColorBrewer::display.brewer.all()
#reds <- RColorBrewer::brewer.pal(9, "Reds")

nrow(over.expressed.H)
p1 <- ggplot(over.expressed.H, aes(x =  reorder(ID, -p.adjust), y = -log10(p.adjust), fill = Count)) + 
  geom_bar(stat = "identity")+ #fill = pair.col[2]) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_label(aes(y = 1, label=GeneRatio),  nudge_y = -0.5, fill = "white", alpha = 0.7, size = 2, label.size = NA)+
  scale_fill_continuous(low = red_blue_20[1], high = red_blue_20[20])+
  scale_y_continuous(expand = expansion(add = c(0, .05))) +
  coord_flip() + theme_bw()+ theme(legend.position = "none")+  labs(x = "")
p1
#pdf("01.TC-GM.overexpressed.H.pdf", width = 7, height = 3, useDingbats = F)
#p1
#dev.off()

nrow(under.expressed.H)
p2 <- ggplot(under.expressed.H, aes(x =  reorder(ID, -p.adjust), y = -log10(p.adjust), fill = Count)) + 
  geom_bar(stat = "identity")+ #fill = pair.col[2]) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_label(aes(y = 1.1, label=GeneRatio),  nudge_y = -0.5, fill = "white", alpha = 0.7, size = 2, label.size = NA)+
  scale_fill_continuous(low = red_blue_20[1], high = red_blue_20[20])+
  scale_y_continuous(expand = expansion(add = c(0, .05))) +
  coord_flip() + theme_bw()+ theme(legend.position = "none")+  labs(x = "")
#pdf("01.TC-GM.underexpressed.H.pdf", width = 7, height = 3, useDingbats = F)
#p2
#dev.off()

nrow(over.expressed.GO)
#over.expressed.GO.1 <- over.expressed.GO[order(over.expressed.GO$p.adjust, decreasing = F)[1:20],]
p3 <- ggplot(over.expressed.GO, aes(x =  reorder(ID, -p.adjust), y = -log10(p.adjust), fill = Count)) + 
  geom_bar(stat = "identity")+ #fill = pair.col[2]) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_label(aes(y = 1.1, label=GeneRatio),  nudge_y = -0.5, fill = "white", alpha = 0.7, size = 2, label.size = NA)+
  scale_fill_continuous(low = red_blue_20[1], high = red_blue_20[20])+
  scale_y_continuous(expand = expansion(add = c(0, .05))) +
  coord_flip() + theme_bw()+ theme(legend.position = "none")+  labs(x = "")
p3

pdf("01.TC-GM.overexpressed.GO.pdf", width = 14, height = 1.2, useDingbats = F)
p3
dev.off()

nrow(under.expressed.GO)
under.expressed.GO.1 <- under.expressed.GO[order(under.expressed.GO$p.adjust, decreasing = F)[1:10],]
p4 <- ggplot(under.expressed.GO.1, aes(x =  reorder(ID, -p.adjust), y = -log10(p.adjust), fill = Count)) + 
  geom_bar(stat = "identity")+ #fill = pair.col[2]) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_label(aes(y = 0.8, label=GeneRatio),  nudge_y = -0.5, fill = "white", alpha = 0.7, size = 2, label.size = NA)+
  scale_fill_continuous(low = red_blue_20[1], high = red_blue_20[20])+
  scale_y_continuous(expand = expansion(add = c(0, .05))) +
  coord_flip() + theme_bw()+ theme(legend.position = "none")+  labs(x = "")
p4
pdf("01.TC-GM.underexpressed.GO.pdf", width = 15, height = 2.5, useDingbats = F)
p4
dev.off()


deg.average <- aggregate(TC.GM.DEG$logFC, list(TC.GM.DEG$gene), mean)
deg.rank = deg.average[,2]
names(deg.rank) <- deg.average[,1]
barplot(sort(deg.rank, decreasing = T))

fgsea.res.H.A <- fgsea(pathways = h.fgsea.db, stats = deg.rank, nPermSimple = 10000, minSize = 10)
fgsea.res.C5.A <- fgsea(pathways = c5.fgsea.db, stats = deg.rank, nPermSimple = 10000, minSize = 10)#Runs fgseaMultilevel
data.table::fwrite(fgsea.res.H.A, file = "01.TC-GM.fgsea-H.csv", sep=",")
data.table::fwrite(fgsea.res.C5.A, file = "01.TC-GM.fgsea-GO.csv", sep=",")

sum(fgsea.res.H.A$padj <0.05 & fgsea.res.H.A$NES>0, na.rm = T)
sum(fgsea.res.H.A$padj <0.05 & fgsea.res.H.A$NES<0, na.rm = T)

topUp <- fgsea.res.H.A %>% filter(NES > 0, padj <0.05) %>%  top_n(n =20, wt= NES)
topDown <- fgsea.res.H.A %>% filter(NES < 0, padj <0.05) %>% top_n(n = 20, wt= -NES)
topPathways <- bind_rows(topUp, topDown) %>% arrange(-NES)
dev.off()

pdf("01.TC-GM.fgsea-H.pdf", width = 8, height = 4, useDingbats = F)
plotGseaTable(h.fgsea.db[topPathways$pathway], 
              deg.rank, 
              fgsea.res.H.A, 
              gseaParam = 0.5)
dev.off()


sum(fgsea.res.C5.A$padj <0.05 & fgsea.res.C5.A$ES>0, na.rm = T)
sum(fgsea.res.C5.A$padj <0.05 & fgsea.res.C5.A$ES<0, na.rm = T)
#\
topUp <- fgsea.res.C5.A %>% filter(NES > 0, padj <0.05) %>%  top_n(n =20, wt= NES)
topDown <- fgsea.res.C5.A %>% filter(NES < 0, padj <0.05) %>% top_n(n = 20, wt= -NES)
topPathways <- bind_rows(topUp, topDown) %>% arrange(-NES)
topPathways$pathway

dev.off()
pdf("01.TC-GM.fgsea-GO.pdf", width = 14, height = 8, useDingbats = F)
plotGseaTable(c5.fgsea.db[topPathways$pathway], 
              deg.rank, 
              fgsea.res.C5.A, 
              gseaParam = 0.5, render = T)
dev.off()



#######################################################################################################################################
###2. IE - GM
IE.list <- grep("IE", unique(y.ann), value = T)
IE.list

IE.GM.DEG <- lapply(IE.list, function(id){
  sample.ID <- unlist(strsplit(id, "\\."))[2]
  message(id, ":", sample.ID)
  #qlf1 <- glmQLFTest(fit, contrast = makeContrasts(id-GM, levels=design))
  qlf1 <- exactTest(y, pair=c("GM", id))
  topTags(qlf1)
  deg <- qlf1$table
  deg$gene <- rownames(deg)
  deg$p.adjust = p.adjust(deg$PValue, method = "fdr")
  deg$sample = sample.ID
  write.csv(deg, file = paste0("02.", sample.ID,".IE-GM.DEG.csv"), row.names = F, quote = F)
  over.expressed <- rownames(deg)[deg$p.adjust<p.cutoff & deg$logFC > 0]; length(over.expressed)
  under.expressed <- rownames(deg)[deg$p.adjust<p.cutoff & deg$logFC < 0]; length(under.expressed)
  p.x <- deg$p.adjust[deg$p.adjust != 0]
  min.p <- min(p.x)
  deg$p.adjust[deg$p.adjust ==0] = min.p/10
  
  p <- ggplot(deg, aes(x = logFC, y = -log10(p.adjust), label = gene)) + 
    geom_point(size = 0.1, color = ifelse(deg$p.adjust <p.cutoff, "black", "gray")) + 
    geom_text_repel(data = deg[deg$gene %in% c(over.expressed, under.expressed),], 
                    min.segment.length = Inf, size = 1, seed = 123)+
    geom_hline(yintercept = -log10(p.cutoff), linetype = "dashed")+
    theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank())+
    scale_y_continuous(expand = expansion(mult = c(0, .1)))
  labs(title = paste0(sample.ID, ": IE v.s. GM"))
  
  pdf(paste0("02.", sample.ID,".IE-GM.DEG.pdf"), width = 4, height = 4, useDingbats = F)
  print(p)
  dev.off()
  deg$sample = sample.ID
  return(deg)
  
})
length(IE.GM.DEG)
IE.GM.DEG <- do.call(rbind, IE.GM.DEG)

ie.up.gene <- table(IE.GM.DEG$gene, IE.GM.DEG$p.adjust <p.cutoff & IE.GM.DEG$logFC> 0)
ie.down.gene <- table(IE.GM.DEG$gene, IE.GM.DEG$p.adjust <p.cutoff & IE.GM.DEG$logFC< 0)
ie.up.gene <- ie.up.gene[ie.up.gene[,2] >= 4,]
ie.down.gene <- ie.down.gene[ie.down.gene[,2] >= 4,]
nrow(ie.up.gene)
nrow(ie.down.gene)

length(intersect(rownames(ie.up.gene), rownames(up.gene)))
length(intersect(rownames(ie.down.gene), rownames(down.gene)))


ie.gene.ann <- data.frame(DEG = c(rep("up", nrow((ie.up.gene))), rep("down", nrow(ie.down.gene))))
rownames(ie.gene.ann) <- c(rownames(ie.up.gene), rownames(ie.down.gene))

roi.cluster <- hclust(as.dist(1 - cor(cpm[rownames(ie.gene.ann),], method = "spearman")), method = "average")
gene.cluster <- hclust(as.dist(1 - cor(t(cpm[rownames(ie.gene.ann),]), method = "spearman")), method = "average")


pdf("02.DEG.heatmap.pdf", width = 12, height = 12, useDingbats = F)
pheatmap(
  cpm[rownames(ie.gene.ann),],
  cluster_cols = roi.cluster,
  #cluster_cols = F,
  cluster_rows = gene.cluster,
  #cluster_rows = F,
  scale = "row",
  annotation_row = ie.gene.ann, 
  annotation_col = roi.ann,
  annotation_colors = annotation.col,
  color = red_blue_50,
  show_rownames = T, show_colnames = F)
dev.off()



write.table(ie.gene.ann, file = "02.IE-GM.all4SampleConsistentDEG.tsv", quote = F, col.names = F)

over.expressed.GO = data.frame(enricher(gene = rownames(ie.up.gene), TERM2GENE = c5, universe = row.names(cpm)))
over.expressed.H = data.frame(enricher(gene = rownames(ie.up.gene), TERM2GENE = h, universe = row.names(cpm)))

###
under.expressed.GO = data.frame(enricher(gene = rownames(ie.down.gene), TERM2GENE = c5, universe = row.names(cpm)))
under.expressed.H = data.frame(enricher(gene = rownames(ie.down.gene), TERM2GENE = h, universe = row.names(cpm)))


write.csv(over.expressed.GO, file = "02.IE-GM.overexpressed.GO.csv", row.names = F)
write.csv(over.expressed.H, file = "02.IE-GM.overexpressed.H.csv", row.names = F)
write.csv(under.expressed.GO, file = "02.IE-GM.underexpressed.GO.csv", row.names = F)
write.csv(over.expressed.H, file = "02.IE-GM.underexpressed.H.csv", row.names = F)

over.expressed.H$ID = gsub("^HALLMARK_", "", over.expressed.H$ID)
under.expressed.H$ID = gsub("^HALLMARK_", "", under.expressed.H$ID)
over.expressed.GO$ID = gsub("^GO\\w\\w\\_", "", over.expressed.GO$ID)
under.expressed.GO$ID = gsub("^GO\\w\\w\\_", "", under.expressed.GO$ID)
#RColorBrewer::display.brewer.all()
#reds <- RColorBrewer::brewer.pal(9, "Reds")

nrow(over.expressed.H)
p1 <- ggplot(over.expressed.H, aes(x =  reorder(ID, -p.adjust), y = -log10(p.adjust), fill = Count)) + 
  geom_bar(stat = "identity")+ #fill = pair.col[2]) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_label(aes(y = 1, label=GeneRatio),  nudge_y = -0.5, fill = "white", alpha = 0.7, size = 2, label.size = NA)+
  scale_fill_continuous(low = red_blue_20[1], high = red_blue_20[20])+
  scale_y_continuous(expand = expansion(add = c(0, .05))) +
  coord_flip() + theme_bw()+ theme(legend.position = "none")+  labs(x = "")
p1
#pdf("01.TC-GM.overexpressed.H.pdf", width = 7, height = 3, useDingbats = F)
#p1
#dev.off()

nrow(under.expressed.H)
p2 <- ggplot(under.expressed.H, aes(x =  reorder(ID, -p.adjust), y = -log10(p.adjust), fill = Count)) + 
  geom_bar(stat = "identity")+ #fill = pair.col[2]) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_label(aes(y = 1.1, label=GeneRatio),  nudge_y = -0.5, fill = "white", alpha = 0.7, size = 2, label.size = NA)+
  scale_fill_continuous(low = red_blue_20[1], high = red_blue_20[20])+
  scale_y_continuous(expand = expansion(add = c(0, .05))) +
  coord_flip() + theme_bw()+ theme(legend.position = "none")+  labs(x = "")
#pdf("01.TC-GM.underexpressed.H.pdf", width = 7, height = 3, useDingbats = F)
#p2
#dev.off()

nrow(over.expressed.GO)
#over.expressed.GO.1 <- over.expressed.GO[order(over.expressed.GO$p.adjust, decreasing = F)[1:20],]
p3 <- ggplot(over.expressed.GO, aes(x =  reorder(ID, -p.adjust), y = -log10(p.adjust), fill = Count)) + 
  geom_bar(stat = "identity")+ #fill = pair.col[2]) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_label(aes(y = 1.1, label=GeneRatio),  nudge_y = -0.5, fill = "white", alpha = 0.7, size = 2, label.size = NA)+
  scale_fill_continuous(low = red_blue_20[1], high = red_blue_20[20])+
  scale_y_continuous(expand = expansion(add = c(0, .05))) +
  coord_flip() + theme_bw()+ theme(legend.position = "none")+  labs(x = "")
p3

pdf("02.IE-GM.overexpressed.GO.pdf", width = 14, height = 1.2, useDingbats = F)
p3
dev.off()

nrow(under.expressed.GO)
under.expressed.GO.1 <- under.expressed.GO[order(under.expressed.GO$p.adjust, decreasing = F)[1:10],]
p4 <- ggplot(under.expressed.GO.1, aes(x =  reorder(ID, -p.adjust), y = -log10(p.adjust), fill = Count)) + 
  geom_bar(stat = "identity")+ #fill = pair.col[2]) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_label(aes(y = 0.8, label=GeneRatio),  nudge_y = -0.5, fill = "white", alpha = 0.7, size = 2, label.size = NA)+
  scale_fill_continuous(low = red_blue_20[1], high = red_blue_20[20])+
  scale_y_continuous(expand = expansion(add = c(0, .05))) +
  coord_flip() + theme_bw()+ theme(legend.position = "none")+  labs(x = "")
p4
pdf("02.IE-GM.underexpressed.GO.pdf", width = 15, height = 2.5, useDingbats = F)
p4
dev.off()


deg.average <- aggregate(IE.GM.DEG$logFC, list(IE.GM.DEG$gene), mean)
deg.rank = deg.average[,2]
names(deg.rank) <- deg.average[,1]
barplot(sort(deg.rank, decreasing = T))

fgsea.res.H.A <- fgsea(pathways = h.fgsea.db, stats = deg.rank, nPermSimple = 10000, minSize = 10)
fgsea.res.C5.A <- fgsea(pathways = c5.fgsea.db, stats = deg.rank, nPermSimple = 10000, minSize = 10)#Runs fgseaMultilevel
data.table::fwrite(fgsea.res.H.A, file = "02.IE-GM.fgsea-H.csv", sep=",")
data.table::fwrite(fgsea.res.C5.A, file = "02.IE-GM.fgsea-GO.csv", sep=",")

sum(fgsea.res.H.A$padj <0.05 & fgsea.res.H.A$NES>0, na.rm = T)
sum(fgsea.res.H.A$padj <0.05 & fgsea.res.H.A$NES<0, na.rm = T)

topUp <- fgsea.res.H.A %>% filter(NES > 0, padj <0.05) %>%  top_n(n =20, wt= NES)
topDown <- fgsea.res.H.A %>% filter(NES < 0, padj <0.05) %>% top_n(n = 20, wt= -NES)
topPathways <- bind_rows(topUp, topDown) %>% arrange(-NES)
dev.off()

pdf("02.IE-GM.fgsea-H.pdf", width = 8, height = 4, useDingbats = F)
plotGseaTable(h.fgsea.db[topPathways$pathway], 
              deg.rank, 
              fgsea.res.H.A, 
              gseaParam = 0.5)
dev.off()


sum(fgsea.res.C5.A$padj <0.05 & fgsea.res.C5.A$ES>0, na.rm = T)
sum(fgsea.res.C5.A$padj <0.05 & fgsea.res.C5.A$ES<0, na.rm = T)
#\
topUp <- fgsea.res.C5.A %>% filter(NES > 0, padj <0.05) %>%  top_n(n =20, wt= NES)
topDown <- fgsea.res.C5.A %>% filter(NES < 0, padj <0.05) %>% top_n(n = 20, wt= -NES)
topPathways <- bind_rows(topUp, topDown) %>% arrange(-NES)
topPathways$pathway

dev.off()
pdf("02.IE-GM.fgsea-GO.pdf", width = 14, height = 8, useDingbats = F)
plotGseaTable(c5.fgsea.db[topPathways$pathway], 
              deg.rank, 
              fgsea.res.C5.A, 
              gseaParam = 0.5, render = T)
dev.off()


#######################################################################################################################################
###3. TC - IE
xx <- as.matrix(table(y.sample.ID, y.ann.0))
xx
xx <- xx[which(xx[,3] >0 & xx[,4] >0),]
xx

TC.IE.DEG <- lapply(rownames(xx), function(sample.ID){
  #sample.ID <- unlist(strsplit(id, "\\."))[2]
  #message(id, ":", sample.ID)
  #qlf1 <- glmQLFTest(fit, contrast = makeContrasts(id-GM, levels=design))
  
  qlf1 <- exactTest(y, pair=c(paste0("IE.", sample.ID), paste0("TC.", sample.ID)))
  topTags(qlf1)
  deg <- qlf1$table
  deg$sample <- sample.ID
  deg$gene <- rownames(deg)
  deg$p.adjust = p.adjust(deg$PValue, method = "fdr")
  write.csv(deg, file = paste0("03.", sample.ID,".TC-IE.DEG.csv"), row.names = F, quote = F)
  over.expressed <- rownames(deg)[deg$p.adjust<p.cutoff & deg$logFC > 0]; length(over.expressed)
  under.expressed <- rownames(deg)[deg$p.adjust<p.cutoff & deg$logFC < 0]; length(under.expressed)
  p.x <- deg$p.adjust[deg$p.adjust != 0]
  min.p <- min(p.x)
  deg$p.adjust[deg$p.adjust ==0] = min.p/10
  
  p <- ggplot(deg, aes(x = logFC, y = -log10(p.adjust), label = gene)) + 
    geom_point(size = 0.1, color = ifelse(deg$p.adjust <p.cutoff, "black", "gray")) + 
    geom_text_repel(data = deg[deg$gene %in% c(over.expressed, under.expressed),], 
                    min.segment.length = Inf, size = 1, seed = 123)+
    geom_hline(yintercept = -log10(p.cutoff), linetype = "dashed")+
    theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank())+
    scale_y_continuous(expand = expansion(mult = c(0, .1)))
  labs(title = paste0(sample.ID, ": TC v.s. IE"))
  
  pdf(paste0("03.", sample.ID,".TC-IE.DEG.pdf"), width = 4, height = 4, useDingbats = F)
  print(p)
  dev.off()
  deg$sample = sample.ID
  return(deg)
  
})

length(TC.IE.DEG)
TC.IE.DEG <- do.call(rbind, TC.IE.DEG)

tc.ie.up.gene <- table(TC.IE.DEG$gene, TC.IE.DEG$p.adjust <p.cutoff & TC.IE.DEG$logFC> 0)
tc.ie.down.gene <- table(TC.IE.DEG$gene, TC.IE.DEG$p.adjust <p.cutoff & TC.IE.DEG$logFC< 0)
tc.ie.up.gene <- tc.ie.up.gene[tc.ie.up.gene[,2] >= 4,]
tc.ie.down.gene <- tc.ie.down.gene[tc.ie.down.gene[,2] >= 4,]
nrow(tc.ie.up.gene)
nrow(tc.ie.down.gene)


tc.ie.gene.ann <- data.frame(DEG = c(rep("up", nrow((tc.ie.up.gene))), rep("down", nrow(tc.ie.down.gene))))
rownames(tc.ie.gene.ann) <- c(rownames(tc.ie.up.gene), rownames(tc.ie.down.gene))

roi.cluster <- hclust(as.dist(1 - cor(cpm[rownames(tc.ie.gene.ann),], method = "spearman")), method = "average")
gene.cluster <- hclust(as.dist(1 - cor(t(cpm[rownames(tc.ie.gene.ann),]), method = "spearman")), method = "average")


pdf("03.DEG.heatmap.pdf", width = 12, height = 12, useDingbats = F)
pheatmap(
  cpm[rownames(tc.ie.gene.ann),],
  cluster_cols = roi.cluster,
  #cluster_cols = F,
  cluster_rows = gene.cluster,
  #cluster_rows = F,
  scale = "row",
  annotation_row = tc.ie.gene.ann, 
  annotation_col = roi.ann,
  annotation_colors = annotation.col,
  color = red_blue_50,
  show_rownames = T, show_colnames = F)
dev.off()



write.table(ie.gene.ann, file = "03.TC-IE.all4SampleConsistentDEG.tsv", quote = F, col.names = F)

over.expressed.GO = data.frame(enricher(gene = rownames(tc.ie.up.gene), TERM2GENE = c5, universe = row.names(cpm)))
over.expressed.H = data.frame(enricher(gene = rownames(tc.ie.up.gene), TERM2GENE = h, universe = row.names(cpm)))

###
under.expressed.GO = data.frame(enricher(gene = rownames(tc.ie.down.gene), TERM2GENE = c5, universe = row.names(cpm)))
under.expressed.H = data.frame(enricher(gene = rownames(tc.ie.down.gene), TERM2GENE = h, universe = row.names(cpm)))


write.csv(over.expressed.GO, file = "03.TC-IE.overexpressed.GO.csv", row.names = F)
write.csv(over.expressed.H, file = "03.TC-IE.overexpressed.H.csv", row.names = F)
write.csv(under.expressed.GO, file = "03.TC-IE.underexpressed.GO.csv", row.names = F)
write.csv(over.expressed.H, file = "03.TC-IE.underexpressed.H.csv", row.names = F)

over.expressed.H$ID = gsub("^HALLMARK_", "", over.expressed.H$ID)
under.expressed.H$ID = gsub("^HALLMARK_", "", under.expressed.H$ID)
over.expressed.GO$ID = gsub("^GO\\w\\w\\_", "", over.expressed.GO$ID)
under.expressed.GO$ID = gsub("^GO\\w\\w\\_", "", under.expressed.GO$ID)
#RColorBrewer::display.brewer.all()
#reds <- RColorBrewer::brewer.pal(9, "Reds")

nrow(over.expressed.H)
p1 <- ggplot(over.expressed.H, aes(x =  reorder(ID, -p.adjust), y = -log10(p.adjust), fill = Count)) + 
  geom_bar(stat = "identity")+ #fill = pair.col[2]) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_label(aes(y = 1, label=GeneRatio),  nudge_y = -0.5, fill = "white", alpha = 0.7, size = 2, label.size = NA)+
  scale_fill_continuous(low = red_blue_20[1], high = red_blue_20[20])+
  scale_y_continuous(expand = expansion(add = c(0, .05))) +
  coord_flip() + theme_bw()+ theme(legend.position = "none")+  labs(x = "")
p1
#pdf("03.TC-IE.overexpressed.H.pdf", width = 7, height = 3, useDingbats = F)
#p1
#dev.off()

nrow(under.expressed.H)
p2 <- ggplot(under.expressed.H, aes(x =  reorder(ID, -p.adjust), y = -log10(p.adjust), fill = Count)) + 
  geom_bar(stat = "identity")+ #fill = pair.col[2]) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_label(aes(y = 1.1, label=GeneRatio),  nudge_y = -0.5, fill = "white", alpha = 0.7, size = 2, label.size = NA)+
  scale_fill_continuous(low = red_blue_20[1], high = red_blue_20[20])+
  scale_y_continuous(expand = expansion(add = c(0, .05))) +
  coord_flip() + theme_bw()+ theme(legend.position = "none")+  labs(x = "")
#pdf("03.TC-IE.underexpressed.H.pdf", width = 7, height = 3, useDingbats = F)
#p2
#dev.off()

nrow(over.expressed.GO)
#over.expressed.GO.1 <- over.expressed.GO[order(over.expressed.GO$p.adjust, decreasing = F)[1:20],]
p3 <- ggplot(over.expressed.GO, aes(x =  reorder(ID, -p.adjust), y = -log10(p.adjust), fill = Count)) + 
  geom_bar(stat = "identity")+ #fill = pair.col[2]) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_label(aes(y = 1.1, label=GeneRatio),  nudge_y = -0.5, fill = "white", alpha = 0.7, size = 2, label.size = NA)+
  scale_fill_continuous(low = red_blue_20[1], high = red_blue_20[20])+
  scale_y_continuous(expand = expansion(add = c(0, .05))) +
  coord_flip() + theme_bw()+ theme(legend.position = "none")+  labs(x = "")
p3

pdf("03.TC-IE.overexpressed.GO.pdf", width = 14, height = 1.2, useDingbats = F)
p3
dev.off()

nrow(under.expressed.GO)
under.expressed.GO.1 <- under.expressed.GO[order(under.expressed.GO$p.adjust, decreasing = F)[1:10],]
p4 <- ggplot(under.expressed.GO.1, aes(x =  reorder(ID, -p.adjust), y = -log10(p.adjust), fill = Count)) + 
  geom_bar(stat = "identity")+ #fill = pair.col[2]) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_label(aes(y = 0.8, label=GeneRatio),  nudge_y = -0.5, fill = "white", alpha = 0.7, size = 2, label.size = NA)+
  scale_fill_continuous(low = red_blue_20[1], high = red_blue_20[20])+
  scale_y_continuous(expand = expansion(add = c(0, .05))) +
  coord_flip() + theme_bw()+ theme(legend.position = "none")+  labs(x = "")
p4
#pdf("03.TC-IE.underexpressed.GO.pdf", width = 15, height = 2.5, useDingbats = F)
#p4
#dev.off()


deg.average <- aggregate(TC.IE.DEG$logFC, list(TC.IE.DEG$gene), mean)
deg.rank = deg.average[,2]
names(deg.rank) <- deg.average[,1]
barplot(sort(deg.rank, decreasing = T))

fgsea.res.H.A <- fgsea(pathways = h.fgsea.db, stats = deg.rank, nPermSimple = 10000, minSize = 10)
fgsea.res.C5.A <- fgsea(pathways = c5.fgsea.db, stats = deg.rank, nPermSimple = 10000, minSize = 10)#Runs fgseaMultilevel
data.table::fwrite(fgsea.res.H.A, file = "03.TC-IE.fgsea-H.csv", sep=",")
data.table::fwrite(fgsea.res.C5.A, file = "03.TC-IE.fgsea-GO.csv", sep=",")

sum(fgsea.res.H.A$padj <0.05 & fgsea.res.H.A$NES>0, na.rm = T)
sum(fgsea.res.H.A$padj <0.05 & fgsea.res.H.A$NES<0, na.rm = T)

topUp <- fgsea.res.H.A %>% filter(NES > 0, padj <0.05) %>%  top_n(n =20, wt= NES)
topDown <- fgsea.res.H.A %>% filter(NES < 0, padj <0.05) %>% top_n(n = 20, wt= -NES)
topPathways <- bind_rows(topUp, topDown) %>% arrange(-NES)
dev.off()

pdf("03.TC-IE.fgsea-H.pdf", width = 8, height = 4, useDingbats = F)
plotGseaTable(h.fgsea.db[topPathways$pathway], 
              deg.rank, 
              fgsea.res.H.A, 
              gseaParam = 0.5)
dev.off()


sum(fgsea.res.C5.A$padj <0.05 & fgsea.res.C5.A$ES>0, na.rm = T)
sum(fgsea.res.C5.A$padj <0.05 & fgsea.res.C5.A$ES<0, na.rm = T)
#\
topUp <- fgsea.res.C5.A %>% filter(NES > 0, padj <0.05) %>%  top_n(n =20, wt= NES)
topDown <- fgsea.res.C5.A %>% filter(NES < 0, padj <0.05) %>% top_n(n = 20, wt= -NES)
topPathways <- bind_rows(topUp, topDown) %>% arrange(-NES)
topPathways$pathway

dev.off()
pdf("03.TC-IE.fgsea-GO.pdf", width = 14, height = 8, useDingbats = F)
plotGseaTable(c5.fgsea.db[topPathways$pathway], 
              deg.rank, 
              fgsea.res.C5.A, 
              gseaParam = 0.5, render = T)
dev.off()










####
gene.ann <- data.frame(gene = unique(c(rownames(up.gene), rownames(down.gene), 
                                       rownames(ie.up.gene), rownames(ie.down.gene),
                                       rownames(tc.ie.up.gene), rownames(tc.ie.down.gene))))

gene.ann$IE <- ifelse(gene.ann$gene %in% rownames(ie.up.gene), "up", 
                      ifelse(gene.ann$gene %in% rownames(ie.down.gene), "down", "non-DEG"))
gene.ann$TCvsIE <- ifelse(gene.ann$gene %in% rownames(tc.ie.up.gene), "up", 
                      ifelse(gene.ann$gene %in% rownames(tc.ie.down.gene), "down", "non-DEG"))
gene.ann$TC <- ifelse(gene.ann$gene %in% rownames(up.gene), "up", 
                      ifelse(gene.ann$gene %in% rownames(down.gene), "down", "non-DEG"))
rownames(gene.ann) <- gene.ann$gene

roi.ann <- data.frame(location = meta$annotation1[match(colnames(cpm), meta$ID)] )
rownames(roi.ann) <- colnames(cpm)
roi.ann$sample <- sapply(rownames(roi.ann), function(x) unlist(strsplit(x, "\\-"))[1])
roi.ann$grade <- tumor.grade[roi.ann$sample]
table(roi.ann)

deg.col <- c(pair.col[c(2,6)], "#808080")
names(deg.col) <- c("down", "up", "non-DEG")
grade.col <- c("#888888", "#000000")
names(grade.col) <- c("low", "high")
annotation.col <- list(location = loc1.col, 
                       grade = grade.col,
                       sample = sample.colors,
                       TC = deg.col,
                       IE = deg.col,
                       TCvsIE = deg.col)

roi.cluster <- hclust(as.dist(1 - cor(cpm[rownames(gene.ann),], method = "spearman")), method = "average")
gene.cluster <- hclust(as.dist(1 - cor(t(cpm[rownames(gene.ann),]), method = "spearman")), method = "average")

pdf("04.all.DEGs.heatmap.pdf", width = 10, height = 12, useDingbats = F)
pheatmap(
  cpm[rownames(gene.ann),],
  cluster_cols = roi.cluster,
  cluster_rows = gene.cluster,
  scale = "row",
  annotation_row = gene.ann[,2:4], 
  annotation_col = roi.ann,
  annotation_colors = annotation.col,
  color = red_blue_50,
  fontsize =6,
  show_rownames = T, show_colnames = F)
dev.off()

write.csv(gene.ann, file = "04.all.DEGs.csv", quote = F)
### Venn diagram?
#library(ggvenn)
library(eulerr)
deg.list.over <- list("TC-normal" = rownames(gene.ann)[gene.ann$TC == "up"], 
                      "IE-normal" = rownames(gene.ann)[gene.ann$IE == "up"], 
                      "TC-IE" = rownames(gene.ann)[gene.ann$TCvsIE == "up"])
#ggvenn(deg.list.over)
p1 <- plot(euler(deg.list.over, shape = "ellipse"), quantities = TRUE, fills = pair.col[c(1,3,5)])

deg.list.under <- list("TC-normal" = rownames(gene.ann)[gene.ann$TC == "down"], 
                       "IE-normal" = rownames(gene.ann)[gene.ann$IE == "down"], 
                       "TC-IE" = rownames(gene.ann)[gene.ann$TCvsIE == "down"])
p2<-plot(euler(deg.list.under, shape = "ellipse"), quantities = TRUE, fills = pair.col[c(1,3,5)])

pdf("04.all.DEGs.veen.pdf", width = 10, height = 4)
cowplot::plot_grid(p1, p2)
dev.off()


#######################################################################################################################################

###5. Tumor core of low grade v.s. high grade
###6. Tumor infiltration edge of low grade v.s. high grade


