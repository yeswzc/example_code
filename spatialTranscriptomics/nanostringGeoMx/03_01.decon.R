library(SpatialDecon)
rm(list = ls())
setwd("/Users/wuz6/Documents/Project/08.spatialTranscriptome/01.RNA/03.deconvolution/")


load("../../00.data/signatures/nature2016-science2017.IDHmut.rda")
ls()
head(cell.type.mean)
markers1 <- data.frame(markers1)
head(markers1)

ggplot(markers1, aes(x = avg_log2FC))+
  geom_density()+
  facet_wrap(~cluster)+
  theme_bw()

table(markers1$avg_log2FC > 6, markers1$cluster)
sum(duplicated(markers1$gene))

markers1[markers1$gene %in% markers1$gene[duplicated(markers1$gene)],]


markers1 <- markers1[! markers1$gene %in% markers1$gene[duplicated(markers1$gene)],]
nrow(markers1)
rownames(cell.type.mean) <- gsub("_", "-", rownames(cell.type.mean))
idh.signature <- cell.type.mean[rownames(cell.type.mean) %in% markers1$gene, ]
nrow(idh.signature)
nrow(markers1)

idh.signature.ann <- data.frame(cell.type = markers1$cluster)
rownames(idh.signature.ann) <- markers1$gene

pdf("../../00.data/signatures/01.IDHmut.signature.zc.pdf", width = 7, height = 5, useDingbats = F)
pheatmap(idh.signature, scale = "row",
         show_rownames = F,
         annotation_row = idh.signature.ann)
dev.off()

rm(markers, markers1, cell.type.mean)

###decon
load("../01.QCfirst/CTA.allGlioma.seuratObj.rda")
head(seurat@meta.data)
seurat
seurat <- subset(seurat, tumor %in% c("IDH_A", "IDH_O"))
seurat

per.observation.mean.neg = seurat@assays$RNA@data["Negative Probe", ]
bg = data.matrix(sweep(seurat@assays$RNA@data * 0, 2, per.observation.mean.neg, "+"))
dim(bg)

###
#idh.signature <- idh.signature[rownames(idh.signature) %in% rownames(seurat), ]

head(idh.signature)
colnames(idh.signature) <- gsub(" ", ".", colnames(idh.signature))
colnames(idh.signature) <- gsub("/", ".", colnames(idh.signature))

res1 <- spatialdecon(norm = data.matrix(seurat@assays$RNA@data), 
                    bg = bg,
                    X = idh.signature)

res1$prop_of_all[1:4,1:4]

### 2 signatures2
load("../../00.data/signatures/01.idh.mtx.spatialDeconMtx.rda")
res2 <- spatialdecon(norm = data.matrix(seurat@assays$RNA@data), 
                     bg = bg,
                     X = idh.mtx)
res2$prop_of_all[1:4,1:4]


pdf("01.IDH.compare.2signature.pdf", width = 10, height = 10, useDingbats = F)
par(mfrow = c(2,2))
plot(res1$prop_of_all[1,], res2$prop_of_all[1,],main = colnames(idh.signature)[1], xlab = "mean top200", ylab = "spatialDecon")
plot(res1$prop_of_all[2,], res2$prop_of_all[2,],main = colnames(idh.signature)[2], xlab = "mean top200", ylab = "spatialDecon")
plot(res1$prop_of_all[3,], res2$prop_of_all[3,],main = colnames(idh.signature)[3], xlab = "mean top200", ylab = "spatialDecon")
plot(res1$prop_of_all[4,], res2$prop_of_all[4,],main = colnames(idh.signature)[4], xlab = "mean top200", ylab = "spatialDecon")
dev.off()

save(res1, file = "01.IDH.topSigture.spatialDecon.res.rda")
save(res2, file = "01.IDH.spatialDeconSignature.spatialDecon.res.rda")

x <- data.frame(t(res1$prop_of_nontumor))
head(x[,1:4])
write.csv(x, file = "01.IDH.decon.csv", quote = F)
head(seurat)
identical(rownames(x), rownames(seurat@meta.data))

x$location <- seurat@meta.data$location

p1 <- ggplot(x, aes(x = location, y = malignant)) +
  geom_boxplot()+
  geom_jitter()+
  theme_bw()
p2 <- ggplot(x, aes(x = location, y = oligodendrocytes)) +
  geom_boxplot()+
  geom_jitter()+
  theme_bw()
p3<- ggplot(x, aes(x = location, y = microglia.macrophage)) +
  geom_boxplot()+
  geom_jitter()+
  theme_bw()
p4 <- ggplot(x, aes(x = location, y = T.cells)) +
  geom_boxplot()+
  geom_jitter()+
  theme_bw()

pdf("01.IDH.decon.pdf", width = 12, height = 10, useDingbats = F)
cowplot::plot_grid(p1,p2,p3,p4, nrow = 2)
dev.off()



####GBM
load("../../00.data/signatures/Neftel.Cell.GBM.rda")
sum(duplicated(markers1$gene))

head(cell.type.mean)
rownames(cell.type.mean) <- gsub("_", "-", rownames(cell.type.mean))

gbm.signature <- cell.type.mean[markers1$gene,]
gbm.signature.ann <- data.frame(type = markers1$cluster)
rownames(gbm.signature.ann) <- markers1$gene

pdf("../../00.data/signatures/Neftel.Cell.GBM.signature.pdf", width = 7, height = 5, useDingbats = F)
pheatmap(gbm.signature, scale = "row", 
         annotation_row = gbm.signature.ann,
         show_rownames = F)
dev.off()

###Cell reports
load("../../00.data/signatures/cellReports2017.GBM.rda")
sum(duplicated(markers1$gene))
markers1 <- markers1[! markers1$gene %in% markers1$gene[duplicated(markers1$gene)],]
head(cell.type.mean)
rownames(cell.type.mean) <- gsub("_", "-", rownames(cell.type.mean))
markers1$gene[!markers1$gene %in% rownames(cell.type.mean)]
markers1 <- markers1[markers1$gene %in% rownames(cell.type.mean), ]

gbm.signature.cellRep <- cell.type.mean[markers1$gene,]
gbm.signature.cellRep.ann <- data.frame(type = markers1$cluster)
rownames(gbm.signature.cellRep.ann) <- markers1$gene

pdf("../../00.data/signatures/cellReports2017.GBM.pdf", width = 7, height = 5, useDingbats = F)
pheatmap(gbm.signature.cellRep, scale = "row", 
         cluster_rows = hclust(dist(gbm.signature.cellRep)),
         cluster_cols = hclust(dist(t(gbm.signature.cellRep))),
         annotation_row = gbm.signature.cellRep.ann,
         show_rownames = F)
dev.off()

###
load("../../00.data/signatures/KelvinNG2021.rda")
sum(! markers1$gene %in% rownames(cell.type.mean))
sum(!rownames(cell.type.mean) %in% markers$gene)
cell.type.mean <- cell.type.mean[rownames(cell.type.mean) %in% markers1$gene, ]

sum(duplicated(markers1$gene))
nrow(markers1)
markers1 <- markers1[! markers1$gene %in% markers1$gene[duplicated(markers1$gene)],]
sum(!markers1$gene %in% rownames(cell.type.mean))

head(markers1)
head(cell.type.mean)
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
GG <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                      "hgnc_symbol"),
            values=rownames(markers),mart= mart)
markers1$ens <- markers1$gene
markers1$gene <- GG$hgnc_symbol[match(markers1$ens, GG$ensembl_gene_id)]
markers1 <- markers1[!is.na(markers1$gene),]
markers1 <- markers1[! markers1$gene %in% markers1$gene[duplicated(markers1$gene)],]

#rownames(cell.type.mean) <- gsub("_", "-", rownames(cell.type.mean))
cell.type.mean <- cell.type.mean[rownames(cell.type.mean) %in% GG$ensembl_gene_id & 
                                   rownames(cell.type.mean) %in% markers1$ens, ]
nrow(cell.type.mean)
rownames(cell.type.mean) <- GG$hgnc_symbol[match(rownames(cell.type.mean), GG$ensembl_gene_id)]
head(cell.type.mean)

sum(!markers1$gene %in% rownames(cell.type.mean))
sum(!rownames(cell.type.mean) %in% markers1$gene)
sum(duplicated(rownames(cell.type.mean)))
sum(duplicated(markers1$gene))
length(intersect(rownames(cell.type.mean), markers1$gene))

gbm.signature.kelvin <- cell.type.mean[markers1$gene,]
gbm.signature.kelvin.ann <- data.frame(type = markers1$cluster)
rownames(gbm.signature.kelvin.ann) <- markers1$gene

pdf("../../00.data/signatures/KelvinNG2021.pdf", width = 7, height = 5, useDingbats = F)
pheatmap(gbm.signature.kelvin, scale = "row", 
         annotation_row = gbm.signature.kelvin.ann,
         show_rownames = F)
dev.off()

###
load("../01.QCfirst/CTA.allGlioma.seuratObj.rda")

seurat <- subset(seurat, tumor == "GBM")
seurat

per.observation.mean.neg = seurat@assays$RNA@data["Negative Probe", ]
bg = data.matrix(sweep(seurat@assays$RNA@data * 0, 2, per.observation.mean.neg, "+"))
dim(bg)

res1 <- spatialdecon(norm = data.matrix(seurat@assays$RNA@data), 
                     bg = bg,
                     X = gbm.signature)

x <- data.frame(t(res$prop_of_all))
head(x)
identical(rownames(x), rownames(seurat@meta.data))
x$location <- seurat@meta.data$location

head(x)
x1 <- reshape2::melt(x, id.vars = c("location"))
head(x1)

pdf("02.GBM.decon.pdf", width = 12, height = 10, useDingbats = F)
ggplot(x1, aes(x = location, y = value))+
  geom_boxplot()+
  geom_jitter()+
  facet_wrap(~variable)+
  theme_bw()
dev.off()
save(res, file = "02.GBM.decon.rda")
write.csv(x, file = "02.GBM.decon.cellRatio.csv", quote = F)

###
res3 <- spatialdecon(norm = data.matrix(seurat@assays$RNA@data), 
                     bg = bg,
                     X = gbm.signature.cellRep)

x <- data.frame(t(res3$prop_of_all))
head(x)
identical(rownames(x), rownames(seurat@meta.data))
x$location <- seurat@meta.data$location

head(x)
x1 <- reshape2::melt(x, id.vars = c("location"))
head(x1)

pdf("02.GBM.decon.cellRep.pdf", width = 12, height = 10, useDingbats = F)
ggplot(x1, aes(x = location, y = value))+
  geom_boxplot()+
  geom_jitter()+
  facet_wrap(~variable)+
  theme_bw()
dev.off()

save(res, file = "02.GBM.decon.cellRep.rda")
write.csv(x, file = "02.GBM.decon.cellRep.cellRatio.csv", quote = F)

head(x)
plot(density(rowSums(x[,c(1,2,4)])))


###
res2 <- spatialdecon(norm = data.matrix(seurat@assays$RNA@data), 
                     bg = bg,
                     X = gbm.signature.kelvin)

x <- data.frame(t(res2$prop_of_all))
head(x)
identical(rownames(x), rownames(seurat@meta.data))
x$location <- seurat@meta.data$location

head(x)
x1 <- reshape2::melt(x, id.vars = c("location"))
head(x1)

pdf("02.GBM.decon.kelvinSig.pdf", width = 12, height = 10, useDingbats = F)
ggplot(x1, aes(x = location, y = value))+
  geom_boxplot()+
  geom_jitter()+
  facet_wrap(~variable)+
  theme_bw()
dev.off()

save(res, file = "02.GBM.decon.kelvinSig.rda")
write.csv(x, file = "02.GBM.decon.kelvinSig.cellRatio.csv", quote = F)



#####################################################################
####WTA
load("../01.QCfirst/WTA/WTA.seuratObj.rda")
head(seurat[[]])
unique(seurat@meta.data$sample)
unique(seurat@meta.data$tumor)

wta.gbm <- subset(seurat, subset = tumor == "GBM")
wta.gbm
grep("Neg", rownames(wta.gbm), value = T)

per.observation.mean.neg = wta.gbm@assays$RNA@data["NegProbe-WTX", ]
bg = data.matrix(sweep(wta.gbm@assays$RNA@data * 0, 2, per.observation.mean.neg, "+"))

res <- spatialdecon(norm = data.matrix(wta.gbm@assays$RNA@data), 
                    bg = bg,
                    X = gbm.signature)

x <- data.frame(t(res$prop_of_all))
head(x)
identical(rownames(x), rownames(wta.gbm@meta.data))
x$location <- wta.gbm@meta.data$location

head(x)
x1 <- reshape2::melt(x, id.vars = c("location"))
head(x1)

pdf("03.wta.GBM.decon.pdf", width = 12, height = 10, useDingbats = F)
ggplot(x1, aes(x = location, y = value))+
  geom_boxplot()+
  geom_jitter()+
  facet_wrap(~variable)+
  theme_bw()
dev.off()

save(res, file = "03.wta.GBM.decon.rda")
write.csv(x, file = "03.wta.GBM.decon.cellRatio.csv", quote = F)

###
unique(seurat@meta.data$tumor)
wta.idh <- subset(seurat, subset = tumor == "IDH_A")
wta.idh
grep("Neg", rownames(wta.idh), value = T)

per.observation.mean.neg = wta.idh@assays$RNA@data["NegProbe-WTX", ]
bg = data.matrix(sweep(wta.idh@assays$RNA@data * 0, 2, per.observation.mean.neg, "+"))

res <- spatialdecon(norm = data.matrix(wta.idh@assays$RNA@data), 
                    bg = bg,
                    X = gbm.signature)

x <- data.frame(t(res$prop_of_all))
head(x)
identical(rownames(x), rownames(wta.idh@meta.data))
x$location <- wta.idh@meta.data$location

head(x)
x1 <- reshape2::melt(x, id.vars = c("location"))
head(x1)

pdf("04.wta.IDH.decon.pdf", width = 12, height = 10, useDingbats = F)
ggplot(x1, aes(x = location, y = value))+
  geom_boxplot()+
  geom_jitter()+
  facet_wrap(~variable)+
  theme_bw()
dev.off()

save(res, file = "04.wta.IDH.decon.rda")
write.csv(x, file = "04.wta.IDH.decon.cellRatio.csv", quote = F)



