library(Seurat)
library(ggplot2)
library(clustree)
library(paletteer)
library(ComplexHeatmap)
library(dplyr)
library(clusterProfiler)
rm(list=ls())
setwd("/Users/wuz6/Documents/Project/lij36/")

RColorBrewer::display.brewer.pal(12, "Paired")
sample.color <- RColorBrewer::brewer.pal(12, "Paired")[c(2,4,8,10)]
names(sample.color) <- c("11H", "6H", "11L", "6L")

#BiocManager::install("monocle")



###
if(0){ #runed
##############################################################################
##############################################################################
all.samples <- list.files("data/")
all.samples
#all.samples <- grep("png|VDJ", all.samples, invert = T, value = T)
all.samples <- grep("vdj", all.samples, invert = T, value = T)
all.samples

###merge single cell RNA
sc.data.list = lapply(all.samples, function(id){
  cat("Reading", id, "\n")
  sc.data = Read10X(data.dir = file.path("data", id, "filtered_feature_bc_matrix/"))
  sc <- CreateSeuratObject(counts = sc.data, min.cells = 3, min.features = 500, project = id); rm(sc.data);
  #sc
  sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^mt-")
  sc[["sampleID"]] = id;
  #sc <- add_clonotype(paste0(id, "-VDJ/"), sc, "T")
  sc <- subset(sc, subset = nFeature_RNA > 500 & percent.mt < 20)
  sc <- RenameCells(sc, paste0(id, "_", colnames(sc)))
  sc <- SCTransform(sc)
  sc <- RunPCA(sc)
  return(sc);
})
colnames(sc.data.list[[1]])[1:4]

#doublet
#
#
##
cat("Combine all samples\n")
all.samples
names(sc.data.list) = all.samples
scRNA = merge(sc.data.list[[1]], y = sc.data.list[-1]) #,add.cell.ids = all.samples)

p00 <- VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p01 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 1)

pdf("01.scQC.pdf", width = 7, height = 5, useDingbats = F)
p00;
p01
dev.off()


scRNA <- NormalizeData(scRNA, assay = "RNA")
scRNA <- ScaleData(scRNA, features = rownames(scRNA), assay = "RNA")
scRNA <- FindVariableFeatures(scRNA, assay = "RNA")

x = data.frame(table(scRNA@meta.data$orig.ident))

p0 <- ggplot(x, aes(x = Var1, y = Freq, fill = Var1))+
  geom_bar(stat = "identity")+
  geom_text(aes(x=Var1, y = 3000, label = Freq))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  #coord_polar()+
  scale_fill_manual(values = sample.color)+
  scale_x_discrete(limits = c("6L", "6H", "11L", "11H"))+
  labs(x="", fill = "", y = "count")
p0

scRNA <- RunPCA(scRNA, assay = "RNA")
ElbowPlot(scRNA, ndims = 40)
scRNA <- RunUMAP(scRNA, assay = "RNA", dims = 1:20)

head(scRNA@reductions$umap@cell.embeddings)

p1 = DimPlot(scRNA, reduction = "umap") + scale_color_manual(values = sample.color) + 
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2", title = "Merge")
p1



###harmony
library(harmony)
scRNA <- RunHarmony(scRNA, "orig.ident")
scRNA <- RunUMAP(scRNA, assay = "RNA", dims = 1:20, reduction = "harmony", reduction.name = "harmonyUMAP")

p2 = DimPlot(scRNA, reduction = "harmonyUMAP") + scale_color_manual(values = sample.color) + 
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2")+
  labs(title = "Harmony")
p2


#Please specify the anchor.features to be used. The expected workflow for integratinge assays produced by SCTransform is
#SelectIntegrationFeatures -> PrepSCTIntegration -> FindIntegrationAnchors.
sel.features = SelectIntegrationFeatures(object.list = sc.data.list)
sc.data.list = PrepSCTIntegration(sc.data.list, anchor.features = sel.features)
data.anchors = FindIntegrationAnchors(sc.data.list, normalization.method = "SCT", anchor.features = sel.features)
data.integrated <- IntegrateData(data.anchors, normalization.method = "SCT")

data.integrated <- SCTransform(data.integrated)
data.integrated <- RunPCA(data.integrated)
data.integrated <- RunUMAP(data.integrated, dims = 1:20, reduction.name = "cca")

head(data.integrated@reductions$cca)
p3 = DimPlot(data.integrated, reduction = "cca") + scale_color_manual(values = sample.color) + 
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2", title = "Seurat integration")
p3

pdf("01.umap_sample_batch.pdf", width = 15, height = 15, useDingbats = F)
cowplot::plot_grid(p0, p1, p2, p3)
dev.off()
#data.integrated@meta.data = data.integrated@meta.data[,1:7]

scRNA@reductions$integrationUMAP = data.integrated@reductions$cca
rownames(scRNA@reductions$integrationUMAP@cell.embeddings) = rownames(scRNA@reductions$umap@cell.embeddings)
tail(scRNA[[]])
###
######################################################################
###add clonotype
tcr <- read.csv("data/vdj_t/filtered_contig_annotations.csv")
tail(tcr$barcode)
# Remove the -1 at the end of each barcode.
# Subsets so only the first line of each barcode is kept,
# as each entry for given barcode will have same clonotype.
tcr <- tcr[!duplicated(tcr$barcode), ]
head(tcr)
unique(tcr$donor)
unique(tcr$origin)
#
tcr$barcode <- gsub("\\-\\d$", "-1", tcr$barcode)
tcr$barcode2 <-paste0(tcr$origin, "_", tcr$barcode, "_", tcr$barcode)
tail(tcr$barcode2)

# Only keep the barcode and clonotype columns.
# We'll get additional clonotype info from the clonotype table.
tcr <- tcr[,c("barcode2", "raw_clonotype_id")]
names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
# Clonotype-centric info.
clono <- read.csv("data/vdj_t/clonotypes.csv")
# Slap the AA sequences onto our original table by clonotype_id.
tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
names(tcr)[names(tcr) == "cdr3s_aa"] <- "cdr3s_aa"
# Reorder so barcodes are first column and set them as rownames.
tcr <- tcr[, c(2,1,3)]
head(tcr)
tail(tcr)
rownames(tcr) <- tcr[,1]
tcr[,1] <- NULL
colnames(tcr) <- paste("T", colnames(tcr), sep="_")
sum(!rownames(tcr) %in% colnames(scRNA)) #993
sum(!colnames(scRNA) %in% rownames(tcr)) #6792
rownames(tcr)[!rownames(tcr) %in% colnames(scRNA)][1:4]
###sc filtered that though have TCR, but the cell is not valid

# Add to the Seurat object's metadata.
scRNA <- AddMetaData(object=scRNA, metadata=tcr)
#data.integrated <- AddMetaData(data.integrated, metadata = tcr)
save(scRNA, data.integrated, file = "01.seurat.obj.rda")
#table(scRNA@meta.data$orig.ident)


}
load("01.seurat.obj.rda")
###clonotype
head(scRNA@meta.data)

DimPlot(scRNA, reduction = "umap", group.by = "T_clonotype_id") + 
  #scale_color_manual(values = sample.color) + 
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2", title = "Merge")


keep.clono <- names(which(table(scRNA@meta.data$T_clonotype_id) > 10))
keep.clono

clono.col = as.character(paletteer_d("ggthemes::Tableau_20"))[1:length(keep.clono)]
names(clono.col) = keep.clono
clono.col

scRNA <- subset(scRNA, subset = T_clonotype_id %in% keep.clono)

xx = cbind(scRNA@reductions$umap@cell.embeddings, 
           scRNA@reductions$harmonyUMAP@cell.embeddings,
           scRNA@reductions$integrationUMAP@cell.embeddings,
           scRNA@meta.data)

colnames(xx)
xx$sampleID <- factor(xx$sampleID, levels = c("6L", "6H", "11L", "11H"))
xx$T_clonotype_id <- factor(xx$T_clonotype_id, levels = paste0("clonotype", 1:length(keep.clono)))



###
p1 = ggplot(xx, aes(x =umap_1, y = umap_2, color = T_clonotype_id)) +
  geom_point(size = .5)+
  scale_colour_manual(values = clono.col)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2", title = "")+
  labs(title = "Merge")
p1
p2 = ggplot(xx, aes(x =harmonyUMAP_1, y = harmonyUMAP_2, color = T_clonotype_id)) +
  geom_point(size = .5)+
  scale_colour_manual(values = clono.col)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2", title = "")+
  labs(title = "Harmony")
p3 = ggplot(xx, aes(x =integrationUMAP_1, y = integrationUMAP_2, color = T_clonotype_id)) +
  geom_point(size = .5)+
  scale_colour_manual(values = clono.col)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2", title = "")+
  labs(title = "Seurat integration")

p1
p2
p3

pdf("01.umap.clonotype.pdf", width = 10, height = 10, useDingbats = F)
cowplot::plot_grid(p1, p2, p3, nrow = 2)
p1 + facet_wrap(~sampleID)
p2 + facet_wrap(~sampleID)
p3 + facet_wrap(~sampleID)
dev.off()
#scRNA <- data.integrated
####scR
scRNA <- FindNeighbors(scRNA, dims = 1:20, assay = "SCT")
scRNA <- FindClusters(scRNA, resolution = 0.5, graph.name = "RNA_snn")#seq(0.4, 1.2, 0.1))


xx = cbind(scRNA@reductions$umap@cell.embeddings, 
           scRNA@reductions$harmonyUMAP@cell.embeddings,
           scRNA@reductions$integrationUMAP@cell.embeddings,
           scRNA@meta.data)

colnames(xx)
xx$sampleID <- factor(xx$sampleID, levels = c("6L", "6H", "11L", "11H"))
xx$T_clonotype_id <- factor(xx$T_clonotype_id, levels = paste0("clonotype", 1:length(keep.clono)))



p1 = ggplot(xx, aes(x =umap_1, y = umap_2, color = seurat_clusters)) +
  geom_point(size = .5)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2", title = "")+
  labs(title = "Merge")

p2 = ggplot(xx, aes(x =harmonyUMAP_1, y = harmonyUMAP_2, color = seurat_clusters)) +
  geom_point(size = .5)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2", title = "")+
  labs(title = "Harmony")


p3 = ggplot(xx, aes(x =integrationUMAP_1, y = integrationUMAP_2, color = seurat_clusters)) +
  geom_point(size = .5)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2", title = "")+
  labs(title = "Seurat integration")

cowplot::plot_grid(p1,p2,p3)


#save for scVelo
#library(SeuratDisk)
#scRNA[["RNA"]] <- as(object = scRNA[["RNA"]], Class = "Assay")
#SaveH5Seurat(scRNA, filename = "01.scTCR2024.h5Seurat", overwrite = T)
#Convert("01.scTCR2024.h5Seurat", dest = "h5ad", overwrite = T)

############################################################################################################
y <- table(scRNA@meta.data$T_clonotype_id, scRNA@meta.data$orig.ident)
y
y1 <- rowSums(y)
library(vegan)
y[1:4,1:4]
###TCR diversity of each sample
#Shannon diversity 
d <- data.frame(apply(y, 2, diversity))
d$x <- rownames(d)
colnames(d) <- c("diversity", "sample")
head(d)



p <- ggplot(d, aes(x = reorder(sample, -diversity), y = diversity, fill = sample))+
  geom_histogram(stat="identity")+
  theme_bw()+
  scale_fill_manual(values = sample.color)+
  theme(aspect.ratio = 1, panel.grid = element_blank())+
  labs(x = "")+
  coord_cartesian(expand = F)
p
pdf("01.scRNA.colontypeDiversity.pdf", width = 5, height = 4, useDingbats = F)
p
dev.off()

#      11H       11L        6H        6L
#0.5217312 1.0433631 0.5427320 1.6155731
y <- y/y1
y <- reshape2::melt(y)
head(y)
y$N <- y1[match(y$Var1, names(y1))]

p1 <- ggplot(y, aes(x = reorder(Var1, -N), y = value, fill = Var2))+
  geom_bar(stat = "identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = sample.color)+
  theme(panel.grid = element_blank())+
  labs(x = "", y = "ratio", fill = "")
  #scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1.0))
p1

library(ggbreak)
y2 = data.frame(table(scRNA@meta.data$T_clonotype_id))

p2 = ggplot(y2, aes(x = reorder(Var1, -Freq), y = Freq))+
  geom_histogram(stat = "identity", fill = "gray", color = NA)+
  #scale_y_break(c(2600, 8000))+
  #scale_y_break(c(670, 2000))+
  theme_bw()+ theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  geom_text(inherit.aes = F, data = y, aes(x = Var1, y = 100, label =N, angle = 90))+
  labs(x= "", y = "count")+
  scale_y_log10(expand = c(0,0))
p2

egg::ggarrange(p2, p1, ncol = 1, nrow = 2, heights = c(5, 10))


y3 = data.frame(table(scRNA@meta.data$orig.ident, scRNA@meta.data$T_clonotype_id))
y3$Var2 = factor(y3$Var2, levels = paste0("clonotype", 1:length(keep.clono)))

p3 = ggplot(y3, aes(x = Var2, y = Freq, fill = Var2))+
  geom_histogram(stat = "identity")+
  scale_fill_manual(values = clono.col)+
  theme_bw()+ 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 12))+
  labs(x= "", y = "count")+
  scale_y_log10(expand = c(0,0))+
  facet_wrap(~Var1, ncol = 1, strip.position = "right")

p3

pdf("01.clontype.sample.freq.pdf", width = 8, height = 8, useDingbats = F)
egg::ggarrange(p2, p1, ncol = 1, nrow = 2, heights = c(5, 40))
dev.off()

colnames(y) <- c("clono", "sample", "ratio", "N")
write.csv(y, file = "01.clontype.sample.freq.csv", quote = F, row.names = F)
###
###

head(scRNA[[]])
table(scRNA@meta.data$T_clonotype_id)
plot(density(table(scRNA@meta.data$T_clonotype_id)))
table(scRNA@meta.data$T_clonotype_id, scRNA@meta.data$orig.ident)
}




###
xx = cbind(scRNA@reductions$umap@cell.embeddings, 
           scRNA@reductions$harmonyUMAP@cell.embeddings,
           scRNA@reductions$integrationUMAP@cell.embeddings,
           scRNA@meta.data)

colnames(xx)
xx$sampleID <- factor(xx$sampleID, levels = c("6L", "6H", "11L", "11H"))
xx$T_clonotype_id <- factor(xx$T_clonotype_id, levels = paste0("clonotype", 1:length(keep.clono)))

colnames(xx)
p1 = ggplot(xx, aes(x =umap_1, y = umap_2, color = sampleID)) +
  geom_point(size = .2)+
  scale_color_manual(values = sample.color)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2", title = "")+
  labs(title = "Merge")
p2 = ggplot(xx, aes(x =harmonyUMAP_1, y = harmonyUMAP_2, color = sampleID)) +
  geom_point(size = .2)+
  scale_color_manual(values = sample.color)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2", title = "")+
  labs(title = "harmony")
p3 = ggplot(xx, aes(x =integrationUMAP_1, y = integrationUMAP_2, color = sampleID)) +
  geom_point(size = .2)+
  scale_color_manual(values = sample.color)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2", title = "")+
  labs(title = "Seurat integration")

pdf("01.umap-clonotype.pdf", width = 12, height = 12, useDingbats = F)
p1 + facet_wrap(~T_clonotype_id)
p2 + facet_wrap(~T_clonotype_id)
p3 + facet_wrap(~T_clonotype_id)
dev.off()
####
#plot cell states





###DEG comparison 2024
#DEG clonotype1
head(scRNA[[]])
scRNA1 = JoinLayers(scRNA, assay = "RNA")
sc.subset = subset(scRNA1, subset = T_clonotype_id == "clonotype1")
Idents(sc.subset) = sc.subset@meta.data$orig.ident
#sc.subset = JoinLayers(sc.subset)
x1 = AverageExpression(sc.subset, assays = "RNA", slot = "data")[[1]]
head(x1)
#
clonotype1.deg.11LH = FindMarkers(sc.subset, ident.1 = "11H", ident.2 = "11L", logfc.threshold = 0, assay = "RNA")
clonotype1.deg.11LH = cbind(clonotype1.deg.11LH, x1[rownames(clonotype1.deg.11LH), c("g11L", "g11H")])
head(clonotype1.deg.11LH)
clonotype1.deg.11LH$deg = ifelse(clonotype1.deg.11LH$p_val>=0.005, "ns",
                                 ifelse(clonotype1.deg.11LH$avg_log2FC >0, "high high", "high low"))
#
clonotype1.deg.6LH = FindMarkers(sc.subset, ident.1 = "6H", ident.2 = "6L", logfc.threshold = 0, assay = "RNA")
clonotype1.deg.6LH = cbind(clonotype1.deg.6LH, x1[rownames(clonotype1.deg.6LH), c("g6L", "g6H")])
clonotype1.deg.6LH$deg = ifelse(clonotype1.deg.6LH$p_val>=0.005, "ns",
                                 ifelse(clonotype1.deg.6LH$avg_log2FC >0, "high high", "high low"))


pair.col = RColorBrewer::brewer.pal(12, "Paired")
p1 = ggplot(clonotype1.deg.11LH, aes(x = avg_log2FC, y = - log10(p_val_adj), color = deg)) +
  geom_point(size = 1)+
  theme_bw()+ theme(panel.grid = element_blank())+
  labs(title = "11L vs 11H (clonotype1)")+
  scale_color_manual(values = c(pair.col[c(6,2)], "gray"))
p1

p2 = ggplot(clonotype1.deg.6LH, aes(x = avg_log2FC, y = - log10(p_val_adj), color = deg)) +
  geom_point(size = 1)+
  theme_bw()+theme(panel.grid = element_blank())+
  labs(title = "6L vs 6H (clonotype1)")+
  scale_color_manual(values = c(pair.col[c(6,2)], "gray"))

###
sc.subset = subset(scRNA1, subset = T_clonotype_id == "clonotype2")
Idents(sc.subset) = sc.subset@meta.data$orig.ident
x1 = AverageExpression(sc.subset, assays = "RNA")[[1]]
head(x1)
#
clonotype2.deg.11LH = FindMarkers(sc.subset, ident.1 = "11H", ident.2 = "11L", logfc.threshold = 0, assay = "RNA")
clonotype2.deg.11LH = cbind(clonotype2.deg.11LH, x1[rownames(clonotype2.deg.11LH), c("g11L", "g11H")])
clonotype2.deg.11LH$deg = ifelse(clonotype2.deg.11LH$p_val>=0.005, "ns",
                                 ifelse(clonotype2.deg.11LH$avg_log2FC >0, "high high", "high low"))
#
clonotype2.deg.6LH = FindMarkers(sc.subset, ident.1 = "6H", ident.2 = "6L", logfc.threshold = 0, assay = "RNA")
clonotype2.deg.6LH = cbind(clonotype2.deg.6LH, x1[rownames(clonotype2.deg.6LH), c("g6L", "g6H")])
clonotype2.deg.6LH$deg = ifelse(clonotype2.deg.6LH$p_val>=0.005, "ns",
                                 ifelse(clonotype2.deg.6LH$avg_log2FC >0, "high high", "high low"))

p3 = ggplot(clonotype2.deg.11LH, aes(x = avg_log2FC, y = - log10(p_val_adj),color = deg)) +
  geom_point(size = 1)+
  theme_bw()+ theme(panel.grid = element_blank(), aspect.ratio = 1)+
  labs(title = "11L vs 11H (clonotype2)")+
  scale_color_manual(values = c(pair.col[c(6,2)], "gray"))

p4 = ggplot(clonotype2.deg.6LH, aes(x = avg_log2FC, y = - log10(p_val_adj),color = deg)) +
  geom_point(size = 1)+
  theme_bw()+ theme(panel.grid = element_blank(), aspect.ratio = 1)+
  labs(title = "6L vs 6H (clonotype2)")+
  scale_color_manual(values = c(pair.col[c(6,2)], "gray"))


cowplot::plot_grid(p1, p2, p3, p4, nrow = 2)



identical(rownames(clonotype1.deg.11LH), rownames(clonotype2.deg.11LH))
use.genes1 = intersect(rownames(clonotype1.deg.11LH), rownames(clonotype2.deg.11LH))
plot(clonotype1.deg.11LH[use.genes1,]$avg_log2FC, clonotype2.deg.11LH[use.genes1,]$avg_log2FC)
d11HL = cbind(clonotype1.deg.11LH[use.genes1,], clonotype2.deg.11LH[use.genes1,])
head(d11HL)
colnames(d11HL)[1:8] = paste0("clonotype1:", colnames(d11HL)[1:8])
table(d11HL$`clonotype1:deg`, d11HL$deg)

d11HL$col = ifelse(d11HL$`clonotype1:deg` == "high high" & d11HL$deg == "high high", "consistent high",
                   ifelse(d11HL$`clonotype1:deg` == "high low" & d11HL$deg == "high low", "consistent low",
                          ifelse(d11HL$`clonotype1:deg` == "ns" & d11HL$deg == "ns", "consistent ns", "inconsistent") ))
table(d11HL$col)

p5 = ggplot(d11HL, aes(x = `clonotype1:avg_log2FC`, y = avg_log2FC, col = col))+
  geom_point(size = .1)+
  theme_bw()+
  theme(panel.grid = element_blank(), aspect.ratio = 1)+
  scale_color_manual(values = c(pair.col[c(6,2)], "gray",  pair.col[10]))+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  labs(x = "log2FC (clonotype1)", y = "log2FC (clonotype2)", col = "11H vs 11L")
p5
#
use.genes2 = intersect(rownames(clonotype1.deg.6LH), rownames(clonotype2.deg.6LH))
plot(clonotype1.deg.6LH[use.genes2,]$avg_log2FC, clonotype2.deg.6LH[use.genes2,]$avg_log2FC)
d6HL = cbind(clonotype1.deg.6LH[use.genes2,], clonotype2.deg.6LH[use.genes2,])
head(d6HL)
colnames(d6HL)[1:8] = paste0("clonotype1:", colnames(d6HL)[1:8])
table(d6HL$`clonotype1:deg`, d6HL$deg)

d6HL$col = ifelse(d6HL$`clonotype1:deg` == "high high" & d6HL$deg == "high high", "consistent high",
                   ifelse(d6HL$`clonotype1:deg` == "high low" & d6HL$deg == "high low", "consistent low",
                          ifelse(d6HL$`clonotype1:deg` == "ns" & d6HL$deg == "ns", "consistent ns", "inconsistent") ))
table(d11HL$col)

p6 = ggplot(d6HL, aes(x = `clonotype1:avg_log2FC`, y = avg_log2FC, col = col))+
  geom_point(size = .1)+
  theme_bw()+
  theme(panel.grid = element_blank(), aspect.ratio = 1)+
  scale_color_manual(values = c(pair.col[c(6,2)], "gray",  pair.col[10]))+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  labs(x = "log2FC (clonotype1)", y = "log2FC (clonotype2)", col = "6H vs 6L")

#

###how much is the overlap?
library(ggVennDiagram)
#avidity L/H, 6/11 simulation
#affinity clonotype1/2
#
#over.expression.11LH, high has high expression
head(clonotype1.deg.11LH)
over.expression.clonotype1.11H = rownames(clonotype1.deg.11LH)[clonotype1.deg.11LH$deg == "high high"]
under.expression.clonotype1.11H = rownames(clonotype1.deg.11LH)[clonotype1.deg.11LH$deg == "high low"]
head(clonotype2.deg.11LH)
over.expression.clonotype2.11H = rownames(clonotype2.deg.11LH)[clonotype2.deg.11LH$deg == "high high"]
under.expression.clonotype2.11H = rownames(clonotype2.deg.11LH)[clonotype2.deg.11LH$deg == "high low"]
length(over.expression.clonotype2.11H)
length(under.expression.clonotype2.11H)



p1 = ggvenn(list("Clonotype 1" = over.expression.clonotype1.11H, 
                 "Clonotype 2" = over.expression.clonotype2.11H),
            fill_color = pair.col[c(1,3)],set_name_size = 4)+ 
  labs(title = "11H overexpressed compared with 11L")

p2 = ggvenn(list("Clonotype 1" = under.expression.clonotype1.11H, 
                 "Clonotype 2" = under.expression.clonotype2.11H),
            fill_color = pair.col[c(2,4)],set_name_size = 4)+ 
  labs(title = "11H under compared with 11L")

#
#overexpression.6LH: high has high expression
head(clonotype1.deg.6LH)
overexpression.clonotype1.6H = rownames(clonotype1.deg.6LH)[clonotype1.deg.6LH$deg == "high high"]
underexpression.clonotype1.6H = rownames(clonotype1.deg.6LH)[clonotype1.deg.6LH$deg == "high low"]
head(clonotype2.deg.6LH)
overexpression.clonotype2.6H = rownames(clonotype2.deg.6LH)[clonotype2.deg.6LH$deg == "high high"]
underexpression.clonotype2.6H = rownames(clonotype2.deg.6LH)[clonotype2.deg.6LH$deg == "high low"]

p3 = ggvenn( list("Clonotype 1" = overexpression.clonotype1.6H, 
                       "Clonotype 2" = overexpression.clonotype2.6H), 
             fill_color = pair.col[c(1,3)], set_name_size = 4)+ 
  labs(title = "6H overexpressed compared with 6L")+ theme(legend.position = "none")
p3

p4 = ggvenn(list("Clonotype 1" = underexpression.clonotype1.6H, 
                       "Clonotype 2" = underexpression.clonotype2.6H),
            fill_color = pair.col[c(2,4)], set_name_size = 4)+ 
  labs(title = "6H underexpressed compared with 6L")
p4

pdf("02.deg.clonotype12.HL.pdf", width = 10, height = 15, useDingbats = F)
cowplot::plot_grid(p1,p2,p3,p4, p5, p6, byrow = T, ncol = 2)
dev.off()






#####
grep("TNFRSF", rownames(scRNA), ignore.case = T, value = T)

grep("CD52", rownames(scRNA), ignore.case = T, value = T)


co.inhibitors = c("Pdcd1", "Cd160", "Lag3", "Tigit", "Btla")
co.stimulators = c("Cd28", "Lair1" , "Icos", "Cd27", "Cd2", "Cd226")
p.features = c(co.inhibitors, co.stimulators)
grep(paste0(p.features, collapse = "|"), rownames(scRNA), ignore.case = T, value = T)
p.features[!p.features %in% rownames(scRNA)]
p.features = p.features[p.features %in% rownames(scRNA)]


Idents(scRNA) = scRNA$orig.ident
pdf("02.inhibitor-stimulatators.dotplot.pdf", width = 4, height = 4, useDingbats = F)
DotPlot(scRNA, features = p.features, cols = pair.col[c(2,6)])+ coord_flip() +
  scale_y_discrete(limits  = c("6L", "6H", "11L", "11H"))+
  theme_bw() + #theme(panel.grid = element_line(linetype = "dashed"))
  theme(panel.grid = element_blank())+
  labs(x = "", y = "")
dev.off()

p.features

real = AggregateExpression(scRNA[p.features,])$RNA
table(scRNA$orig.ident)

#randomly make 20 pseudobulk, 200 cells each time
set.seed(123, kind = "L'Ecuyer-CMRG")
sel.cell.list = lapply(1:20, function(k){
    sel.cells = sapply(unique(scRNA$orig.ident), function(x){
          sample(colnames(scRNA)[scRNA$orig.ident == x], 500)
    })
    sel.cells
})

library(edgeR)
pseudobulks = lapply(1:length(sel.cell.list), function(x){
      use.cells = as.vector(sel.cell.list[[x]])
      b = AverageExpression(scRNA[,use.cells])$RNA
      y <- DGEList(counts=b)
      y <- y[,keep.lib.sizes=FALSE]
      y <- normLibSizes(y)
      b = cpm(y)[p.features,]
      colnames(b) = paste0(colnames(b), "-", x)
      b
})
pseudobulks[[1]]

pseudobulks.df = data.matrix(do.call(cbind, pseudobulks))

pseudobulks.df = t(apply(pseudobulks.df, 1, scale))
colnames(pseudobulks.df) = colnames(do.call(cbind, pseudobulks))


library(ComplexHeatmap)

pseudobulks.df.meta = data.frame(row.names = colnames(pseudobulks.df), 
                                 sample = substr(colnames(pseudobulks.df), 2, 4))
pseudobulks.df.meta$sample = gsub("-", "", pseudobulks.df.meta$sample)

sample.cluster = cluster_within_group(pseudobulks.df, pseudobulks.df.meta$sample)



gene.meta = data.frame(row.names = c(co.inhibitors, co.stimulators), group = c(rep("co-inhibitors", length(co.inhibitors)),
                                                                           rep("co-stimulators", length(co.stimulators))) )
gene.meta = gene.meta[rownames(pseudobulks.df),,drop =F]

cluster.gene = cluster_within_group(t(pseudobulks.df), gene.meta$group)

gene.meta

pair.col = RColorBrewer::brewer.pal(12, "Paired")

top.ann = HeatmapAnnotation(sample = pseudobulks.df.meta$sample, col = list(sample = sample.color))
row.ann = rowAnnotation(group = gene.meta$group, col = list(group = c("co-inhibitors" = pair.col[8], "co-stimulators" = pair.col[4])))


h = Heatmap(pseudobulks.df, name = "z-score",
        cluster_columns = sample.cluster, cluster_rows = cluster.gene,
        top_annotation = top.ann, left_annotation = row.ann,
        show_column_names = F, split = 2, column_split = 4,
        show_row_dend = F, show_column_dend = F)
h

pdf("03.groups.InhibitorStimulator.heatmap.pdf", width = 10, height = 5, useDingbats = F)
draw(h)
dev.off()

###clonotype 1-2
set.seed(123, kind = "L'Ecuyer-CMRG")
sel.cell.list = lapply(1:20, function(k){
  sel.cells = sapply(unique(scRNA$orig.ident), function(x){
    sample(colnames(scRNA)[scRNA$orig.ident == x & scRNA@meta.data$T_clonotype_id == "clonotype1"], 200)
  })
  sel.cells
})

pseudobulks = lapply(1:length(sel.cell.list), function(x){
  use.cells = as.vector(sel.cell.list[[x]])
  b = AverageExpression(scRNA[, use.cells])$RNA
  y <- DGEList(counts=b)
  y <- y[,keep.lib.sizes=FALSE]
  y <- normLibSizes(y)
  b = cpm(y)[p.features,]
  colnames(b) = paste0(colnames(b), "-", x)
  b
})

pseudobulks.df = data.matrix(do.call(cbind, pseudobulks))
pseudobulks.df = t(apply(pseudobulks.df, 1, scale))
colnames(pseudobulks.df) = colnames(do.call(cbind, pseudobulks))
#
pseudobulks.df.meta = data.frame(row.names = colnames(pseudobulks.df), 
                                 sample = substr(colnames(pseudobulks.df), 2, 4))
pseudobulks.df.meta$sample = gsub("-", "", pseudobulks.df.meta$sample)

sample.cluster = cluster_within_group(pseudobulks.df, pseudobulks.df.meta$sample)

co.inhibitors = c("Pdcd1", "Cd160", "Lag3", "Tigit", "Btla")
co.stimulators = c("Cd28","Cd69","Lair1" , "Icos", "Cd27", "Cd2", "Cd226")

gene.meta = data.frame(row.names = c(co.inhibitors, co.stimulators), group = c(rep("co-inhibitors", length(co.inhibitors)),
                                                                               rep("co-stimulators", length(co.stimulators))) )
gene.meta = gene.meta[rownames(pseudobulks.df),,drop =F]

cluster.gene = cluster_within_group(t(pseudobulks.df), gene.meta$group)

gene.meta

pair.col = RColorBrewer::brewer.pal(12, "Paired")

top.ann = HeatmapAnnotation(sample = pseudobulks.df.meta$sample, col = list(sample = sample.color))
row.ann = rowAnnotation(group = gene.meta$group, col = list(group = c("co-inhibitors" = pair.col[8], "co-stimulators" = pair.col[4])))


h = Heatmap(pseudobulks.df, name = "z-score",
            cluster_columns = sample.cluster, cluster_rows = cluster.gene,
            top_annotation = top.ann, left_annotation = row.ann,
            show_column_names = F, split = 2, column_split = 4,
            show_row_dend = F, show_column_dend = F)
h

pdf("03.groups.clonetype1.InhibitorStimulator.heatmap.pdf", width = 10, height = 5, useDingbats = F)
draw(h)
dev.off()

###
###clonotype 2
set.seed(123, kind = "L'Ecuyer-CMRG")
sel.cell.list = lapply(1:20, function(k){
  sel.cells = sapply(unique(scRNA$orig.ident), function(x){
    sample(colnames(scRNA)[scRNA$orig.ident == x & scRNA@meta.data$T_clonotype_id == "clonotype2"], 200, replace = T)
  })
  sel.cells
})

pseudobulks = lapply(1:length(sel.cell.list), function(x){
  use.cells = as.vector(sel.cell.list[[x]])
  b = AverageExpression(scRNA[, use.cells])$RNA
  y <- DGEList(counts=b)
  y <- y[,keep.lib.sizes=FALSE]
  y <- normLibSizes(y)
  b = cpm(y)[p.features,]
  colnames(b) = paste0(colnames(b), "-", x)
  b
})

pseudobulks.df = data.matrix(do.call(cbind, pseudobulks))
pseudobulks.df = t(apply(pseudobulks.df, 1, scale))
colnames(pseudobulks.df) = colnames(do.call(cbind, pseudobulks))
#
pseudobulks.df.meta = data.frame(row.names = colnames(pseudobulks.df), 
                                 sample = substr(colnames(pseudobulks.df), 2, 4))
pseudobulks.df.meta$sample = gsub("-", "", pseudobulks.df.meta$sample)

sample.cluster = cluster_within_group(pseudobulks.df, pseudobulks.df.meta$sample)

co.inhibitors = c("Pdcd1", "Cd160", "Lag3", "Tigit", "Btla")
co.stimulators = c("Cd28","Cd69","Lair1" , "Icos", "Cd27", "Cd2", "Cd226")

gene.meta = data.frame(row.names = c(co.inhibitors, co.stimulators), group = c(rep("co-inhibitors", length(co.inhibitors)),
                                                                               rep("co-stimulators", length(co.stimulators))) )
gene.meta = gene.meta[rownames(pseudobulks.df),,drop =F]

cluster.gene = cluster_within_group(t(pseudobulks.df), gene.meta$group)

gene.meta

pair.col = RColorBrewer::brewer.pal(12, "Paired")

top.ann = HeatmapAnnotation(sample = pseudobulks.df.meta$sample, col = list(sample = sample.color))
row.ann = rowAnnotation(group = gene.meta$group, col = list(group = c("co-inhibitors" = pair.col[8], "co-stimulators" = pair.col[4])))


h = Heatmap(pseudobulks.df, name = "z-score",
            cluster_columns = sample.cluster, cluster_rows = cluster.gene,
            top_annotation = top.ann, left_annotation = row.ann,
            show_column_names = F, split = 2, column_split = 4,
            show_row_dend = F, show_column_dend = F)
h

pdf("03.groups.clonetype2.InhibitorStimulator.heatmap.pdf", width = 10, height = 5, useDingbats = F)
draw(h)
dev.off()






###DEG
jianping.features1 <- c("Cd2", "Ifngr1", "Pdcd1", "Lag3", "Cd28", "Cd69", "Ly6a", "Ly6e", "Gramd3", "Cd8b1", "Cd8a", "Cd6", "Tigit")
#jianping.features2 <- c("TCR", "Cd8", "Cd3", "Tim3", "Klrg1", "IL-4r")
grep(paste0(unique(c(jianping.features1, jianping.features2)), collapse = "|"), rownames(scRNA), ignore.case = T, value = T)

grep(paste0(c("CD6", "TIGIT"), collapse = "|"), rownames(scRNA), ignore.case = T, value = T)
grep("tcrb", rownames(scRNA), ignore.case = T, value = T)



#############DEG analysis

library(fgsea)
mouse.c5 <- msigdbr::msigdbr(species = "mouse", category = "C5") %>%
  dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
mouse.c5.split <- split(x = mouse.c5$gene_symbol, f = mouse.c5$gs_name)


head(scRNA)
scRNA@meta.data$HLgroup <- ifelse(grepl("H",scRNA@meta.data$orig.ident), "H", "L")
###
scRNA.c1 <- subset(scRNA, subset = T_clonotype_id == "clonotype1")
table(scRNA.c1@meta.data$orig.ident, scRNA.c1@meta.data$T_clonotype_id)
Idents(scRNA.c1) <- scRNA.c1@meta.data$HLgroup

#library(ggpubr)
#p <- VlnPlot(scRNA.c1, features = jianping.features1, cols = as.character(sample.color[c(3,1)]), group.by = "HLgroup") & 
#  theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank())
#p

#p <- DotPlot(scRNA.c1, features = jianping.features1, group.by = "HLgroup") + RotatedAxis() &
#  theme_bw() + theme(panel.grid = element_blank())
#pdf("02.clontype1.Vlnplt.pdf", width = 10, height = 10, useDingbats = F)  
#pdf("02.clontype1.Dotplot.pdf", width = 8, height = 3, useDingbats = F)
#p
#dev.off()

deg1 <- FindMarkers(scRNA.c1, logfc.threshold=0, test.use = "wilcox", ident.1 = "H", ident.2 = "L")
x <-AverageExpression(scRNA.c1)[[1]]
deg1 <- cbind(deg1, x[rownames(deg1),])
head(deg1)

write.csv(deg1, file = "02.clontype1.HvsL.DEG.csv")
deg1.rank <- deg1$avg_log2FC
names(deg1.rank) <- rownames(deg1)
fgsea.1 <- fgsea(pathways = mouse.c5.split, deg1.rank)
head(fgsea.1[order(pval), ], n = 10)
x <- do.call(cbind, fgsea.1)
write.csv(x, file = "02.clonotype1.HvsL.GO_fgsea.csv", row.names = F)


###
scRNA.c1 <- subset(scRNA, subset = T_clonotype_id == "clonotype2") #, "clonotype2"))
table(scRNA.c1@meta.data$orig.ident, scRNA.c1@meta.data$T_clonotype_id)
#Identify DEG between H v.s. L groups
Idents(scRNA.c1) <- scRNA.c1@meta.data$HLgroup

head(scRNA.c1[[]])
p <- VlnPlot(scRNA.c1, features = jianping.features1, cols = as.character(sample.color[c(3,1)]), group.by = "HLgroup") & 
  theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank())

p <- DotPlot(scRNA.c1, features = jianping.features1, group.by = "HLgroup") + RotatedAxis() &
  theme_bw() + theme(panel.grid = element_blank())
p

#pdf("02.clontype2.Vlnplt.pdf", width = 10, height = 10, useDingbats = F)  
pdf("02.clontype2.Dotplot.pdf", width = 8, height = 3, useDingbats = F)
p
dev.off()

deg2 <- FindMarkers(scRNA.c1, logfc.threshold=0, test.use = "wilcox", ident.1 = "H", ident.2 = "L")
x <-AverageExpression(scRNA.c1)[[1]]
deg2 <- cbind(deg2, x[rownames(deg2),])
head(deg2)
write.csv(deg2, file = "02.clontype2.HvsL.DEG.csv")
###
deg2.rank <- deg2$avg_log2FC
names(deg2.rank) <- rownames(deg2)
fgsea.2 <- fgsea(pathways = mouse.c5.split, deg2.rank)
head(fgsea.2[order(pval), ], n = 10)
x <- do.call(cbind, fgsea.2)
write.csv(x, file = "02.clonotype2.HvsL.GO_fgsea.csv", row.names = F)



deg.signif1 <- deg1[deg1$p_val_adj<0.05 & abs(deg1$avg_log2FC)>0.25,]
deg.hilight1 <- deg.signif1[rownames(deg.signif1) %in% c(jianping.features1, jianping.features2),]
deg.hilight1$name <- rownames(deg.hilight1)
deg.hilight1

deg.signif2 <- deg2[deg2$p_val_adj<0.05 & abs(deg2$avg_log2FC)>0.25,]
deg.hilight2 <- deg.signif2[rownames(deg.signif2) %in% c(jianping.features1, jianping.features2),]
deg.hilight2$name <- rownames(deg.hilight2)
deg.hilight2

library(ggrepel)

p <- ggplot(deg1, aes(x=avg_log2FC, y = -log10(p_val_adj)))+
  annotate("rect", xmin = 0.25, xmax = Inf, ymin = -log10(0.05), ymax = Inf, fill = sample.color[3])+
  annotate("text", x = 2, y = 200, label = "H overexpression")+
  annotate("rect", xmin = -0.25, xmax = -Inf, ymin = -log10(0.05), ymax = Inf, fill = sample.color[1])+
  annotate("text", x = -2, y = 200, label = "L overexpression")+
  geom_point(color = "gray2")+
  geom_text_repel(inherit.aes = F, data = deg.hilight1, 
                  aes(x = avg_log2FC, y = -log10(p_val_adj), label = name), 
                  fontface = "bold", box.padding = 0.5, max.overlaps = Inf, color = "black")+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank())+
  labs(x = "") #+ coord_cartesian(expand = F)
p  

pdf("02.clonotype1.HvsL.DEG.volcano.pdf", width = 8, height = 7, useDingbats = F)
p
dev.off()

###
deg.hilight2
p <- ggplot(deg2, aes(x=avg_log2FC, y = -log10(p_val_adj)))+
  annotate("rect", xmin = 0.25, xmax = Inf, ymin = -log10(0.05), ymax = Inf, fill = sample.color[3])+
  annotate("text", x = 2, y = 120, label = "H overexpression")+
  annotate("rect", xmin = -0.25, xmax = -Inf, ymin = -log10(0.05), ymax = Inf, fill = sample.color[1])+
  annotate("text", x = -2, y = 120, label = "L overexpression")+
  geom_point(color = "gray2")+
  geom_text_repel(inherit.aes = F, data = deg.hilight2, 
                  aes(x = avg_log2FC, y = -log10(p_val_adj), label = name), 
                  fontface = "bold", box.padding = 0.5, max.overlaps = Inf, color = "black")+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank())+
  labs(x = "") #+ coord_cartesian(expand = F)
p  

pdf("02.clonotype2.HvsL.DEG.volcano.pdf", width = 8, height = 7, useDingbats = F)
p
dev.off()



###overlaped DEG
xx <- intersect(rownames(deg.signif1), rownames(deg.signif2))
sum(deg.signif1[xx,"avg_log2FC"] * deg.signif2[xx, "avg_log2FC"] < 0 )

xx
rm(scRNA.c1)
scRNA.c12 <- subset(scRNA, subset = T_clonotype_id %in% c("clonotype1", "clonotype2"))
Idents(scRNA.c12) <- paste0(scRNA.c12@meta.data$HLgroup, "-", scRNA.c12@meta.data$T_clonotype_id)
table(scRNA.c12@meta.data$orig.ident, scRNA.c12@meta.data$T_clonotype_id)

make.heatmap <- function(seurat, markers){
  ht_opt$message = FALSE
  #seurat <- sc.subset.subset
  #markers <- rownames(mk1)
  mat<- seurat[["RNA"]]@data[rownames(markers), ] %>% as.matrix()
  
  row.ann <- rowAnnotation(DEG= ifelse(markers$avg_log2FC > 0, "high", "low"),
                           col = list(DEG=c("high" = "red", "low" = "blue")))
  
  ## scale the rows
  mat<- t(scale(t(mat)))
  #cluster_anno<- Idents(seurat)
  cluster_anno <- HeatmapAnnotation(HL = seurat[[]]$HLgroup, clonotype = seurat[[]]$T_clonotype_id,
                                    col = list(HL = c("H"= as.character(sample.color[3]), "L" = as.character(sample.color[1]))))
  col_fun = circlize::colorRamp2(c(-3, 0, 3), c("#FF00FF", "black", "#FFFF00"))
  #plot the heatmap
  
  h <- Heatmap(mat, name = "Expression", 
               column_split = Idents(seurat),#factor(cluster_anno),
               row_split = factor(ifelse(markers$avg_log2FC > 0, "high", "low")),
               cluster_columns = TRUE,
               show_column_dend = FALSE,
               cluster_column_slices = TRUE,
               column_title_gp = gpar(fontsize = 8),
               column_gap = unit(0.5, "mm"),
               cluster_rows = TRUE,
               show_row_dend = FALSE,
               col = col_fun,
               row_names_gp = gpar(fontsize = 4),
               column_title_rot = 90,
               left_annotation = row.ann,
               top_annotation = cluster_anno, #HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
               show_column_names = FALSE,
               use_raster = TRUE,
               raster_quality = 4)
  h
}


h <- make.heatmap(scRNA.c12, deg1[xx,])

png("02.OverlapDEG.png", width = 12, height = 8, units = "in", res = 200)
draw(h)
dev.off()

scRNA.c1 <- subset(scRNA, subset = T_clonotype_id == "clonotype2")
table(scRNA.c1@meta.data$orig.ident, scRNA.c1@meta.data$T_clonotype_id)
Idents(scRNA.c1) <- scRNA.c1@meta.data$HLgroup
h <- make.heatmap(scRNA.c1, deg1[xx,])
png("02.overlapedDEG.in.clontype2.heatmap.png", width = 12, height = 8, units = "in", res = 200)
draw(h)
dev.off()

###
jianping.features1 %in% rownames(scRNA)
jianping.features2 %in% rownames(scRNA)

deg["Ly6a",]

###

#RidgePlot(scRNA.c12, features = jianping.features1, ncol = 2)

#DotPlot(scRNA.c12, features = jianping.features1) + RotatedAxis()

#VlnPlot(scRNA.c12, features = jianping.features1, cols = as.character(sample.color[c(3,1)]), split.by = "HLgroup", split.plot = T) & 
#  theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank())





###
library(fgsea)
mouse.c5 <- msigdbr::msigdbr(species = "mouse", category = "C5") %>%
  dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
mouse.c5.split <- split(x = mouse.c5$gene_symbol, f = mouse.c5$gs_name)

deg.rank <- deg$avg_log2FC
names(deg.rank) <- rownames(deg)
res.fgsea <- fgsea(pathways = mouse.c5.split, deg.rank)
head(res.fgsea[order(pval), ], n = 10)

save(res.fgsea, file = "02.HvsL.DEG.GO-_fgsea.rda")

x <- do.call(cbind, res.fgsea)
head(x)
write.csv(x, file = "02.HvsL.DEG.GO-_fgsea.csv", row.names = F)

###
topPathwaysUp <- res.fgsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- res.fgsea[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(mouse.c5.split[topPathways], deg.rank, res.fgsea, gseaParam=0.5)

###permutation test to find p-value and logFC
#deg.p.log2fc.permutation <- lapply(1:1000, function(k){
#    message(k)
#    Idents(scRNA.main) <- sample(scRNA.main@meta.data$HLgroup)
#    deg <- FindMarkers(scRNA.main, logfc.threshold=0, test.use = "wilcox", ident.1 = "H", ident.2 = "L", verbose = F)
#    x1 <- deg$p_val[order(deg$p_val, decreasing = F)][round(nrow(deg)*0.05)]
#    x2 <- abs(deg$avg_log2FC[order(abs(deg$avg_log2FC), decreasing = T)][round(nrow(deg)*0.05)])
#    c(x1, x2)
#})

nrow(deg)
sum(deg$p_val_adj<0.05 & abs(deg$avg_log2FC)>0.25)


####Not Done

library(monocle)
#remotes::install_github('satijalab/seurat-wrappers', force = T)
#library(SeuratWrappers)

sc.subset <- subset(scRNA, subset = T_clonotype_id %in% c("clonotype4", "clonotype2", "clonotype3", "clonotype4"))
scRNA
sc.subset
#cds <- importCDS(scRNA)
#Extract data, phenotype data, and feature data from the SeuratObject
#data <- as(as.matrix(scRNA@assays$RNA@data), 'sparseMatrix')
fData <- data.frame(gene_short_name = row.names(as(as.matrix(sc.subset@assays$RNA@data), 'sparseMatrix')), 
                    row.names = row.names(as(as.matrix(sc.subset@assays$RNA@data), 'sparseMatrix')))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
monocle_cds <- newCellDataSet(sc.subset@assays$RNA@data,
                              phenoData = new('AnnotatedDataFrame', data = sc.subset@meta.data),
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
monocle_cds

#estimateSizeFactors() and estimateDispersions() will only work, and are only needed, 
#if you are working with a CellDataSet with a negbinomial() or negbinomial.size() expression family.
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

#HSMM <- detectGenes(HSMM, min_expr = 0.1)

#Trajectory step 1: choose genes that define a cell's progress
diff_test_res <- differentialGeneTest(monocle_cds, fullModelFormulaStr = "~cluster", cores = 2)
#ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:1000]
tail(ordering_genes)
##we need to set DEGs in the HSMM object, the next several functions will depend on them.
monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
plot_ordering_genes(monocle_cds)

#Trajectory step 2: reduce data dimensionality
monocle_cds <- reduceDimension(monocle_cds, max_components = 2, method = 'DDRTree')

#Trajectory step 3: order cells along the trajectory
monocle_cds <- orderCells(monocle_cds)
pdf("test.monocle.pdf", width = 10, height = 10, useDingbats = F)
plot_cell_trajectory(monocle_cds, color_by = "orig.ident")
plot_cell_trajectory(monocle_cds, color_by = "cluster")
dev.off()


q(save = "no")


###
p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "T_clonotype_id", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)






##save data for scVelo
library(SeuratDisk)
SaveH5Seurat(scRNA, filename = "scTCR.h5Seurat")
Convert("scTCR.h5Seurat", dest = "h5ad")
#clonotype1 and 2
head(scRNA[[]])
sc.subset <- subset(scRNA, subset = T_clonotype_id == "clonotype1")
SaveH5Seurat(sc.subset, filename = "scTCRclontype1.h5Seurat")
Convert("scTCRclontype1.h5Seurat", dest = "h5ad")
###
sc.subset <- subset(scRNA, subset = T_clonotype_id == "clonotype2")
SaveH5Seurat(sc.subset, filename = "scTCRclontype2.h5Seurat")
Convert("scTCRclontype2.h5Seurat", dest = "h5ad")

###DEG
##################################################################################################################
head(scRNA@meta.data)

mouse.c5 <- msigdbr::msigdbr(species = "mouse", category = "C5") %>%
  dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
mouse.c5.list = split(x = mouse.c5$gene_symbol, f = mouse.c5$gs_name)

DefaultAssay(scRNA) = "RNA"

lapply(c("clonotype1", "clonotype2"), function(kk){
  #kk <- "clonotype1"
  sc.subset <- subset(scRNA, subset = T_clonotype_id %in% kk) # "clonotype2"))
  Idents(sc.subset) <- factor(sc.subset@meta.data$orig.ident)
  mk1 <- FindMarkers(sc.subset, ident.1 = "6L", ident.2 = "6H", logfc.threshold = 0.25, min.pct = 0.2)
  sc.subset.subset <- subset(sc.subset, subset = orig.ident %in% c("6L", "6H"))
  x <- AverageExpression(sc.subset.subset, group.by = "ident", assays = "RNA", slot = 'data')[[1]]
  mk1 <- cbind(mk1, x[rownames(mk1),])
  write.csv(mk1, file = paste0("03.", kk, ".6HvsL.markers.csv"), quote = F)
  h <- make.heatmap(sc.subset.subset, markers = mk1)
  png(paste0("03.", kk, ".6HvsL.png"), width = 10, height = 10, units = "in", res = 150)
  #print(DoHeatmap(sc.subset.subset, features = rownames(mk1), assay = "RNA")+ NoLegend())
  draw(h)
  dev.off()
  go <- NULL
  go <- enricher(gene = rownames(mk1)[mk1$avg_log2FC > 0], TERM2GENE = mouse.c5, pvalueCutoff = 1, qvalueCutoff = 1)@result
  #go <- go[go$p.adjust <0.05,]
  if(nrow(go) > 0) write.csv(go, file = paste0("03.", kk, ".6HvsL.logFCpositive.GO.csv"))
  go <- NULL
  go <- enricher(gene = rownames(mk1)[mk1$avg_log2FC < 0], TERM2GENE = mouse.c5, pvalueCutoff = 1, qvalueCutoff = 1)@result
  #go <- go[go$p.adjust <0.05,]
  if(nrow(go) > 0) write.csv(go, file = paste0("03.", kk, ".6HvsL.logFCnegative.GO.csv"))
  
  ###
  mk2 <- FindMarkers(sc.subset, ident.1 = "11L", ident.2 = "11H", logfc.threshold = 0.25, min.pct = 0.2)
  sc.subset.subset <- subset(sc.subset, subset = orig.ident %in% c("11L", "11H"))
  x <- AverageExpression(sc.subset.subset, group.by = "ident", assays = "RNA", slot = 'data')[[1]]
  mk2 <- cbind(mk2, x[rownames(mk2),])
  write.csv(mk2, file = paste0("03.", kk, ".11HvsL.markers.csv"), quote = F)
  h <- make.heatmap(sc.subset.subset, markers = mk2)
  png(paste0("03.",kk, ".11HvsL.png"), width = 10, height = 10, units = "in", res = 150)
  #print(DoHeatmap(sc.subset.subset, features = rownames(mk2), assay = "RNA")+ NoLegend())
  draw(h)
  dev.off()
  go <- NULL
  go <- enricher(gene = rownames(mk2)[mk2$avg_log2FC > 0], TERM2GENE = mouse.c5, pvalueCutoff = 1, qvalueCutoff = 1)@result
  #go <- go[go$p.adjust <0.05,]
  if(nrow(go) > 0) write.csv(go, file = paste0("03.", kk, ".11HvsL.logFCpositive.GO.csv"))
  go <- NULL
  go <- enricher(gene = rownames(mk2)[mk2$avg_log2FC < 0], TERM2GENE = mouse.c5, pvalueCutoff = 1, qvalueCutoff = 1)@result
  #go <- go[go$p.adjust <0.05,]
  if(nrow(go) > 0) write.csv(go, file = paste0("03.", kk, ".11HvsL.logFCnegative.GO.csv"))
  
  ###
  #mk3 <- FindMarkers(sc.subset, ident.1 = "6L", ident.2 = "11L", logfc.threshold = 0.25, min.pct = 0.2)
  #sc.subset.subset <- subset(sc.subset, subset = orig.ident %in% c("6L", "11L"))
  #x <- AverageExpression(sc.subset.subset, group.by = "ident", assays = "RNA", slot = 'data')[[1]]
  #mk3 <- cbind(mk3, x[rownames(mk3),])
  #write.csv(mk3, file = paste0("03.", kk, ".L6v11.markers.csv"), quote = F)
  #h <- make.heatmap(sc.subset.subset, markers = mk3)
  #png(paste0("03.", kk, ".L6vs11.png"), width = 10, height = 10, units = "in", res = 150)
  #draw(h)
  #print(DoHeatmap(sc.subset.subset, features = rownames(mk2), assay = "RNA")+ NoLegend())
  #dev.off()
  #go <- NULL
  #go <- enricher(gene = rownames(mk3)[mk3$avg_log2FC > 0], TERM2GENE = mouse.c5)@result
  #go <- go[go$p.adjust <0.05,]
  #if(nrow(go) > 0) write.csv(go, file = paste0("03.", kk, ".L6vs11.logFCpositive.GO.csv"))
  #go <- NULL
  #go <- enricher(gene = rownames(mk3)[mk3$avg_log2FC < 0], TERM2GENE = mouse.c5)@result
  #go <- go[go$p.adjust <0.05,]
  #if(nrow(go) > 0) write.csv(go, file = paste0("03.", kk, ".L6vs11.logFCnegative.GO.csv"))
  ###
  #mk4 <- FindMarkers(sc.subset, ident.1 = "6H", ident.2 = "11H", logfc.threshold = 0.25, min.pct = 0.2)
  #sc.subset.subset <- subset(sc.subset, subset = orig.ident %in% c("6H", "11H"))
  #x <- AverageExpression(sc.subset.subset, group.by = "ident", assays = "RNA", slot = 'data')[[1]]
  #mk4 <- cbind(mk4, x[rownames(mk4),])
  #write.csv(mk4, file = paste0("03.", kk, ".H6vs11.markers.csv"), quote = F)
  #h <- make.heatmap(sc.subset.subset, markers = mk4)
  #png(paste0("03.", kk, ".H6vs11.png"), width = 10, height = 10, units = "in", res = 150)
  #draw(h)
  #print(DoHeatmap(sc.subset.subset, features = rownames(mk3), assay = "RNA")+ NoLegend())
  #dev.off()
  #go <- NULL
  #go <- enricher(gene = rownames(mk4)[mk4$avg_log2FC > 0], TERM2GENE = mouse.c5)@result
  #go <- go[go$p.adjust <0.05,]
  #if(nrow(go) > 0) write.csv(go, file = paste0("03.", kk, ".H6vs11.logFCpositive.GO.csv"))
  #go <- NULL
  #go <- enricher(gene = rownames(mk4)[mk4$avg_log2FC < 0], TERM2GENE = mouse.c5)@result
  #go <- go[go$p.adjust <0.05,]
  #if(nrow(go) > 0) write.csv(go, file = paste0("03.", kk, ".H6vs11.logFCnegative.GO.csv"))
  ###
})


###
HL.6h = read.csv("03.clonotype1.6HvsL.logFCnegative.GO.csv", head =T, row.names = 1)
HL.6l = read.csv("03.clonotype1.6HvsL.logFCpositive.GO.csv", head =T, row.names = 1)

HL.11h = read.csv("03.clonotype1.11HvsL.logFCnegative.GO.csv", head =T, row.names = 1)
HL.11l = read.csv("03.clonotype1.11HvsL.logFCpositive.GO.csv", head =T, row.names = 1)


ggplot(HL.6h[order(HL.6h$p.adjust, decreasing = F)[1:50],], aes(x = reorder(ID, -p.adjust), y = -log10(p.adjust) )) +
  geom_histogram(stat = "identity")+
  theme_bw()+ theme(panel.grid = element_blank()) +
  coord_flip()


ggplot(HL.11h[order(HL.11h$p.adjust, decreasing = F)[1:10],], aes(x = reorder(ID, -p.adjust), y = -log10(p.adjust) )) +
  geom_histogram(stat = "identity")+
  theme_bw()+ theme(panel.grid = element_blank()) +
  coord_flip()


nrow(HL.6n)
nrow(HL.6p)

keep.go = unique(c(HL.6h$ID[HL.6h$p.adjust<0.01],
                   HL.6l$ID[HL.6l$p.adjust<0.01],
                   HL.11h$ID[HL.11h$p.adjust<0.01],
                   HL.11l$ID[HL.11l$p.adjust<0.01]))


length(keep.go)


#HL.6p = HL.6p[order(HL.6p$p.adjust, decreasing = F)[1:10],]
#HL.6n = HL.6n[order(HL.6n$p.adjust, decreasing = F)[1:10],]
#HL.11p = HL.11p[order(HL.11p$p.adjust, decreasing = F)[1:10],]
#HL.11n = HL.11n[order(HL.11n$p.adjust, decreasing = F)[1:10],]


data = data.frame(ID = keep.go,
                  "6L" = HL.6l$p.adjust[match(keep.go, HL.6l$ID)],
                  "6H" = HL.6h$p.adjust[match(keep.go, HL.6h$ID)],
                  "11L" = HL.11l$p.adjust[match(keep.go, HL.11l$ID)],
                  "11H" = HL.11h$p.adjust[match(keep.go, HL.11h$ID)])

head(data)

write.csv(data, file = "01.HL116.clontype1.go.csv", row.names = F)



d1 = reshape2::melt(data, id.vars = c("ID"))
head(d1)
d1$variable = gsub("^X", "", d1$variable)

use.go = data.frame(readxl::read_xlsx("use.go.xlsx"))[,1]
use.go

custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

ggplot(d1[d1$ID%in% use.go,], aes(x = variable, y = ID, fill = -log10(value))) +
  geom_tile()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_gradient2(limits = c(2,20), low= custom_magma[1:111],  
                       mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 15,na.value = "white")


#####d1#

data = rbind(HL.6n, HL.6p, HL.11n, HL.11p)
head(data, 2)

data$yy = ifelse(grepl("L", data$group), log10(data$p.adjust), - log10(data$p.adjust))
data$time = ifelse(grepl("11", data$group), "11", "6")

table(data$yy>0)


ggplot(data, aes(x = reorder(ID, yy), y = yy, fill = group)) +
  geom_histogram(stat = "identity")+
  theme_bw()+ theme(panel.grid = element_blank()) +
  facet_wrap(~time, scales = "free")+
  scale_fill_manual(values = sample.color)+
  coord_flip()


###VIM

sc.subset <- subset(scRNA, subset = T_clonotype_id %in% "clonotype1") 

library(ggpubr)
x1 = data.frame(VIM = sc.subset@assays$RNA$data["Vim",], Group = sc.subset@meta.data$orig.ident )

my_comparisons <- list( c("11H", "11L"), c("6H", "6L") ) 

p1 = ggplot(x1, aes(x = Group, y = VIM, fill = Group))+
  geom_jitter(color = "gray", size = .2)+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom='errorbar', linetype=1, width=0.2)+
  scale_fill_manual(values = sample.color)+
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 1)+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+
  labs(title = "clonotype1")
  
sc.subset <- subset(scRNA, subset = T_clonotype_id %in% "clonotype2") 

x2 = data.frame(VIM = sc.subset@assays$RNA$data["Vim",], Group = sc.subset@meta.data$orig.ident )


p2 = ggplot(x2, aes(x = Group, y = VIM, fill = Group))+
  geom_jitter(color = "gray", size = .2)+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom='errorbar', linetype=1, width=0.2)+
  scale_fill_manual(values = sample.color)+
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 1)+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+
  labs(title = "clonotype2")

p2
pdf("Vim.expression.pdf", width = 10, height = 4, useDingbats = F)
cowplot::plot_grid(p1, p2)
dev.off()

ggboxplot(ToothGrowth, x = "dose", y = "len",
          color = "dose", palette = "jco")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value
