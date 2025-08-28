library(Seurat)
library(ggplot2)
library(clustree)
library(paletteer)
library(ComplexHeatmap)
library(dplyr)
library(clusterProfiler)
library(CytoTRACE)
rm(list=ls())
dev.off()
setwd("/Users/wuz6/Documents/Project/lij36/2025/")

RColorBrewer::display.brewer.pal(12, "Paired")
sample.color <- RColorBrewer::brewer.pal(12, "Paired")[c(2,4,8,10)]
names(sample.color) <- c("11H", "6H", "11L", "6L")

#BiocManager::install("monocle")



###
if(0){ #runed
##############################################################################
##############################################################################
all.samples <- list.files("../data/")
all.samples
#all.samples <- grep("png|VDJ", all.samples, invert = T, value = T)
all.samples <- grep("vdj", all.samples, invert = T, value = T)
all.samples

###merge single cell RNA
if(0){
  sc.data.list = lapply(all.samples, function(id){
    cat("Reading", id, "\n")
    sc = Read10X(data.dir = file.path("data", id, "filtered_feature_bc_matrix/"))
    sc <- CreateSeuratObject(counts = sc, min.cells = 3, min.features = 500, project = id);
    sc <- RenameCells(sc, paste0(id, "_", colnames(sc)))
    sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^mt-")
    #sc <- add_clonotype(paste0(id, "-VDJ/"), sc, "T")
    sc <- subset(sc, subset = nFeature_RNA > 500 & percent.mt < 20)
    sce <- scDblFinder(sc[['RNA']]$counts)
    #sce <- findDoubletClusters(sce)
    sce.dbl <- data.frame(sce@colData); rm(sce)
    sc <- AddMetaData(sc, sce.dbl)
    #sc <- subset(sc, scDblFinder.class == "singlet")
    sc$orig.ident <- id
    return(sc);
  })
  colnames(sc.data.list[[1]])[1:4]
}

load("sc.rna.filterDoublet.rda")
scRNA = merge(sc.data.list[[1]], y = sc.data.list[-1]) #,add.cell.ids = all.samples)
head(scRNA@meta.data)
scRNA <- subset(scRNA, subset = scDblFinder.class == "singlet")
rm(sc.data.list)

##change cell name to normal name
colnames(scRNA) <- sapply(colnames(scRNA), function(x){ paste0(strsplit(x, "_")[[1]][1:2], collapse = "_") })

###
p00 <- VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p01 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 1)

p00
p01

pdf("01.scQCpass.pdf", width = 7, height = 5, useDingbats = F)
p00;
p01
dev.off()


scRNA <- NormalizeData(scRNA, assay = "RNA")
scRNA <- ScaleData(scRNA, features = rownames(scRNA), assay = "RNA")
scRNA <- FindVariableFeatures(scRNA, assay = "RNA")

x = data.frame(table(scRNA@meta.data$orig.ident))

p0 <- ggplot(x, aes(x = Var1, y = Freq, fill = Var1))+
  geom_bar(stat = "identity")+
  geom_text(aes(x=Var1, y = Freq+100, label = Freq))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  #coord_polar()+
  scale_fill_manual(values = sample.color)+
  scale_x_discrete(limits = c("6L", "6H", "11L", "11H"))+
  labs(x="", fill = "", y = "count")
p0

scRNA <- RunPCA(scRNA, assay = "RNA")


pdf("01.pc_elbowplot.pdf", width = 7, height = 4, useDingbats = F)
ElbowPlot(scRNA, ndims = 40)
dev.off()

scRNA <- RunUMAP(scRNA, assay = "RNA", dims = 1:6, reduction.name = "umap_pc6")
scRNA <- RunUMAP(scRNA, assay = "RNA", dims = 1:9, reduction.name = "umap_pc9")
scRNA <- RunUMAP(scRNA, assay = "RNA", dims = 1:15, reduction.name = "umap_pc15")

scRNA <- RunTSNE(scRNA, assay = "RNA", dims = 1:15, reduction.name = "tsne_pc15")

p1 <- DimPlot(scRNA, reduction = "umap_pc6", label = F, group.by = "orig.ident") + scale_color_manual(values = sample.color)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "umap1", y = "umap2", title = "pc 1:6")

p2 <- DimPlot(scRNA, reduction = "umap_pc9", label = F, group.by = "orig.ident") + scale_color_manual(values = sample.color)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank())+
  labs(x = "umap1", y = "umap2", title = "pc 1:9")

p3 <- DimPlot(scRNA, reduction = "umap_pc15", label = F, group.by = "orig.ident") + scale_color_manual(values = sample.color)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "umap1", y = "umap2", title = "pc 1:15")

pdf("01.pc_umap.nPC.pdf", width = 24, height = 8, useDingbats = F)
cowplot::plot_grid(p1,p2,p3, nrow = 1)
dev.off()


d <- cbind(scRNA@meta.data, scRNA@reductions$umap_pc15@cell.embeddings, scRNA@reductions$tsne_pc15@cell.embeddings)

p3 <- ggplot(d, aes(x = umappc15_1, y = umappc15_2, color = orig.ident))+
  #DimPlot(scRNA, reduction = "umap_pc15", label = F, group.by = "orig.ident") + 
  geom_point(size = .1)+
  scale_color_manual(values = sample.color)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "umap1", y = "umap2", title = "pc 1:15")+
  facet_wrap(~orig.ident)
p3

pdf("01.umap-sample.pdf", width = 10, height = 10, useDingbats = F)
p3
dev.off()


p1 <- DimPlot(scRNA, reduction = "tsne_pc15", label = F, group.by = "orig.ident") + scale_color_manual(values = sample.color)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "tsne1", y = "tsne2", title = "pc 1:15")
p1
p2 <- ggplot(d, aes(x = tSNE_1, y = tSNE_2, color = orig.ident))+
  #DimPlot(scRNA, reduction = "umap_pc15", label = F, group.by = "orig.ident") + 
  geom_point(size = .1)+
  scale_color_manual(values = sample.color)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "tsne1", y = "tsne2", title = "pc 1:15")+
  facet_wrap(~orig.ident)
p2

pdf("01.tsne-sample.pdf", width = 10, height = 10, useDingbats = F)
p1
p2
dev.off()


rm(p1, p2, p3)
scRNA <- FindNeighbors(scRNA, reduction = "pca", assay = "RNA", dims = 1:15)
scRNA <- FindClusters(scRNA, resolution = 0.4, algorithm = 4, cluster.name = "pca_cluster")


p1 <- DimPlot(scRNA, reduction = "umap_pc15", label = T, group.by = "pca_cluster") + 
  scale_color_brewer(palette = "Paired") +
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2", title = "")
p1

p2 <- DimPlot(scRNA, reduction = "tsne_pc15", label = T, group.by = "pca_cluster") + 
  scale_color_brewer(palette = "Paired") +
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "tsne1", y = "tsne2", title = "")

pdf("01.cluster-umap-tsne.pdf", width = 16, height = 8, useDingbats = F)
cowplot::plot_grid(p1, p2)
dev.off()

scRNA <- JoinLayers(scRNA)
table(Idents(scRNA))

mk <- FindAllMarkers(scRNA, assay = "RNA")
head(mk)
write.csv(mk, file = "01.RNA.pc15.cluster.res04.markers.csv", quote = F, row.names = F)

###GO
library(clusterProfiler)
mouse.gene.set <- msigdbr::msigdbr(species = "mouse", category = "C5")  %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
table(mk$cluster)

lapply(as.character(unique(mk$cluster)), function(k){
  x = ''
  x <- enricher(gene = mk$gene[mk$cluster==k], TERM2GENE = mouse.gene.set)
  write.csv(x, file = paste0("01.RNA.pc15.cluster",k, ".GO.csv"), row.names = F)
})


###harmony
library(harmony)
DefaultAssay(scRNA) <- "RNA"
scRNA <- RunHarmony(scRNA, "orig.ident")
ElbowPlot(scRNA, reduction = "harmony", ndims = 50)

scRNA <- FindNeighbors(scRNA, reduction = "harmony", assay = "RNA")
scRNA <- RunUMAP(scRNA, assay = "RNA", dims = 1:6, reduction = "harmony", reduction.name = "harmonyUMAP6") 
scRNA <- RunUMAP(scRNA, assay = "RNA", dims = 1:9, reduction = "harmony", reduction.name = "harmonyUMAP9")
scRNA <- RunUMAP(scRNA, assay = "RNA", dims = 1:15, reduction = "harmony", reduction.name = "harmonyUMAP15")
scRNA <- RunTSNE(scRNA, assay = "RNA", dims = 1:15, reduction = "harmony", reduction.name = "harmonyTSNE")

scRNA@reductions$harmonyUMAP6
p1 = DimPlot(scRNA, reduction = "harmonyUMAP6", group.by = "orig.ident") + scale_color_manual(values = sample.color) + 
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank())+
  labs(x = "umap1", y = "umap2", title = "Harmony PC 1:6")
    


p2 = DimPlot(scRNA, reduction = "harmonyUMAP9", group.by = "orig.ident") + scale_color_manual(values = sample.color) + 
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2", title = "Harmony PC 1:9")

p3 = DimPlot(scRNA, reduction = "harmonyUMAP15", group.by = "orig.ident") + scale_color_manual(values = sample.color) + 
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2", title = "Harmony PC 1:15")
 
pdf("02.harmoneyUMAP.pcComparison.pdf", width = 15, height = 5, useDingbats = F)
cowplot::plot_grid(p1, p2, p3, nrow = 1)
dev.off()

rm(p1, p2, p3)


DimPlot(scRNA, reduction = "harmonyUMAP15", group.by = "pca_cluster") + 
  scale_color_brewer(palette = "Paired")+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2", title = "Harmony PC 1:15")

###
scRNA <- FindNeighbors(scRNA, reduction = "harmony", assay = "RNA")
scRNA <- FindClusters(scRNA, resolution = 0.4, algorithm = 4, cluster.name = "harmony_cluster")

DimPlot(scRNA, reduction = "harmonyUMAP15", label = T, group.by = "harmony_cluster") +
  scale_color_brewer(palette = "Paired")

DimPlot(scRNA, reduction = "harmonyTSNE", label = T, group.by = "harmony_cluster") +
  scale_color_brewer(palette = "Paired")



###Funcitonal clusters?
scRNA@meta.data$cluster <- as.character(scRNA@meta.data$pca_cluster)
DimPlot(scRNA, reduction = "umap_pc15", label = T, group.by = "cluster") + scale_color_brewer(palette = "Paired")

scRNA@meta.data$cluster[scRNA@meta.data$cluster %in% c(4,8)] <- "umapR"
scRNA@meta.data$cluster[scRNA@meta.data$cluster %in% c(7)] <- "umapL"
scRNA@meta.data$cluster[!scRNA@meta.data$cluster %in% c("umapR","umapL")] <- "others"

p <- DimPlot(scRNA, reduction = "umap_pc15", label = T, group.by = "cluster") + 
  scale_color_brewer(palette = "Paired")+
  labs(x= "umap1", y = "umap2")
pdf("05.funciton.clusters.umap.pdf", width = 5, height = 5, useDingbats = F)
p
dev.off()

Idents(scRNA) <- factor(scRNA@meta.data$cluster)

markers <- FindAllMarkers(scRNA)
head(markers)
write.csv(markers, file = "05.funciton.clusters.markers.csv", row.names = F, quote = F)



###
######################################################################
###add clonotype
tcr <- read.csv("../data/vdj_t/filtered_contig_annotations.csv")
tail(tcr, 2)

# Remove the -1 at the end of each barcode.
# Subsets so only the first line of each barcode is kept,
# as each entry for given barcode will have same clonotype.
tcr <- tcr[!duplicated(tcr$barcode), ]
head(tcr)
unique(tcr$donor)
unique(tcr$origin)
tcr$barcode[1:4]
tcr$origin[1:4]
#
tcr$barcode <- gsub("\\-\\d$", "-1", tcr$barcode)
tcr$barcode2 <-paste0(tcr$origin, "_", tcr$barcode)
tail(tcr$barcode2)

# Only keep the barcode and clonotype columns.
# We'll get additional clonotype info from the clonotype table.
tcr <- tcr[,c("barcode2", "raw_clonotype_id")]
names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
# Clonotype-centric info.
clono <- read.csv("../data/vdj_t/clonotypes.csv")
head(clono)

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
sum(!rownames(tcr) %in% colnames(scRNA)) #3370
sum(!colnames(scRNA) %in% rownames(tcr)) #6064
sum(rownames(tcr) %in% colnames(scRNA)) #13845
sum(colnames(scRNA) %in% rownames(tcr)) #13845

rownames(tcr)[!rownames(tcr) %in% colnames(scRNA)][1:4]
###sc filtered that though have TCR, but the cell is not valid

# Add to the Seurat object's metadata.
scRNA <- AddMetaData(object=scRNA, metadata=tcr)

#data.integrated <- AddMetaData(data.integrated, metadata = tcr)
save(scRNA, file = "01.seurat.obj.rda")
#table(scRNA@meta.data$orig.ident)


}

load("01.seurat.obj.rda")
tableau.col <- as.character(paletteer::paletteer_d("ggthemes::Tableau_10"))


##
grep("HPS", rownames(scRNA), ignore.case = T, value = T)

p1 <- FeaturePlot(scRNA, features = c("Gzmb", "Ifng", "Prf1", "Tnf"), reduction = "umap_pc15") & theme(aspect.ratio = 1)
#CXCR5 not found
p2 <- FeaturePlot(scRNA, features = c("Pdcd1", "Tox", "Tcf7", "Lag3", "Ctla4", "Tigit"), reduction = "umap_pc15") & theme(aspect.ratio = 1)
p2
#应激应答
hspgenes <- grep("HPS", rownames(scRNA), ignore.case = T, value = T)
p3 <- FeaturePlot(scRNA, features = c(hspgenes, "Brca1", "Brca2", "Atm"), reduction = "umap_pc15") & theme(aspect.ratio = 1)
p3
#代谢
p4 <- FeaturePlot(scRNA, features = c("Srebf1", "Hmgcr", "Cpt1a", "Acat1"), reduction = "umap_pc15")& theme(aspect.ratio = 1)
p4
ggsave(p1, file = "02.ImmuneActivity.pdf", width = 10, height = 10, useDingbats = F)
ggsave(p2, file = "02.Stemness.pdf", width = 10, height = 10, useDingbats = F)
ggsave(p3, file = "02.StressResponse.pdf", width = 10, height = 10, useDingbats = F)
ggsave(p4, file = "02.Metabolism.pdf", width = 10, height = 10, useDingbats = F)

##CytoTRACE
GeneCounts <- as.matrix(scRNA[['RNA']]$counts)
iOrd <- rowSums(GeneCounts>0)
sum(iOrd>100)
GeneCounts <- GeneCounts[iOrd>100,] # only keep genes expressing in more than 100 cell
CytoTRACE.score <- CytoTRACE(GeneCounts, enableFast = TRUE, ncores = 5, subsamplesize = 1000)
#
plotCytoTRACE(CytoTRACE.score)
scRNA$CytoTRACE <- CytoTRACE.score$CytoTRACE[colnames(scRNA)]

temp <- RColorBrewer::brewer.pal(11, "Spectral")
temp[6] <- "gold"
rbPal <- colorRampPalette(temp)


data <- cbind(scRNA@reductions$umap_pc15@cell.embeddings, scRNA[[]])

p1 <- ggplot2::ggplot(data, ggplot2::aes(x=data[,1], y=data[,2])) +
  ggplot2::geom_point(ggplot2::aes(colour = CytoTRACE), size = 0.5)+
  ggplot2::scale_colour_gradientn(name = "Predicted\norder",
                                  colours = rev(rbPal(50)),
                                  guide = ggplot2::guide_colourbar(ticks.colour = "black",
                                                                   ticks.linewidth = 1,
                                                                   frame.colour = "black"),
                                  breaks = seq(0, 1, 0.2),
                                  labels=c("0.0 (More diff.)", 0.2,0.4, 0.6, 0.8, "1.0 (Less diff.)"))+
  ggplot2::labs(x = colnames(data)[1], y = colnames(data)[2], title = "CytoTRACE")+
  ggpubr::theme_pubr()+
  ggplot2::theme(
    aspect.ratio = 1,
    legend.text = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 14),
    plot.title = ggplot2::element_text(size = 21, hjust = 0.5),
    axis.title.x = ggplot2::element_text(size = 18),
    axis.title.y = ggplot2::element_text(size = 18),
    #axis.text=ggplot2::element_text(size=16),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position="left",
    plot.margin = ggplot2::unit(c(0.5,1,0.5,1), "cm")
    
  )+
  labs(x = "umap", y = "")
p1
ggsave("02.CytoTRACE.umap_pc15.pdf", plot = p1, width = 10, height = 10, useDingbats = F)




#plot cell states
cell.state <- read.csv("03.cellState.prediction.csv", head =T, row.names = 1)
cell.state <- cell.state[,c("predicted.celltype"), drop = F]
head(cell.state)

scRNA <- AddMetaData(scRNA, cell.state)

p <- DimPlot(scRNA, reduction = "umap_pc15", group.by = "predicted.celltype") + 
  scale_color_manual(values = tableau.col[c(1,3,5,7,9)]) + 
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2")
pdf("03.cellState.prediction.umap_pc15.pdf", width = 5, height = 5, useDingbats = F)
p
dev.off()

###clonotype
#ElbowPlot(scRNA, reduction = "harmony")
#scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:10)
#scRNA <- FindClusters(scRNA, resolution = 0.8)

head(scRNA@meta.data)

DimPlot(scRNA, reduction = "umap_pc15", group.by = "T_clonotype_id") + 
  #scale_color_manual(values = sample.color) + 
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2")


keep.clono <- names(which(table(scRNA@meta.data$T_clonotype_id) > 10))
keep.clono

clono.col = as.character(paletteer_d("ggthemes::Tableau_20"))[1:length(keep.clono)]
names(clono.col) = keep.clono
clono.col

scRNA <- subset(scRNA, subset = T_clonotype_id %in% keep.clono)

xx = cbind(scRNA@reductions$umap_pc15@cell.embeddings, 
           scRNA@reductions$harmonyUMAP15@cell.embeddings,
           scRNA@meta.data)

colnames(xx)

xx$sampleID <- factor(xx$orig.ident, levels = c("6L", "6H", "11L", "11H"))
xx$T_clonotype_id <- factor(xx$T_clonotype_id, levels = paste0("clonotype", 1:length(keep.clono)))



###
p1 = ggplot(xx, aes(x =umappc15_1, y = umappc15_2, color = T_clonotype_id)) +
  geom_point(size = .5)+
  scale_colour_manual(values = clono.col)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2", title = "")+
  labs(title = "Merge")
p1

p2 = ggplot(xx, aes(x =harmonyUMAP15_1, y = harmonyUMAP15_2, color = T_clonotype_id)) +
  geom_point(size = .5)+
  scale_colour_manual(values = clono.col)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(x = "umap1", y = "umap2", title = "")+
  labs(title = "Harmony")
p1
p2

pdf("04.umap.clonotype.pdf", width = 10, height = 10, useDingbats = F)
cowplot::plot_grid(p1, p2, nrow = 2)
p1 + facet_wrap(~sampleID)
p2 + facet_wrap(~sampleID)
dev.off()
#scRNA <- data.integrated
####scR

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

pdf("04.scRNA.colontypeDiversity.pdf", width = 5, height = 4, useDingbats = F)
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
  geom_text(aes(x = Var1, y = 10, label =Freq, angle = 90))+
  labs(x= "", y = "count")+
  scale_y_log10(expand = c(0,0))
p2

egg::ggarrange(p2, p1, ncol = 1, nrow = 2, heights = c(5, 10))

pdf("04.clontype.sample.freq.pdf", width = 8, height = 8, useDingbats = F)
egg::ggarrange(p2, p1, ncol = 1, nrow = 2, heights = c(5, 40))
dev.off()


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


y4 = data.frame(table(scRNA@meta.data$orig.ident, scRNA@meta.data$T_clonotype_id)/rowSums(table(scRNA@meta.data$orig.ident, scRNA@meta.data$T_clonotype_id)))

p4 = ggplot(y4, aes(x = "", y = Freq, fill = Var2))+
  geom_bar(stat = "identity", width = 1)+
  scale_fill_manual(values = clono.col)+
  theme_bw()+ 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 12))+
  labs(x= "", y = "count")+
  coord_polar("y")+
  facet_wrap(~Var1, strip.position = "right", nrow = 1)

p4


pdf("04.sample-clonotype.freq.pdf", width = 15, height = 5, useDingbats = F)
p3
p4
dev.off()

colnames(y) <- c("clono", "sample", "ratio", "N")
write.csv(y, file = "04.sample-clontype.freq.csv", quote = F, row.names = F)
###
###

head(scRNA[[]])
table(scRNA@meta.data$T_clonotype_id)
plot(density(table(scRNA@meta.data$T_clonotype_id)))
table(scRNA@meta.data$T_clonotype_id, scRNA@meta.data$orig.ident)


####
#DEG clonotype1
head(scRNA[[]])
scRNA1 = JoinLayers(scRNA, assay = "RNA")

##
scRNA1@meta.data$avidity_plus_clonotype <- paste0(scRNA1@meta.data$orig.ident, "_", scRNA1@meta.data$T_clonotype_id)
unique(scRNA1[[]]$avidity_plus_clonotype)
Idents(scRNA1) <- scRNA1@meta.data$avidity_plus_clonotype

###
clonotype1_11Hvsclonotype2_11H = FindMarkers(scRNA1, ident.1 = "11H_clonotype1", ident.2 = "11H_clonotype2", logfc.threshold = 0, assay = "RNA")
clonotype1_11Lvsclonotype2_11L = FindMarkers(scRNA1, ident.1 = "11L_clonotype1", ident.2 = "11L_clonotype2", logfc.threshold = 0, assay = "RNA")
clonotype1_6Hvsclonotype2_6H = FindMarkers(scRNA1, ident.1 = "6H_clonotype1", ident.2 = "6H_clonotype2", logfc.threshold = 0, assay = "RNA")
clonotype1_6Lvsclonotype2_6L = FindMarkers(scRNA1, ident.1 = "6L_clonotype1", ident.2 = "6L_clonotype2", logfc.threshold = 0, assay = "RNA")
#
x1 <- average_expression <- AverageExpression(scRNA1, 
                             assays = "RNA", slot = "data")[[1]][,c("g11H-clonotype1", "g11H-clonotype2", "g11L-clonotype1", "g11L-clonotype2")]
#
colnames(clonotype1_11Hvsclonotype2_11H) <- paste0("11H_clonotype1vs2", c("p_val", "avg_logFC", "pct1", "pct2", "p_val_adj"))
colnames(clonotype1_11Lvsclonotype2_11L) <- paste0("11L_clonotype1vs2", c("p_val", "avg_logFC", "pct1", "pct2", "p_val_adj"))
colnames(clonotype1_6Hvsclonotype2_6H) <- paste0("6H_clonotype1vs2", c("p_val", "avg_logFC", "pct1", "pct2", "p_val_adj"))
colnames(clonotype1_6Lvsclonotype2_6L) <- paste0("6L_clonotype1vs2", c("p_val", "avg_logFC", "pct1", "pct2", "p_val_adj"))

###
clonotype1.deg.11HvsL = FindMarkers(scRNA1, ident.1 = "11H_clonotype1", ident.2 = "11L_clonotype1", logfc.threshold = 0, assay = "RNA")
clonotype1.deg.6HvsL = FindMarkers(scRNA1, ident.1 = "6H_clonotype1", ident.2 = "6L_clonotype1", logfc.threshold = 0, assay = "RNA")
colnames(clonotype1.deg.11HvsL) <- paste0("clonotype1_11Hvs11L", c("p_val", "avg_logFC", "pct1", "pct2", "p_val_adj"))
colnames(clonotype1.deg.6HvsL) <- paste0("clonotype1_6Hvs6L", c("p_val", "avg_logFC", "pct1", "pct2", "p_val_adj"))
###
clonotype2.deg.11LH = FindMarkers(scRNA1, ident.1 = "11H_clonotype2", ident.2 = "11L_clonotype2", logfc.threshold = 0, assay = "RNA")
clonotype2.deg.6LH = FindMarkers(scRNA1, ident.1 = "6H_clonotype2", ident.2 = "6L_clonotype2", logfc.threshold = 0, assay = "RNA")
colnames(clonotype2.deg.11LH) <- paste0("clonotype2_11Hvs11L", c("p_val", "avg_logFC", "pct1", "pct2", "p_val_adj"))
colnames(clonotype2.deg.6LH) <- paste0("clonotype2_6Hvs6L", c("p_val", "avg_logFC", "pct1", "pct2", "p_val_adj"))
#
deg <- cbind(x1, clonotype1_11Hvsclonotype2_11H[rownames(x1), ], clonotype1_11Lvsclonotype2_11L[rownames(x1), ], 
             clonotype1_6Hvsclonotype2_6H[rownames(x1), ], clonotype1_6Lvsclonotype2_6L[rownames(x1), ],
             clonotype1.deg.11HvsL[rownames(x1), ], clonotype1.deg.6HvsL[rownames(x1), ],
             clonotype2.deg.11LH[rownames(x1), ], clonotype2.deg.6LH[rownames(x1), ])

write.csv(deg, file = "05.deg.clonotype12.HL.csv", quote = F, row.names = T)

##
#avidity11H.clonotype1vs2 <- deg[which(deg$`11H_clonotype1vs2:avg_logFC` > 0.25 & deg$`11H_clonotype1vs2:p_val` < 0.05),]
#avidity11L.clonotype1vs2 <- deg[which(deg$`11L_clonotype1vs2:avg_logFC` >0.25 & deg$`11L_clonotype1vs2:p_val` <0.05),]
#avidity6H.clonotype1vs2 <- deg[which(deg$`6H_clonotype1vs2:avg_logFC` & deg$`6H_clonotype1vs2:p_val` < 0.05),]
#avidity6L.clonotype1vs2 <- deg[which(deg$`6L_clonotype1vs2:avg_logFC` & deg$`6L_clonotype1vs2:p_val` < 0.05),]
library(clusterProfiler)
mouse.gene.set <- msigdbr::msigdbr(species = "mouse", category = "C5")  %>% dplyr::distinct(gs_name, gene_symbol, gs_exact_source) %>% as.data.frame()

#GO terms from HT
if(1){
  go_ids <- list(
    T_CELL_ACTIVATION_STRENGTH = c(
      "GO:0070371", "GO:0042110", "GO:0046631", "GO:0050870", "GO:0004715",
      "GO:0042101", "GO:0050852", "GO:0050851", "GO:0050863", "GO:0019900"
    ),
    
    CALCIUM_SIGNALING_AND_ION_REGULATION = c(
      "GO:0051279", "GO:0050849", "GO:0050848", "GO:0019722", "GO:1905664",
      "GO:0051928", "GO:0006816", "GO:0055074"
    ),
    
    CYTOKINE_NETWORKS = c(
      "GO:0019221", "GO:0001819", "GO:0001817", "GO:0005125", "GO:0005126",
      "GO:0038110", "GO:0035723", "GO:0140888", "GO:0033209"
    ),
    
    CELL_PROLIFERATION_AND_DNA_REPLICATION = c(
      "GO:1902969", "GO:0046641", "GO:0006260", "GO:0006270", "GO:0022616",
      "GO:0017116", "GO:0000082", "GO:0044839", "GO:0051781", "GO:0010564"
    ),
    
    CHROMOSOME_ORGANIZATION = c(
      "GO:0051276", "GO:0006325", "GO:0000793", "GO:0007076", "GO:0098813",
      "GO:0051225", "GO:0010965"
    ),
    
    TERMINAL_EFFECTOR_DIFFERENTIATION = c(
      "GO:0045595", "GO:0046638", "GO:0030217", "GO:0045596", "GO:0046632", "GO:0033077"
    ),
    
    APOPTOSIS_SENSITIVITY = c(
      "GO:0042981", "GO:0070232", "GO:0043065", "GO:0097193", "GO:0097191",
      "GO:0006915", "GO:0043066"
    ),
    
    MEMBRANE_ORGANIZATION_AND_SYNAPSE_FORMATION = c(
      "GO:0045121", "GO:0009898", "GO:0001771", "GO:0044853", "GO:0001765", "GO:0061024"
    ),
    
    CYTOSKELETAL_ORGANIZATION = c(
      "GO:0099513", "GO:0015629", "GO:0061640", "GO:0008017", "GO:0030864", "GO:0007015"
    ),
    
    METABOLISM = c(
      "GO:0006096", "GO:0033539", "GO:0046034", "GO:0006119", "GO:0006099",
      "GO:0045333", "GO:0006629", "GO:0008203", "GO:0090322", "GO:0071456"
    ),
    
    EPIGENETICS_AND_CHROMATIN = c(
      "GO:0006325", "GO:0006338", "GO:0035035", "GO:0042054", "GO:0040029",
      "GO:0006334", "GO:0031507", "GO:0006304"
    ),
    
    TRANSCRIPTIONAL_CONTROL = c(
      "GO:0010468", "GO:0045944", "GO:0001228", "GO:0000976",
      "GO:0010467", "GO:0010628", "GO:0032774"
    ),
    
    STRESS_RESPONSES = c(
      "GO:0006974", "GO:0034599", "GO:0034976", "GO:0140467", "GO:1903894",
      "GO:0000423", "GO:0010508", "GO:0062197"
    ),
    
    CHEMOKINE_AND_MIGRATION = c(
      "GO:0019956", "GO:0008009", "GO:0140131", "GO:0010818", "GO:0048020", "GO:0016477"
    ),
    
    NEGATIVE_REGULATION_AND_FEEDBACK = c(
      "GO:0050868", "GO:0032088", "GO:1902532", "GO:0001818", "GO:0048523"
    )
  )
}

p_cutoff <- 0.01
logFC_cutoff <- 0.1
if(1){
avidity11H.clonotype1vs2 <- deg[which(deg$`11H_clonotype1vs2avg_logFC` > logFC_cutoff & deg$`11H_clonotype1vs2p_val` < p_cutoff),]
avidity11L.clonotype1vs2 <- deg[which(deg$`11L_clonotype1vs2avg_logFC` > logFC_cutoff & deg$`11L_clonotype1vs2p_val` <p_cutoff),]
avidity6H.clonotype1vs2 <- deg[which(deg$`6H_clonotype1vs2avg_logFC`> logFC_cutoff & deg$`6H_clonotype1vs2p_val` < p_cutoff),]
avidity6L.clonotype1vs2 <- deg[which(deg$`6L_clonotype1vs2avg_logFC`>logFC_cutoff & deg$`6L_clonotype1vs2p_val` < p_cutoff),]

nrow(avidity11H.clonotype1vs2)
nrow(avidity11L.clonotype1vs2)
nrow(avidity6H.clonotype1vs2)
nrow(avidity6L.clonotype1vs2)

clonotype1.11H <- deg[which(deg$clonotype1_11Hvs11Lavg_logFC > logFC_cutoff & deg$clonotype1_11Hvs11Lp_val < p_cutoff),]
clonotype1.6H <- deg[which(deg$clonotype1_6Hvs6Lavg_logFC > logFC_cutoff & deg$clonotype1_6Hvs6Lp_val < p_cutoff),]
#view(clonotype1.6H)
clonotype2.11H <- deg[which(deg$clonotype2_11Hvs11Lavg_logFC > logFC_cutoff & deg$clonotype2_11Hvs11Lp_val < p_cutoff),]
clonotype2.6H <- deg[which(deg$clonotype2_6Hvs6Lavg_logFC > logFC_cutoff & deg$clonotype2_6Hvs6Lp_val < p_cutoff),]
nrow(clonotype1.11H)
nrow(clonotype1.6H)
nrow(clonotype2.11H)
nrow(clonotype2.6H)

library(ggVennDiagram)
p1 <- ggVennDiagram(list(clonotype1.11H = rownames(clonotype1.11H), 
                         clonotype2.11H = rownames(clonotype2.11H),
                         clonotype1.6H = rownames(clonotype1.6H),
                         clonotype2.6H = rownames(clonotype2.6H)
                         #avidity11H.clonotype1vs2 = rownames(avidity11H.clonotype1vs2),
                         #avidity11L.clonotype1vs2 = rownames(avidity11L.clonotype1vs2),
                        ),
                    category.names = c("clonotype1\n11Hvs11L", 
                                       "clonotype2\n11Hvs11L",
                                       "clonotype1\n6Hvs6L",
                                       "clonotype2\n6Hvs6L"
                                       #"avidity11H:clonotype1vs2",
                                       #"avidity11L:clonotype1vs2"
                                       ),
                    label_alpha = 0)+
  scale_fill_distiller(palette = "Reds", direction = 1)+
  labs(title = "High avidity overexpressed across clonotype 1/2")
#p1

ggsave(p1, file = paste0("05.sameClono.highAvidityOverexpressed_p",p_cutoff, "logFC", logFC_cutoff, "_venn.pdf"), width = 6, height = 5, useDingbats = F)

p1.1 <- ggVennDiagram(list(#clonotype1.11H = rownames(clonotype1.11H), 
                         #clonotype2.11H = rownames(clonotype2.11H),
                         #clonotype1.6H = rownames(clonotype1.6H),
                         #clonotype2.6H = rownames(clonotype2.6H)),
                         avidity11H.clonotype1vs2 = rownames(avidity11H.clonotype1vs2),
                         avidity11L.clonotype1vs2 = rownames(avidity11L.clonotype1vs2),
                         avidity6H.clonotype1vs2 = rownames(avidity11H.clonotype1vs2),
                         avidity6L.clonotype1vs2 = rownames(avidity11L.clonotype1vs2)),
                    category.names = c(#"clonotype1:11H", 
                                       #"clonotype2:11H",
                                       #"clonotype1:6H", 
                                       #"clonotype2:6H",
                                       "stimu. 11H\nclonotype1vs2",
                                       "stimu. 11L\nclonotype1vs2",
                                       "stimu. 6H\nclonotype1vs2",
                                       "stimu. 6L\nclonotype1vs2"
                                       ),
                    label_alpha = 0)+
  scale_fill_distiller(palette = "Reds", direction = 1)+
  labs(title = "Clonotype1 overexpressed vs Clonotype2 across different stimulations")
#p1.1
ggsave(p1.1, file = paste0("06.sameAvidity_clono1Overexpressed_p", p_cutoff, "logFC", logFC_cutoff,"_venn.pdf"), width = 5, height = 5, useDingbats = F)

avidity11H.clonotype1vs2 <- deg[which(deg$`11H_clonotype1vs2avg_logFC` < -1*logFC_cutoff & deg$`11H_clonotype1vs2p_val` < p_cutoff),]
avidity11L.clonotype1vs2 <- deg[which(deg$`11L_clonotype1vs2avg_logFC` < -1*logFC_cutoff & deg$`11L_clonotype1vs2p_val` < p_cutoff),]
avidity6H.clonotype1vs2 <- deg[which(deg$`6H_clonotype1vs2avg_logFC` < -1*logFC_cutoff & deg$`6H_clonotype1vs2p_val` < p_cutoff ),]
avidity6L.clonotype1vs2 <- deg[which(deg$`6L_clonotype1vs2avg_logFC` < -1*logFC_cutoff & deg$`6L_clonotype1vs2p_val` < p_cutoff),]
nrow(avidity11H.clonotype1vs2)
nrow(avidity11L.clonotype1vs2)
nrow(avidity6H.clonotype1vs2)
nrow(avidity6L.clonotype1vs2)

clonotype1.11L <- deg[which(deg$clonotype1_11Hvs11Lavg_logFC < -1*logFC_cutoff & deg$clonotype1_11Hvs11Lp_val <p_cutoff ),]
clonotype1.6L <- deg[which(deg$clonotype1_6Hvs6Lavg_logFC < -1*logFC_cutoff & deg$clonotype1_6Hvs6Lp_val < p_cutoff),]
clonotype2.11L <- deg[which(deg$clonotype2_11Hvs11Lavg_logFC < -1*logFC_cutoff & deg$clonotype2_11Hvs11Lp_val < p_cutoff),]
clonotype2.6L <- deg[which(deg$clonotype2_6Hvs6Lavg_logFC < -1*logFC_cutoff & deg$clonotype2_6Hvs6Lp_val < p_cutoff),]
nrow(clonotype1.11L)
nrow(clonotype1.6L)
nrow(clonotype2.11L)
nrow(clonotype2.6L)

p2 <- ggVennDiagram(list(clonotype1.11L = rownames(clonotype1.11L), 
                         clonotype2.11L = rownames(clonotype2.11L),
                         clonotype1.6L = rownames(clonotype1.6L),
                         clonotype2.6L = rownames(clonotype2.6L)),
                    category.names = c("clonotype1\n11Lvs11H",
                                       "clonotype2\n11Lvs11H",
                                       "clonotype1\n6Lvs6H", 
                                       "clonotype2\n6Lvs6H"),
                    
            label_alpha = 0)+
  scale_fill_distiller(palette = "Blues", direction = 1)+
  labs(title = "Low avidity overexpressed across clonotype 1/2")
#p2
ggsave(p2, file = paste0("05.sameClono_lowAvidityOverexpressed_p", p_cutoff, "logFC", logFC_cutoff,".venn.pdf"), width = 5, height = 5, useDingbats = F)

p2.1 <- ggVennDiagram(list(
                         #clonotype1.11L = rownames(clonotype1.11L), 
                         #clonotype2.11L = rownames(clonotype2.11L),
                         #clonotype1.6L = rownames(clonotype1.6L),
                         #clonotype2.6L = rownames(clonotype2.6L)),
                         avidity11H.clonotype1vs2 = rownames(avidity11L.clonotype1vs2),
                         avidity11L.clonotype1vs2 = rownames(avidity11H.clonotype1vs2),
                         avidity6H.clonotype1vs2 = rownames(avidity6H.clonotype1vs2),
                         avidity6L.clonotype1vs2 = rownames(avidity6L.clonotype1vs2)),
                    category.names = c(#"clonotype1:11Lvs11H", 
                                       #"clonotype2:11Lvs11H",
                                       #"clonotype1:6Lvs6H", 
                                       #"clonotype2:6Lvs6H"),
                                       "stimu. 11H\nclonotype1vs2",
                                       "stimu. 11L\nclonotype1vs2",
                                       "stimu. 6H\nclonotype1vs2",
                                       "stimu. 6L\nclonotype1vs2"),
                    label_alpha = 0)+
  scale_fill_distiller(palette = "Blues", direction = 1)+
  labs(title = "Clonotype2 overexpressed vs Clonotype1 across different stimulations")
#p2.1

ggsave(p2.1, file = paste0("06.sameAvidity_clono2higher_avidity_p",p_cutoff, "_logFC",logFC_cutoff,"_venn.pdf"), width = 5, height = 5, useDingbats = F)



##Go enrichment


length(unique(unlist(go_ids))) #112
sum(unique(unlist(go_ids)) %in% mouse.gene.set$gs_exact_source)  #90

#1. clonotype1 11H vs 11L
#2. clonotype1 6H vs 6L
#3. clonotype2 11H vs 11L
#4. clonotype2 6H vs 6L
compare.groups <- list("clonotype1_11Hvs11L" = c("clonotype1_11Hvs11Lavg_logFC", "clonotype1_11Hvs11Lp_val"),
                       "clonotype1_6Hvs6L" = c("clonotype1_6Hvs6Lavg_logFC", "clonotype1_6Hvs6Lp_val"),
                       "clonotype2_11Hvs11L" = c("clonotype2_11Hvs11Lavg_logFC", "clonotype2_11Hvs11Lp_val"),
                       "clonotype2_6Hvs6L" = c("clonotype2_6Hvs6Lavg_logFC", "clonotype2_6Hvs6Lp_val"))

up.go <- lapply(1:length(compare.groups), function(k){
  #x = ''
  use.group <- names(compare.groups)[k]
  message(paste0("Processing group: ", use.group))
  use.genes1 <- rownames(deg)[which(deg[, compare.groups[[k]][1]] > logFC_cutoff & deg[, compare.groups[[k]][2]] < p_cutoff)]
  message("N over expressed in high stimulation:", length(use.genes1))
  x1 <- enricher(gene = use.genes1, TERM2GENE = mouse.gene.set, pvalueCutoff = 1)@result
  x1$GOID <- mouse.gene.set$gs_exact_source[match(x1$ID, mouse.gene.set$gs_name)]
  x1$group <- use.group
  x1$direction <- "up"
  #return(x1)
  #
  use.genes2 <- rownames(deg)[which(deg[, compare.groups[[k]][1]] < -1*logFC_cutoff & deg[, compare.groups[[k]][2]] < p_cutoff)]
  message("N under expressed:", length(use.genes2))
  x2 <- enricher(gene = use.genes2, TERM2GENE = mouse.gene.set, pvalueCutoff = 1)@result
  x2$GOID <- mouse.gene.set$gs_exact_source[match(x2$ID, mouse.gene.set$gs_name)]
  x2$group <- use.group
  x2$direction <- "down"
  x <- rbind(x1, x2)
  x
})

sum(up.go[[1]]$p.adjust < 0.05)
sum(up.go[[2]]$p.adjust < 0.05)
sum(up.go[[3]]$p.adjust < 0.05)
sum(up.go[[4]]$p.adjust < 0.05)

up.go.top200.GOID <- lapply(up.go, function(x){
  #x[order(x$p.adjust),][1:100,]$GOID
  x[x$p.adjust < 0.05,]$GOID[1:50]
})
up.go.top200.GOID <- unique(unlist(up.go.top200.GOID))

sum(unique(unlist(go_ids)) %in% up.go.top200.GOID) #18 (top100), #41 (any significant) #top50: 15
sum(unique(unlist(go_ids)) %in% mouse.gene.set$gs_exact_source) #90

###
library(viridis)
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

up.go.df <- do.call(rbind, up.go)
head(up.go.df)
up.go.df <- up.go.df[up.go.df$direction == "up",]
#pivot wider to cluster
up.go.df.wider <- data.frame(tidyr::pivot_wider(up.go.df[, c("GOID", "group", "p.adjust")], names_from = group, 
                                                values_from = p.adjust, values_fn = mean))
rownames(up.go.df.wider) <- up.go.df.wider[,1]
up.go.df.wider <- up.go.df.wider[,-1]
up.go.df.wider <- up.go.df.wider[rownames(up.go.df.wider) %in% up.go.top200.GOID,]
up.go.df.wider[is.na(up.go.df.wider)] <- 1
nrow(up.go.df.wider)
apply(up.go.df.wider, 2, function(x) sum(x<0.05 & !is.na(x)))


#
group.cluster <- hclust(as.dist(1-cor(up.go.df.wider)), method = "average")

GO.cluster <- hclust(dist((-log10(up.go.df.wider))), method = "ward.D2")

up.go.df <- up.go.df[up.go.df$GOID %in% up.go.top200.GOID,]
y.order <- data.frame(GOID= GO.cluster$labels[GO.cluster$order],
                      y = 1:length(GO.cluster$labels[GO.cluster$order]),
                      GOID2 = up.go.df$ID[match(GO.cluster$labels[GO.cluster$order], up.go.df$GOID)])
head(y.order)
up.go.df$y <- y.order$y[match(up.go.df$GOID, y.order$GOID)]

unique(up.go.df$group)

up.go.df$p.adjust2 <- ifelse(up.go.df$p.adjust < 1e-30 & !is.na(up.go.df$p.adjust), 1e-30, up.go.df$p.adjust)

p <- ggplot(up.go.df, aes(x = group, y = y, fill = -log10(p.adjust2))) +
  geom_tile(width = 1) +
  scale_fill_gradientn(colours = custom_magma, na.value = "white", limits = c(-log10(5e-2), 30))+
  #facet_wrap(~direction, ncol = 1, strip.position = "right", scales = "free") +
  #scale_y_discrete(limits = GO.cluster$labels[GO.cluster$order]) +
  scale_y_continuous(breaks = y.order$y, labels = y.order$GOID,
                     sec.axis = sec_axis(~., breaks = y.order$y, labels = y.order$GOID2),
                     expand = expansion(mult = c(0, 0.01)))+
  scale_x_discrete(limits = c("clonotype1_11Hvs11L", "clonotype2_11Hvs11L","clonotype1_6Hvs6L", "clonotype2_6Hvs6L"))+ 
                   #guide = guide_axis(n.dodge = 2)) + #limits = group.cluster$labels[group.cluster$order], 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 3), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  labs(x = "", y = "",fill = "-log10(p.adjust)")

p

pdf(paste0("05.sameClono_differentAvidityDEG_p",p_cutoff, "_logFC",logFC_cutoff,"_GOtop50.enrichment.pdf"), width = 8, height = 8, useDingbats = F)
print(p)
dev.off()

up.go.df.wider$yorder <- y.order$y[match(rownames(up.go.df.wider), y.order$GOID)]
write.csv(up.go.df.wider, file = paste0("05.sameClono_differentAvidityDEG_p", p_cutoff, "_logFC", logFC_cutoff, "_GOtop50.enrichment.csv"), quote = F, row.names = T)

###GO compare clonotype 1 vs 2 across different stimulation
rm(up.go, up.go.df, up.go.df.wider, up.go.top200.GOID, y.order)
compare.groups <- list("11H_clonotype1vs2" = c("11H_clonotype1vs2avg_logFC", "11H_clonotype1vs2p_val"),
                       "11L_clonotype1vs2" = c("11L_clonotype1vs2avg_logFC", "11L_clonotype1vs2p_val"),
                       "6H_clonotype1vs2" = c("6H_clonotype1vs2avg_logFC", "6H_clonotype1vs2p_val"),
                       "6L_clonotype1vs2" = c("6L_clonotype1vs2avg_logFC", "6L_clonotype1vs2p_val"))
up.go <- lapply(1:length(compare.groups), function(k){
  
  #x = ''
  use.group <- names(compare.groups)[k]
  message(paste0("Processing group: ", use.group))
  use.genes1 <- rownames(deg)[which(deg[, compare.groups[[k]][1]] > logFC_cutoff & deg[, compare.groups[[k]][2]] < p_cutoff)]
  message("N over expressed in high stimulation:", length(use.genes1))
  x1 <- enricher(gene = use.genes1, TERM2GENE = mouse.gene.set, pvalueCutoff = 1)@result
  x1$GOID <- mouse.gene.set$gs_exact_source[match(x1$ID, mouse.gene.set$gs_name)]
  x1$group <- use.group
  x1$direction <- "up"
  #return(x1)
  #
  use.genes2 <- rownames(deg)[which(deg[, compare.groups[[k]][1]] < -1*logFC_cutoff & deg[, compare.groups[[k]][2]] < p_cutoff)]
  message("N under expressed:", length(use.genes2))
  x2 <- enricher(gene = use.genes2, TERM2GENE = mouse.gene.set, pvalueCutoff = 1)@result
  x2$GOID <- mouse.gene.set$gs_exact_source[match(x2$ID, mouse.gene.set$gs_name)]
  x2$group <- use.group
  x2$direction <- "down"
  x <- rbind(x1, x2)
  x
})

sum(up.go[[1]]$p.adjust < 0.05)
sum(up.go[[2]]$p.adjust < 0.05)
sum(up.go[[3]]$p.adjust < 0.05)
sum(up.go[[4]]$p.adjust < 0.05)

up.go.top200.GOID <- lapply(up.go, function(x){
  #x[order(x$p.adjust),][1:100,]$GOID
  x[x$p.adjust < 0.05,]$GOID[1:50]
})
up.go.top200.GOID <- unique(unlist(up.go.top200.GOID))

sum(unique(unlist(go_ids)) %in% up.go.top200.GOID) #18 (top100), #41 (any significant) #top50: 15
sum(unique(unlist(go_ids)) %in% mouse.gene.set$gs_exact_source) #90

###
library(viridis)
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

up.go.df <- do.call(rbind, up.go)
head(up.go.df)
up.go.df <- up.go.df[up.go.df$direction == "up",]
#pivot wider to cluster
up.go.df.wider <- data.frame(tidyr::pivot_wider(up.go.df[, c("GOID", "group", "p.adjust")], names_from = group, 
                                                values_from = p.adjust, values_fn = mean))
rownames(up.go.df.wider) <- up.go.df.wider[,1]
up.go.df.wider <- up.go.df.wider[,-1]
up.go.df.wider <- up.go.df.wider[rownames(up.go.df.wider) %in% up.go.top200.GOID,]
up.go.df.wider[is.na(up.go.df.wider)] <- 1
nrow(up.go.df.wider)
apply(up.go.df.wider, 2, function(x) sum(x<0.05 & !is.na(x)))


#
group.cluster <- hclust(as.dist(1-cor(up.go.df.wider)), method = "average")

GO.cluster <- hclust(dist((-log10(up.go.df.wider))), method = "ward.D2")

up.go.df <- up.go.df[up.go.df$GOID %in% up.go.top200.GOID,]
y.order <- data.frame(GOID= GO.cluster$labels[GO.cluster$order],
                      y = 1:length(GO.cluster$labels[GO.cluster$order]),
                      GOID2 = up.go.df$ID[match(GO.cluster$labels[GO.cluster$order], up.go.df$GOID)])
head(y.order)
up.go.df$y <- y.order$y[match(up.go.df$GOID, y.order$GOID)]

unique(up.go.df$group)

up.go.df$p.adjust2 <- ifelse(up.go.df$p.adjust < 1e-30 & !is.na(up.go.df$p.adjust), 1e-30, up.go.df$p.adjust)

p <- ggplot(up.go.df, aes(x = group, y = y, fill = -log10(p.adjust2))) +
  geom_tile(width = 1) +
  scale_fill_gradientn(colours = custom_magma, na.value = "white", limits = c(-log10(5e-2), 30))+
  #facet_wrap(~direction, ncol = 1, strip.position = "right", scales = "free") +
  #scale_y_discrete(limits = GO.cluster$labels[GO.cluster$order]) +
  scale_y_continuous(breaks = y.order$y, labels = y.order$GOID,
                     sec.axis = sec_axis(~., breaks = y.order$y, labels = y.order$GOID2),
                     expand = expansion(mult = c(0, 0.01)))+
  scale_x_discrete(limits = c("11H_clonotype1vs2", "11L_clonotype1vs2","6H_clonotype1vs2", "6L_clonotype1vs2"))+ 
  #guide = guide_axis(n.dodge = 2)) + #limits = group.cluster$labels[group.cluster$order], 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 3), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  labs(x = "", y = "",fill = "-log10(p.adjust)")

#p

pdf(paste0("06.sameStimulation_clono1_2DEG_p",p_cutoff, "_logFC",logFC_cutoff,"_GOtop50.enrichment.pdf"), width = 8, height = 8, useDingbats = F)
print(p)
dev.off()

up.go.df.wider$yorder <- y.order$y[match(rownames(up.go.df.wider), y.order$GOID)]
write.csv(up.go.df.wider, file = paste0("06.sameStimulation_clono1_2DEG_p", p_cutoff, "_logFC", logFC_cutoff, "_GOtop50.enrichment.csv"), quote = F, row.names = T)


}  

"GO:0070371" %in% up.go.df$GOID
"GO:0042110" %in% up.go.df$GOID


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

mean.experssion <- data.frame(AverageExpression(scRNA[p.features,])$RNA)
mean.experssion$gene <- rownames(mean.experssion)
head(mean.experssion)

mean.experssion <- reshape2::melt(mean.experssion)
head(mean.experssion)


feature.ratio <- lapply( unique(scRNA@meta.data$orig.ident) , function(cluster){
  message(cluster)
  subset_cells <- scRNA[,scRNA$orig.ident == cluster]
  r <- rowMeans(subset_cells[['RNA']]$counts[p.features,] >0)
  print(r)
  r
})
feature.ratio <- data.frame(do.call(cbind, feature.ratio))
colnames(feature.ratio) <- unique(scRNA@meta.data$orig.ident)
head(feature.ratio)
feature.ratio$gene <- rownames(feature.ratio)

feature.ratio <- reshape2::melt(feature.ratio)
head(feature.ratio)
feature.ratio$x <- paste0(feature.ratio$gene, "-", feature.ratio$variable)

head(mean.experssion)
mean.experssion$variable <- gsub("^g", "", mean.experssion$variable)
mean.experssion$x <- paste0(mean.experssion$gene, "-", mean.experssion$variable)


identical(mean.experssion$x, feature.ratio$x)

nrow(mean.experssion)
nrow(feature.ratio)

mean.experssion$ratio <- feature.ratio$value[match(mean.experssion$x, feature.ratio$x)]
head(mean.experssion)

p <- ggplot(mean.experssion, aes(x = variable, y = gene, size = ratio, fill = value))+
  geom_point(shape = 21)+
  theme_bw()+
  scale_size(range = c(5,10))+
  scale_fill_gradient2(low = pair.col[2], mid = "white", high = pair.col[6], midpoint = 1)+
  scale_y_discrete(limits = p.features)+
  scale_x_discrete(limits = c("11H", "6H", "11L", "6L"))+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth =0.5))+
  labs(x ="", y = "")
p
pdf("06.inhibitorMarkers.dotPlot.pdf", width = 3, height = 4, useDingbats = F)
p
dev.off()



end here 2025-02-19

###
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


###DEG between clonotype1 - 2



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
