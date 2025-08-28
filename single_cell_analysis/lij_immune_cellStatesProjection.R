#remotes::install_github("carmonalab/STACAS")
#remotes::install_github("carmonalab/ProjecTILs")
setwd("/Users/wuz6/Documents/Project/lij36/2025/")

library(ProjecTILs)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(plotly)
library(scRepertoire)
library(paletteer)

rm(list)
ref <- load.reference.map()
cell.states.col = c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000", "#87f6a5", "#e812dd")
names(cell.states.col) = levels(ref$functional.cluster)



if(0){
#ref <- readRDS("ref_TILAtlas_mouse_v1.rds")
data(query_example_seurat)

query.projected <- Run.ProjecTILs(query_example_seurat, ref=ref)

query.projected

head(query.projected@meta.data)
rm(list=ls())

####https://carmonalab.github.io/ProjecTILs_CaseStudies/Xiong19_TCR.html
library(gridExtra)
library(ggplot2)
library(plotly)
library(ProjecTILs)
library(scRepertoire)

projectID <- "Xiong_TIL"
libIDtoSampleID <- c("Mouse 1", "Mouse 2", "Mouse 3", "Mouse 4")
names(libIDtoSampleID) <- 4:7

options(Seurat.object.assay.version = "v3") #
exp_mat <- Read10X("./E-MTAB-7919/")
querydata <- CreateSeuratObject(counts = exp_mat, project = projectID, min.cells = 3,
                                min.features = 50)
#querydata[["RNA"]] <- as(object = querydata[["RNA"]], Class = "Assay")

querydata$Sample <- substring(colnames(querydata), 18)
table(querydata$Sample)
querydata$SampleLabel <- factor(querydata$Sample, levels = c(4:7), labels = libIDtoSampleID)
table(querydata$SampleLabel)

###VDJ
libIDtoSampleID_VDJ <- 4:7
names(libIDtoSampleID_VDJ) <- 35:38

tcr_dir = "E-MTAB-7918/"
vdj.list <- list()
for (i in 1:length(libIDtoSampleID_VDJ)) {
  s <- names(libIDtoSampleID_VDJ)[i]
  vdj.list[[i]] <- read.csv(sprintf("%s/filtered_contig_annotations_%s.csv", tcr_dir,
                                    s), as.is = T)
  
  # Rename barcodes to match scRNA-seq suffixes
  vdj.list[[i]]$barcode <- sub("\\d$", "", vdj.list[[i]]$barcode)
  vdj.list[[i]]$barcode <- paste0(vdj.list[[i]]$barcode, libIDtoSampleID_VDJ[i])
  vdj.list[[i]]$raw_clonotype_id <- paste0(vdj.list[[i]]$raw_clonotype_id, "-",
                                           libIDtoSampleID_VDJ[i])
  vdj.list[[i]]$SampleLabel <- libIDtoSampleID_VDJ[i]
  
}

combined <- combineTCR(vdj.list, samples = libIDtoSampleID_VDJ, ID = names(libIDtoSampleID_VDJ),
                       cells = "T-AB", removeNA = T, removeMulti = T)

for (i in seq_along(combined)) {
  combined[[i]] <- stripBarcode(combined[[i]], column = 1, connector = "_", num_connects = 3)
}

querydata <- combineExpression(combined, querydata, cloneCall = "gene", group.by = "none")

ref <- load.reference.map()

query.projected <- make.projection(querydata, ref = ref, ncores = 2)

p1 <- plot.projection(ref)
p2 <- plot.projection(ref, query.projected, linesize = 0.5, pointsize = 0.5)
grid.arrange(p1, p2, ncol = 2)

head(query.projected@meta.data)

query.projected <- cellstate.predict(ref = ref, query = query.projected)

query.projected@meta.data$functional.cluster[1:4]
plot.statepred.composition(ref, query.projected, metric = "Percent")


}


####Our data
sample.names <- list.files("../data/", pattern = "[6|11]")
projectID <- "avidity"
options(Seurat.object.assay.version = "v3") #
load("sc.rna.filterDoublet.rda")
scRNA = merge(sc.data.list[[1]], y = sc.data.list[-1]) #,add.cell.ids = all.samples)
head(scRNA@meta.data)
scRNA <- subset(scRNA, subset = scDblFinder.class == "singlet")
rm(sc.data.list)

colnames(scRNA) <- sapply(colnames(scRNA), function(x){ paste0(strsplit(x, "_")[[1]][1:2], collapse = "_") })
colnames(scRNA)[1:4]

####
vdj = read.csv("../data/vdj_t/filtered_contig_annotations.csv", as.is = T)
#vdj.list[[i]]$barcode <- sub("\\d$", "", vdj.list[[i]]$barcode)
#vdj.list[[i]]$barcode <- paste0(vdj.list[[i]]$barcode, libIDtoSampleID_VDJ[i])
vdj$raw_clonotype_id <- paste0(vdj$origin, "-", vdj$raw_clonotype_id)
vdj$SampleLabel <- vdj$origin
vdj.list = split(vdj, vdj$SampleLabel)

combined <- combineTCR(vdj.list, samples = sample.names, cells = "T-AB", removeNA = F, removeMulti = T)
combined[[1]]$barcode[1:4]

#for (i in seq_along(combined)) {
#  combined[[i]] <- stripBarcode(combined[[i]], column = 1, connector = "_", num_connects = 3)
#}
#######
querydata <- combineExpression(combined, scRNA, cloneCall = "gene", group.by = "none")


p0 <- plot.projection(ref) & theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 1)

keep.ref.celltype <- grep("CD8", unique(ref[[]]$functional.cluster), value = T)

ref1 <- subset(ref, subset = functional.cluster %in% keep.ref.celltype)
VariableFeatures(ref1)
###use preselected featires
ref1 <- RunPCA(ref1, assay = "integrated", reduction.name = "integrated.pca")
ElbowPlot(ref1, reduction = "integrated.pca", ndims = 50)
ref1 <- FindNeighbors(ref1, reduction = "integrated.pca")
ref1@reductions$umap.old<- ref1@reductions$umap

###
ref1 <- RunUMAP(ref1, assay = "integrated", #reduction.name = "ref1.umap", 
                dims = 1:15, nn.name = NULL, graph = NULL, features = NULL, return.model =T)
DimPlot(ref1, reduction = "umap")

#ref1@reductions$umap <- ref1@reductions$ref1.umap

p1 <- plot.projection(ref1) & labs(title = "Reference CD8") & theme_bw()+ theme(panel.grid = element_blank(), aspect.ratio = 1) 
p1


###
scRNA.v3 <- as(object = scRNA[["RNA"]], Class = "Assay")
scRNA.v3 <- CreateSeuratObject(scRNA.v3)
scRNA.v3 <- NormalizeData(scRNA.v3)
scRNA.v3 <- ScaleData(scRNA.v3)

query.projected <- make.projection(scRNA.v3, ref = ref1, ncores = 2)
DimPlot(query.projected)

p2 <- plot.projection(ref1, query.projected, linesize = 0.2, pointsize = 0.5) & theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 1)
p2


query.projected <- cellstate.predict(ref = ref1, query = query.projected)
query.projected@meta.data$functional.cluster = factor(query.projected@meta.data$functional.cluster, 
                                                      levels = levels(ref@meta.data$functional.cluster)[levels(ref@meta.data$functional.cluster) %in% unique(query.projected@meta.data$functional.cluster)])
quantile(query.projected@meta.data$functional.cluster.conf)
plot(density(query.projected@meta.data$functional.cluster.conf))
abline(v = 0.53)
query.projected@meta.data$functional.cluster.highconf = query.projected@meta.data$functional.cluster
query.projected@meta.data$functional.cluster.highconf[query.projected@meta.data$functional.cluster.conf < 0.52] = NA

head(query.projected[[]])

x1 <- cbind(ref1@reductions$umap@cell.embeddings, ref1@meta.data)
x2 <- cbind(query.projected@reductions$umap@cell.embeddings, query.projected@meta.data)

p3 <- ggplot(x1, aes(x = umap_1, y = umap_2))+
  geom_point(color = "gray", size = 1)+
  geom_point(inherit.aes = F, data = x2, aes(x = UMAP_1, y = UMAP_2, color = functional.cluster.highconf), size =1)+
  scale_color_manual(values = cell.states.col)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank())
p3


pdf("03.function.cells.paperMethod.pdf", width = 10, height = 10)
cowplot::plot_grid(p0, p1, p2,p3,  nrow = 2)
dev.off()
##Something is wrong with the plot, seems the query umap did not update from ref1 new umap.


###transfer using seurat
anchors <- FindTransferAnchors(reference = ref1, query = scRNA.v3, dims = 1:15, reduction = "rpca", reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = ref1$functional.cluster, dims = 1:15)

scRNA.v3 <- AddMetaData(scRNA.v3, metadata = predictions)
scRNA.v3 <- MapQuery(anchorset = anchors, reference = ref1, query = scRNA.v3,
                           refdata = list(celltype = "functional.cluster"), reference.reduction = "pca", reduction.model = "umap")


plot(density(scRNA.v3@meta.data$prediction.score.max))
abline(v = 0.5)
quantile(scRNA.v3@meta.data$prediction.score.max, 0.1)
FeaturePlot(scRNA.v3, features = "prediction.score.max")
scRNA.v3@meta.data$predicted.celltype[scRNA.v3@meta.data$prediction.score.max< 0.5] <- 'NA'

table(scRNA.v3@meta.data$predicted.celltype)

p1 <- DimPlot(ref1, reduction = "umap", group.by = "functional.cluster", label = TRUE, label.size = 3,
              repel = TRUE) &  ggtitle("Reference CD8") & 
  scale_color_manual(values = cell.states.col) &
  theme_bw() &
  theme(aspect.ratio = 1, panel.grid = element_blank()) &
  guides(colour = guide_legend(override.aes = list(size=3))) 
p1

p2 <- DimPlot(scRNA.v3, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) &
  labs(title = "Query transferred labels", x = "ref umap1", y = "ref umap2") &
  scale_color_manual(values = cell.states.col)&
  theme_bw() &
  theme(aspect.ratio = 1, panel.grid = element_blank()) &
  guides(colour = guide_legend(override.aes = list(size=3))) 
p2

p1 + p2



###
###
x1 = data.frame(ref1@reductions$umap@cell.embeddings)
x1$state = ref1@meta.data$functional.cluster

x2 = data.frame(scRNA.v3@reductions$ref.umap@cell.embeddings)
x2$state = scRNA.v3@meta.data$predicted.celltype

p3 = ggplot(x1, aes(x = umap_1, y = umap_2, color = state)) + 
  geom_point(inherit.aes = F,data = x2, aes(x = refUMAP_1, y = refUMAP_2), color = "gray",size = .5, shape = 2)+
  geom_point(size = .1)+
  theme_bw()+
  scale_color_manual(values = cell.states.col)+
  theme(aspect.ratio = 1, panel.grid = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=3))) 

p4 = ggplot(x2, aes(x = refUMAP_1, y = refUMAP_2, color = state)) + 
  geom_point(inherit.aes = F,data = x1, aes(x = umap_1, y = umap_2), color = "gray",size = .1)+
  geom_point(size = .5, shape = 2)+
  theme_bw()+
  scale_color_manual(values = cell.states.col)+
  theme(aspect.ratio = 1, panel.grid = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  labs(x = "umap1", y = "umap2")

cowplot::plot_grid(p3, p4)

pdf("03.function.cell.seuratTransfer.pdf", width = 10, height = 10, useDingbats = F)
grid.arrange(p1, p2,p3,p4, nrow = 2)
dev.off()


markers = c("CD4", "TCF7", "CCR7", "PDCD1","HAVCR2", "GZMB", "GZMK", "CX3CR1", "CXCR5", "FOXP3")
grep(paste0(markers, collapse = "|"), rownames(query.projected), ignore.case = T, value = T)

use.markers = c("Tcf7","Ccr7","Pdcd1", "Havcr2", "Gzmb", "Gzmk", "Cx3cr1", "Foxp3")

pdf("03.umap.cellStates.projectionMarkers.pdf", width = 10, height = 10, useDingbats = F)
FeaturePlot(scRNA.v3, features = use.markers)
dev.off()

mm = scRNA.v3@meta.data
head(mm)

write.csv(mm, file = "03.cellState.prediction.csv")

table(mm$predicted.celltype)


mm1 <- mm[mm$predicted.celltype != "NA", ]
table(mm1$predicted.celltype)

x = data.frame(table(mm1$orig.ident, as.character(mm1$predicted.celltype)))

x$Var1 = factor(x$Var1, levels = c("6L", "6H", "11L", "11H"))

p = ggplot(x, aes(x = Var1, y = Freq, fill = Var2))+
  geom_bar(stat="identity", width = 1)+
  theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank(), axis.ticks = element_blank())+
  #coord_polar("y", start = 0)+
  scale_fill_manual(values = cell.states.col[names(cell.states.col) %in% unique(as.character(x$Var2))])+
  labs(fill = "predicted cell states", x = "", y = "")+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))+
  scale_x_discrete(expand = c(0,0))
p

pdf("03.cellStates.projection.ratio.pdf", width = 5, height = 5, useDingbats = F)
p
dev.off()

write.csv(x, file = "03.cellStates.projection.ratio.csv", row.names = F)



identical(colnames(scRNA), colnames(scRNA.v3))

scRNA@meta.data$func <- scRNA.v3@meta.data$predicted.celltype[match(rownames(scRNA@meta.data), rownames(scRNA.v3@meta.data))]

p1 <- DimPlot(scRNA, reduction = "umap_pc15",  group.by = "func", pt.size = .05)+
  scale_color_manual(values = cell.states.col)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "umap1", y = "umap2", title = "")
p1

p2 <- DimPlot(scRNA, reduction = "tsne_pc15",  group.by = "func", pt.size = .05)+
  scale_color_manual(values = cell.states.col)+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "tsne1", y = "tsne2", title = "")
p2

pdf("03.function.cells.pdf", width = 10, height = 5, useDingbats = F)
cowplot::plot_grid(p1, p2)
dev.off()

VariableFeatures(ref)



###clonotype based.....
load("../01.seurat.obj.rda")
colnames(scRNA) = paste0(sapply(rownames(scRNA@meta.data), function(x) strsplit(x, "\\-1")[[1]][1]), "-1")
keep.clono <- names(which(table(scRNA@meta.data$T_clonotype_id) > 10))
keep.clono
clono.col = as.character(paletteer_d("ggthemes::Tableau_20"))[1:length(keep.clono)]
names(clono.col) = keep.clono
clono.col
scRNA <- subset(scRNA, subset = T_clonotype_id %in% keep.clono)

sum(!rownames(scRNA@meta.data) %in% rownames(mm))
rownames(scRNA@meta.data)[!rownames(scRNA@meta.data) %in% rownames(mm)]

scRNA = AddMetaData(scRNA, mm)

xx = cbind(scRNA@reductions$umap@cell.embeddings, scRNA@reductions$integrationUMAP@cell.embeddings, scRNA@meta.data)
xx$functional.cluster2 = xx$functional.cluster.highconf
#xx$functional.cluster2[xx$functional.cluster.conf <0.52] = NA
xx$sampleID <- factor(xx$sampleID, levels = c("6L", "6H", "11L", "11H"))
xx$T_clonotype_id <- factor(xx$T_clonotype_id, levels = paste0("clonotype", 1:length(keep.clono)))


sample.color <- RColorBrewer::brewer.pal(12, "Paired")[c(2,4,8,10)]
names(sample.color) <- c("11H", "6H", "11L", "6L")

p.sample = DimPlot(scRNA, group.by = "orig.ident") + 
  scale_color_manual(values = sample.color) + 
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 1)+
  labs(title = "")+
  guides(colour = guide_legend(override.aes = list(size=5))) 
p.sample

p.clonotype = DimPlot(scRNA, group.by = "T_clonotype_id") + 
  scale_color_paletteer_d("ggthemes::Tableau_20")+
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 1)+
  labs(title = "")+
  guides(colour = guide_legend(override.aes = list(size=5))) 
p.clonotype


p1 = ggplot(xx, aes(x = umap_1, y = umap_2, col = functional.cluster.highconf))+
  geom_point(size = .1)+
  scale_color_manual(values = cell.states.col) + 
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 1)+
  facet_wrap(~sampleID, nrow = 1)+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  labs(x = "umap1", y = "umap2", title = "Merge")
p1

p2 = ggplot(xx, aes(x = cca_1, y = cca_2, col = functional.cluster.highconf))+
  geom_point(size = .1)+
  scale_color_manual(values = cell.states.col) + 
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 1)+
  facet_wrap(~sampleID, nrow = 1)+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  labs(x = "umap1", y = "umap2", title = "Seurat integration")
p2



p3 = ggplot(xx[xx$T_clonotype_id %in% c("clonotype1", "clonotype2", "clonotype3", "clonotype4"),], 
            aes(x = umap_1, y = umap_2, col = functional.cluster.highconf))+
  geom_point(size = .1)+
  scale_color_manual(values = cell.states.col) + 
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 1)+
  facet_grid(T_clonotype_id ~ sampleID, switch = "y")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  labs(x = "umap1", y = "umap2", title = "Merge")
p3

p4 = ggplot(xx[xx$T_clonotype_id %in% c("clonotype1", "clonotype2", "clonotype3", "clonotype4"),], 
            aes(x = cca_1, y = cca_2, col = functional.cluster.highconf))+
  geom_point(size = .1)+
  scale_color_manual(values = cell.states.col) + 
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 1)+
  facet_grid(T_clonotype_id ~ sampleID, switch = "y")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  labs(x = "umap1", y = "umap2", title = "Seurat integration")
p4

pdf("03.cellState.UMAP.pdf", width = 16, height = 8, useDingbats = F)
p1
p2
p3
p4
dev.off()

FeaturePlot()

xx$orig.ident <- factor(xx$orig.ident, levels = c("6L", "6H", "11L", "11H"))

xx2 = data.frame(table(xx$orig.ident, xx$functional.cluster.highconf))
#View(table(xx$orig.ident, xx$functional.cluster)/rowSums(table(xx$orig.ident, xx$functional.cluster)))
xx2 = data.frame(table(xx$orig.ident, xx$functional.cluster.highconf)/rowSums(table(xx$orig.ident, xx$functional.cluster.highconf)))
head(xx2)
write.csv(xx2, file = "03.cellStates.freq.csv")

xx3 = xx[xx$T_clonotype_id == "clonotype1",]
xx3 = data.frame(table(xx3$orig.ident, xx3$functional.cluster)/rowSums(table(xx3$orig.ident, xx3$functional.cluster)))
xx4 = xx[xx$T_clonotype_id == "clonotype2",]
xx4 = data.frame(table(xx4$orig.ident, xx4$functional.cluster)/rowSums(table(xx4$orig.ident, xx4$functional.cluster)))

ggplot(xx2, aes(x= Var1, y = Freq, fill = Var2))+
  geom_histogram(stat="identity")+
  scale_color_manual(values = cell.states.col) + 
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 1)

ggplot(xx2, aes(x= Var1, y = Freq, fill = Var2))+
  geom_bar(stat="identity", width = 0.8)+
  scale_color_manual(values = cell.states.col) + 
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 1)+
  coord_polar("y", start=0)


xx2$Var1 <- factor(xx2$Var1, levels = c("6L", "6H", "11L", "11H"))
xx3$Var1 <- factor(xx3$Var1, levels = c("6L", "6H", "11L", "11H"))

p1 = ggplot(xx2, aes(x= Var1, y = Freq, fill = Var2))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = cell.states.col) + 
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 1)+
  facet_wrap(~Var2, scales = "free")+
  scale_x_discrete(limits = c("6L", "6H", "11L", "11H"))+
  labs(x="", fill = "")
p1

p2 = ggplot(xx3, aes(x= Var1, y = Freq, fill = Var2))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = cell.states.col) + 
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 1)+
  facet_wrap(~Var2, scales = "free")+
  scale_x_discrete(limits = c("6L", "6H", "11L", "11H"))+
  labs(x="", fill = "", title = "clonotype1")

p3 = ggplot(xx4, aes(x= Var1, y = Freq, fill = Var2))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = cell.states.col) + 
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 1)+
  facet_wrap(~Var2, scales = "free")+
  scale_x_discrete(limits = c("6L", "6H", "11L", "11H"))+
  labs(x="", fill = "", title = "clonotype2")

pdf("03.cellStates.freq.histogram.pdf", width = 20, height = 10, useDingbats = F)
cowplot::plot_grid(p1, p2, p3)
dev.off()


table(xx$orig.ident, xx$functional.cluster)

sample.size = 2000
use.states = unique(xx$functional.cluster)
use.states = use.states[!is.na(use.states)]
use.states
set.seed(123, kind = "L'Ecuyer-CMRG")

x.freq = lapply(1:100, function(k){
    ff = t(sapply(unique(xx$orig.ident), function(s.name){
        xx1 = xx[sample(rownames(xx)[xx$orig.ident ==s.name], sample.size),]
        y = table(xx1$functional.cluster)[use.states]
        y[is.na(y)] = 0
        names(y) = use.states
        y
    }))
    
    data.frame(ff/rowSums(ff))
})
x.freq[[1]]
rowSums(x.freq[[1]])

x.freq = lapply(x.freq, function(x){x = data.frame(x); x$sample = rownames(x);x})

x.freq = data.frame(do.call(rbind, x.freq))
head(x.freq)

x.freq.df = reshape2::melt(x.freq, id.vars = c("sample"))

head(x.freq.df)

x.freq.mean = aggregate(x.freq.df$value, list(x.freq.df$sample, x.freq.df$variable), mean)
x.freq.sd = aggregate(x.freq.df$value, list(x.freq.df$sample, x.freq.df$variable), sd)
x.freq.mean$sd = x.freq.sd$x

p = ggplot(x.freq.mean, aes(x = Group.1, y = x, fill = Group.2))+
  geom_histogram(stat = "identity")+
  geom_errorbar( aes(x=Group.1, ymin=x, ymax=x+sd), width=0.3, colour="black", alpha=0.9, size=1.3)+
  scale_fill_manual(values = cell.states.col) + 
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 1)+
  facet_wrap(~Group.2, scales = "free")+
  scale_x_discrete(limits = c("11H", "11L", "6H", "6L"))+
  labs(x = "", y = "", fill = "")
p  

#pdf("03.cellStates.freq.histogram.pdf", width = 10, height = 10, useDingbats = F)
#p
#dev.off()

