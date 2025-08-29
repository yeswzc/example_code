setwd("/Users/wuz6/Documents/Project/08.spatialTranscriptome/01.RNA/01.QCfirst/")
#R/4.1.2
rm(list=ls())
library(ggplot2)
library(ggrepel)
library(ggridges)
#library(dendextend)
library(Seurat)
library(harmony)
library(pheatmap)
library(ComplexHeatmap)
#library(umap)
#library(sva)
library(clusterProfiler);
library(dplyr)
library(paletteer)
library(edgeR)

#
pair.col = RColorBrewer::brewer.pal(12, "Paired")
red_blue_30 <- rev(as.character(paletteer::paletteer_c("ggthemes::Classic Red-Blue", 30)))

my.thm <- theme_bw()+ theme(aspect.ratio = 1, panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

loc.shape <- c("TC" =15, "IE" = 2, "GM" = 1, "MVP" = 11, "PN" = 4, "BV" = 3)

load("CTA.allGlioma.seuratObj.rda")
load("CTA.colors.shape.rda")

#########################################################################################################
#########################################################################################################
#########################################################################################################
length(VariableFeatures(seurat))
sub.seurat <- subset(seurat, subset = tumor == "IDH_A")
sample.A.colors = RColorBrewer::brewer.pal(12, "Paired")[c(1:4,12,5:10)]
names(sample.A.colors) <- unique(sub.seurat@meta.data$sample)
sub.seurat <- subset(seurat, subset = tumor == "IDH_O")
sample.O.colors = as.character(paletteer::paletteer_d("ggsci::default_nejm"))
names(sample.O.colors) <- unique(sub.seurat@meta.data$sample)
sub.seurat <- subset(seurat, tumor == "GBM" | tumor == "GBMped")
sample.G.colors = c(as.character(paletteer::paletteer_d("ggsci::category10_d3")), cancer.col['GBMped'])
names(sample.G.colors) <- unique(sub.seurat@meta.data$sample)

sample.colors <- c(sample.A.colors, sample.O.colors, sample.G.colors)

#########################################################################################################
tmp <- seurat[-which(rownames(seurat) == "Negative Probe"),]
tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 300)
x <- HVFInfo(tmp)
head(x)
plot(density(x$variance.standardized))
q <- quantile(x$variance.standardized, 0.75)
q
sum(x$variance.standardized >= q)
VariableFeatures(seurat) <- rownames(x)[x$variance.standardized >= q];rm(tmp)
#VariableFeatures(seurat) <- VariableFeatures(tmp) ;rm(tmp)
length(VariableFeatures(seurat))

#seurat <- SCTransform(seurat)
seurat <- RunPCA(seurat, verbose = F)
x <- seurat@reductions$pca@stdev^2/sum(seurat@reductions$pca@stdev^2)
x
sum(x[1:8])
sum(!unique(seurat$sample) %in% names(sample.colors))

p1 <- DimPlot(seurat, reduction = "pca", group.by = "sample", cols = sample.colors, pt.size = 1,
              shape.by = "location")+ my.thm + labs(title = "")#+ loc.shape
p1
p2 <- DimPlot(seurat, reduction = "pca", group.by = "location", cols = location.col, pt.size = 1)+ my.thm + labs(title = "")
p2
p3 <- DimPlot(seurat, reduction = "pca", group.by = "batch",pt.size = 1)+ my.thm + labs(title = "")
p4 <- DimPlot(seurat, reduction = "pca", group.by = "grade",pt.size = 1)+ my.thm + labs(title = "")
p5 <- DimPlot(seurat, reduction = "pca", group.by = "sex", cols = pair.col[c(5,2)], pt.size = 1)+ my.thm + labs(title = "")
dev.off()

pdf("all.glioma/01.PCA.pdf", width = 20, height = 14, useDingbats = F)
#pdf("AIDH/01.AIDH.PCA.pdf", width = 15, height = 10, useDingbats = F)#
#pdf("GBM/01.GBM.PCA.pdf", width = 15, height = 10, useDingbats = F)
cowplot::plot_grid(p1,p2, p3, p4, p5, nrow = 2)
dev.off()



table(seurat@meta.data$location)
pdf("all.glioma/01.pca.elbow.pdf", width = 5, height = 4, useDingbats = F)
ElbowPlot(seurat) + theme(aspect.ratio = 1)
dev.off()

loc.shape <- c("TC" =15, "IE" = 2, "GM" = 1, "MVP" = 11, "PN" = 4, "BV" = 3)
loc.shape <- loc.shape[names(loc.shape) %in% seurat@meta.data$location]
loc.shape
loc.shape <- scale_shape_manual(values = loc.shape)

seurat <- RunUMAP(seurat, dims = 1:10, verbose = F, n.neighbors = 15L)
###Plot UMAP
p1 <- DimPlot(seurat, reduction = "umap", group.by = "sample", cols = sample.colors, 
              pt.size = 1,shape.by = "location")+ my.thm + labs(x="", y = "", title = "") + loc.shape
p1
p2 <- DimPlot(seurat, reduction = "umap", group.by = "location", cols = location.col, pt.size = 1)+ my.thm + labs(x="", y = "", title = "")
p2
#p3 <- DimPlot(seurat, reduction = "umap", group.by = "Tumor", cols = cancer.col, pt.size = .5)+ my.thm + labs(title = "")
p3 <- DimPlot(seurat, reduction = "umap", group.by = "tumor", cols = cancer.col, pt.size = 1, shape.by = "location")+ my.thm + labs(title = "")+loc.shape
p4 <- DimPlot(seurat, reduction = "umap", group.by = "grade", cols = grade.col, pt.size = 1, shape.by = "location")+ my.thm + labs(x="", y = "", title = "")+loc.shape
p5 <- DimPlot(seurat, reduction = "umap", group.by = "sex", cols = pair.col[c(5,2)], pt.size = 1, shape.by = "location")+ my.thm + labs(x="", y = "", title = "")+loc.shape

pdf("all.glioma/01.umap.pdf", width = 20, height = 15, useDingbats = F)
#pdf("AIDH/01.AIDH.umap.pdf", width = 15, height = 10, useDingbats = F)
#pdf("GBM/01.GBM.umap.pdf", width = 15, height = 10, useDingbats = F)
cowplot::plot_grid(p1,p2,p3, p4, p5, nrow = 2)
dev.off()

if(0){
  res <- cbind(seurat@meta.data, seurat@reductions$umap@cell.embeddings)
  
  library(plotly)
  
  p1 <- plot_ly(x = res$UMAP_1, y = res$UMAP_2, text = paste0(rownames(res), ":", res$location),
                color = res$sample, colors = sample.colors, symbol = res$location, 
                symbols = c("TC" =15, "IE" = 2, "GM" = 1, "MVP" = 11, "PN" = 4, "BV" = 3),
                type = "scatter", mode = "markers",
                marker = list(size = 4), legendgroup= "sample")
  
  
  p1
  
  p1 <- p1 %>%  
    #  add_annotations( text="location:", xref="paper", yref="paper",
    #                   x=1.02, xanchor="left",
    #                   y=0.9, yanchor="bottom",   
    #                   legendtitle=TRUE, showarrow=FALSE ) %>%
    #add_annotations( text="sample:", xref="paper", yref="paper",
    #                 x=1.02, xanchor="left",
    #                 y=0.7, yanchor="bottom",   
    #                 legendtitle=TRUE, showarrow=FALSE ) %>%
    layout(showlegend = TRUE, 
           xaxis = list(title = "umap 1", zeroline = FALSE),
           yaxis = list(title = "umap 2", zeroline = FALSE))
  
  p1
  
  
}

seurat <- RunTSNE(seurat, dims = 1:10, verbose = F)
###Plot TSNE
p1 <- DimPlot(seurat, reduction = "tsne", group.by = "sample", cols = sample.colors, pt.size = 1,shape.by = "location")+ my.thm + labs(x="", y = "", title = "")+loc.shape
p1
p2 <- DimPlot(seurat, reduction = "tsne", group.by = "location", cols = location.col, pt.size = 1)+ my.thm + labs(x="", y = "", title = "")
p2
#p3 <- DimPlot(seurat, reduction = "tsne", group.by = "Tumor", cols = cancer.col, pt.size = .5)+ my.thm + labs(title = "")
p3 <- DimPlot(seurat, reduction = "tsne", group.by = "tumor", cols = cancer.col, pt.size = 1, shape.by = "location")+ my.thm + labs(title = "")+loc.shape
p4 <- DimPlot(seurat, reduction = "tsne", group.by = "grade", cols = grade.col, pt.size = 1, shape.by = "location")+ my.thm + labs(x="", y = "", title = "")+loc.shape
p5 <- DimPlot(seurat, reduction = "tsne", group.by = "sex", cols = pair.col[c(5,2)], pt.size = 1, shape.by = "location")+ my.thm + labs(x="", y = "", title = "")+loc.shape

pdf("all.glioma/01.tsne.pdf", width = 20, height = 15, useDingbats = F)
#pdf("AIDH/01.AIDH.tsne.pdf", width = 15, height = 10, useDingbats = F)
#pdf("GBM/01.GBM.tsne.pdf", width = 15, height = 10, useDingbats = F)
cowplot::plot_grid(p1,p2,p3, p4, p5, nrow = 2)
dev.off()


### Hierachical clustering
seurat@assays$RNA@data[1:4, 1:4]
identical(colnames(seurat@assays$RNA@scale.data), rownames(seurat@meta.data))
head(seurat@meta.data)
data.variable.ann <- data.frame(seurat@meta.data[,c(5:10)])
head(data.variable.ann)

set.seed(123)
#seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 800)
#uses.genes.idx <- sample(1:nrow(seurat), 500)



###################
m <- t(scale(t(data.matrix(seurat@assays$RNA@data[VariableFeatures(seurat),]))))
m[1:4,1:4]
dim(m)
#gene.cluster <- hclust(dist(seurat@assays$RNA@scale.data[VariableFeatures(seurat),], method = "euclidean"), method = "ward.D2")
#roi.cluster <- hclust(dist(t(as.matrix(seurat@assays$RNA@scale.data[VariableFeatures(seurat),])), method = "euclidean"), method = "ward.D2")
#gene.cluster <- hclust(as.dist(1 - cor(t(m), method = "pearson")), method = "average")
#roi.cluster <- hclust(as.dist(1 - cor(m, method = "pearson")), method = "average")

dev.off()
plot(density(m))

dim(m)

set.seed(123)
#sample.colors <- sample(as.character(paletteer::paletteer_d("ggsci::default_igv")), 30)
#names(sample.colors) <- sort(unique(seurat@meta.data$sample))

sample.ann <- HeatmapAnnotation(batch = seurat@meta.data$batch,
                                sex = seurat@meta.data$sex,
                                grade = seurat@meta.data$grade,
                                location = seurat@meta.data$location,
                                tumor = seurat@meta.data$tumor,
                                sample = seurat@meta.data$sample,
                                col = list(sample = sample.colors[names(sample.colors) %in% seurat@meta.data$sample],
                                           tumor = cancer.col,
                                           location = location.col[names(location.col) %in% seurat@meta.data$location],
                                           sex = sex.col[names(sex.col) %in% seurat@meta.data$sex],
                                           grade = grade.col[names(grade.col) %in% seurat@meta.data$grade],
                                           batch = batch.col[names(batch.col) %in% seurat@meta.data$batch]))

m <- t(scale(t(data.matrix(seurat@assays$RNA@data[VariableFeatures(seurat),]))))
m1 <- m; m1[m1 > 4] =4; m1[m1< -4] = -4

gene.km <- kmeans(m, centers = 10, nstart = 500)$cluster
gene.km[1:4]
write.csv(gene.km, file = "all.glioma/01.variable.Q75genes.cluster.csv", quote = F)
gene.cluster <- cluster_within_group(t(m), gene.km)

roi.km <- kmeans(t(m), centers = 10, nstart = 500)$cluster
roi.cluster <- cluster_within_group(m, roi.km)

gene.to.mark <- read.table("all.glioma/genes.to.mark.txt", head =F)[,1]
gene.to.mark
gene.to.mark = rowAnnotation(foo = anno_mark(at = which(rownames(m) %in% gene.to.mark), 
                                   labels = rownames(m)[which(rownames(m) %in% gene.to.mark)]))

dev.off()


identical(names(roi.km), rownames(seurat@meta.data))
Idents(seurat) <- roi.km

#mk <- FindAllMarkers(seurat, only.pos = T, logfc.threshold = 0.5)
#mk <- mk[rownames(mk) %in% rownames(m),]
#nrow(mk)
#table(mk$cluster)
#write.csv(mk, file = "all.glioma/01.Q75.ROI.kmeans.markers.csv", quote = F, row.names = F)

h <- Heatmap(m1,  name = "z-score",
             top_annotation = sample.ann,
             right_annotation = gene.to.mark,
             col =  red_blue_30,
             border = T,
             cluster_rows = gene.cluster, 
             cluster_columns = roi.cluster,
             row_split = 10, column_split = 10,
             show_row_dend = T, show_column_dend = T,
             show_row_names = F, show_column_names = F)

png("all.glioma/01.variable.Q75genes.clusterHeatmap.png", width = 15, height = 15, units = "in", res = 100)
draw(h)
dev.off()



h <- Heatmap(m1,  name = "z-score",
             top_annotation = sample.ann,
             right_annotation = gene.to.mark,
             col =  red_blue_30,
             border = T,
             cluster_rows = gene.cluster, 
             cluster_columns = roi.cluster,
             #row_split = 4, column_split = 2,
             show_row_dend = F, show_column_dend = F,
             show_row_names = F, show_column_names = F)

png("all.glioma/01.variable.Q75genes.clusterHeatmap.v2.png", width = 15, height = 15, units = "in", res = 100)
draw(h)
dev.off()


h <- Heatmap(m1,  name = "z-score",
             top_annotation = sample.ann,
             right_annotation = gene.to.mark,
             col =  red_blue_30,
             border = T,
             cluster_rows = gene.cluster, 
             cluster_columns = roi.cluster,
             row_split = 10, column_split = 10,
             show_row_dend = F, show_column_dend = F,
             show_row_names = F, show_column_names = F)

png("all.glioma/01.variable.Q75genes.clusterHeatmap.v3.png", width = 15, height = 15, units = "in", res = 100)
draw(h)
dev.off()

table(gene.km)

c5.db = msigdbr::msigdbr(species = "human", category = "C5") %>% dplyr::distinct(gene_symbol, gs_name) %>% as.data.frame()

e1 <- enricher(gene = names(gene.km[gene.km==1]), TERM2GENE = c5.db)@result #universe = names(gene.km), 
write.csv(e1, file = "all.glioma/02.gene.cluster.k1.go.csv", row.names = F)

sum(e1$p.adjust <0.05)
e1$ID[e1$p.adjust <0.05][1:10]
e1$ID <- gsub("\\_", " ", e1$ID); e1$ID <- gsub("^GO\\w\\w", " ", e1$ID); e1$ID <- gsub("^HP\\s", " ", e1$ID)

p1 <- ggplot(e1[order(e1$p.adjust, decreasing = F)[1:20],], aes(x = reorder(ID, -p.adjust), y = -log10(p.adjust))) +
  geom_histogram(stat = "identity")+
  scale_x_discrete(expand = expansion(mult = c(0,0)), 
                   labels = function(x) stringr::str_wrap(x, width = 40))+
  theme_bw()+
  theme(aspect.ratio = 1)+
  coord_flip()+
  labs(x = "", main = "k1")


e2 <- enricher(gene = names(gene.km[gene.km==7]), TERM2GENE = c5.db)@result
write.csv(e2, file = "all.glioma/02.gene.cluster.k7.go.csv", row.names = F)
sum(e2$p.adjust <0.05)
e2$ID[e2$p.adjust <0.05][1:4]
e2$ID <- gsub("\\_", " ", e2$ID); e2$ID <- gsub("^GO\\w\\w", " ", e2$ID); e2$ID <- gsub("^HP\\s", " ", e2$ID)
p2 <- ggplot(e2[order(e2$p.adjust, decreasing = F)[1:20],], aes(x = reorder(ID, -p.adjust), y = -log10(p.adjust))) +
  geom_histogram(stat = "identity")+
  scale_x_discrete(expand = expansion(mult = c(0,0)), 
                   labels = function(x) stringr::str_wrap(x, width = 40))+
  theme_bw()+
  theme(aspect.ratio = 1)+
  coord_flip()+
  labs(x = "", main = "k7")

##cell cycle
e3 <- enricher(gene = names(gene.km[!gene.km %in% c(1,7)]), TERM2GENE = c5.db)@result #universe = names(gene.km),
write.csv(e3, file = "all.glioma/02.gene.cluster.kNon17.go.csv", row.names = F)
sum(e3$p.adjust <0.05)
e3$ID[e4$p.adjust <0.05][1:10]
e3$ID <- gsub("\\_", " ", e3$ID); e3$ID <- gsub("^GO\\w\\w", " ", e3$ID); e3$ID <- gsub("^HP\\s", " ", e3$ID)
p3 <- ggplot(e3[order(e3$p.adjust, decreasing = F)[1:20],], aes(x = reorder(ID, -p.adjust), y = -log10(p.adjust))) +
  geom_histogram(stat = "identity")+
  scale_x_discrete(expand = expansion(mult = c(0,0)), 
                   labels = function(x) stringr::str_wrap(x, width = 40))+
  theme_bw()+
  theme(aspect.ratio = 1)+
  coord_flip()+
  labs(x = "", main = "k.non1.non7")


pdf("all.glioma/02.gene.cluster.GO.pdf", width = 20, height = 10, useDingbats = F)
cowplot::plot_grid(p2,p3)
dev.off()


###Check the overlap of the k4 genes and public signatures
x <- data.frame(readxl::read_xlsx("../../00.data/Ravi.Spatial.xlsx", skip = 1))
table(x$Regional.Programm)
gbm.spatial.programs <- list("Radial.Glia" = x$Gene[x$Regional.Programm == "Radial.Glia"],
                             "Reactive.Hypoxia" = x$Gene[x$Regional.Programm == "Reactive.Hypoxia"],
                             "Reactive.Immune" = x$Gene[x$Regional.Programm == "Reactive.Immune"],
                             "Spatial.NPC" = x$Gene[x$Regional.Programm == "Regional.NPC"],
                             "Spatial.OPC" = x$Gene[x$Regional.Programm == "Regional.OPC"])
x <- data.frame(readxl::read_xlsx("../../00.data/Richards.CureatedSignatureList.xlsx"))
colnames(x)
gbm.rna.subtype <- list("Verhaak.Proneural" = x$Verhaak_CancerCell_2010_Proneural[!is.na(x$Verhaak_CancerCell_2010_Proneural)],
                        "Verhaak.Neural" = x$Verhaak_CancerCell_2010_Neural[!is.na(x$Verhaak_CancerCell_2010_Neural)],
                        "Verhaak.Mesenchymal" = x$Verhaak_CancerCell_2010_Mesenchymal[!is.na(x$Verhaak_CancerCell_2010_Mesenchymal)],
                        "Verhaak.Classical" = x$Verhaak_CancerCell_2010_Classical[!is.na(x$Verhaak_CancerCell_2010_Classical)])

x <- data.frame(readxl::read_xlsx("../../00.data/GBM.cellStates.signature.xlsx", skip = 4))
head(x)
gbm.cell.states.signature <- list("MES" = c(x$MES1, x$MES2),
                                  "AC" = x$AC[!is.na(x$AC)],
                                  "OPC" = x$OPC[!is.na(x$OPC)],
                                  "NPC" = c(x$NPC1, x$NPC2),
                                  "G1S" = x$G1.S[!is.na(x$G1.S)],
                                  "G2M" = x$G2.M[!is.na(x$G2.M)])

table(gene.km)

res <- lapply(1:length(gbm.spatial.programs), function(k){
  x1 <- names(gbm.spatial.programs)[k]
  r <- lapply(unique(gene.km), function(j){
    G <- names(which(gene.km == j))
    #x2 <- length(intersect(G, gbm.spatial.programs[[k]]))/length(union(G, gbm.spatial.programs[[k]])) #Jac
    x2 <- mean(G %in% gbm.spatial.programs[[k]])
    c(x1, j, x2)
  })
  r <- do.call(rbind, r)
  r
})
res <- data.frame(do.call(rbind, res))
res
###

res <- lapply(1:length(gbm.rna.subtype), function(k){
  x1 <- names(gbm.rna.subtype)[k]
  r <- lapply(unique(gene.km), function(j){
    G <- names(which(gene.km == j))
    #x2 <- length(intersect(G, gbm.rna.subtype[[k]]))/length(union(G, gbm.rna.subtype[[k]])) #Jac
    x2 <- mean(G %in% gbm.rna.subtype[[k]])
    c(x1, j, x2)
  })
  r <- do.call(rbind, r)
  r
})
res <- data.frame(do.call(rbind, res))
res
###

res <- lapply(1:length(gbm.cell.states.signature), function(k){
  x1 <- names(gbm.cell.states.signature)[k]
  r <- lapply(unique(gene.km), function(j){
    G <- names(which(gene.km == j))
    #x2 <- length(intersect(G, gbm.cell.states.signature[[k]]))/length(union(G, gbm.cell.states.signature[[k]])) #jac
    x2 <- mean(G %in% gbm.cell.states.signature[[k]])
    c(x1, j, x2)
  })
  r <- do.call(rbind, r)
  r
})
res <- data.frame(do.call(rbind, res))
res
###

###It may not make sense to cluster all 3 tumors together because CTA cann't tell the 3 classes.

################shared variable genes across samples
if(1){
  
  xx <- lapply(unique(seurat@meta.data$sample), function(sp){
    x <- subset(seurat, sample == sp)
    x <- x[-which(rownames(x) == "Negative Probe"),]
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 500)
    #x <- HVFInfo(x)
    #q <- quantile(x$variance.standardized, 0.85)
    #rownames(x)[x$variance.standardized >= q]
    VariableFeatures(x)
    
  })
  xx.freq <- table(unlist(xx))
  #xx.freq
  plot(density(xx.freq)) 
  q <- quantile(xx.freq, 0.3)
  q
  sum(xx.freq > q) #85% quantile
  length(names(xx.freq[xx.freq > q]))
  VariableFeatures(seurat) <- names(xx.freq[xx.freq > q])
  length(VariableFeatures(seurat))
  
  #m[1:4,1:4]
  #plot(density(as.vector(m)))
  #gene.cluster <- hclust(as.dist(1 - cor(t(m), method = "pearson")), method = "average")
  #roi.cluster <- hclust(as.dist(1 - cor(m, method = "pearson")), method = "average")
  
  sample.ann <- HeatmapAnnotation(batch = seurat@meta.data$batch,
                                  sex = seurat@meta.data$sex,
                                  grade = seurat@meta.data$grade,
                                  location = seurat@meta.data$location,
                                  tumor = seurat@meta.data$tumor,
                                  sample = seurat@meta.data$sample,
                                  col = list(sample = sample.colors[names(sample.colors) %in% seurat@meta.data$sample],
                                             tumor = cancer.col,
                                             location = location.col[names(location.col) %in% seurat@meta.data$location],
                                             sex = sex.col[names(sex.col) %in% seurat@meta.data$sex],
                                             grade = grade.col[names(grade.col) %in% seurat@meta.data$grade],
                                             batch = batch.col[names(batch.col) %in% seurat@meta.data$batch]))
  
  m <- t(scale(t(data.matrix(seurat@assays$RNA@data[VariableFeatures(seurat),]))))
  gene.km <- kmeans(m, centers = 4, nstart = 500)$cluster
  gene.cluster <- cluster_within_group(t(m), gene.km)
  roi.km <- kmeans(t(m), centers = 2, nstart = 500)$cluster
  roi.cluster <- cluster_within_group(m, roi.km)
  
  m[m>4] = 4; m[m< -4] = -4
  h <- Heatmap(m,  name = "z-score",
               top_annotation = sample.ann,
               col =  red_blue_30,
               border = T,
               cluster_rows = gene.cluster, 
               cluster_columns = roi.cluster,
               row_split = 4, column_split = 2,
               #row_km = 4,
               #column_km = 2,
               show_row_dend = T, show_column_dend = T,
               show_row_names = F, show_column_names = F)
  
  png("all.glioma/01.shared.perSample.variable.genes.clusterHeatmap.png", width = 15, height = 15, units = "in", res = 100)
  draw(h)
  dev.off()
}





dev.off()
####Clustering for each tumor type
xx <- lapply(unique(seurat@meta.data$tumor), function(Tumor){
  message(Tumor)
  #Tumor <- "IDH_A"
  #dev.off()
  if(Tumor == "IDH_A"){
    sample.col <- sample.A.colors
  }else if(Tumor == "IDH_O"){
    sample.col <- sample.O.colors
  }else if(Tumor == "GBM"){
    sample.col <- sample.G.colors
  }else if(Tumor == "GBMped"){
    sample.col <- sample.G.colors
  }
  
  
  sub.seurat <- subset(seurat, tumor == Tumor)
  
  if(0 & length(unique(sub.seurat@meta.data$sample) > 1)){
    message("Feature selection by sample")
    out.prefix <- paste0("all.glioma/03.", Tumor, ".shared.variableGenes.")
    xx <- lapply(unique(sub.seurat$sample), function(sp){
      x <- subset(sub.seurat, sample == sp)
      x <- FindVariableFeatures(x, nfeatures = 500)
      VariableFeatures(x)
    })
    
    xx.freq <- table(unlist(xx))
    q <- quantile(xx.freq, 0.3)
    x <- names(xx.freq[xx.freq >= q])
    message("N variable genes = ", length(x))
    write.table(x, file = paste0(out.prefix, "csv"), quote = F, col.names = F, row.names = F)
    VariableFeatures(sub.seurat) <- x
  }else{
    message("Feature selection from all.")
    out.prefix <- paste0("all.glioma/03.", Tumor, ".09variableGenes.")
    tmp <- sub.seurat[-which(rownames(sub.seurat) == "Negative Probe"),]
    tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 200)
    xxx <- HVFInfo(tmp)
    #plot(density(xxx$variance.standardized))
    q <- quantile(xxx$variance.standardized, 0.9)
    #sum(xxx$variance.standardized > q)
    v.g <- rownames(xxx)[xxx$variance.standardized > q]
    #VariableFeatures(sub.seurat) <- VariableFeatures(tmp);rm(tmp)
    VariableFeatures(sub.seurat) <- v.g;rm(tmp, v.g)
    
  }
  loc.shape <- c("TC" =15, "IE" = 2, "GM" = 1, "MVP" = 11, "PN" = 4, "BV" = 3)
  loc.shape <- loc.shape[names(loc.shape) %in% sub.seurat@meta.data$location]
  loc.shape
  loc.shape <- scale_shape_manual(values = loc.shape)
  sub.seurat <- RunPCA(sub.seurat, verbose = F)
  sub.seurat <- RunUMAP(sub.seurat, dims = 1:5, verbose = F, n.neighbors = 15L)
  ###Plot UMAP
  p1 <- DimPlot(sub.seurat, reduction = "umap", group.by = "sample", cols = sample.col, 
                pt.size = 1,shape.by = "location")+ my.thm + labs(x="", y = "", title = "") + loc.shape
  #p1
  p2 <- DimPlot(sub.seurat, reduction = "umap", group.by = "location", 
                cols = location.col[names(location.col) %in% sub.seurat@meta.data$location], pt.size = 1)+ 
    my.thm + labs(x="", y = "", title = "")
  #p2
  p4 <- DimPlot(sub.seurat, reduction = "umap", group.by = "grade", cols = grade.col, pt.size = 1, shape.by = "location")+ my.thm + labs(x="", y = "", title = "")+loc.shape
  p5 <- DimPlot(sub.seurat, reduction = "umap", group.by = "sex", cols = pair.col[c(5,2)], pt.size = 1, shape.by = "location")+ my.thm + labs(x="", y = "", title = "")+loc.shape
  
  pdf(paste0(out.prefix, "umap.pdf"), width = 20, height = 15, useDingbats = F)
  print(cowplot::plot_grid(p1,p2,p4, p5, nrow = 2))
  dev.off()
  
  
  m <- t(scale(t(data.matrix(sub.seurat@assays$RNA@data[VariableFeatures(sub.seurat),]))))
  
  gene.cluster <- hclust(as.dist(1 - cor(t(m), method = "pearson")), method = "average")
  roi.cluster <- hclust(as.dist(1 - cor(m, method = "pearson")), method = "average")
  
  #gene.km <- kmeans(m, centers = 4, nstart = 500)$cluster
  #gene.cluster <- cluster_within_group(t(m), gene.km)
  #roi.km <- kmeans(t(m), centers = 6, nstart = 500)$cluster
  #roi.cluster <- cluster_within_group(m, roi.km)
  
  #data.variable.ann <- data.frame(batch = sub.seurat@meta.data$batch,
  #                                sex = sub.seurat@meta.data$sex,
  #                                grade = sub.seurat@meta.data$grade,
  #                                location = sub.seurat@meta.data$location,
  #                                sample = sub.seurat@meta.data$sample
  #)
  
  sample.ann <- HeatmapAnnotation(batch = sub.seurat@meta.data$batch,
                                  sex = sub.seurat@meta.data$sex,
                                  grade = sub.seurat@meta.data$grade,
                                  location = sub.seurat@meta.data$location,
                                  sample = sub.seurat@meta.data$sample,
                                  col = list(sample = sample.col[names(sample.col) %in% sub.seurat@meta.data$sample],
                                             location = location.col[names(location.col) %in% sub.seurat@meta.data$location],
                                             sex = sex.col[names(sex.col) %in% sub.seurat@meta.data$sex],
                                             grade = grade.col[names(grade.col) %in% sub.seurat@meta.data$grade],
                                             batch = batch.col[names(batch.col) %in% sub.seurat@meta.data$batch]))
  
  
  m[m>4] =4; m[m< -4] = -4;
  h <- Heatmap(m,  name = "z-score",
               top_annotation = sample.ann,
               col =  red_blue_30,
               border = T,
               cluster_rows = gene.cluster, cluster_columns = roi.cluster,
               #row_split = 4, column_split = 6,
               show_row_dend = T, show_column_dend = T,
               show_row_names = F, show_column_names = F)
  png( paste0(out.prefix,"cHeatmap.png"), width = 8, height = 8, units = "in", res = 100)
  draw(h)
  dev.off()
  VariableFeatures(sub.seurat) 
  
})



xx[[1]]
variable.gene.overlap <- lapply(1:length(xx), function(k){
  x <- ifelse(unique(unlist(xx)) %in% xx[[k]], 1, 0)
  x
})
variable.gene.overlap <- data.frame(do.call(cbind, variable.gene.overlap))
rownames(variable.gene.overlap) <- unique(unlist(xx))
colnames(variable.gene.overlap) <- unique(seurat@meta.data$tumor)
head(variable.gene.overlap)

library(UpSetR)
#print(head(genes.upset))
write.csv(variable.gene.overlap, file = "variableGenes/03.variableGnesUpset.csv")
pdf(file = "variableGenes/03.variableGnesUpset.pdf", width = 15, height = 8, useDingbats = F)
UpSetR::upset(variable.gene.overlap, sets = colnames(variable.gene.overlap), nintersects = NA)
dev.off()

###
sub.seurat <- subset(seurat, tumor == "GBM" | tumor == "GBMped")
unique(sub.seurat@meta.data$sample)
sample.G.colors = c(as.character(paletteer::paletteer_d("ggsci::category10_d3")), cancer.col['GBMped'])
names(sample.G.colors) <- unique(sub.seurat@meta.data$sample)
sample.G.colors

tmp <- sub.seurat[-which(rownames(sub.seurat) == "Negative Probe"),]
tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 800)
VariableFeatures(sub.seurat) <- VariableFeatures(tmp);rm(tmp)
sub.seurat <- RunPCA(sub.seurat, verbose = F)
#x <- sub.seurat@reductions$pca@stdev^2/sum(sub.seurat@reductions$pca@stdev^2)
#x
#sum(x[1:8])

table(sub.seurat@meta.data$location)
ElbowPlot(sub.seurat)

sub.seurat <- RunUMAP(sub.seurat, dims = 1:5, verbose = F, n.neighbors = 15L)
###Plot UMAP
p1 <- DimPlot(sub.seurat, reduction = "umap", group.by = "sample", cols = sample.G.colors, pt.size = 1,shape.by = "location")+ my.thm + labs(x="", y = "", title = "")+loc.shape
p1
p2 <- DimPlot(sub.seurat, reduction = "umap", group.by = "location", cols = location.col, pt.size = 1)+ my.thm + labs(x="", y = "", title = "")
p2
#p3 <- DimPlot(seurat, reduction = "umap", group.by = "Tumor", cols = cancer.col, pt.size = .5)+ my.thm + labs(title = "")
p3 <- DimPlot(sub.seurat, reduction = "umap", group.by = "batch", pt.size = 1, shape.by = "location")+ my.thm + labs(title = "")+loc.shape
p4 <- DimPlot(sub.seurat, reduction = "umap", group.by = "grade", cols = grade.col, pt.size = 1, shape.by = "location")+ my.thm + labs(x="", y = "", title = "")+loc.shape
p5 <- DimPlot(sub.seurat, reduction = "umap", group.by = "sex", cols = pair.col[c(5,2)], pt.size = 1, shape.by = "location")+ my.thm + labs(x="", y = "", title = "")+loc.shape

pdf("all.glioma/03.GBM.pGBM.800Genes.umap.pdf", width = 15, height = 10, useDingbats = F)
cowplot::plot_grid(p1,p2,p3, p4, p5, nrow = 2)
dev.off()







###end here
exit!
  
  
  
  gene.ann <- data.frame(gene = unique(c(gbm.over, gbm.under, idh.a.over, idh.a.under)),
                         class = "GBM")
gene.ann$class[gene.ann$gene %in% gbm.under] = "IDH"
gene.ann$class[gene.ann$gene %in% idh.a.over] = "AIDH"
gene.ann$class[gene.ann$gene %in% idh.a.under] = "OIDH"
rownames(gene.ann) <- gene.ann$gene
gene.ann <- gene.ann[rownames(gene.ann) %in% rownames(seurat),2, drop = F]
nrow(gene.ann)

pdf("all.glioma/01.tcgaDEG.clusterHeatmap.pdf", width = 10, height = 15, useDingbats = F)
pheatmap(seurat@assays$RNA@scale.data[rownames(gene.ann),], 
         #cluster_rows = gene.cluster, cluster_cols = roi.cluster,
         #scale = "row",
         annotation_col = data.variable.ann,
         annotation_row = gene.ann,
         annotation_colors = list(sample= sample.colors,
                                  location = location.col,
                                  tumor = cancer.col,
                                  sex = sex.col,
                                  grade = grade.col,
                                  batch = batch.col),
         show_rownames = F, show_colnames = F, 
         color = red_blue_20)

dev.off()


###DEG IDH/GBM
head(seurat[[]])
table(seurat@meta.data$tumor)

seurat@meta.data$IDH <- ifelse(seurat@meta.data$tumor == "GBM", "WT", "Mutant")
seurat.core <- subset(seurat, subset = location == "TC")

seurat.core
Idents(seurat.core) <- seurat.core@meta.data$tumor
Idents(seurat.core) <- seurat.core@meta.data$IDH
table(Idents(seurat.core))

res <- FindAllMarkers(seurat.core, test.use = "t", only.pos = T)
head(res)

###Verhaak subtypes
cgp_gene_sets <- msigdbr::msigdbr(species = "human", category = "C2", subcategory = "CGP") %>% dplyr::distinct(gene_symbol, gs_name) %>% as.data.frame()

#grep("VERHAAK_GLIOBLASTOMA", cgp_gene_sets$gs_name, value = T)
#c5.db = msigdbr::msigdbr(species = "human", category = "C5") %>% dplyr::distinct(gene_symbol, gs_name) %>% as.data.frame()
head(cgp_gene_sets)

VERHAAK_GLIOBLASTOMA <- cgp_gene_sets[grep("VERHAAK_GLIOBLASTOMA", cgp_gene_sets$gs_name),]
head(VERHAAK_GLIOBLASTOMA)
rownames(VERHAAK_GLIOBLASTOMA) <- VERHAAK_GLIOBLASTOMA$gene_symbol
VERHAAK_GLIOBLASTOMA$gs_name <- gsub("VERHAAK_GLIOBLASTOMA_", "", VERHAAK_GLIOBLASTOMA$gs_name)
VERHAAK_GLIOBLASTOMA <- VERHAAK_GLIOBLASTOMA[rownames(VERHAAK_GLIOBLASTOMA) %in% rownames(seurat),1, drop = F]
head(VERHAAK_GLIOBLASTOMA)
table(VERHAAK_GLIOBLASTOMA$gs_name)
colnames(VERHAAK_GLIOBLASTOMA) <- "subtype"
# CLASSICAL MESENCHYMAL      NEURAL   PRONEURAL 
#30          51           4          28
head(seurat.core)
verhaak.subtype.col <- pair.col[c(1,3,5,7)]
names(verhaak.subtype.col) <- unique(VERHAAK_GLIOBLASTOMA$subtype)

col.ann <- seurat.core@meta.data[,c(5,7)]
head(col.ann)
pdf("all.glioma/01.spatial.VerhaakSubtype.pdf", width = 10, height = 10, useDingbats = F)
pheatmap(seurat.core@assays$RNA@scale.data[rownames(VERHAAK_GLIOBLASTOMA),], 
         annotation_row = VERHAAK_GLIOBLASTOMA,
         annotation_col = col.ann,
         annotation_colors = list(tumor = cancer.col, sample = sample.colors, subtype = verhaak.subtype.col),
         show_rownames = F, show_colnames = F,
         color = red_blue_20
)
dev.off()

#######
VizDimLoadings(seurat, dims = 1:2, reduction = "pca")
seurat@reductions$pca@feature.loadings[1:4,1:4]

pc1.top.features <- rownames(seurat@reductions$pca@feature.loadings)[order(seurat@reductions$pca@feature.loadings[,1], decreasing = T )[1:30]]
pc2.top.features <- rownames(seurat@reductions$pca@feature.loadings)[order(seurat@reductions$pca@feature.loadings[,2], decreasing = T )[1:30]]


pc1.top.GO <- enricher(pc1.top.features, TERM2GENE = c5.db)
View(pc1.top.GO@result)
#write.csv(pc1.top.GO@result, file = "AIDH/03.PC1.topLoading.GO.csv", row.names = F)
pc1.top.GO <- pc1.top.GO@result[order(pc1.top.GO@result$p.adjust, decreasing = F)[1:20],]
pc1.top.GO$ID <- gsub("^GO\\w\\w\\_", "", pc1.top.GO$ID)
pc1.top.GO$ID <- gsub("\\_", " ", pc1.top.GO$ID)

p1 <- ggplot(pc1.top.GO, aes(x = reorder(ID, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", fill = "#0072B5", color = "black")+
  coord_flip()+ theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank())+
  #scale_x_discrete(label = function(x) stringr::str_trunc(x, 30))+
  #scale_x_discrete(labels = scales::label_wrap(10))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x,35), expand = expansion(mult = .03)) +
  scale_y_discrete(expand = expansion(mult = c(0, .05))) +
  labs(x= "")
p1    
#pdf("AIDH/03.PC1.topLoading.GO.top20.pdf", width = 10, height = 10, useDingbats = F)
p1
dev.off()

#
pc2.top.GO <- enricher(pc2.top.features, TERM2GENE = c5.db)
View(pc2.top.GO@result)
#write.csv(pc2.top.GO@result, file = "AIDH/03.PC2.topLoading.GO.csv", row.names = F)
pc2.top.GO <- pc2.top.GO@result[order(pc2.top.GO@result$p.adjust, decreasing = F)[1:20],]
pc2.top.GO$ID <- gsub("^GO\\w\\w\\_", "", pc2.top.GO$ID)
pc2.top.GO$ID <- gsub("\\_", " ", pc2.top.GO$ID)
pc2.top.GO$GeneRatio <- as.numeric(as.character(pc2.top.GO$GeneRatio))

p2 <- ggplot(pc2.top.GO, aes(x = reorder(ID, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", fill = "#0072B5", color = "black")+
  coord_flip()+ theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank())+
  #scale_x_discrete(label = function(x) stringr::str_trunc(x, 30))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x,35), expand = expansion(mult = .03)) +
  scale_y_discrete(expand = expansion(mult = c(0, .05))) +
  #scale_x_discrete(labels = scales::label_wrap(10))+
  labs(x= "")
p2
#pdf("AIDH/03.PC2.topLoading.GO.top20.pdf", width = 10, height = 10, useDingbats = F)
p2
dev.off()

pc1.low <- rownames(seurat@reductions$pca@feature.loadings)[order(seurat@reductions$pca@feature.loadings[,1], decreasing = F )[1:30]]
pc1.low.GO <- enricher(pc1.low, TERM2GENE = c5.db)
View(pc1.low.GO@result)

pc2.low <- rownames(seurat@reductions$pca@feature.loadings)[order(seurat@reductions$pca@feature.loadings[,2], decreasing = F )[1:30]]
pc2.low.GO <- enricher(pc2.low, TERM2GENE = c5.db)
View(pc2.low.GO@result)

###
VariableFeatures(seurat)
data <- seurat@assays$RNA@scale.data
dim(data)
data[1:4,1:4]


roi.ann <- data.frame(sample = seurat@meta.data$sample,
                      tumor = seurat@meta.data$tumor,
                      sex = seurat@meta.data$sex,
                      grade = seurat@meta.data$grade,
                      location = seurat@meta.data$location)
rownames(roi.ann) <- rownames(seurat@meta.data)
head(roi.ann)
gene.cluster <- hclust(dist(data, method = "euclidean"), method = "ward.D2")
roi.cluster <- hclust(dist(t(data), method = "euclidean"), method = "ward.D2")

dev.off()
#pdf("AIDH//04.AIDH.variableGenes.heatmap.pdf", width =20, height = 20, useDingbats = F)
#pdf("all.glioma/04.glioma.variableGenes.heatmap.pdf", width =30, height = 30, useDingbats = F)
pheatmap(data, show_rownames = F, show_colnames = F, color = red_blue_20,
         cluster_rows = gene.cluster, cluster_cols = roi.cluster,
         annotation_col = roi.ann, 
         annotation_colors = list(sample = sample.colors, 
                                  tumor = cancer.col,
                                  grade = grade.col, sex = sex.col,
                                  location = location.col),
         cellwidth = 1, cellheight = 1, border_color = NA)
dev.off()

###

#seurat <- FindNeighbors(seurat, dims = 1:10)
#seurat <- FindClusters(seurat, resolution = 0.5)
head(seurat[[]])
table(seurat@meta.data$sex)
seurat.harmony <- RunHarmony(seurat, group.by.vars = c("sample")) #Sex
seurat.harmony <- RunUMAP(seurat.harmony, reduction = "harmony", dims = 1:10, verbose = F, n.neighbors = 20)

seurat.harmony@reductions$harmony[1:4,1:4]

p1 <- DimPlot(seurat.harmony, reduction = "harmony", group.by = "sample", cols = sample.colors, pt.size = 1,shape.by = "location")+ my.thm + labs(title = "")+loc.shape
p2 <- DimPlot(seurat.harmony, reduction = "harmony", group.by = "location", cols = location.col, pt.size = 1)+ my.thm + labs(title = "")
p3 <- DimPlot(seurat.harmony, reduction = "harmony", group.by = "grade", cols = grade.col, pt.size = 1)+ my.thm + labs(title = "")
p4 <- DimPlot(seurat.harmony, reduction = "harmony", group.by = "sex", cols = pair.col[c(5,2)], pt.size = 1)+ my.thm + labs(title = "")

#pdf("all.glioma/01.harmony.PCA.pdf", width = 20, height = 20, useDingbats = F)
#pdf("AIDH/01.AIDH.harmony.PCA.pdf", width = 10, height = 10, useDingbats = F)
#pdf("GBM/01.GBM.harmony.PCA.pdf", width = 10, height = 10, useDingbats = F)
cowplot::plot_grid(p1,p2, p3, p4, nrow = 2)
dev.off()

p1 <- DimPlot(seurat.harmony, reduction = "umap", group.by = "sample", cols = sample.colors, pt.size = 1,shape.by = "location")+ my.thm + labs(title = "")+loc.shape
p2 <- DimPlot(seurat.harmony, reduction = "umap", group.by = "location", cols = location.col, pt.size = 1)+ my.thm + labs(title = "")
p3 <- DimPlot(seurat.harmony, reduction = "umap", group.by = "grade", cols = grade.col, pt.size = 1)+ my.thm + labs(title = "")
p4 <- DimPlot(seurat.harmony, reduction = "umap", group.by = "sex", cols = pair.col[c(5,2)], pt.size = 1)+ my.thm + labs(title = "")

#pdf("all.glioma/01.harmony.umap.pdf", width = 20, height = 20, useDingbats = F)
#pdf("AIDH/01.AIDH.harmony.umap.pdf", width = 10, height = 10, useDingbats = F)
#pdf("GBM/01.GBM.harmony.umap.pdf", width = 15, height = 10, useDingbats = F)
cowplot::plot_grid(p1,p2, p3, p4, nrow = 2)
dev.off()

###harmony batch correction




##cluster and plot, one by one
unique(seurat@meta.data$sample)
unique(seurat@meta.data$sample0)
dev.off()
all.data.files <- data.files.G
all.sample.names <- sample.names.G

loc.shape <- c("TC" =15, "IE" = 2, "GM" = 4, "MVP" = 8, "PN" = 10, "BV" = 13)
loc.shape <- loc.shape[names(loc.shape) %in% seurat@meta.data$location]
loc.shape

out.prefix <- "all.glioma/02.glima"
#out.prefix <- "AIDH/02.AIDH."
#out.prefix <- "GBM/02.GBM."
dev.off()
all.spatial.plots <- lapply(1:length(unique(seurat@meta.data$sample)), function(i){
  #i = 3
  sample.id <- unique(seurat@meta.data$sample)[i]
  message("Runnning ", sample.id)
  sub.seurat <- subset(seurat, sample == sample.id)
  sub.seurat@meta.data$ROI <- sapply(rownames(sub.seurat@meta.data), function(x) unlist(strsplit(x, "\\-"))[2])
  
  sub.seurat <- RunPCA(sub.seurat, npcs = 20, verbose = F)
  n.k = 20 #default
  n.k = 15
  #if(ncol(sub.seurat) < 30) n.k <- 15
  sub.seurat <- RunUMAP(sub.seurat, dims = 1:10, n.neighbors = n.k, verbose = F)
  sub.seurat <- FindNeighbors(sub.seurat, dims = 1:10, k.param = n.k)
  sub.seurat <- FindClusters(sub.seurat, resolution = 0.6)
  sub.seurat@meta.data$cluster <- sub.seurat@meta.data$RNA_snn_res.0.6
  
  deg <- FindAllMarkers(sub.seurat, test.use = "t", only.pos = T)
  deg <- deg[deg$p_val_adj <0.01,]
  if(nrow(deg) > 5){
    roi.ann <- sub.seurat@meta.data[,c("location","cluster"), drop =F]
    keep.genes <- unique(deg$gene)
    gene.ann <- data.frame(cluser = as.character(deg$cluster[match(keep.genes, deg$gene)]))
    rownames(gene.ann) <- keep.genes
    multi.over.genes <- unique(names(which(table(deg$gene)>1)))
    gene.ann$cluser[rownames(gene.ann) %in% multi.over.genes] <- "Multi"
    
    pdf(file = paste0(out.prefix, sample.id, ".clusterDEG.heatmap.pdf"), width = 10, height = 10, useDingbats = F)
    print(pheatmap(sub.seurat@assays$RNA@data[unique(deg$gene),], 
                   annotation_col = roi.ann, 
                   annotation_row = gene.ann,
                   #annotation_colors = list(cluster = pair.col),
                   scale = "row",
                   show_colnames = F, show_rownames = F,
                   color = red_blue_20))
    dev.off()
    write.csv(deg, file = paste0(out.prefix, sample.id, ".DEG.csv"))
    x <- lapply(unique(deg$cluster), function(k){
      #k =2
      g <- unlist(sapply(rownames(deg)[deg$cluster == k], function(x) strsplit(x, "\\/")))
      if(length(g) < 5) return(1)
      e0 <- enricher(gene = g, TERM2GENE = c5.db)
      write.csv(e0@result, file = paste0(out.prefix, "", sample.id, ".cluster", k, ".GO.csv"), row.names = F)
      1;
    })
  }
  sub.loc.col <- location.col[names(location.col) %in% sub.seurat@meta.data$location]
  sub.loc.shape = loc.shape
  sub.loc.shape <- loc.shape[names(loc.shape) %in% sub.seurat@meta.data$location]
  message(sub.loc.shape)
  sub.loc.shape <- scale_shape_manual(values = sub.loc.shape)
  message(sub.loc.shape$limits)
  p1 <- DimPlot(sub.seurat, reduction = "pca", group.by = "location", cols = sub.loc.col, shape.by = "location", pt.size = 1)+ my.thm + 
    labs(title = "PCA", x= "", y = "")+sub.loc.shape
  p2 <- DimPlot(sub.seurat, reduction = "umap", group.by = "location", cols = sub.loc.col, shape.by = "location",pt.size = 1)+ my.thm + 
    labs(title = "UMAP", x= "", y = "")+sub.loc.shape
  p3 <- DimPlot(sub.seurat, reduction = "umap", group.by = "cluster", shape.by = "location", pt.size = 1)+ my.thm + 
    labs(title = "UMAP", x= "", y = "") + scale_colour_brewer(palette = "Dark2")+sub.loc.shape
  if(length(unique(sub.seurat@meta.data$Batch)) > 1){
    p1.1 <- DimPlot(sub.seurat, reduction = "pca", group.by = "batch", shape.by = "location", pt.size = 1)+ my.thm + labs(title = "", x= "", y = "")+sub.loc.shape
    p2.1 <- DimPlot(sub.seurat, reduction = "umap", group.by = "batch", shape.by = "location",pt.size = 1)+ my.thm + labs(title = "", x= "", y = "")+sub.loc.shape
    
    pdf(paste0(out.prefix, sample.id, ".batch.pca.umap.pdf"), width = 10,height = 5, useDingbats = F)
    print(cowplot::plot_grid(p1.1, p2.1, nrow = 1))
    dev.off()
  }
  sample.id.all.batch <- as.character(unique(sub.seurat@meta.data$sample0))
  
  spatial.plots <- lapply(1:length(sample.id.all.batch), function(k){
    #k = 1
    sample.id.x <- sample.id.all.batch[k]
    message("Spatial plot for ", sample.id.x)
    #idx <- which(all.sample.names == sample.id.x)
    #idx
    #data.file <- all.data.files[idx]
    #message(sample.id.x, ":", data.file)
    
    #xy <- read_nanostring_RNA_XY(data.file)
    #xy$cluster <- sub.seurat@meta.data$cluster[match(xy$ROILabel, sub.seurat@meta.data$ROI)]
    #xy$Location <- sub.seurat@meta.data$Location[match(xy$ROILabel, sub.seurat@meta.data$ROI)]
    xy <- sub.seurat@meta.data
    xy <- xy[xy$sample0 == sample.id.x,]
    #head(xy)
    p <- ggplot(xy, aes(x = spatialY, y = -spatialX, col = cluster, shape = location))+geom_point() + 
      my.thm+labs(x = "", y = "", title = sample.id.x)+scale_colour_brewer(palette = "Dark2")+sub.loc.shape
    p
  })
  
  pdf(file = paste0(out.prefix, sample.id, ".pca.umap.cluster.pdf"), width = 10, height = 10, useDingbats = F)
  print(cowplot::plot_grid(plotlist = append(list(p1,p2,p3), spatial.plots)))
  dev.off()
  
  spatial.plots
  
})

all.spatial.plots.1 <- do.call(c, all.spatial.plots)

pdf(paste0(out.prefix, "all.spatial.cluster.pdf"), width = 15, height = 10, useDingbats = F)
cowplot::plot_grid(plotlist = all.spatial.plots.1)
dev.off()


### common entiched GO
all.go.files <- list.files("AIDH/", pattern = "02.*.GO.csv$", full.names = T)
all.go <- lapply(all.go.files, function(f){
  xx <- unlist(strsplit(f, split = "\\."))
  sample.id <- xx[3]
  k <- gsub("cluster", "", xx[4])
  message(sample.id, "-", k)
  x <- read.csv(f, head =T)
  x$cluster <- paste0(sample.id, "-", k)
  x <- x[x$p.adjust <0.05 & x$Count >=10, ]
  x
})
unlist(lapply(all.go, length))
head(all.go[[3]])
head(all.go[[2]])

all.go <- data.frame(do.call(rbind, all.go))
x <- data.frame(table(all.go$ID))
quantile(x$Freq)
quantile(x$Freq, 0.85)

p <- ggplot(x, aes(x = Freq)) + 
  #geom_line(stat = "density")+
  geom_density()+
  theme_bw()+
  geom_vline(xintercept = 5, linetype = "dashed")+
  theme(aspect.ratio = 1, panel.grid = element_blank())+
  labs(x = "Frequency", y = "density")
p

pdf("AIDH/05.clusters.GO.frequency.pdf", width = 5, height = 4, useDingbats = F)
p
dev.off()

#4 is 85% quantile
sum(x$Freq >=5)
keep.go <- as.character(x$Var1[x$Freq >= 5])

all.go.keep <- all.go[all.go$ID %in% keep.go,]

y.order <- sort(table(all.go.keep$ID))
p <- ggplot(all.go.keep, aes(x= cluster, y = ID, fill = -log10(p.adjust)))+
  geom_tile()+
  theme_bw()+
  scale_y_discrete(limits = names(y.order),
                   label = function(x) stringr::str_trunc(x, 50))+
  theme(axis.text.x = element_text(angle = 90), panel.grid = element_blank())+
  scale_fill_viridis_c(option = "plasma")

pdf("AIDH/06.cluster.GO.share.pdf", width = 15, height = 20, useDingbats = F)
p
dev.off()

keep.go.list <- lapply(keep.go, function(id){
  genes <- all.go$geneID[which(all.go$ID == id)]
  genes <- unlist(strsplit(genes, split = "../"))
  genes
})
names(keep.go.list) <- keep.go

keep.go.scores <- UCell::ScoreSignatures_UCell(data, keep.go.list)
colnames(keep.go.scores) <- gsub("_UCell", "", colnames(keep.go.scores))
keep.go.scores[1:4,1:4]

identical(colnames(data), rownames(seurat@meta.data))
colnames(seurat@meta.data)
roi.ann <- seurat@meta.data[,c(5,6,8,9)]
pheatmap(t(keep.go.scores), scale = "row", show_colnames = F, color = red_blue_20,
         annotation_col = roi.ann,
         annotation_colors = list(sample = sample.colors,
                                  sex = sex.col,
                                  location = location.col),
         fontsize_row = 5,
         cutree_rows = 5)

pc <- prcomp(keep.go.scores)
pc$x[1:4,1:4]
res <- data.frame(pc$x[,1:5])
res$sample <- seurat@meta.data$sample
ggplot(res, aes(x = PC1, y = PC2, col = sample))+
  geom_point()+
  scale_color_manual(values = sample.colors)+
  theme_bw()

###TC only
seurat <- subset(seurat, subset = Location == "TC")
seurat
unique(seurat$sample)

variable.features <- lapply(unique(seurat$sample), function(k){
  sub.data <- subset(seurat, subset = sample == k)
  sub.data <- FindVariableFeatures(sub.data, nfeatures = 100)
  VariableFeatures(sub.data)
  #x <- HVFInfo(sub.data)#[VariableFeatures(sub.data), ]
  #x$gene <- rownames(x)
  #x$sample <- k
  #x
})

dim(variable.features[[1]])
variable.features <- unlist(variable.features)
plot(density(table(variable.features)))

quantile(table(variable.features))
sum(table(variable.features)>4)

most.vst <- names(which(table(variable.features) >4))

data <- seurat@assays$RNA@data[most.vst,]
head(seurat@meta.data)
data.ann <- seurat@meta.data[,c(5,8,9)]
head(data.ann)

pheatmap(data, scale = "row",
         annotation_col = data.ann, 
         show_colnames = F,
         annotation_colors = list(sample = sample.colors, Grade = grade.col, Sex = sex.col))


sample.here <- unique(seurat@meta.data$sample)[1]
sub.data <- subset(seurat, subset = sample == sample.here)
identical(colnames(sub.data), rownames(sub.data@meta.data))
sub.data <- cbind(sub.data@meta.data, sub.data@assays$RNA@data[most.vst[1],])
colnames(sub.data)[ncol(sub.data)] = most.vst[1]
head(sub.data)

ggplot(sub.data, aes_string(x = "spatialX", y = "spatialY", col = most.vst[1]))+
  geom_point()+
  my.thm+
  scale_color_viridis_c()





###Analyze sample, one by one

i = 1
sample.id <- unique(seurat@meta.data$sample)[i]
message("Runnning ", sample.id)
sub.seurat <- subset(seurat, sample == sample.id)
sub.seurat@meta.data$ROI <- sapply(rownames(sub.seurat@meta.data), function(x) unlist(strsplit(x, "\\-"))[2])
sub.seurat <- RunPCA(sub.seurat, npcs = 20, verbose = F)
#n.k = 20 #default
n.k = 15
#if(ncol(sub.seurat) < 30) n.k <- 15
sub.seurat <- RunUMAP(sub.seurat, dims = 1:10, n.neighbors = n.k, verbose = F)
sub.seurat <- FindNeighbors(sub.seurat, dims = 1:10, k.param = n.k)
sub.seurat <- FindClusters(sub.seurat, resolution = 0.8)
sub.seurat@meta.data$cluster <- sub.seurat@meta.data$RNA_snn_res.0.8

sub.loc.col <- location.col[names(location.col) %in% sub.seurat@meta.data$location]
sub.loc.shape = loc.shape
sub.loc.shape <- loc.shape[names(loc.shape) %in% sub.seurat@meta.data$location]
message(sub.loc.shape)
sub.loc.shape <- scale_shape_manual(values = sub.loc.shape)
message(sub.loc.shape$limits)
p1 <- DimPlot(sub.seurat, reduction = "pca", group.by = "location", cols = sub.loc.col, shape.by = "location", pt.size = 1)+ my.thm + 
  labs(title = "PCA", x= "", y = "")+sub.loc.shape
p2 <- DimPlot(sub.seurat, reduction = "umap", group.by = "location", cols = sub.loc.col, shape.by = "location",pt.size = 1)+ my.thm + 
  labs(title = "UMAP", x= "", y = "")+sub.loc.shape
p3 <- DimPlot(sub.seurat, reduction = "umap", group.by = "cluster", shape.by = "location", pt.size = 1)+ my.thm + 
  labs(title = "UMAP", x= "", y = "") + scale_colour_brewer(palette = "Dark2")+sub.loc.shape
cowplot::plot_grid(p1,p2,p3, nrow = 1)

#deg <- FindMarkers(sub.seurat, ident.1 = 0, ident.2 = 2)
deg <- FindAllMarkers(sub.seurat, test.use = "t", only.pos = T)
head(deg)
nrow(deg)
table(deg$cluster)
e0 <- enricher(gene = unlist(sapply(rownames(deg)[deg$cluster == 0], function(x) strsplit(x, "\\/"))), TERM2GENE = c5.db)
e1 <- enricher(gene = unlist(sapply(rownames(deg)[deg$cluster == 1], function(x) strsplit(x, "\\/"))), TERM2GENE = c5.db)
e2 <- enricher(gene = unlist(sapply(rownames(deg)[deg$cluster == 2], function(x) strsplit(x, "\\/"))), TERM2GENE = c5.db)
e0@result$ID[1:10]
e1@result$ID[1:10]
e2@result$ID[1:10]
write.csv(deg, file = paste0("05.", sample.id, ".DEG.csv"))
write.csv(e0@result, file = paste0("AIDH/05.", sample.id, ".cluster0.GO.csv"), row.names = F)
write.csv(e1@result, file = paste0("AIDH/05.", sample.id, ".cluster1.GO.csv"), row.names = F)
write.csv(e2@result, file = paste0("AIDH/05.", sample.id, ".cluster2.GO.csv"), row.names = F)

###
i = 2
sample.id <- unique(seurat@meta.data$sample)[i]
message("Runnning ", sample.id)
sub.seurat <- subset(seurat, sample == sample.id)
sub.seurat@meta.data$ROI <- sapply(rownames(sub.seurat@meta.data), function(x) unlist(strsplit(x, "\\-"))[2])
sub.seurat <- RunPCA(sub.seurat, npcs = 20, verbose = F)
#n.k = 20 #default
n.k = 15
#if(ncol(sub.seurat) < 30) n.k <- 15
sub.seurat <- RunUMAP(sub.seurat, dims = 1:10, n.neighbors = n.k, verbose = F)
sub.seurat <- FindNeighbors(sub.seurat, dims = 1:10, k.param = n.k)
sub.seurat <- FindClusters(sub.seurat, resolution = 0.8)
sub.seurat@meta.data$cluster <- sub.seurat@meta.data$RNA_snn_res.0.8

sub.loc.col <- location.col[names(location.col) %in% sub.seurat@meta.data$location]
sub.loc.shape = loc.shape
sub.loc.shape <- loc.shape[names(loc.shape) %in% sub.seurat@meta.data$location]
message(sub.loc.shape)
sub.loc.shape <- scale_shape_manual(values = sub.loc.shape)
message(sub.loc.shape$limits)
p1 <- DimPlot(sub.seurat, reduction = "pca", group.by = "location", cols = sub.loc.col, shape.by = "location", pt.size = 1)+ my.thm + 
  labs(title = "PCA", x= "", y = "")+sub.loc.shape
p2 <- DimPlot(sub.seurat, reduction = "umap", group.by = "location", cols = sub.loc.col, shape.by = "location",pt.size = 1)+ my.thm + 
  labs(title = "UMAP", x= "", y = "")+sub.loc.shape
p3 <- DimPlot(sub.seurat, reduction = "umap", group.by = "cluster", shape.by = "location", pt.size = 1)+ my.thm + 
  labs(title = "UMAP", x= "", y = "") + scale_colour_brewer(palette = "Dark2")+sub.loc.shape
cowplot::plot_grid(p1,p2,p3, nrow = 1)

#deg <- FindMarkers(sub.seurat, ident.1 = 0, ident.2 = 2)
deg <- FindAllMarkers(sub.seurat, test.use = "t", only.pos = T)
head(deg)
nrow(deg)
table(deg$cluster)
e0 <- enricher(gene = unlist(sapply(rownames(deg)[deg$cluster == 0], function(x) strsplit(x, "\\/"))), TERM2GENE = c5.db)
e1 <- enricher(gene = unlist(sapply(rownames(deg)[deg$cluster == 1], function(x) strsplit(x, "\\/"))), TERM2GENE = c5.db)
e2 <- enricher(gene = unlist(sapply(rownames(deg)[deg$cluster == 2], function(x) strsplit(x, "\\/"))), TERM2GENE = c5.db)
e0@result$ID[1:10]
e1@result$ID[1:10]
e2@result$ID[1:10]
write.csv(deg, file = paste0("05.", sample.id, ".DEG.csv"))
#write.csv(e0@result, file = paste0("AIDH/05.", sample.id, ".cluster0.GO.csv"), row.names = F)
write.csv(e1@result, file = paste0("AIDH/05.", sample.id, ".cluster1.GO.csv"), row.names = F)
write.csv(e2@result, file = paste0("AIDH/05.", sample.id, ".cluster2.GO.csv"), row.names = F)



###
i = 3
sample.id <- unique(seurat@meta.data$sample)[i]
message("Runnning ", sample.id)
sub.seurat <- subset(seurat, sample == sample.id)
sub.seurat@meta.data$ROI <- sapply(rownames(sub.seurat@meta.data), function(x) unlist(strsplit(x, "\\-"))[2])
sub.seurat <- RunPCA(sub.seurat, npcs = 20, verbose = F)
#n.k = 20 #default
n.k = 15
#if(ncol(sub.seurat) < 30) n.k <- 15
sub.seurat <- RunUMAP(sub.seurat, dims = 1:10, n.neighbors = n.k, verbose = F)
sub.seurat <- FindNeighbors(sub.seurat, dims = 1:10, k.param = n.k)
sub.seurat <- FindClusters(sub.seurat, resolution = 0.8)
sub.seurat@meta.data$cluster <- sub.seurat@meta.data$RNA_snn_res.0.8

sub.loc.col <- location.col[names(location.col) %in% sub.seurat@meta.data$location]
sub.loc.shape = loc.shape
sub.loc.shape <- loc.shape[names(loc.shape) %in% sub.seurat@meta.data$location]
message(sub.loc.shape)
sub.loc.shape <- scale_shape_manual(values = sub.loc.shape)
message(sub.loc.shape$limits)
p1 <- DimPlot(sub.seurat, reduction = "pca", group.by = "location", cols = sub.loc.col, shape.by = "location", pt.size = 1)+ my.thm + 
  labs(title = "PCA", x= "", y = "")+sub.loc.shape
p2 <- DimPlot(sub.seurat, reduction = "umap", group.by = "location", cols = sub.loc.col, shape.by = "location",pt.size = 1)+ my.thm + 
  labs(title = "UMAP", x= "", y = "")+sub.loc.shape
p3 <- DimPlot(sub.seurat, reduction = "umap", group.by = "cluster", shape.by = "location", pt.size = 1)+ my.thm + 
  labs(title = "UMAP", x= "", y = "") + scale_colour_brewer(palette = "Dark2")+sub.loc.shape
cowplot::plot_grid(p1,p2,p3, nrow = 1)

#deg <- FindMarkers(sub.seurat, ident.1 = 0, ident.2 = 2)
deg <- FindAllMarkers(sub.seurat, test.use = "t", only.pos = T)
head(deg)
nrow(deg)
head(sub.seurat[[]])
roi.ann <- sub.seurat@meta.data[,16, drop =F]
head(roi.ann)
gene.ann <- data.frame(cluser = deg$cluster)
rownames(gene.ann) <- deg$gene

pheatmap(sub.seurat@assays$RNA@data[unique(deg$gene),], annotation_col = roi.ann, scale = "row",
         show_colnames = F, show_rownames = F,
         color = red_blue_20)

table(deg$cluster)
e0 <- enricher(gene = unlist(sapply(rownames(deg)[deg$cluster == 0], function(x) strsplit(x, "\\/"))), TERM2GENE = c5.db)
e1 <- enricher(gene = unlist(sapply(rownames(deg)[deg$cluster == 1], function(x) strsplit(x, "\\/"))), TERM2GENE = c5.db)
e2 <- enricher(gene = unlist(sapply(rownames(deg)[deg$cluster == 2], function(x) strsplit(x, "\\/"))), TERM2GENE = c5.db)
e3 <- enricher(gene = unlist(sapply(rownames(deg)[deg$cluster == 3], function(x) strsplit(x, "\\/"))), TERM2GENE = c5.db)

e0@result$ID[1:10]
e1@result$ID[1:10]
e2@result$ID[1:10]
e3@result$ID[1:10]

write.csv(deg, file = paste0("05.", sample.id, ".DEG.csv"))
write.csv(e0@result, file = paste0("AIDH/05.", sample.id, ".cluster0.GO.csv"), row.names = F)
write.csv(e1@result, file = paste0("AIDH/05.", sample.id, ".cluster1.GO.csv"), row.names = F)
write.csv(e2@result, file = paste0("AIDH/05.", sample.id, ".cluster2.GO.csv"), row.names = F)
write.csv(e3@result, file = paste0("AIDH/05.", sample.id, ".cluster3.GO.csv"), row.names = F)

###
