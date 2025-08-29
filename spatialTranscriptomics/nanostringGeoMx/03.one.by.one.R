library(Seurat)
setwd("/Users/wuz6/Documents/Project/08.spatialTranscriptome/01.RNA/01.QCfirst/")
rm(list=ls())
red_blue_30 <- rev(as.character(paletteer::paletteer_c("ggthemes::Classic Red-Blue", 30)))
pair.col = RColorBrewer::brewer.pal(12, "Paired")
dark2.col = RColorBrewer::brewer.pal(8, "Dark2")

library(ComplexHeatmap)
library(pheatmap)
library(clusterProfiler)


c5.db = msigdbr::msigdbr(species = "human", category = "C5") %>% dplyr::distinct(gene_symbol, gs_name) %>% as.data.frame()


my.thm <- theme_bw()+ theme(aspect.ratio = 1, panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
#
location.col = c("#808080", pair.col[c(6,10,2,12)], "#000000")
names(location.col) = c("GM", "BV", "MVP", "IE", "PN", "TC")
location.col

grade.col <- c("#BDBDBD", "#616161", "#212121")
names(grade.col) <- c("II", "III", "IV")
sex.col <- pair.col[c(5,2)]
names(sex.col) <- c("F", "M")
#


################################################################################################
load("CTA.allGlioma.seuratObj.rda")
load("CTA.colors.shape.rda")

###
if(0){
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
sample.colors
}
##cluster and plot, one by one
unique(seurat@meta.data$sample)
unique(seurat@meta.data$sample0)
dev.off()

all.sample.names <- names(sample.colors)
all.sample.names

loc.shape <- c("TC" =15, "IE" = 2, "GM" = 4, "MVP" = 8, "PN" = 10, "BV" = 13)
loc.shape <- loc.shape[names(loc.shape) %in% seurat@meta.data$location]
loc.shape

library(clustree)
dir.create("one.by.one")
out.prefix <- "one.by.one/01."


y.n.col <- c(pair.col[4], "#808080"); names(y.n.col) <- c("Y", "N")
#
sample.info <- data.frame(readxl::read_xlsx("../../00.data/Brain Samples Selected for Spatial Profiling.xlsx"))  
sample.info$Seurat.res[match(unique(seurat@meta.data$sample), sample.info$ID)]

dev.off();
all.spatial.plots <- lapply(1:length(unique(seurat@meta.data$sample)), function(i){
  #i = 1
  sample.id <- unique(seurat@meta.data$sample)[i]
  message("Runnning ", sample.id)
  
  sub.seurat <- subset(seurat, sample == sample.id)
  sub.seurat <- FindVariableFeatures(sub.seurat, nfeatures = 300)
  sub.seurat@meta.data$ROI <- sapply(rownames(sub.seurat@meta.data), function(x) unlist(strsplit(x, "\\_"))[2])
  #if()
  sub.seurat <- RunPCA(sub.seurat, npcs = 20, features = VariableFeatures(sub.seurat), verbose = F)
  #n.k = 20 #default
  n.k = 15
  #if(ncol(sub.seurat) < 30) n.k <- 15
  sub.seurat <- RunUMAP(sub.seurat, dims = 1:10, n.neighbors = n.k, verbose = F)
  sub.seurat <- FindNeighbors(sub.seurat, dims = 1:10, k.param = n.k)
  sub.seurat <- FindClusters(sub.seurat, resolution = seq(0.2, 1.0, 0.1))
  if(0){
    pdf(file = paste0(out.prefix, sample.id, ".clustree.pdf"), width = 10, height = 10, useDingbats = F)
    print(clustree(sub.seurat, prefix = "RNA_snn_res."))
    dev.off()
    return(1)
  }
  
  
  select.res <- paste0("RNA_snn_res.", sample.info$Seurat.res[which(sample.info$ID == sample.id)])
  message(select.res)
  sub.seurat@meta.data$cluster <- sub.seurat@meta.data[[select.res]] ##!!!!!!!
  Idents(sub.seurat) <-  sub.seurat@meta.data$cluster
  x <- sub.seurat@meta.data
  write.csv(x, file = paste0(out.prefix, sample.id, ".cluster.csv"))
  cluster.col <- dark2.col[1:length(unique(sub.seurat@meta.data$cluster))]
  names(cluster.col) <- sort(unique(sub.seurat@meta.data$cluster))
  
  sub.seurat.tc <- subset(sub.seurat, location == "TC")
  deg <- FindAllMarkers(sub.seurat.tc, test.use = "t", only.pos = T)
  deg <- deg[deg$p_val_adj <0.01,]
  
  
  #roi.ann <- sub.seurat.tc@meta.data[,c("location","cluster"), drop =F]
  
          
  roi.ann <- sub.seurat.tc@meta.data[,c("cluster"), drop =F]
  roi.top.ann <- HeatmapAnnotation(cluster = sub.seurat.tc@meta.data$cluster,col = list(cluster = cluster.col))
  if(1){
    #sub.seurat <- FindVariableFeatures(sub.seurat, nfeatures = 100)
    #dev.off()
    m <- t(scale(t(data.matrix(sub.seurat.tc@assays$RNA@data[VariableFeatures(sub.seurat),]))))
    m[m>4] <-4
    m[m< -4] <- -4
    
    h <- Heatmap(m, name = " ",
                                 top_annotation = roi.top.ann,
                 row_dend_reorder = T, column_dend_reorder = T,
                 show_column_names = F, show_row_names = F, col = red_blue_30)
    
    pdf(file = paste0(out.prefix, sample.id, ".300.heatmap.pdf"), width = 10, height = 10, useDingbats = F)
    draw(h)
    dev.off()

  }
  #return(1)
  
  if(nrow(deg) > 5){
    
    keep.genes <- unique(deg$gene)
    gene.ann <- data.frame(cluster = as.character(deg$cluster[match(keep.genes, deg$gene)]))
    rownames(gene.ann) <- keep.genes
    multi.over.genes <- unique(names(which(table(deg$gene)>1)))
    gene.ann$cluster[rownames(gene.ann) %in% multi.over.genes] <- "Multi"
    if("Multi" %in% gene.ann$cluster) cluster.col['Multi'] <- "#000000"#dark2.col[length(cluster.col) + 1]
    
    gene.ann <- rowAnnotation(cluster = gene.ann$cluster, col = list(cluster = cluster.col))
    m <- t(scale(t(data.matrix(sub.seurat.tc@assays$RNA@data[unique(deg$gene),]))))
    m[m>4] <-4
    m[m< -4] <- -4
    h <- Heatmap(m, name = " ",
                 top_annotation = roi.top.ann,
                 left_annotation = gene.ann,
                 row_dend_reorder = T, column_dend_reorder = T,
                 show_column_names = F, show_row_names = F, col = red_blue_30)
    

    pdf(file = paste0(out.prefix, sample.id, ".clusterDEG.heatmap.pdf"), width = 10, height = 10, useDingbats = F)
    draw(h)
    dev.off()
    if(nrow(deg) > 2) write.csv(deg, file = paste0(out.prefix, sample.id, ".DEG.csv"))
    x <- lapply(unique(deg$cluster), function(k){
      #k =2
      g <- unlist(sapply(rownames(deg)[deg$cluster == k], function(x) strsplit(x, "\\/")))
      if(length(g) < 10) return(1)
      e0 <- enricher(gene = g, TERM2GENE = c5.db)@result
      e0 <- e0[e0$p.adjust < 0.05,]
      write.csv(e0, file = paste0(out.prefix, "", sample.id, ".cluster", k, ".GO.csv"), row.names = F)
      1;
    })
  }
  ###
  
  sub.loc.col <- location.col[names(location.col) %in% sub.seurat@meta.data$location]
  sub.loc.shape <- loc.shape
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
    p <- ggplot(xy, aes(x = spatialX, y = -spatialY, col = cluster, shape = location))+geom_point() + 
      my.thm+labs(x = "", y = "", title = sample.id.x)+scale_colour_brewer(palette = "Dark2")+sub.loc.shape
    p
  })
  
  pdf(file = paste0(out.prefix, sample.id, ".pca.umap.cluster.pdf"), width = 10, height = 10, useDingbats = F)
  print(cowplot::plot_grid(plotlist = append(list(p1,p2,p3), spatial.plots)))
  dev.off()
  
  spatial.plots 
  
})
###
all.spatial.plots.1 <- do.call(c, all.spatial.plots)

pdf(paste0(out.prefix, "all.spatial.cluster.pdf"), width = 15, height = 10, useDingbats = F)
cowplot::plot_grid(plotlist = all.spatial.plots.1)
dev.off()



########################################################################

###ENDed here






####common DEG?
all.deg.files <- list.files("one.by.one/", pattern = "*.DEG.csv", full.names = T)
all.deg.files
all.deg <- lapply(all.deg.files, function(f){
  xx <- unlist(strsplit(f, split = "\\."))
  sample.id <- xx[4]
  #k <- gsub("cluster", "", xx[4])
  message(sample.id)
  x <- read.csv(f, head =T, row.names = 1)
  x$cluster <- paste0(sample.id, "_", x$cluster)
  x <- x[x$p_val_adj <0.01 & x$avg_log2FC > 0.2, ]
  if(nrow(x) < 1) return(NULL)
  x
})
tail(all.deg[[1]])

all.deg.genes <- unique(unlist(lapply(1:length(all.deg), function(x){all.deg[[x]]$gene})))
length(all.deg.genes)

all.deg.upset <- lapply(all.deg, function(x){ return(ifelse(all.deg.genes %in% x$gene, 1, 0)) })
all.deg.upset <- data.frame(do.call(cbind, all.deg.upset))
colnames(all.deg.upset) <- sapply(all.deg, function(x) unique(x$cluster))

#pdf("test.pdf")
UpSetR::upset(all.deg.upset)
dev.off()

all.deg.df <- data.frame(do.call(rbind, all.deg))
m <- table(all.deg.df$gene, all.deg.df$cluster)
m <- t(m) %*% m
pheatmap(m)

##jaccard index
all.deg.jac <- lapply(unique(all.deg.df$cluster), function(x){
  x1 <- all.deg.df[all.deg.df$cluster == x,]
  some.jac <- lapply(unique(all.deg.df$cluster), function(y){
      y1 <- all.deg.df[all.deg.df$cluster == y,]
      overlap <- intersect(x1$gene, y1$gene)
      all.genes <- unique(c(x1$gene, y1$gene))
      jac <- length(overlap)/length(all.genes)
      #jac
      return(c(x1$cluster[1], y1$cluster[1], jac))
  })
  some.jac <- do.call(rbind, some.jac)
  some.jac
})

all.deg.jac <- data.frame(do.call(rbind, all.deg.jac))
colnames(all.deg.jac) <- c("C1", "C2", "Jaccard")
head(all.deg.jac)

m <- reshape2::dcast(all.deg.jac, C1~C2)
m <- apply(m[,-1], 2, as.numeric)
rownames(m) <- colnames(m)
m[1:4,1:4]

m.cluster <- hclust(as.dist(1-m), method = "ward.D2")
pdf("one.by.one/02.cluster.DEG.jac.pdf", width = 10, height = 9, useDingbats = F)
pheatmap(m, cluster_rows = m.cluster, cluster_cols = m.cluster, color = red_blue_30)
dev.off();dev.off();dev.off();

identical(rownames(m), colnames(m))
idx <- grep("^G", colnames(m))
idx
m1 <- m[idx, idx]
m1.cluster <- hclust(as.dist(1-m1), method = "ward.D2")
pdf("one.by.one/02.cluster.DEG.GBM.jac.pdf", width = 10, height = 9, useDingbats = F)
pheatmap(m1, cluster_rows = m1.cluster, cluster_cols = m1.cluster, color = red_blue_30)
dev.off();dev.off();dev.off();

idx <- grep("^A", colnames(m))
idx
m1 <- m[idx, idx]
m1.cluster <- hclust(as.dist(1-m1), method = "ward.D2")
pdf("one.by.one/02.cluster.DEG.IDHA.jac.pdf", width = 10, height = 9, useDingbats = F)
pheatmap(m1, cluster_rows = m1.cluster, cluster_cols = m1.cluster, color = red_blue_30)
dev.off();dev.off();dev.off();

idx <- grep("^O", colnames(m))
idx
m1 <- m[idx, idx]
m1.cluster <- hclust(as.dist(1-m1), method = "ward.D2")
pdf("one.by.one/02.cluster.DEG.IDHO.jac.pdf", width = 10, height = 9, useDingbats = F)
pheatmap(m1, cluster_rows = m1.cluster, cluster_cols = m1.cluster, color = red_blue_30)
dev.off();dev.off();dev.off();

### common entiched GO/ or genes?
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
