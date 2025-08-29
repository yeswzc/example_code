rm(list = ls())
setwd("/Users/wuz6/Documents/Project/08.spatialTranscriptome/01.RNA/01.edgeR.QCfirst/")

#library(dbscan)
library(Seurat)
library(NMF)
library(ggplot2)
source("../../src/read_nanostring.R")

###
meta = data.frame(readxl::read_xlsx("../../00.data/GeoMx DSP ROI Groups Annotations.xlsx", sheet = "ROIs_Primary Annotations"), check.names = F)
head(meta)
meta <- reshape2::melt(meta, id.vars = c("roi"))
head(meta)
colnames(meta) = c("ROI", "sample", "annotation")
meta$ID = paste0(meta$sample, "-", meta$ROI)
unique(meta$sample)
meta$ID = gsub("\\-1", "", meta$ID)
meta$ID
table(meta$annotation)
meta$annotation[meta$annotation == "NG"] = "GM"
meta$annotation1 <- gsub("\\-\\w", "", meta$annotation)

pair.col <- RColorBrewer::brewer.pal(12, "Paired")
x.col = c("#808080", pair.col[6])
names(x.col) = c("GM", "BV")
cancer.color <- gplots::col2hex(c('chocolate', 'chocolate2', 'coral', 'coral2', 'brown', 'brown2'))
names(cancer.color) <- paste0("TC-", LETTERS[1:6])
cancer.color
ie.color <- gplots::col2hex(c('cyan', 'cyan2','cadetblue', 'cadetblue2', 'deepskyblue', 'deepskyblue2', 
                              'darkslategray1', 'darkslategray3', 'cornflowerblue'))
names(ie.color) <- paste0("IE-", LETTERS[1:9])
ie.color
loc.col <- c(x.col, cancer.color, ie.color) 
loc.col

###
data.files = list.files("../../00.data/CTA/raw//", pattern = "C*xlsx$", full.names = T)
data.files
sample.names = c(1399, 1905, 2023, 2751, 3984, 4461, 604, 6614, "8147M", "8147O")

sample.colors = RColorBrewer::brewer.pal(12, "Paired")[c(1:10)]
names(sample.colors) = sample.names

###
data <- read.csv("00.IDHastrocytoma.CTA.cpm.filtered.csv", head =T, row.names = 1, check.names = F)
data[1:4,1:4]

###


all.p.list <- lapply(sample.names, function(id){
  set.seed(123)
  #id <- sample.names[1]
  cta.file <- grep(id, data.files, value = T)
  message(id, ": ", cta.file)
  sub.data <- data[,grep(id, colnames(data))]
  message("Sub data nrows: ", dim(sub.data)[1], ", ncols: ", dim(sub.data)[2])

  
  res <- nmf(sub.data, rank = 10, method = "brunet", seed= 'nndsvd')
  h <- coef(res)
  rownames(h) <- paste0("nmf", 1:10)
  #w <- basis(res)
  #colnames(w) <- colnames(h)
  #s <- featureScore(res)
  #s <- extractFeatures(res)

  spatial.coord <- read_nanostring_RNA_XY(cta.file)
  spatial.coord$ID <- paste0(id, "-", spatial.coord$ROILabel)
  spatial.coord$location <- meta$annotation1[match(spatial.coord$ID, meta$ID)]
  spatial.coord <- spatial.coord[match(colnames(h), spatial.coord$ID),]

  res1 <- lapply(2:4, function(k){
     km <- kmeans(t(h), centers = k)
     message(max(km$cluster))
      cluster <- data.frame(ROI = colnames(sub.data), cluster = factor(km$cluster))
     #cluster$roi <- unlist(sapply(cluster$ROI, function(x) unlist(strsplit(x, "\\-"))[2]))
     rownames(cluster) <- unlist(sapply(cluster$ROI, function(x) unlist(strsplit(x, "\\-"))[2]))
     res <- spatial.coord;
     res$cluster <- factor(km$cluster)
     res$k = k
     #p <- ggplot(res, aes(x = ROICoordinateX/1e3, y = ROICoordinateY/1e3, color = cluster, shape = location)) +
     #   geom_point(size = 1)+ theme_bw() + 
      # theme(panel.grid = element_blank(), aspect.ratio = 1, legend.position = "right")+
      # scale_color_brewer(palette = "Dark2") + 
      # labs(x = "", y = "", title = paste0(id, " ", "K = ", k), color = "", shape = "")
     #return(p)
     res
    })
  #return(p.list)
  res <- do.call(rbind, res1)
  res$sample <- id
  res
  
})

res <- do.call(rbind, all.p.list)
res$k = paste0("k = ", res$k)
p <- ggplot(res, aes(x = ROICoordinateY, y = -ROICoordinateX, color = cluster, shape = location)) +
  geom_point(alpha = 0.6)+ theme_bw() + 
  facet_wrap(vars(sample, k), scales = "free", ncol = 6)+
  theme(panel.grid = element_blank(), aspect.ratio = 1, legend.position = "right",
        axis.text = element_blank(), axis.ticks = element_blank())+
  scale_color_brewer(palette = "Set1") + 
  scale_shape_manual(values=c("TC" = 15, "IE" = 16, "GM" = 4, "BV" = 2))+
  labs(x = "", y = "", color = "", shape = "")

pdf("04.ROI.cluster.pdf", width = 15, height = 15, useDingbats = F)
#cowplot::plot_grid(plotlist = all.p, ncol = 4, label_size = 7)
p
dev.off()

