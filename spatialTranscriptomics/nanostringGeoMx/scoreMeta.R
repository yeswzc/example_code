rm(list = ls())
setwd("/Users/wuz6/Documents/Project/08.spatialTranscriptome/01.RNA/01.edgeR.QCfirst/")

#library(dbscan)
#library(NMF)
library(ggplot2)
library(Seurat)
source("../../src/read_nanostring.R")
red_blue_20 = c("#124984","#1f63a8","#2f79b5","#3f8ec0","#5fa5cd","#87beda","#a7d0e4",
                "#c5dfec","#dbeaf2","#edf2f5","#f9f0eb","#fbe3d4","#fbd0b9","#f7b799",
                "#f09c7b","#e17860","#d25849","#c13639","#ae172a","#8a0b25");


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

programs <- read.csv("08.meta.programs.csv", head =T)

all.scores.umap <- lapply(sample.names, function(id){
  set.seed(123)
  #id <- sample.names[1]
  cta.file <- grep(id, data.files, value = T)
  spatial.xy <- read_nanostring_RNA_XY(cta.file)
  spatial.xy[,1] <- paste0(id, "-", spatial.xy[,1])
  
  message(id, ": ", cta.file)
  sub.data <- data[,grep(id, colnames(data))]
  message("Sub data nrows: ", dim(sub.data)[1], ", ncols: ", dim(sub.data)[2])
  seurat <- CreateSeuratObject(counts = sub.data)
  meta.program.scores <- lapply(1:ncol(programs), function(k){
    genes <- programs[,k]
    genes <- genes[!is.na(genes) & genes != ""]
    scores <- AddModuleScore(seurat, features = genes, name = "sig", ctrl = 50)
    score <- rowMeans(scores[[]][,grep("sig", colnames(scores[[]]))])
    score;
    
  })
  meta.program.scores <- do.call(cbind, meta.program.scores)
  meta.program.scores <- apply(meta.program.scores, 2, as.numeric)
  colnames(meta.program.scores) <- paste0("MP", 1:ncol(programs))
  pc <- prcomp(t(sub.data))
  set.seed(123)
  um <- umap::umap(pc$x[,1:10])$layout
  colnames(um) <- c("umap1", "umap2")
  res <- data.frame(cbind(um, meta.program.scores))
  res$sample = id
  res;
  res <- cbind(res, spatial.xy[match(rownames(res), spatial.xy[,1]),])
  res;
  
 })

all.scores.umap <- do.call(rbind, all.scores.umap)
all.scores.umap$ann <- meta$annotation1[match(rownames(all.scores.umap), meta$ID)]
head(all.scores.umap)

p.list <- lapply(unique(all.scores.umap$sample), function(id){
  p.data <- all.scores.umap[all.scores.umap$sample == id,]
  ps <- lapply(paste0("MP", 1:4), function(mp){
    mid.point <- median(p.data[,mp])
    p <- ggplot(p.data, aes_string(x= "ROICoordinateY", y = "-ROICoordinateX", col = mp, shape = "ann")) +
      #facet_grid(~variable, scales = "free")+
      geom_point(size =0.5)+
      scale_color_gradient2(low = red_blue_20[1], high = red_blue_20[20], mid = red_blue_20[10], midpoint = mid.point)+
      theme_bw()+ theme(aspect.ratio = 1, legend.position = "none", 
                        axis.text = element_blank(), axis.ticks = element_blank(),
                        axis.title = element_text(size =2),
                        panel.grid = element_blank())+
      scale_shape_manual(values=c("TC" = 15, "IE" = 16, "GM" = 4, "BV" = 2))+
      labs(x = "", y = "", color = "", shape = "", title = paste0(id, "-", mp))
    p
  })
  #ps <- do.call(c, ps)
  ps
})
p.list <- do.call(c, p.list)
p.list[[2]]

#cowplot::plot_grid(plotlist = p.list, ncol = 4)
dev.off()

pdf("08.sample.spatial.metaScore.pdf", width = 10, height = 10, useDingbats = F)
cowplot::plot_grid(plotlist = p.list, ncol = 8)
dev.off()


p.list <- lapply(unique(all.scores.umap$sample), function(id){
  p.data <- all.scores.umap[all.scores.umap$sample == id,]
  ps <- lapply(paste0("MP", 1:4), function(mp){
    mid.point <- median(p.data[,mp])
    p <- ggplot(p.data, aes_string(x= "umap1", y = "umap2", col = mp, shape = "ann")) +
      #facet_grid(~variable, scales = "free")+
      geom_point(size =0.5)+
      scale_color_gradient2(low = red_blue_20[1], high = red_blue_20[20], mid = red_blue_20[10], midpoint = mid.point)+
      theme_bw()+ theme(aspect.ratio = 1, legend.position = "none", 
                        axis.text = element_blank(), axis.ticks = element_blank(),
                        axis.title = element_text(size =2),
                        panel.grid = element_blank())+
      scale_shape_manual(values=c("TC" = 15, "IE" = 16, "GM" = 4, "BV" = 2))+
      labs(x = "", y = "", shape = "", title = paste0(id, "-", mp))
    p
  })
  #ps <- do.call(c, ps)
  ps
})
p.list <- do.call(c, p.list)
p.list[[1]]

pdf("08.sample.umap.metaScore.pdf", width = 10, height = 10, useDingbats = F)
cowplot::plot_grid(plotlist = p.list, ncol = 8)
dev.off()
