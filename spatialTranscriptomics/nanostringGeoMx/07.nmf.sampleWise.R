###NMF deconvolve/cluster
rm(list=ls())
library(NMF)
library(ggplot2)
library(pheatmap)
library(umap)
library(clusterProfiler);
library(dplyr)
source("../../src/read_nanostring.R")
library(fgsea)

library(Seurat)

c5 = msigdbr::msigdbr(species = "human", category = "C5")
db = c5 %>% dplyr::distinct(gene_symbol, gs_name) %>% as.data.frame()
db.H = msigdbr::msigdbr(species = "human", category = "H") %>% dplyr::distinct(gene_symbol, gs_name) %>% as.data.frame()
db.KEGG = msigdbr::msigdbr(species = "human", category = "C2", subcategory = "KEGG") %>% dplyr::distinct(gene_symbol, gs_name) %>% as.data.frame()
db.REACTOME = msigdbr::msigdbr(species = "human", category = "C2", subcategory = "REACTOME") %>% dplyr::distinct(gene_symbol, gs_name) %>% as.data.frame()

c5.fgsea.db = split(x = c5$gene_symbol, f = c5$gs_name)

red_blue_20 = c("#124984","#1f63a8","#2f79b5","#3f8ec0","#5fa5cd","#87beda","#a7d0e4",
                "#c5dfec","#dbeaf2","#edf2f5","#f9f0eb","#fbe3d4","#fbd0b9","#f7b799",
                "#f09c7b","#e17860","#d25849","#c13639","#ae172a","#8a0b25");

data.files = list.files("../../00.data/CTA/raw//", pattern = "C*xlsx$", full.names = T)

sample.names = c(1399, 1905, 2023, 2751, 3984, 4461, 604, 6614, "8147M", "8147O")
sample.colors = RColorBrewer::brewer.pal(12, "Paired")[c(1:10)]
names(sample.colors) = sample.names


seurat.score <- function(x, gene.matrix){
  #message("...")
  seurat <- CreateSeuratObject(counts = x)
  #print(seurat)
  scores <- lapply(1:ncol(gene.matrix), function(k){
    genes <- gene.matrix[,k]
    s <- AddModuleScore(seurat, features = genes, name = "sig", ctrl = 50)
    score <- rowMeans(s[[]][,grep("sig", colnames(s[[]]))])
    score;
  })
  scores <- do.call(cbind, scores)
  scores
}

###
setwd("/Users/wuz6/Documents/Project/08.spatialTranscriptome/01.RNA/01.edgeR.QCfirst/")

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

data <- read.csv("00.IDHastrocytoma.CTA.cpm.filtered.csv", head =T, row.names = 1, check.names = F)
data[1:4,1:4]
###
#Keep only TC
data.ann <- meta$annotation1[match(colnames(data), meta$ID)]
data <- data[,data.ann == "TC"]
dim(data)
data.samples <- sapply(colnames(data), function(x) unlist(strsplit(x, "\\-"))[1])

#xx <- nmf(data, 10)
#row as genes, columns as samples/ROIs
#w <- basis(xx)
#h <- coef(xx)

k = 1
### selecting r/k
id <- sample.names[k]
id




xx <- lapply(sample.names, function(id){

test.r = 2:5
estim.r <- nmf(data[, data.samples == id], test.r, nrun=5, seed=12345, method = "lee")
#estim.r$measures
pdf(paste0("07.", id, ".NMF.rK.test.pdf"), width = 10, height = 10, useDingbats = F)
print(plot(estim.r))
dev.off()

### Selecting r based on when cophenetic starts to drop fast
coph <- estim.r$measures$cophenetic
coph.diff <- sapply(2:length(coph), function(k){
   diff <- coph[k-1] - coph[k]
})
R = which.max(coph.diff)+1

###
message(id, ":r = ", R)
#overfitting?
#consensusmap(estim.r, annCol=data, labCol=NA, labRow=NA)
### Selecting method, skip
###Run with r =4, method = lee
res <- nmf(data[, data.samples == id], R, method = "lee", seed = 12345, nrun = 5)
#row as genes, columns as samples/ROIs
w <- basis(res)
#h <- coef(res)
#w[1:2,1:2]
#h[1:2,1:2]

um <- data.frame(umap::umap(prcomp(t(data[, data.samples == id]))$x[,1:10])$layout)
colnames(um) <- c("umap1", "umap2")
#um <- umap::umap(t(h))$layout
#um <- data.frame(um)
#colnames(um) <- c("umap1", "umap2")
um$sample <- sapply(rownames(um), function(x) unlist(strsplit(x, "\\-"))[1])
#um <- cbind(um1, um)
#head(um)

f <- extractFeatures(res, 50L)
f <- lapply(f, function(x) rownames(res)[x])
f <- do.call(cbind, f)
colnames(f) <- paste0("NMF", 1:ncol(f))
write.csv(f, file = paste0("07.", id, ".features.csv"), row.names = F, quote = F,)

nmf.scores <- seurat.score(data[, data.samples == id], f)
colnames(nmf.scores) <- paste0("NMF", 1:R)
identical(rownames(um), rownames(nmf.scores))
cta.file <- grep(id, data.files, value = T)
spatial.coord <- read_nanostring_RNA_XY(cta.file)

um <- cbind(um, nmf.scores)
um$ROI <- sapply(rownames(um), function(x) unlist(strsplit(x, "\\-"))[2])

um <- cbind(um, spatial.coord[match(um$ROI, spatial.coord$ROILabel), ])

p.list <- lapply(1:R, function(k){
  mid <- median(um[, paste0("NMF", k)])
  p <- ggplot(um, aes_string(x = 'umap1', y = 'umap2', col = paste0("NMF", k))) + 
    geom_point(size =1)+
    scale_color_gradient2(low = red_blue_20[1], mid = red_blue_20[10], high = red_blue_20[20], midpoint = mid)+
    theme_bw() + theme(aspect.ratio = 1, legend.key.size = unit(0.5, 'cm'))
  p
})



pdf( paste0("07.", id, ".NMFscore.umap.pdf"), width = 10, height = 10, useDingbats = F)
print(cowplot::plot_grid(plotlist = p.list))
dev.off()

p.list <- lapply(1:R, function(k){
  mid <- median(um[, paste0("NMF", k)])
  p <- ggplot(um, aes_string(x = 'ROICoordinateX', y = 'ROICoordinateY', col = paste0("NMF", k))) + 
    geom_point(size =1)+
    scale_color_gradient2(low = red_blue_20[1], mid = red_blue_20[10], high = red_blue_20[20], midpoint = mid)+
    theme_bw() + 
    theme(aspect.ratio = 1, axis.text = element_blank(), 
          axis.ticks = element_blank(),
          legend.key.size = unit(0.5, 'cm'))+
    labs(x = "", y = "")
  p
})

pdf( paste0("07.", id, ".NMFscore.xy.pdf"), width = 5, height = 4, useDingbats = F)
print(cowplot::plot_grid(plotlist = p.list))
dev.off()




})

meta.files <- list.files("./", pattern = "07.*.features.csv")
.. <- lapply(meta.files, function(file){
  f <- read.csv(file, head =T, stringsAsFactors = F)
  id <- unlist(strsplit(file,"\\."))[2]
  message(id)
  .. <- lapply(1:ncol(f), function(k){
    #go <- enricher(gene = f[,k], TERM2GENE = db, minGSSize = 5)@result #universe = row.names(data)),
    #go$ID = gsub("GO\\w\\w_", "", go$ID)
    #write.csv(go, file = paste0("07.", id, ".NMF.K", k, ".GO.csv"), row.names = F)
    #h <- enricher(gene = f[,k], TERM2GENE = db.H, minGSSize = 5)@result #universe = row.names(data)),
    #write.csv(h, file = paste0("07.", id, ".NMF.K", k, ".H.csv"), row.names = F)
    kegg <- enricher(gene = f[,k], TERM2GENE = db.KEGG, minGSSize = 5)@result
    write.csv(kegg, file = paste0("07.", id, ".NMF.K", k, ".KEGG.csv"), row.names = F)
    react <- enricher(gene = f[,k], TERM2GENE = db.REACTOME, minGSSize = 5)@result
    write.csv(react, file = paste0("07.", id, ".NMF.K", k, ".REACTOME.csv"), row.names = F)
    
    #gene.rank = w[,k]
    #names(gene.rank) <- rownames(w)
    #barplot(sort(gene.rank, decreasing = T))
    #fgsea.res.C5.A <- fgsea(pathways = c5.fgsea.db, stats = gene.rank, nPermSimple = 10000, minSize = 10)#Runs fgseaMultilevel
    #data.table::fwrite(fgsea.res.C5.A, file = paste0("07.", id, ".NMF.K", k, ".GOfgsea.csv"), sep=",")
    })
})
rm(..)


################################################################################
################################################################################
################################################################################
##Run R = 4 for each samples
xx <- lapply(sample.names, function(id){
  id = sample.names[1]
  res <- nmf(data[, data.samples == id], 4, method = "lee", seed = 12345, nrun = 5)
  
  f <- extractFeatures(res, 50L)
  f <- lapply(f, function(x) rownames(res)[x])
  f <- do.call(cbind, f)
  colnames(f) <- paste0("NMF", 1:ncol(f))
  write.csv(f, file = paste0("08.", id, ".k4.features.csv"), row.names = F, quote = F,)
})

###Identify shared meta program
meta.files <- list.files("./", pattern = "08.*.k4.features.csv")
meta.files

meta.programs <- lapply(meta.files, function(file.name){
   id <- unlist(strsplit(file.name, "\\."))[2]
   f <- read.csv(file.name, head =T, stringsAsFactors = F)#, nrow = 20)
   colnames(f) <- paste0(id, "-", colnames(f))
   f
})
#head(meta.programs[[1]])
meta.programs <- do.call(cbind, meta.programs)
dim(meta.programs)

meta.jac <- lapply(1:(ncol(meta.programs)), function(i){
   x1 <- colnames(meta.programs)[i]
   jac <- lapply(1:ncol(meta.programs), function(j){
      x2 <- colnames(meta.programs[j])
      l1 <- length(intersect(meta.programs[,i], meta.programs[,j]))
      l2 <- length(unique(c(meta.programs[,i], meta.programs[,j])))
      j <- l1/l2
      c(x1, x2, j)
   })
   jac <- do.call(rbind, jac)
   jac
})
meta.jac <- data.frame(do.call(rbind, meta.jac))
head(meta.jac)

meta.jac[,3] <- as.numeric(meta.jac[,3])
meta.jac.m <- reshape2::dcast(meta.jac, X1~X2)
rownames(meta.jac.m) <- meta.jac.m[,1]
meta.jac.m <- meta.jac.m[,-1]
meta.jac.m[1:4,1:4]

c <- hclust(as.dist(1-meta.jac.m), method = "average")

dev.off()
pdf("08.meta.programs.Jaccard.pdf", width = 10, height = 10, useDingbats = F)
pheatmap(meta.jac.m, 
         cluster_rows = c, cluster_cols = c, 
         color = red_blue_20, 
         cutree_cols =16, cutree_rows = 16,
         cellwidth = 10, cellheight = 10,
         main = "Jaccard index")

pheatmap(meta.jac.m, 
         cluster_rows = c, cluster_cols = c, 
         color = red_blue_20, 
         #cutree_cols =16, cutree_rows = 16,
         cellwidth = 10, cellheight = 10,
         main = "Jaccard index")

dev.off()

m <- meta.jac.m[c$order, c$order]
write.csv(m, file = "08.meta.jaccard.csv", quote = F)

library(dendextend)
x <- sort(cutree(c, k = 16))
keep = which(table(x) >= 3)
keep

meta.programs.overlap <- lapply(1:length(keep), function(k){
   cut <- keep[k]
   keep.programs <- names(x[x == cut])
   genes <- as.vector(as.matrix(meta.programs[,keep.programs]))
   tb <- table(genes)
   names(tb)[tb >=2]
   #sort(unique(genes))
})

length(meta.programs.overlap)

meta.genes <- unique(unlist(meta.programs.overlap))
dim(data)
length(data.samples)

data.ann <- data.frame(sample = data.samples)
rownames(data.ann) <- names(data.samples)

gene.cluster <- hclust(as.dist(1 - cor(t(data[meta.genes,]), method = "pearson")), method = "average")
sample.cluster <- hclust(as.dist(1 - cor(data[meta.genes,], method = "pearson")), method = "average")

duplicated <- table(unlist(meta.programs.overlap))
duplicated <- names(which(duplicated >=2))
meta.programs.overlap <- lapply(1:length(meta.programs.overlap), function(k){
  G <- meta.programs.overlap[[k]]
  G <- G[!G %in% duplicated]
  G
})

gene.ann <- lapply(1:length(meta.programs.overlap), function(k){
  G <- meta.programs.overlap[[k]]
  #G <- G[!G %in% duplicated]
  #message(G)
  go <- enricher(gene = G, TERM2GENE = db, minGSSize = 5)@result #universe = row.names(data)),
  go$ID = gsub("GO\\w\\w_", "", go$ID)
  write.csv(go, file = paste0("08.meta.K", k, ".GO.csv"), row.names = F)
  h <- enricher(gene = G, TERM2GENE = db.H, minGSSize = 5)@result #universe = row.names(data)),
  write.csv(h, file = paste0("08.meta.K", k, ".H.csv"), row.names = F)
  kegg <- enricher(gene = G, TERM2GENE = db.KEGG, minGSSize = 5)@result
  write.csv(kegg, file = paste0("08.meta.K", k, ".KEGG.csv"), row.names = F)
  react <- enricher(gene = G, TERM2GENE = db.REACTOME, minGSSize = 5)@result
  write.csv(react, file = paste0("08.meta.K", k, ".REACTOME.csv"), row.names = F)
  
   ifelse(meta.genes %in% G, "Y", "N")
})


x <- do.call(cbind, meta.programs.overlap)
write.csv(x, file = "08.meta.programs.csv", quote = F, row.names = F)

gene.ann <- data.frame(do.call(cbind, gene.ann))
row.names(gene.ann) <- meta.genes
colnames(gene.ann) <- paste0("meta", 1:ncol(gene.ann))

in.col <- c("Y" = "green", "N" = "black")

gene.order <- lapply(meta.programs.overlap, function(G){
  h <- hclust(as.dist(1 - cor(t(data[G,]), method = "pearson")), method = "average")
  h$label[h$order]
})
gene.order <- unlist(gene.order)

gene.cluster <- hclust(as.dist(1 - cor(t(data[meta.genes,]), method = "pearson")), method = "average")
sample.cluster <- hclust(as.dist(1 - cor(data[meta.genes,], method = "pearson")), method = "average")

#gene.cluster <- hclust(dist(data[meta.genes,], method = "euclidean"), method = "ward.D2")
sample.cluster <- hclust(dist(t(data[meta.genes,]), method = "euclidean"), method = "ward.D2")
dev.off();dev.off();dev.off();dev.off();
dev.off()
pdf("08.meta.programs.profiles.pdf", width = 10, height = 10, useDingbats = F)
pheatmap(data[gene.order,], 
         annotation_row = gene.ann,
         scale = "row",
         col = red_blue_20, show_colnames = F, 
         cluster_rows = F,
         #cluster_rows = gene.cluster,
         cluster_cols = sample.cluster,
         annotation_col = data.ann, 
         annotation_colors = list(sample = sample.colors, 
                                  meta1 = in.col,meta2 = in.col,
                                  meta3 = in.col,meta4 = in.col),
         cutree_rows = 4,
         cutree_cols = 5,
         fontsize = 6)

dev.off()

#### The heatmap of meta programs displayed sample specific patterm, meta 2-4 are highly similar to each other
###Meta 1 is high in 1399, 8147, 4461
###Meta 2 is high in 1905, 3984, 2751, and 2023, some of 604
###Meta 3 is high in some of 1905, 3984, 2751, 2023, high in some of 1399, 8147, 4461
###Meta 4 is high in some of 1905, 3984, 2751, some of 1399, 8147, 4461, 604, a few of 2023


