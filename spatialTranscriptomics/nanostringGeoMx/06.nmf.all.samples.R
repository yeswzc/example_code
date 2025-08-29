###NMF deconvolve/cluster
rm(list=ls())
library(NMF)
library(ggplot2)
library(pheatmap)
library(umap)
library(clusterProfiler);
library(dplyr)

db = msigdbr::msigdbr(species = "human", category = "C5") %>% dplyr::distinct(gene_symbol, gs_name) %>% as.data.frame()


red_blue_20 = c("#124984","#1f63a8","#2f79b5","#3f8ec0","#5fa5cd","#87beda","#a7d0e4",
                "#c5dfec","#dbeaf2","#edf2f5","#f9f0eb","#fbe3d4","#fbd0b9","#f7b799",
                "#f09c7b","#e17860","#d25849","#c13639","#ae172a","#8a0b25");

sample.names = c(1399, 1905, 2023, 2751, 3984, 4461, 604, 6614, "8147M", "8147O")
sample.colors = RColorBrewer::brewer.pal(12, "Paired")[c(1:10)]
names(sample.colors) = sample.names

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

#xx <- nmf(data, 10)
#row as genes, columns as samples/ROIs
#w <- basis(xx)
#h <- coef(xx)


### selecting r/k
estim.r <- nmf(data, 3:12, nrun=5, seed=12345)
pdf("06.NMF.rK.test.pdf", width = 10, height = 10, useDingbats = F)
plot(estim.r)
dev.off()
#overfitting?
#consensusmap(estim.r, annCol=data, labCol=NA, labRow=NA)
### Selecting method
res.multi.method <- nmf(data, 4, list('brunet', 'lee', 'ns'), seed=12345, .options='t')
compare(res.multi.method)

#method   seed rng    metric rank sparseness.basis sparseness.coef silhouette.coef silhouette.basis residuals
#brunet brunet random   1        KL    4       0.04352269      0.03784516       0.7082316        0.4814079  2361.569
#lee       lee random   1 euclidean    4       0.03797525      0.04143375       0.8313189        0.5288883 23726.876
#nsNMF   nsNMF random   1        KL    4       0.05061990      0.07996412       0.6479233        0.5299533  2362.913
#niter    cpu cpu.all nrun
#brunet  2000 44.277  44.277    1
#lee     1450 19.237  19.237    1
#nsNMF   2000 47.198  47.198    1

rm(estim.r, res.multi.method)


###Run with r =4, method = lee
res <- nmf(data, 4, method = "lee", seed = 12345, nrun = 5)
#row as genes, columns as samples/ROIs
w <- basis(res)
h <- coef(res)
w[1:4,1:4]
h[1:4,1:4]

um <- umap::umap(t(h))$layout
um <- data.frame(um)
colnames(um) <- c("umap1", "umap2")
um$sample <- sapply(rownames(um), function(x) unlist(strsplit(x, "\\-"))[1])
head(um)

pdf("06.NMF.all.sample.umap.pdf", width = 5, height = 4, useDingbats = F)
ggplot(um, aes(x = umap1, y = umap2, col = sample)) + 
  geom_point()+
  scale_color_manual(values = sample.colors)+
  theme_bw() + theme(aspect.ratio = 1)
dev.off()


f <- extractFeatures(res, 30L)
f <- lapply(f, function(x) rownames(res)[x])
f <- do.call("rbind", f)
head(f)
