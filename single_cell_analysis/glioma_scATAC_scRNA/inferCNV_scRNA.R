setwd("/Users/wuz6/Project/01.single_cell_ATAC/IDH.glioma/scRNA/CNV/")
rm(list=ls())
source("/Users/wuz6/Project/01.single_cell_ATAC/IDHwt.GBM/scRNA/CNA/inferCNV.src.R") #
load("/Users/wuz6/Project/01.single_cell_ATAC/IDHwt.GBM/scRNA/CNA/hg38.cellranger3.0.0.gene_annotation.rda")




### CEH
load("CEH.Matrix.rda")
max(suva_expr)
min(suva_expr)

#suva_expr = center_expression(suva_expr)
max(suva_expr)
min(suva_expr)
plot(density(rowMeans(suva_expr)))

chr.bin.genes = make_cnv_bins(gene_annotation,keep.genes = rownames(suva_expr), chromosomes = 1:22)


suva_expr.cnv.raw = cnv.bin(t(suva_expr), chr.bin.genes )

head(suva_expr.cnv.raw[,1:4])
max(suva_expr.cnv.raw); min(suva_expr.cnv.raw)

#center the data, sequencing depth
suva_expr.cnv.scaled = scale(suva_expr.cnv.raw, center = T, scale = T)
max(suva_expr.cnv.scaled);min(suva_expr.cnv.scaled)

plot(density(as.vector(suva_expr.cnv.scaled)))
suva_expr.cnv.scaled[suva_expr.cnv.scaled>3] = 3
suva_expr.cnv.scaled[suva_expr.cnv.scaled < -3] = -3

###CEH
cnv.cluster = make.initial.cnvHeatmap(suva_expr.cnv.scaled, meta, output = "CEH.1-1pearson.cnv.png", k.cut = 3, 
                                      cluster.method = "average", cluster.distance = "pearson")
dim(suva_expr.cnv.scaled)

findsplitK(cnv.cluster, n.cut = 3, k.init = 7, N = 50)
k = 15
cnv.cluster.cut = cutree(cnv.cluster, k)
sort(table(cnv.cluster.cut), decreasing = T)[1:4]
cnv.cluster.cut[cnv.cluster.cut>2] = 3
cnv.cluster = make.initial.cnvHeatmap(suva_expr.cnv.scaled, meta, output = "CEH.1-2pearson.cnv.png", k.cut = 10, 
                                      cluster.method = "average", cluster.distance = "pearson", 
                                      cluster.cut = cnv.cluster.cut)

cnv.cluster.cut[cnv.cluster.cut>2] = 2
table(cnv.cluster.cut)
cnv.final = subtractControl(suva_expr.cnv.scaled, cnv.cluster.cut, k.control = 2, choice = "mean")
head(cnv.final[,1:4])
cnv.cluster2 = make.initial.cnvHeatmap(cnv.final, meta, output="CEH.1-3.cnvCtl.png", k.cut = 3, 
                                       cluster.method = "average", cluster.distance = "pearson",
                                       cluster.cut = cnv.cluster.cut)
#cannt't tell for now
findsplitK(cnv.cluster2, k.init = 20)
k2 = get_subdendrograms(as.dendrogram(cnv.cluster2), k = 2)
length(k2[[1]])
x = k2[[1]]
length(x)
findsplitK(k2[[1]], N = 20)
sub.cut = cutree(k2[[1]], k = 5)
table(sub.cut)
sub.cut[sub.cut>2] = 2

sub.cut2 = cutree(k2[[2]], k = 2)
sub.cut2[sub.cut2 > 0] = 2

cnv.cluster.cut2 = c(sub.cut, sub.cut2)

cnv.cluster.cut2 = cnv.cluster.cut2[colnames(cnv.final)]
cnv.cluster2 = make.initial.cnvHeatmap(cnv.final, meta, output="CEH.1-4.cnvCtl.png", k.cut = 6, 
                                       cluster.method = "average", cluster.distance = "pearson",
                                       cluster.cut = cnv.cluster.cut2)
table(cnv.cluster.cut2)
table(cnv.cluster.cut)

candi.malignant = names(which(cnv.cluster.cut2 ==1))
candi.malignant[1:10]

#write.table(candi.malignant, file = "CEH.candi.malignant.txt", row.names = F, col.names = F, quote = F)
####







############################################################################################################
rm(chr.bin.genes, cnv.cluster, cnv.cluster2, cnv.final, meta, suva_expr, suva_expr.cnv.raw, suva_expr.cnv.scaled, 
   cnv.cluster.cut, cnv.cluster.cut2, k)
load("FNF.Matrix.rda")
max(suva_expr)
min(suva_expr)

plot(density(rowMeans(suva_expr)))
chr.bin.genes = make_cnv_bins(gene_annotation,keep.genes = rownames(suva_expr), chromosomes = 1:22)
suva_expr.cnv.raw = cnv.bin(t(suva_expr), chr.bin.genes )

head(suva_expr.cnv.raw[,1:4])
max(suva_expr.cnv.raw); min(suva_expr.cnv.raw)

#center the data, sequencing depth
suva_expr.cnv.scaled = scale(suva_expr.cnv.raw, center = T, scale = T)
max(suva_expr.cnv.scaled);min(suva_expr.cnv.scaled)

plot(density(as.vector(suva_expr.cnv.scaled)))
suva_expr.cnv.scaled[suva_expr.cnv.scaled>3] = 3
suva_expr.cnv.scaled[suva_expr.cnv.scaled < -3] = -3

###FNF
cnv.cluster = make.initial.cnvHeatmap(suva_expr.cnv.scaled, meta, output = "FNF.1-1pearson.cnv.png", k.cut = 3, 
                                      cluster.method = "average", cluster.distance = "pearson")
dim(suva_expr.cnv.scaled)

findsplitK(cnv.cluster, N = 20, n.cut = 3, k.init = 3)
k = 28
cnv.cluster.cut = cutree(cnv.cluster, k)
sort(table(cnv.cluster.cut), decreasing = T)[1:4]
cnv.cluster.cut[cnv.cluster.cut>10] = 10
cnv.cluster.cut[cnv.cluster.cut!=2] = 10
cnv.cluster = make.initial.cnvHeatmap(suva_expr.cnv.scaled, meta, output = "FNF.1-2pearson.cnv.png", k.cut = 10, 
                                      cluster.method = "average", cluster.distance = "pearson", 
                                      cluster.cut = cnv.cluster.cut)
table(cnv.cluster.cut)

#k2 is candidately cancer
cnv.cluster.cut[cnv.cluster.cut > 2] = 1
table(cnv.cluster.cut)
cnv.final = subtractControl(suva_expr.cnv.scaled, cnv.cluster.cut, k.control = 1, choice = "mean")
head(cnv.final[,1:4])
cnv.cluster2 = make.initial.cnvHeatmap(cnv.final, meta, output="FNF.1-3.cnvCtl.png", k.cut = 3, 
                                       cluster.method = "average", cluster.distance = "pearson",
                                       cluster.cut = cnv.cluster.cut)


library(dendextend)
k2 = get_subdendrograms(as.dendrogram(cnv.cluster2), k = 2)
nleaves(k2[[1]])
nleaves(k2[[2]])
findsplitK(k2[[1]], N = 80, k.init = 5, n.cut = 3)
sub.cut = cutree(k2[[1]], k = 42)
sort(table(sub.cut), decreasing = T)[1:4]
sub.cut[sub.cut>3] = 4

sub.cut2 = cutree(k2[[2]], k = 2)
sub.cut2[sub.cut2 > 0] = 6

cnv.cluster.cut2 = c(sub.cut, sub.cut2)
cnv.cluster.cut2 = cnv.cluster.cut2[colnames(cnv.final)];
table(cnv.cluster.cut2)
cnv.cluster.cut2[cnv.cluster.cut2 != 2] = 1
cnv.cluster2 = make.initial.cnvHeatmap(cnv.final, meta, output="FNF.1-4.cnvCtl.png", k.cut = 6, 
                                       cluster.method = "average", cluster.distance = "pearson",
                                       cluster.cut = cnv.cluster.cut2)

table(cnv.cluster.cut2)
table(cnv.cluster.cut)

candi.malignant = names(which(cnv.cluster.cut2 ==2))
candi.malignant[1:10]

write.table(candi.malignant, file = "FNF.candi.malignant.txt", row.names = F, col.names = F, quote = F)
####





############################################################################################################
rm(chr.bin.genes, cnv.cluster, cnv.cluster2, cnv.final, meta, suva_expr, suva_expr.cnv.raw, suva_expr.cnv.scaled, 
   cnv.cluster.cut, cnv.cluster.cut2, k)
load("GPU.Matrix.rda")
max(suva_expr)
min(suva_expr)

plot(density(rowMeans(suva_expr)))
chr.bin.genes = make_cnv_bins(gene_annotation,keep.genes = rownames(suva_expr), chromosomes = 1:22)
suva_expr.cnv.raw = cnv.bin(t(suva_expr), chr.bin.genes )

head(suva_expr.cnv.raw[,1:4])
max(suva_expr.cnv.raw); min(suva_expr.cnv.raw)

#center the data, sequencing depth
suva_expr.cnv.scaled = scale(suva_expr.cnv.raw, center = T, scale = T)
max(suva_expr.cnv.scaled);min(suva_expr.cnv.scaled)

plot(density(as.vector(suva_expr.cnv.scaled)))
suva_expr.cnv.scaled[suva_expr.cnv.scaled>3] = 3
suva_expr.cnv.scaled[suva_expr.cnv.scaled < -3] = -3

###GPU
cnv.cluster = make.initial.cnvHeatmap(suva_expr.cnv.scaled, meta, output = "GPU.1-1pearson.cnv.png", k.cut = 3, 
                                      cluster.method = "average", cluster.distance = "pearson")
dim(suva_expr.cnv.scaled)

#findsplitK(cnv.cluster,)
cnv.cluster.cut = cutree(cnv.cluster, 3)
sort(table(cnv.cluster.cut), decreasing = T)[1:4]
#cnv.cluster.cut[cnv.cluster.cut>5] = 
cnv.cluster = make.initial.cnvHeatmap(suva_expr.cnv.scaled, meta, output = "GPU.1-2pearson.cnv.png", k.cut = 10, 
                                      cluster.method = "average", cluster.distance = "pearson", 
                                      cluster.cut = cnv.cluster.cut)
#k2 is candidately cancer
cnv.cluster.cut[cnv.cluster.cut>2] = 2
table(cnv.cluster.cut)
cnv.final = subtractControl(suva_expr.cnv.scaled, cnv.cluster.cut, k.control = 2, choice = "mean")
head(cnv.final[,1:4])
cnv.cluster2 = make.initial.cnvHeatmap(cnv.final, meta, output="GPU.1-3.cnvCtl.png", k.cut = 3, 
                                       cluster.method = "average", cluster.distance = "pearson",
                                       cluster.cut = cnv.cluster.cut)

cnv.cluster2 = make.initial.cnvHeatmap(cnv.final, meta, output="GPU.1-3.cnvCtl.png", k.cut = 3, 
                                       cluster.method = "average", cluster.distance = "pearson")
#cannt't tell for now
findsplitKrev(cnv.cluster2, N = 640, n.cut = 1, k.max = 297)

cnv.cluster.cut2 = cutree(cnv.cluster2, 232)
sort(table(cnv.cluster.cut2), decreasing = T)[1:4]
cnv.cluster.cut2[cnv.cluster.cut2>2] =  2
cnv.cluster = make.initial.cnvHeatmap(cnv.final, meta, output = "GPU.1-4pearson.cnv.png", k.cut = 10, 
                                      cluster.method = "average", cluster.distance = "pearson", 
                                      cluster.cut = cnv.cluster.cut2)


table(cnv.cluster.cut)

#some cells looks to be non-malignant cells were 
candi.malignant = names(which(cnv.cluster.cut2 ==1))
candi.malignant[1:10]
write.table(candi.malignant, file = "GPU.candi.malignant.txt", row.names = F, col.names = F, quote = F)






############################################################################################################
rm(chr.bin.genes, cnv.cluster, cnv.cluster2, cnv.final, meta, suva_expr, suva_expr.cnv.raw, suva_expr.cnv.scaled, 
   cnv.cluster.cut, cnv.cluster.cut2, k)
rm(candi.malignant)
load("JAS.Matrix.rda")
max(suva_expr)
min(suva_expr)

plot(density(rowMeans(suva_expr)))
chr.bin.genes = make_cnv_bins(gene_annotation,keep.genes = rownames(suva_expr), chromosomes = 1:22)
suva_expr.cnv.raw = cnv.bin(t(suva_expr), chr.bin.genes )

head(suva_expr.cnv.raw[,1:4])
max(suva_expr.cnv.raw); min(suva_expr.cnv.raw)

#center the data, sequencing depth
suva_expr.cnv.scaled = scale(suva_expr.cnv.raw, center = T, scale = T)
max(suva_expr.cnv.scaled);min(suva_expr.cnv.scaled)

plot(density(as.vector(suva_expr.cnv.scaled)))
suva_expr.cnv.scaled[suva_expr.cnv.scaled>3] = 3
suva_expr.cnv.scaled[suva_expr.cnv.scaled < -3] = -3

###JAS
cnv.cluster = make.initial.cnvHeatmap(suva_expr.cnv.scaled, meta, output = "JAS.1-1pearson.cnv.png", k.cut = 3, 
                                      cluster.method = "average", cluster.distance = "pearson")
dim(suva_expr.cnv.scaled)

findsplitK(cnv.cluster, N = 50)
k = 19
cnv.cluster.cut = cutree(cnv.cluster, k)
sort(table(cnv.cluster.cut), decreasing = T)[1:4]
cnv.cluster.cut[cnv.cluster.cut>2] = 3
cnv.cluster = make.initial.cnvHeatmap(suva_expr.cnv.scaled, meta, output = "JAS.1-2pearson.cnv.png", k.cut = 10, 
                                      cluster.method = "average", cluster.distance = "pearson", 
                                      cluster.cut = cnv.cluster.cut)
#k2 is candidately cancer
cnv.cluster.cut[cnv.cluster.cut>2] = 2
table(cnv.cluster.cut)
cnv.final = subtractControl(suva_expr.cnv.scaled, cnv.cluster.cut, k.control = 2, choice = "mean")
head(cnv.final[,1:4])
cnv.cluster2 = make.initial.cnvHeatmap(cnv.final, meta, output="JAS.1-3.cnvCtl.png", k.cut = 3, 
                                       cluster.method = "average", cluster.distance = "pearson",
                                       cluster.cut = cnv.cluster.cut)


findsplitK(cnv.cluster2, N = 10, k.init = 3)
cnv.cluster.cut2 = cutree(cnv.cluster2, 2)
table(cnv.cluster.cut2)
cnv.cluster2 = make.initial.cnvHeatmap(cnv.final, meta, output="JAS.1-4.cnvCtl.png", k.cut = 6, 
                                       cluster.method = "average", cluster.distance = "pearson",
                                       cluster.cut = cnv.cluster.cut2)

table(cnv.cluster.cut2)
table(cnv.cluster.cut)

candi.malignant = names(which(cnv.cluster.cut2 ==1))
candi.malignant[1:10]
write.table(candi.malignant, file = "JAS.candi.malignant.txt", row.names = F, col.names = F, quote = F)
####





############################################################################################################
rm(chr.bin.genes, cnv.cluster, cnv.cluster2, cnv.final, meta, suva_expr, suva_expr.cnv.raw, suva_expr.cnv.scaled, 
   cnv.cluster.cut, cnv.cluster.cut2, k)
load("JWS.Matrix.rda")
max(suva_expr)
min(suva_expr)

plot(density(rowMeans(suva_expr)))
chr.bin.genes = make_cnv_bins(gene_annotation,keep.genes = rownames(suva_expr), chromosomes = 1:22)
suva_expr.cnv.raw = cnv.bin(t(suva_expr), chr.bin.genes )

head(suva_expr.cnv.raw[,1:4])
max(suva_expr.cnv.raw); min(suva_expr.cnv.raw)

#center the data, sequencing depth
suva_expr.cnv.scaled = scale(suva_expr.cnv.raw, center = T, scale = T)
max(suva_expr.cnv.scaled);min(suva_expr.cnv.scaled)

plot(density(as.vector(suva_expr.cnv.scaled)))
suva_expr.cnv.scaled[suva_expr.cnv.scaled>3] = 3
suva_expr.cnv.scaled[suva_expr.cnv.scaled < -3] = -3

###JWS
cnv.cluster = make.initial.cnvHeatmap(suva_expr.cnv.scaled, meta, output = "JWS.1-1pearson.cnv.png", k.cut = 3, 
                                      cluster.method = "average", cluster.distance = "pearson")
dim(suva_expr.cnv.scaled)

findsplitK(cnv.cluster, n.cut = 3, k.init = 3, N = 10)
k = 7
cnv.cluster.cut = cutree(cnv.cluster, k)
sort(table(cnv.cluster.cut), decreasing = T)[1:4]
#cnv.cluster.cut[cnv.cluster.cut>5] = 5
cnv.cluster = make.initial.cnvHeatmap(suva_expr.cnv.scaled, meta, output = "JWS.1-2pearson.cnv.png", k.cut = 10, 
                                      cluster.method = "average", cluster.distance = "pearson", 
                                      cluster.cut = cnv.cluster.cut)
#k3 is candidately cancer
cnv.cluster.cut[cnv.cluster.cut != 3] = 1
table(cnv.cluster.cut)
cnv.final = subtractControl(suva_expr.cnv.scaled, cnv.cluster.cut, k.control = 1, choice = "mean")
head(cnv.final[,1:4])
cnv.cluster2 = make.initial.cnvHeatmap(cnv.final, meta, output="JWS.1-3.cnvCtl.png", k.cut = 3, 
                                       cluster.method = "average", cluster.distance = "pearson",
                                       cluster.cut = cnv.cluster.cut)

findsplitKrev(cnv.cluster2, N = 50)
library(dendextend)
k2 = get_subdendrograms(as.dendrogram(cnv.cluster2), k = 2)
nleaves(k2[[1]])
nleaves(k2[[2]])

findsplitKrev(k2[[2]], N = 30, k.max = 90, n.cut = 1)
sub.cut = cutree(k2[[2]], k = 29)
sort(table(sub.cut), decreasing = T)[1:4]
sub.cut[sub.cut>3] = 4

sub.cut2 = cutree(k2[[1]], k = 2)
sub.cut2[sub.cut2 > 0] = 6

cnv.cluster.cut2 = c(sub.cut, sub.cut2)
cnv.cluster.cut2 = cnv.cluster.cut2[colnames(cnv.final)];
table(cnv.cluster.cut2)
#cnv.cluster.cut2[cnv.cluster.cut2 != 1] = 2

#findsplitKrev(cnv.cluster2, N = 30, k.max = 300)
#k = 125
#cnv.cluster.cut2 = cutree(cnv.cluster2, k)
#sort(table(cnv.cluster.cut2), decreasing = T)[1:10]
#cnv.cluster.cut2[!cnv.cluster.cut2 %in% c(1,2,14)] = 10

cnv.cluster2 = make.initial.cnvHeatmap(cnv.final, meta, output="JWS.1-4.cnvCtl.png", k.cut = 6, 
                                       cluster.method = "average", cluster.distance = "pearson",
                                       cluster.cut = cnv.cluster.cut2)

table(cnv.cluster.cut2)
table(cnv.cluster.cut)

candi.malignant = names(which(cnv.cluster.cut2 ==2))
candi.malignant[1:10]
##not very good
write.table(candi.malignant, file = "JWS.candi.malignant.txt", row.names = F, col.names = F, quote = F)
####







############################################################################################################
rm(chr.bin.genes, cnv.cluster, cnv.cluster2, cnv.final, meta, suva_expr, suva_expr.cnv.raw, suva_expr.cnv.scaled, 
   cnv.cluster.cut, cnv.cluster.cut2, k)
load("RPO.Matrix.rda")
max(suva_expr)
min(suva_expr)

plot(density(rowMeans(suva_expr)))
chr.bin.genes = make_cnv_bins(gene_annotation,keep.genes = rownames(suva_expr), chromosomes = 1:22)
suva_expr.cnv.raw = cnv.bin(t(suva_expr), chr.bin.genes )

head(suva_expr.cnv.raw[,1:4])
max(suva_expr.cnv.raw); min(suva_expr.cnv.raw)

#center the data, sequencing depth
suva_expr.cnv.scaled = scale(suva_expr.cnv.raw, center = T, scale = T)
max(suva_expr.cnv.scaled);min(suva_expr.cnv.scaled)

plot(density(as.vector(suva_expr.cnv.scaled)))
suva_expr.cnv.scaled[suva_expr.cnv.scaled>3] = 3
suva_expr.cnv.scaled[suva_expr.cnv.scaled < -3] = -3

###RPO
cnv.cluster = make.initial.cnvHeatmap(suva_expr.cnv.scaled, meta, output = "RPO.1-1pearson.cnv.png", k.cut = 3, 
                                      cluster.method = "average", cluster.distance = "pearson")
dim(suva_expr.cnv.scaled)

findsplitK(cnv.cluster, N = 20, k.init = 10, n.cut = 3)
k = 30
cnv.cluster.cut = cutree(cnv.cluster, k)
sort(table(cnv.cluster.cut), decreasing = T)[1:4]
cnv.cluster.cut[cnv.cluster.cut>4] = 5
cnv.cluster.cut[cnv.cluster.cut < 3] = 3
cnv.cluster = make.initial.cnvHeatmap(suva_expr.cnv.scaled, meta, output = "RPO.1-2pearson.cnv.png", k.cut = 10, 
                                      cluster.method = "average", cluster.distance = "pearson", 
                                      cluster.cut = cnv.cluster.cut)
#k4 is candidately cancer
cnv.cluster.cut[cnv.cluster.cut != 4] = 1
table(cnv.cluster.cut)
cnv.final = subtractControl(suva_expr.cnv.scaled, cnv.cluster.cut, k.control = 1, choice = "mean")
head(cnv.final[,1:4])
cnv.cluster2 = make.initial.cnvHeatmap(cnv.final, meta, output="RPO.1-3.cnvCtl.png", k.cut = 3, 
                                       cluster.method = "average", cluster.distance = "pearson",
                                       cluster.cut = cnv.cluster.cut)


#findsplitK(cnv.cluster2, N = 50, k.init = 50)
library(dendextend)
k2 = get_subdendrograms(as.dendrogram(cnv.cluster2), k = 2)
nleaves(k2[[1]])
nleaves(k2[[2]])
findsplitK(k2[[2]], N = 10)
sub.cut = cutree(k2[[2]], k = 13)
table(sub.cut)
sub.cut[sub.cut>3] = 4

sub.cut2 = cutree(k2[[1]], k = 2)
sub.cut2[sub.cut2 > 0] = 6

cnv.cluster.cut2 = c(sub.cut, sub.cut2)

cnv.cluster.cut2 = cnv.cluster.cut2[colnames(cnv.final)];
#cnv.cluster.cut2[cnv.cluster.cut2 = 3] = 7
cnv.cluster2 = make.initial.cnvHeatmap(cnv.final, meta, output="RPO.1-4.cnvCtl.png", k.cut = 6, 
                                       cluster.method = "average", cluster.distance = "pearson",
                                       cluster.cut = cnv.cluster.cut2)

table(cnv.cluster.cut2)
table(cnv.cluster.cut)

candi.malignant = names(which(cnv.cluster.cut2 ==2))
candi.malignant[1:10]
write.table(candi.malignant, file = "RPO.candi.malignant.txt", row.names = F, col.names = F, quote = F)

candi.malignant = names(which(cnv.cluster.cut2 !=6))
candi.malignant[1:10]
write.table(candi.malignant, file = "RPO.candi.malignant.v2.txt", row.names = F, col.names = F, quote = F)



####





############################################################################################################
rm(chr.bin.genes, cnv.cluster, cnv.cluster2, cnv.final, meta, suva_expr, suva_expr.cnv.raw, suva_expr.cnv.scaled, 
   cnv.cluster.cut, cnv.cluster.cut2, k)
load("JWS.Matrix.rda")
max(suva_expr)
min(suva_expr)

plot(density(rowMeans(suva_expr)))
chr.bin.genes = make_cnv_bins(gene_annotation,keep.genes = rownames(suva_expr), chromosomes = 1:22)
suva_expr.cnv.raw = cnv.bin(t(suva_expr), chr.bin.genes )

head(suva_expr.cnv.raw[,1:4])
max(suva_expr.cnv.raw); min(suva_expr.cnv.raw)

#center the data, sequencing depth
suva_expr.cnv.scaled = scale(suva_expr.cnv.raw, center = T, scale = T)
max(suva_expr.cnv.scaled);min(suva_expr.cnv.scaled)

plot(density(as.vector(suva_expr.cnv.scaled)))
suva_expr.cnv.scaled[suva_expr.cnv.scaled>3] = 3
suva_expr.cnv.scaled[suva_expr.cnv.scaled < -3] = -3

###JWS
cnv.cluster = make.initial.cnvHeatmap(suva_expr.cnv.scaled, meta, output = "JWS.1-1pearson.cnv.png", k.cut = 3, 
                                      cluster.method = "average", cluster.distance = "pearson")
dim(suva_expr.cnv.scaled)

findsplitK(cnv.cluster, n.cut = 3, k.init = 3, N = 10)
k = 7
cnv.cluster.cut = cutree(cnv.cluster, k)
sort(table(cnv.cluster.cut), decreasing = T)[1:4]
#cnv.cluster.cut[cnv.cluster.cut>5] = 5
cnv.cluster = make.initial.cnvHeatmap(suva_expr.cnv.scaled, meta, output = "JWS.1-2pearson.cnv.png", k.cut = 10, 
                                      cluster.method = "average", cluster.distance = "pearson", 
                                      cluster.cut = cnv.cluster.cut)
#k2 is candidately cancer
cnv.cluster.cut[cnv.cluster.cut != 3] = 1
table(cnv.cluster.cut)
cnv.final = subtractControl(suva_expr.cnv.scaled, cnv.cluster.cut, k.control = 1, choice = "mean")
head(cnv.final[,1:4])
cnv.cluster2 = make.initial.cnvHeatmap(cnv.final, meta, output="JWS.1-3.cnvCtl.png", k.cut = 3, 
                                       cluster.method = "average", cluster.distance = "pearson",
                                       cluster.cut = cnv.cluster.cut)

findsplitK(cnv.cluster2, N = 50, k.init = 3)
library(dendextend)
k2 = get_subdendrograms(as.dendrogram(cnv.cluster2), k = 2)
nleaves(k2[[1]])
nleaves(k2[[2]])
findsplitK(k2[[2]], N = 10)
sub.cut = cutree(k2[[2]], k = 3)
table(sub.cut)
#sub.cut[sub.cut>3] = 4

sub.cut2 = cutree(k2[[1]], k = 2)
sub.cut2[sub.cut2 > 0] = 6

cnv.cluster.cut2 = c(sub.cut, sub.cut2)
cnv.cluster.cut2 = cnv.cluster.cut2[colnames(cnv.final)];
table(cnv.cluster.cut2)
cnv.cluster.cut2[cnv.cluster.cut2 != 1] = 2
cnv.cluster2 = make.initial.cnvHeatmap(cnv.final, meta, output="JWS.1-4.cnvCtl.png", k.cut = 6, 
                                       cluster.method = "average", cluster.distance = "pearson",
                                       cluster.cut = cnv.cluster.cut2)

table(cnv.cluster.cut2)
table(cnv.cluster.cut)

candi.malignant = names(which(cnv.cluster.cut2 ==2))
candi.malignant[1:10]
##not very good
write.table(candi.malignant, file = "JWS.candi.malignant.txt", row.names = F, col.names = F, quote = F)
####





############################################################################################################
rm(chr.bin.genes, cnv.cluster, cnv.cluster2, cnv.final, meta, suva_expr, suva_expr.cnv.raw, suva_expr.cnv.scaled, 
   cnv.cluster.cut, cnv.cluster.cut2, k)
load("SAB.Matrix.rda")
max(suva_expr)
min(suva_expr)

plot(density(rowMeans(suva_expr)))
chr.bin.genes = make_cnv_bins(gene_annotation,keep.genes = rownames(suva_expr), chromosomes = 1:22)
suva_expr.cnv.raw = cnv.bin(t(suva_expr), chr.bin.genes )

head(suva_expr.cnv.raw[,1:4])
max(suva_expr.cnv.raw); min(suva_expr.cnv.raw)

#center the data, sequencing depth
suva_expr.cnv.scaled = scale(suva_expr.cnv.raw, center = T, scale = T)
max(suva_expr.cnv.scaled);min(suva_expr.cnv.scaled)

plot(density(as.vector(suva_expr.cnv.scaled)))
suva_expr.cnv.scaled[suva_expr.cnv.scaled>3] = 3
suva_expr.cnv.scaled[suva_expr.cnv.scaled < -3] = -3

###SAB
cnv.cluster = make.initial.cnvHeatmap(suva_expr.cnv.scaled, meta, output = "SAB.1-1pearson.cnv.png", k.cut = 3, 
                                      cluster.method = "average", cluster.distance = "pearson")
dim(suva_expr.cnv.scaled)

findsplitK(cnv.cluster)
k = 46
cnv.cluster.cut = cutree(cnv.cluster, k)
sort(table(cnv.cluster.cut), decreasing = T)[1:4]
cnv.cluster.cut[cnv.cluster.cut>3] = 3
cnv.cluster = make.initial.cnvHeatmap(suva_expr.cnv.scaled, meta, output = "SAB.1-2pearson.cnv.png", k.cut = 10, 
                                      cluster.method = "average", cluster.distance = "pearson", 
                                      cluster.cut = cnv.cluster.cut)
#k1 is candidately cancer
cnv.cluster.cut[cnv.cluster.cut > 2] = 2
table(cnv.cluster.cut)
cnv.final = subtractControl(suva_expr.cnv.scaled, cnv.cluster.cut, k.control = 2, choice = "mean")
head(cnv.final[,1:4])
cnv.cluster2 = make.initial.cnvHeatmap(cnv.final, meta, output="SAB.1-3.cnvCtl.png", k.cut = 3, 
                                       cluster.method = "average", cluster.distance = "pearson",
                                       cluster.cut = cnv.cluster.cut)

findsplitK(cnv.cluster2, N = 1600)

library(dendextend)
k2 = get_subdendrograms(as.dendrogram(cnv.cluster2), k = 2)
nleaves(k2[[1]])
nleaves(k2[[2]])
findsplitKrev(k2[[1]], N = 1630, n.cut = 1, k.max = 480)
sub.cut = cutree(k2[[1]], k = 448)
sort(table(sub.cut), decreasing = T)[1:4]
sub.cut[sub.cut>2] = 3

sub.cut2 = cutree(k2[[2]], k = 2)
sub.cut2[sub.cut2 > 0] = 6

cnv.cluster.cut2 = c(sub.cut, sub.cut2)
cnv.cluster.cut2 = cnv.cluster.cut2[colnames(cnv.final)];
table(cnv.cluster.cut2)
#cnv.cluster.cut2[cnv.cluster.cut2 != 1] = 2
cnv.cluster2 = make.initial.cnvHeatmap(cnv.final, meta, output="SAB.1-4.cnvCtl.png", k.cut = 6, 
                                       cluster.method = "average", cluster.distance = "pearson",
                                       cluster.cut = cnv.cluster.cut2)

table(cnv.cluster.cut2)
table(cnv.cluster.cut)

candi.malignant = names(which(cnv.cluster.cut2 ==1))
candi.malignant[1:10]
##not very good
write.table(candi.malignant, file = "SAB.candi.malignant.txt", row.names = F, col.names = F, quote = F)
####






