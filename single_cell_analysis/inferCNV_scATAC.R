
rm(list=ls())
getwd()
setwd("/Users/wuz6/Project/01.single_cell_ATAC/IDHwt.GBM/scATAC/CNV/")
#devtools::install_github("jokergoo/ComplexHeatmap") #version 
source("cnv.plot.src.R")
library(dendextend)

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
list.files()
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

# GBM
#?#YNA	6/17/19	Ö	Ö	x	V321	61	Male	DMG, K27	0.5	DMG, K27
#NO#KYL	1/3/20	Ö	Ö	x	W909	29	Male	EPN, MPE	1	EPN, MPE

#RCA	10/11/19	Ö	Ö	x	W327	73	Male	MTGF_GBM	1	GBM, MID ----
#JDD	4/15/19	x	Ö	x	U840	49	Male	MTGF_GBM	1	GBM, RTK I
#DMD	5/13/19	x	Ö	x	V059	64	Male	MTGF_GBM	1	GBM, RTK I
#CCH	8/23/19	Ö	Ö	x	V818	47	Female	MTGF_GBM	1	GBM, RTK I
#ROB	10/30/19	Ö	Ö	x	W413	68	Male	MTGF_GBM	1	GBM, RTK II -------

#######################################################################################################################################



#######################################################################################################################################
# RCA
cnv = readRDS("RCA.CNVraw.rds")
pos = readRDS("RCA.CNVpos.rds")
dim(cnv)
barcodes = readRDS("RCA.barcodes.rds")
colnames(cnv) = gsub("RCA#", "", colnames(cnv))
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3
dir.create("RCA")

#1.raw correlation and CNV
cnv.cluster = make_correlation_heatmap(cnv, "RCA/RCA.1_1.correlation.png", barcodes)
p = make_raw_cnv_plot(cnv, pos, cnv.cluster, "RCA/RCA.1_2.CNVraw.png")



findsplitK(cnv.cluster, 200)
k.split = 10
sort(table(cutree(cnv.cluster, k.split)), decreasing = T)[1:10]
cell.cluster.cut = cutree(cnv.cluster, k.split)
cell.cluster.cut[cell.cluster.cut > 2] = 2
p  = make_leftAnn_cnvplot(cnv, pos, "RCA/RCA.1_3.CNV.clusterAnn.png", cnv.cluster, cell.cluster.cut, plot.dend = T)

# subtract normal cell CNV 
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 2, choice = "median")

#2. test if cnv2 has better clustering
cnv.cluster2 = make_correlation_heatmap(cnv2, "RCA/RCA.2_1.correlation.png", barcodes)
p = make_raw_cnv_plot(cnv2, pos, cnv.cluster2, "RCA/RCA.2_2.cnv2.png")
#resulting dendrogram is hard to splot
library(dendextend)
subdend = get_subdendrograms(as.dendrogram(cnv.cluster2), k = 2)

#subdend[[2]] #nomral

findsplitK(subdend[[1]], N = 100, k.init = 3)
sort( table(cutree(subdend[[1]], k = 487)), decreasing = T)[1:10]
rm(subdend)

dim(cnv)
findsplitKrev(cnv.cluster2, 2850, k.max = 325, n.cut = 1)
k.split = 235
cell.cluster.cut.2 = cutree(cnv.cluster2, k.split)
sort(table(cell.cluster.cut.2), decreasing = T)[1:4]
cell.cluster.cut.2[cell.cluster.cut.2>2] =2
p = make_leftAnn_cnvplot(cnv2, pos, "RCA/RCA.2_3.cnv2.clusterAnn.png", cnv.cluster2, cell.cluster.cut.2, plot.dend = T)

# make final plot, 
# as subclone observed, use dend to order
p = make_final_cnv_plot_Dendreorder(cnv2, pos, "RCA/RCA.4.final.png", cnv.cluster2, cell.cluster.cut.2, k.malignant = 1)
p = make_final_cnv_plot_Dendreorder(cnv2, pos, "RCA/RCA.4_2.final.png", cnv.cluster2, cell.cluster.cut.2, k.malignant = 1, cluster.ann = T)

candidate.malignant = names(which(cell.cluster.cut == 1))
length(candidate.malignant)
write.table(candidate.malignant, "RCA/RCA.4.CNVmalignant.txt", quote = F, row.names = F, col.names = F)

#Not perfect, some 'normal' cells on the top looks like malignant cells; subclone cluster not good



#######################################################################################################################################
rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)
#JDD
cnv = readRDS("JDD.CNVraw.rds")
pos = readRDS("JDD.CNVpos.rds")
barcodes = readRDS("JDD.barcodes.rds")
head(cnv[,1:2])
colnames(cnv) = gsub("JDD#", "", colnames(cnv))
dim(cnv)
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dir.create("JDD")
#1.raw correlation and CNV plots
cnv.cluster = make_correlation_heatmap(cnv, "JDD/JDD.1_1.correlation.png", barcodes)
p = make_raw_cnv_plot(cnv, pos, cnv.cluster, "JDD/JDD.1_2.CNVraw.png")

# find candicate malignant cells
dim(cnv)
findsplitK(cnv.cluster, 2000)
k.split = 36
sort(table(cutree(cnv.cluster, k.split)), decreasing = T)[1:10]
cell.cluster.cut = cutree(cnv.cluster, k.split)
table(cell.cluster.cut)
cell.cluster.cut[cell.cluster.cut>2] = 1

p = make_leftAnn_cnvplot(cnv, pos, "JDD/JDD.1_3.CNVraw.clusterAnn.png",cnv.cluster, cell.cluster.cut)


#2. subtract normal control, suppose cluster 1 is malignant
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 2, choice = "median")
cnv.cluster2 = make_correlation_heatmap(cnv2, "JDD/JDD.2_1.cnv2.correlation.png", barcodes, plot.dend = T) #not able to make the plot with dend
#plot(cnv.cluster2, ann = F,axes = F, labels=FALSE)
p = make_raw_cnv_plot(cnv2, pos, cnv.cluster2, "JDD/JDD.2_2.cnv2.cluster.png")
#plot(cnv.cluster2, ann = F,axes = F, labels=FALSE)

dim(cnv)
sort(table(cutree(cnv.cluster2, 3110)),decreasing = T)[1:2]
table(cell.cluster.cut)
#cancer cell < 4698
findsplitKrev(cnv.cluster2, 4500, k.max = 3000, n.cut = 1)
findsplitKrev(cnv.cluster2, 4360, k.max = 2400, n.cut = 1 )


#2385 is close, but... maybe 
#k.split = 1075 is too small

#k 200 is too small
#k 2622 is high, some non-malignant cells were clustered as malignant
k.split = 1122
cell.cluster.cut.2 = cutree(cnv.cluster2, k.split)
sort(table(cell.cluster.cut.2), decreasing = T)[1:4]
cell.cluster.cut.2[cell.cluster.cut.2>2] =2
p = make_leftAnn_cnvplot(cnv2, pos, "JDD/JDD.2_3.cnv2.clusterAnn.png", cnv.cluster2, cell.cluster.cut.2, plot.dend = T)


#make final plot
p = make_final_cnv_plot_Dendreorder(cnv2, pos, "JDD/JDD.4.final.png", cnv.cluster2, cell.cluster.cut.2, k.malignant = 1)
p = make_final_cnv_plot_Dendreorder(cnv2, pos, "JDD/JDD.4_2.final.png", cnv.cluster2, cell.cluster.cut.2, k.malignant = 1, cluster.ann = T)

table(cell.cluster.cut.2)

candidate.malignant = names(which(cell.cluster.cut.2 == 1))
length(candidate.malignant)
write.table(candidate.malignant, "JDD/JDD.4.CNVmalignant.txt", quote = F, row.names = F, col.names = F)




#######################################################################################################################################
rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)
#DMD
cnv = readRDS("DMD.CNVraw.rds")
pos = readRDS("DMD.CNVpos.rds")
barcodes = readRDS("DMD.barcodes.rds")
head(cnv[,1:2])
colnames(cnv) = gsub("DMD#", "", colnames(cnv))
dim(cnv)
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dir.create("DMD")
#1.raw correlation and CNV plot
cnv.cluster = make_correlation_heatmap(cnv, "DMD/DMD.1_1.correlation.png", barcodes)
p = make_raw_cnv_plot(cnv, pos, cnv.cluster, "DMD/DMD.1_2.CNVraw.png")
# find candicate malignant cells
dim(cnv)
findsplitK(cnv.cluster, 1000)
k.split = 12
sort(table(cutree(cnv.cluster, k.split)), decreasing = T)[1:10]
cell.cluster.cut = cutree(cnv.cluster, k.split)
table(cell.cluster.cut)
cell.cluster.cut[cell.cluster.cut>2] = 3
p = make_leftAnn_cnvplot(cnv, pos, "DMD/DMD.1_3.CNVraw.clusterAnn.png", cnv.cluster, cell.cluster.cut)

#cluster 1 looks like malignant

#2 subtract normal control
cell.cluster.cut[cell.cluster.cut>2] =2
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 2, choice = "median")
cnv.cluster2 = make_correlation_heatmap(cnv2, "DMD/DMD.2_1.cnv2.correlation.png", barcodes)
p = make_raw_cnv_plot(cnv2, pos, cnv.cluster2, "DMD/DMD.2_2.cnv2.cluster.png")

#find malignant cluster
dim(cnv)
findsplitKrev(cnv.cluster2, 1900, k.max = 3000, n.cut = 1)

#2583 is too high, more cells are malignant
k.split = 2181 
sort(table(cutree(cnv.cluster2, k.split)), decreasing = T)[1:10]

cell.cluster.cut.2 = cutree(cnv.cluster2, k.split)
cell.cluster.cut.2[cell.cluster.cut.2>2] = 2
table(cell.cluster.cut.2)

p = make_leftAnn_cnvplot(cnv2, pos, "DMD/DMD.2_3.CNVraw.clusterAnn.png", cnv.cluster2, cell.cluster.cut.2)


#4.make final plot
p = make_final_cnv_plot_Dendreorder(cnv2, pos, "DMD/DMD.4_2.final.png", cnv.cluster2, cell.cluster.cut.2, k.malignant = 1, cluster.ann = T)

p = make_final_cnv_plot_Dendreorder(cnv2, pos, "DMD/DMD.4.final.png", cnv.cluster2, cell.cluster.cut.2, k.malignant = 1)

candidate.malignant = names(which(cell.cluster.cut.2 == 1))
length(candidate.malignant)
write.table(candidate.malignant, "DMD/DMD.4.CNVmalignant.txt", quote = F, row.names = F, col.names = F)





#######################################################################################################################################
rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)
#CCH
cnv = readRDS("CCH.CNVraw.rds")
pos = readRDS("CCH.CNVpos.rds")
barcodes = readRDS("CCH.barcodes.rds")
head(cnv[,1:2])
colnames(cnv) = gsub("CCH#", "", colnames(cnv))
dim(cnv)

barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dir.create("CCH")
#1.raw correlation and CNV plot
cnv.cluster = make_correlation_heatmap(cnv, "CCH/CCH.1_1.correlation.png", barcodes)
p = make_raw_cnv_plot(cnv, pos, cnv.cluster, "CCH/CCH.1_2.CNVraw.png")
# find candicate malignant cells
dim(cnv)

findsplitK(cnv.cluster, 10)
k.split = 4
cell.cluster.cut = cutree(cnv.cluster, k.split)
sort(table(cell.cluster.cut), decreasing = T)[1:4]
cell.cluster.cut[cell.cluster.cut>2] = 2
p = make_leftAnn_cnvplot(cnv, pos, "CCH/CCH.1_3.CNVraw.clusterAnn.png", cnv.cluster, cell.cluster.cut)
table(cell.cluster.cut)
#2. subtract control, looks cluster 1 is malignant
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 2, choice = "median")
cnv.cluster2 = make_correlation_heatmap(cnv2, "CCH/CCH.2_1.correlation.png", barcodes)
#not good #cnv.cluster2 = make_correlation_heatmap(cnv2, "CCH/CCH.2_1.correlation.png", barcodes, cluster.method = "complete")
p = make_raw_cnv_plot(cnv2, pos, cnv.cluster2, "CCH/CCH.2_2.CNV2.png")

#re-cluster is 
dim(cnv)
findsplitKrev(cnv.cluster2, 800, k.max = 50 ,n.cut = 1)
findsplitK(cnv.cluster2, 4, k.init = 10)

k.split = 36
cell.cluster.cut.2 = cutree(cnv.cluster2, k.split)
sort(table(cell.cluster.cut.2), decreasing = T)[1:4]
cell.cluster.cut.2[cell.cluster.cut.2>2] = 2
p = make_leftAnn_cnvplot(cnv2, pos, "CCH/CCH.2_3.CNV2.clusterAnn.png", cnv.cluster2, cell.cluster.cut.2)

#4.make final plot
p = make_final_cnv_plot_Dendreorder(cnv2, pos, "CCH/CCH.4.final.png", cnv.cluster2, cell.cluster.cut.2, k.malignant = 1)

candidate.malignant = names(which(cell.cluster.cut.2 == 1))
length(candidate.malignant)
write.table(candidate.malignant, "CCH/CCH.4.CNVmalignant.txt", quote = F, row.names = F, col.names = F)

#low cell number 


#######################################################################################################################################
rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)
#ROB
cnv = readRDS("ROB.CNVraw.rds")
pos = readRDS("ROB.CNVpos.rds")
barcodes = readRDS("ROB.barcodes.rds")
head(cnv[,1:2])
colnames(cnv) = gsub("ROB#", "", colnames(cnv))
dim(cnv)
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dir.create("ROB")
#1.raw correlation and CNV plot
cnv.cluster = make_correlation_heatmap(cnv, "ROB/ROB.1_1.correlation.png", barcodes, plot.dend = F)
p = make_raw_cnv_plot(cnv, pos, cell.cluster = F, png.file = "ROB/ROB.1_2.CNVraw.png")
p = make_rawOrdered_cnv_plot(cnv, pos, cell.cluster = cnv.cluster, png.file = "ROB/ROB.1_2.CNVraw.png")

plot(cnv.cluster, ann = F,axes = F, labels=FALSE)
# find candicate malignant cells
findsplitK(cnv.cluster, k.init = 5)
k.split = 44
cell.cluster.cut = cutree(cnv.cluster, k.split)
sort(table(cell.cluster.cut), decreasing = T)[1:4]
cell.cluster.cut[cell.cluster.cut>2] = 3
p = make_leftAnn_cnvplot(cnv, pos, "ROB/ROB.1_3.CNVraw.clusterAnn.png", cell.cluster = cnv.cluster, cell.cluster.cut)

cell.cluster.cut[cell.cluster.cut>2] = 2
#2.subtract CNV control
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 2, choice = "median")
#cnv.cluster2 = make_correlation_heatmap(cnv2, "ROB/ROB.2_1.cnv2.correlation.png", barcodes) #failed
cnv.cluster2 = make_correlation_heatmap(cnv2, "ROB/ROB.2_1.cnv2.correlation.png", barcodes, plot.dend = F)
p = make_rawOrdered_cnv_plot(cnv2, pos, cnv.cluster2, "ROB/ROB.2_2.cnv2.png")

plot(cnv.cluster2, ann = F,axes = F, labels=FALSE)

findsplitK(cnv.cluster2, N = 100, k.init = 10, n.cut = 3)
k.split = 87
sort(table(cutree(cnv.cluster2, 20)), decreasing = T)[1:4]

cell.cluster.cut.2 = cutree(cnv.cluster2, k.split)
sort(table(cutree(cnv.cluster2,  k.split)), decreasing = T)[1:4]
cell.cluster.cut.2[cell.cluster.cut.2>2] = 2
p = make_leftAnn_cnvplot(cnv2, pos, "ROB/ROB.2_3.cnv2.clusterAnn.png", cnv.cluster2, cell.cluster.cut.2)
#4.make final plot
#the cluster performance is not good, some cells without CNV were clustered with CNV group
p = make_final_cnv_plot(cnv2, pos, "ROB/ROB.4.final.png", cell.cluster.cut.2, k.malignant = 1)


candidate.malignant = names(which(cell.cluster.cut.2 == 1))
length(candidate.malignant)
write.table(candidate.malignant, "ROB/ROB.4.CNVmalignant.txt", quote = F, row.names = F, col.names = F)



#######################################################################################################################################
rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)
#ROB
cnv = readRDS("YNA.CNVraw.rds")
pos = readRDS("YNA.CNVpos.rds")
barcodes = readRDS("YNA.barcodes.rds")
head(cnv[,1:2])
colnames(cnv) = gsub("YNA#", "", colnames(cnv))
dim(cnv)
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dir.create("YNA")
#1.raw correlation and CNV plot
cnv.cluster = make_correlation_heatmap(cnv, "YNA/YNA.1_1.correlation.png", barcodes)
p = make_raw_cnv_plot(cnv, pos, cnv.cluster, png.file = "YNA/YNA.1_2.CNVraw.png")



#not able to tell
# find candicate malignant cells
ncol(cnv)
findsplitK(cnv.cluster)
k.split = 352
cell.cluster.cut = cutree(cnv.cluster, k.split)
sort(table(cell.cluster.cut), decreasing = T)[1:4]
cell.cluster.cut[cell.cluster.cut>5] = 6
p = make_leftAnn_cnvplot(cnv, pos, "YNA/YNA.1_3.CNVraw.clusterAnn.png", cell.cluster = cnv.cluster, cell.cluster.cut)

#2.subtract CNV control
#cell.cluster.cut[cell.cluster.cut !=4] = 1
table(cell.cluster.cut)
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 5, choice = "median")
#cnv.cluster2 = make_correlation_heatmap(cnv2, "ROB/ROB.2_1.cnv2.correlation.png", barcodes) #failed
cnv.cluster2 = make_correlation_heatmap(cnv2, "YNA/YNA.2_1.cnv2.correlation.png", barcodes)
p = make_raw_cnv_plot(cnv2, pos, cnv.cluster2, "YNA/YNA.2_2.cnv2.png")


findsplitK(cnv.cluster2)
k.split = 87
sort(table(cutree(cnv.cluster2, 20)), decreasing = T)[1:4]

cell.cluster.cut.2 = cutree(cnv.cluster2, k.split)
sort(table(cutree(cnv.cluster2,  k.split)), decreasing = T)[1:4]
cell.cluster.cut.2[cell.cluster.cut.2>2] = 2
p = make_leftAnn_cnvplot(cnv2, pos, "YNA/YNA.2_3.cnv2.clusterAnn.png", cnv.cluster2, cell.cluster.cut.2)
#4.make final plot
#the cluster performance is not good, some cells without CNV were clustered with CNV group
p = make_final_cnv_plot(cnv2, pos, "YNA/YNA.4.final.png", cell.cluster.cut.2, k.malignant = 1)


candidate.malignant = names(which(cell.cluster.cut.2 == 1))
length(candidate.malignant)
#write.table(candidate.malignant, "YNA/YNA.4.CNVmalignant.txt", quote = F, row.names = F, col.names = F)



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################











#DBG	A IDH
#JGO	A IDH, HG
#CEH	A IDH, HG
#JAS	A IDH, HG
#SAB  A IDH

rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)
#DBG
cnv = readRDS("DBG.CNVraw.rds")
dim(cnv)
pos = readRDS("DBG.CNVpos.rds")
barcodes = readRDS("DBG.barcodes.rds")
head(cnv[,1:2])
colnames(cnv) = gsub("DBG#", "", colnames(cnv))
dim(cnv)
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dir.create("DBG")
#1.raw correlation and CNV plot

#methylation CNV: 7p,10q,12q loss, 7q,12qEnd gain,
cnv.cluster = make_correlation_heatmap(cnv, "DBG/DBG.1_1.correlation.png", barcodes)
p = make_raw_cnv_plot(cnv, pos, cell.cluster = cnv.cluster, png.file = "DBG/DBG.1_2.CNVraw.png")

findsplitK(cnv.cluster)

k.split = 85
cell.cluster.cut = cutree(cnv.cluster, k.split)
sort(table(cell.cluster.cut), decreasing = T)[1:4]
cell.cluster.cut[cell.cluster.cut>2] = 3
p = make_leftAnn_cnvplot(cnv, pos, "DBG/DBG.1_3.CNVraw.clusterAnn.png", cell.cluster = cnv.cluster, 
                         cell.cluster.cut)



cell.cluster.cut[cell.cluster.cut>2] = 2
#2.subtract CNV control
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 2, choice = "median")
#cnv.cluster2 = make_correlation_heatmap(cnv2, "ROB/ROB.2_1.cnv2.correlation.png", barcodes) #failed
cnv.cluster2 = make_correlation_heatmap(cnv2, "DBG/DBG.2_1.cnv2.correlation.png", barcodes)
p = make_rawOrdered_cnv_plot(cnv2, pos, cnv.cluster2, "DBG/DBG.2_2.cnv2.png")

#plot(cnv.cluster2, ann = F,axes = F, labels=FALSE)
dim(cnv)
k2 = get_subdendrograms(as.dendrogram(cnv.cluster2), k = 2)
nleaves(k2[[1]])
nleaves(k2[[2]])

findsplitK(k2[[1]], N = 20)
sub.cut1 = cutree(k2[[1]], k = 66)
sort(table(sub.cut1), decreasing = T)[1:4]
sub.cut1[sub.cut1>2] = 2;
table(sub.cut1);

sub.cut2 = cutree(k2[[2]], k = 2)
sub.cut2[sub.cut2>0] = 3

cell.cluster.cut.2 = c(sub.cut1, sub.cut2)[colnames(cnv2)]
p = make_leftAnn_cnvplot(cnv2, pos, "DBG/DBG.2_3.cnv2.clusterAnn.png", cnv.cluster2, cell.cluster.cut.2)
#4.make final plot
#the cluster performance is not good, some cells without CNV were clustered with CNV group
p = make_final_cnv_plot(cnv2, pos, "DBG/DBG.4.final.png", cell.cluster.cut.2, k.malignant = 1)



candidate.malignant = names(which(cell.cluster.cut.2 == 1))
length(candidate.malignant)
write.table(candidate.malignant, "DBG/DBG.4.CNVmalignant.txt", quote = F, row.names = F, col.names = F)




#######################################################################################################################################


rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)
#JGO
cnv = readRDS("JGO.CNVraw.rds")
pos = readRDS("JGO.CNVpos.rds")
barcodes = readRDS("JGO.barcodes.rds")
head(cnv[,1:2])
colnames(cnv) = gsub("JGO#", "", colnames(cnv))
dim(cnv)
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dir.crea te("JGO")
#1.raw correlation and CNV plot
cnv.cluster = make_correlation_heatmap(cnv, "JGO/JGO.1_1.correlation.png", barcodes)
p = make_raw_cnv_plot(cnv, pos, cell.cluster = cnv.cluster, png.file = "JGO/JGO.1_2.CNVraw.png")

findsplitK(cnv.cluster)
k.split = 20
cell.cluster.cut = cutree(cnv.cluster, k.split)
sort(table(cell.cluster.cut), decreasing = T)[1:4]
cell.cluster.cut[cell.cluster.cut>2] = 3
p = make_leftAnn_cnvplot(cnv, pos, "JGO/JGO.1_3.CNVraw.clusterAnn.png", cell.cluster = cnv.cluster, 
                         cell.cluster.cut)



cell.cluster.cut[cell.cluster.cut>2] = 1
#2.subtract CNV control
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 1, choice = "median")

cnv.cluster2 = make_correlation_heatmap(cnv2, "JGO/JGO.2_1.cnv2.correlation.png", barcodes)
#p = make_rawOrdered_cnv_plot(cnv2, pos, cnv.cluster2, "JGO/JGO.2_2.cnv2.png")
p = make_raw_cnv_plot(cnv2, pos, cnv.cluster2, "JGO/JGO.2_2.cnv2.png")

ncol(cnv2)
findsplitKrev(cnv.cluster2, n.cut = 1, N = 2440, k.max = 1000)
kk = 844

cell.cluster.cut.2 = cutree(cnv.cluster2, k = kk)
sort(table(cell.cluster.cut.2), decreasing = T)[1:4]
cell.cluster.cut.2[!cell.cluster.cut.2 %in% c(3)] = 1
p = make_leftAnn_cnvplot(cnv2, pos, "JGO/JGO.2_3.cnv2.clusterAnn.png", cnv.cluster2, cell.cluster.cut.2)

#4.make final plot
#the cluster performance is not good, some cells without CNV were clustered with CNV group
p = make_final_cnv_plot(cnv2, pos, "JGO/JGO.4.final.png", cell.cluster.cut.2, k.malignant = 3)


candidate.malignant = names(which(cell.cluster.cut.2 == 3))
length(candidate.malignant)
write.table(candidate.malignant, "JGO/JGO.4.CNVmalignant.txt", quote = F, row.names = F, col.names = F)

#######################################################################################################################################


rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)
#CEH
cnv = readRDS("CEH.CNVraw.rds")
pos = readRDS("CEH.CNVpos.rds")
barcodes = readRDS("CEH.barcodes.rds")
head(cnv[,1:2])
colnames(cnv) = gsub("CEH#", "", colnames(cnv))
dim(cnv)
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
max(as.vector(cnv))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dir.create("CEH")
#1.raw correlation and CNV plot
#methylation CEH: 4gain
cnv.cluster = make_correlation_heatmap(cnv, "CEH/CEH.1_1.correlation.png", barcodes)
p = make_raw_cnv_plot(cnv, pos, cell.cluster = cnv.cluster, png.file = "CEH/CEH.1_2.CNVraw.png")

ncol(cnv)
findsplitK(cnv.cluster, N = 10)

k.split = 6
cell.cluster.cut = cutree(cnv.cluster, k.split)
sort(table(cell.cluster.cut), decreasing = T)[1:4]
cell.cluster.cut[cell.cluster.cut>2] = 3
p = make_leftAnn_cnvplot(cnv, pos, "CEH/CEH.1_3.CNVraw.clusterAnn.png", cell.cluster = cnv.cluster, 
                         cell.cluster.cut)



cell.cluster.cut[cell.cluster.cut>2] = 1
#2.subtract CNV control
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 1, choice = "median")
cnv.cluster2 = make_correlation_heatmap(cnv2, "CEH/CEH.2_1.cnv2.correlation.png", barcodes)
p = make_leftAnn_cnvplot(cnv2, pos, "CEH/CEH.2_2.cnv2.png", cnv.cluster2, cell.cluster.cut)

findsplitK(cnv.cluster2, N = 10, k.init = 10)
p = make_leftAnn_cnvplot(cnv2, pos, "CEH/CEH.2_3.cnv2.clusterAnn.png", cnv.cluster2, cell.cluster.cut)
#4.make final plot
#the cluster performance is not good, some cells without CNV were clustered with CNV group
p = make_final_cnv_plot(cnv2, pos, "CEH/CEH.4.final.png", cell.cluster.cut, k.malignant = 2)

table(cell.cluster.cut)

candidate.malignant = names(which(cell.cluster.cut == 2))
length(candidate.malignant)
write.table(candidate.malignant, "CEH/CEH.4.CNVmalignant.txt", quote = F, row.names = F, col.names = F)

#######################################################################################################################################


rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)
#JAS
cnv = readRDS("JAS.CNVraw.rds")
pos = readRDS("JAS.CNVpos.rds")
barcodes = readRDS("JAS.barcodes.rds")
head(cnv[,1:2])
colnames(cnv) = gsub("JAS#", "", colnames(cnv))
dim(cnv)
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dir.create("JAS")
#1.raw correlation and CNV plot
cnv.cluster = make_correlation_heatmap(cnv, "JAS/JAS.1_1.correlation.png", barcodes)
p = make_raw_cnv_plot(cnv, pos, cell.cluster = cnv.cluster, png.file = "JAS/JAS.1_2.CNVraw.png")
ncol(cnv)
findsplitK(cnv.cluster, N = 10)
k.split = 26
cell.cluster.cut = cutree(cnv.cluster, k.split)
sort(table(cell.cluster.cut), decreasing = T)[1:4]
cell.cluster.cut[cell.cluster.cut>2] = 3
p = make_leftAnn_cnvplot(cnv, pos, "JAS/JAS.1_3.CNVraw.clusterAnn.png", cell.cluster = cnv.cluster, 
                         cell.cluster.cut)



cell.cluster.cut[cell.cluster.cut>2] = 2
#2.subtract CNV control
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 2, choice = "median")
cnv.cluster2 = make_correlation_heatmap(cnv2, "JAS/JAS.2_1.cnv2.correlation.png", barcodes)
p = make_raw_cnv_plot(cnv2, pos, cell.cluster = cnv.cluster2, "JAS/JAS.2_2.cnv2.png")
  
ncol(cnv2)
findsplitKrev(cnv.cluster2, N = 2950, n.cut = 1, k.max =33)

sort(table(cutree(cnv.cluster2, k = 16)), decreasing = T)[1:4]

cell.cluster.cut.2 = cutree(cnv.cluster2, k = 16)
cell.cluster.cut.2[cell.cluster.cut.2>2] = 3

p = make_leftAnn_cnvplot(cnv2, pos, "JAS/JAS.2_3.cnv2.clusterAnn.png", cnv.cluster2, cell.cluster.cut.2)


#4.make final plot
p = make_final_cnv_plot(cnv2, pos, "JAS/JAS.4.final.png", cell.cluster.cut.2, k.malignant = 1)



candidate.malignant = names(which(cell.cluster.cut.2 == 1))
length(candidate.malignant)
write.table(candidate.malignant, "JAS/JAS.4.CNVmalignant.txt", quote = F, row.names = F, col.names = F)

###SAB
###VBH

rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)

sample.name = "VBH"
cnv = readRDS(paste0(sample.name, ".CNVraw.rds"))
pos = readRDS(paste0(sample.name, ".CNVpos.rds"))
barcodes = readRDS(paste0(sample.name, ".barcodes.rds"))
head(cnv[,1:2])
colnames(cnv) = gsub(paste0(sample.name,"#"), "", colnames(cnv))
dim(cnv)
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dir.create(sample.name)
#1.raw correlation and CNV plot
cnv.cluster = make_correlation_heatmap(cnv, paste0(sample.name,"/", sample.name, ".1_1.correlation.png"), barcodes)
p = make_raw_cnv_plot(cnv, pos, cell.cluster = cnv.cluster, png.file = paste0(sample.name, "/", sample.name, ".1_2.CNVraw.png"))
ncol(cnv)

load("/Users/wuz6/Project/01.single_cell_ATAC/IDH.glioma/scATAC/01.all_cell/vbh.predicted.celltype.rda")
head(cnv[,1:4])
head(vbh)
rownames(vbh) = gsub("VBH#", "", rownames(vbh))
table(vbh[,2])

vbh[,2][!vbh[,2] %in% c(3,9,12, 14,16,19,20)] = 1
cell.col = c("Malignant" = "red",  "Oligodendrocyte" = "green", "Microglia/Macrophage" = "blue")
cluster.col = rainbow(length(unique(vbh$cluster)))
names(cluster.col) = unique(vbh$cluster)
cluster.col


col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#0a3b70", "white", "#760521"))

png("VBH/VBH.1_3.CNV.celltype.png", width = 10, height = 5, units = "in", res = 150)
Heatmap(t(cnv),cluster_columns = F,
        left_annotation = rowAnnotation(ct = vbh[colnames(cnv),4], cluster = vbh[colnames(cnv),2],
                                        col = list(ct = cell.col, cluster = cluster.col)),
        row_dend_gp = gpar(lwd = 0.2),
        show_row_names = F,show_column_names = F,
        cluster_rows = cnv.cluster, show_row_dend = T, 
        column_split = pos$seqnames, cluster_column_slices = F,
        column_title_rot = 90, column_title_side = "bottom",
        border = TRUE, column_gap = unit(0.2, "mm"),
        name = "cnv", col = col_fun,
        use_raster = TRUE, raster_quality = 2,
        heatmap_legend_param = list(direction = "horizontal"))

dev.off()


findsplitK(cnv.cluster)
k.split = 54
cell.cluster.cut = cutree(cnv.cluster, k.split)
sort(table(cell.cluster.cut), decreasing = T)[1:4]
cell.cluster.cut[cell.cluster.cut>2] = 3
p = make_leftAnn_cnvplot(cnv, pos, paste0(sample.name, "/", sample.name,".1_3.CNVraw.clusterAnn.png"), cell.cluster = cnv.cluster, 
                         cell.cluster.cut)



cell.cluster.cut[cell.cluster.cut>2] = 1
#2.subtract CNV control
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 1, choice = "median")
cnv.cluster2 = make_correlation_heatmap(cnv2, paste0(sample.name, "/", sample.name, ".2_1.cnv2.correlation.png"), barcodes)
p = make_raw_cnv_plot(cnv2, pos, cell.cluster = cnv.cluster2, paste0(sample.name, "/", sample.name,".2_2.cnv2.png"))

ncol(cnv2)
findsplitK(cnv.cluster2, N = 500, k.init = 20, k.max = 2000)
findsplitKrev(cnv.cluster2, N = 1590, n.cut = 1, k.max = 2302)

kk = 1862
sort(table(cutree(cnv.cluster2, k = kk)), decreasing = T)[1:4]

cell.cluster.cut.2 = cutree(cnv.cluster2, k = kk)
cell.cluster.cut.2[!cell.cluster.cut.2 %in% c(1,3,10)] = 4
table(cell.cluster.cut.2)
p = make_leftAnn_cnvplot(cnv2, pos, paste0(sample.name, "/", sample.name,".2_3.cnv2.clusterAnn.png"), cnv.cluster2, cell.cluster.cut.2)


png("VBH/VBH.2_3.CNV.celltype.png", width = 10, height = 5, units = "in", res = 150)
Heatmap(t(cnv2),cluster_columns = F,
        left_annotation = rowAnnotation(ct =vbh[colnames(cnv2),4], cluster = vbh[colnames(cnv2),2], 
                                        col = list(ct = cell.col, cluster = cluster.col) ),
        row_dend_gp = gpar(lwd = 0.2),
        show_row_names = F,show_column_names = F,
        cluster_rows = cnv.cluster2, show_row_dend = T, 
        column_split = pos$seqnames, cluster_column_slices = F,
        column_title_rot = 90, column_title_side = "bottom",
        border = TRUE, column_gap = unit(0.2, "mm"),
        name = "cnv", col = col_fun,
        use_raster = TRUE, raster_quality = 2,
        heatmap_legend_param = list(direction = "horizontal"))

dev.off()


#4.make final plot
p = make_final_cnv_plot(cnv2, pos, paste0(sample.name, "/", sample.name,".4.final.png"), cell.cluster.cut.2, k.malignant = 3)



candidate.malignant = names(which(cell.cluster.cut.2 == 3))
length(candidate.malignant)
write.table(candidate.malignant, paste0(sample.name, "/", sample.name,".4.CNVmalignant.txt"), quote = F, row.names = F, col.names = F)



###SAB
rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)

sample.name = "SAB"
cnv = readRDS(paste0(sample.name, ".CNVraw.rds"))
pos = readRDS(paste0(sample.name, ".CNVpos.rds"))
barcodes = readRDS(paste0(sample.name, ".barcodes.rds"))
head(cnv[,1:2])
colnames(cnv) = gsub(paste0(sample.name,"#"), "", colnames(cnv))
dim(cnv)
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dir.create(sample.name)
#1.raw correlation and CNV plot
#methylation: 1p loss, 7q,8q gain, 19 loss...
cnv.cluster = make_correlation_heatmap(cnv, paste0(sample.name,"/", sample.name, ".1_1.correlation.png"), barcodes)
p = make_raw_cnv_plot(cnv, pos, cell.cluster = cnv.cluster, png.file = paste0(sample.name, "/", sample.name, ".1_2.CNVraw.png"))
ncol(cnv)
findsplitK(cnv.cluster)
k.split = 29
cell.cluster.cut = cutree(cnv.cluster, k.split)
sort(table(cell.cluster.cut), decreasing = T)[1:4]
cell.cluster.cut[cell.cluster.cut>2] = 3
p = make_leftAnn_cnvplot(cnv, pos, paste0(sample.name, "/", sample.name,".1_3.CNVraw.clusterAnn.png"), cell.cluster = cnv.cluster, 
                         cell.cluster.cut)



cell.cluster.cut[cell.cluster.cut>2] = 1
#2.subtract CNV control
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 1, choice = "median")
cnv.cluster2 = make_correlation_heatmap(cnv2, paste0(sample.name, "/", sample.name, ".2_1.cnv2.correlation.png"), barcodes)
p = make_raw_cnv_plot(cnv2, pos, cell.cluster = cnv.cluster2, paste0(sample.name, "/", sample.name,".2_2.cnv2.png"))

ncol(cnv2)

findsplitKrev(cnv.cluster2, N = 2199, n.cut = 1, k.max = 1400)

kk = 1356
sort(table(cutree(cnv.cluster2, k = kk)), decreasing = T)[1:4]
cell.cluster.cut.2 = cutree(cnv.cluster2, k = kk)
cell.cluster.cut.2[!cell.cluster.cut.2 %in% c(2, 18)] = 1
table(cell.cluster.cut.2)
p = make_leftAnn_cnvplot(cnv2, pos, paste0(sample.name, "/", sample.name,".2_3.cnv2.clusterAnn.png"), cnv.cluster2, cell.cluster.cut.2)


#4.make final plot
p = make_final_cnv_plot(cnv2, pos, paste0(sample.name, "/", sample.name,".4.final.png"), cell.cluster.cut.2, k.malignant = 2)



candidate.malignant = names(which(cell.cluster.cut.2 == 2))
length(candidate.malignant)
write.table(candidate.malignant, paste0(sample.name, "/", sample.name,".4.CNVmalignant.txt"), quote = F, row.names = F, col.names = F)



### O IDH, 1p19q co-del
###VBH

rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)

sample.name = "VBH"
cnv = readRDS(paste0(sample.name, ".CNVraw.rds"))
pos = readRDS(paste0(sample.name, ".CNVpos.rds"))
barcodes = readRDS(paste0(sample.name, ".barcodes.rds"))
head(cnv[,1:2])
colnames(cnv) = gsub(paste0(sample.name,"#"), "", colnames(cnv))
dim(cnv)
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dir.create(sample.name)
#1.raw correlation and CNV plot
cnv.cluster = make_correlation_heatmap(cnv, paste0(sample.name,"/", sample.name, ".1_1.correlation.png"), barcodes)
p = make_raw_cnv_plot(cnv, pos, cell.cluster = cnv.cluster, png.file = paste0(sample.name, "/", sample.name, ".1_2.CNVraw.png"))
ncol(cnv)
findsplitK(cnv.cluster)
k.split = 54
cell.cluster.cut = cutree(cnv.cluster, k.split)
sort(table(cell.cluster.cut), decreasing = T)[1:4]
cell.cluster.cut[cell.cluster.cut>2] = 3
p = make_leftAnn_cnvplot(cnv, pos, paste0(sample.name, "/", sample.name,".1_3.CNVraw.clusterAnn.png"), cell.cluster = cnv.cluster, 
                         cell.cluster.cut)



cell.cluster.cut[cell.cluster.cut>2] = 1
#2.subtract CNV control
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 1, choice = "median")
cnv.cluster2 = make_correlation_heatmap(cnv2, paste0(sample.name, "/", sample.name, ".2_1.cnv2.correlation.png"), barcodes)
p = make_raw_cnv_plot(cnv2, pos, cell.cluster = cnv.cluster2, paste0(sample.name, "/", sample.name,".2_2.cnv2.png"))

ncol(cnv2)
findsplitK(cnv.cluster2, N = 500, k.init = 20, k.max = 2000)
findsplitKrev(cnv.cluster2, N = 1590, n.cut = 1, k.max = 2302)

kk = 1862
sort(table(cutree(cnv.cluster2, k = kk)), decreasing = T)[1:4]

cell.cluster.cut.2 = cutree(cnv.cluster2, k = kk)
cell.cluster.cut.2[!cell.cluster.cut.2 %in% c(1,3,10)] = 4
table(cell.cluster.cut.2)
p = make_leftAnn_cnvplot(cnv2, pos, paste0(sample.name, "/", sample.name,".2_3.cnv2.clusterAnn.png"), cnv.cluster2, cell.cluster.cut.2)


#4.make final plot
p = make_final_cnv_plot(cnv2, pos, paste0(sample.name, "/", sample.name,".4.final.png"), cell.cluster.cut.2, k.malignant = 3)



candidate.malignant = names(which(cell.cluster.cut.2 == 3))
length(candidate.malignant)
write.table(candidate.malignant, paste0(sample.name, "/", sample.name,".4.CNVmalignant.txt"), quote = F, row.names = F, col.names = F)


###GPU

rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)

sample.name = "GPU"
cnv = readRDS(paste0(sample.name, ".CNVraw.rds"))
pos = readRDS(paste0(sample.name, ".CNVpos.rds"))
barcodes = readRDS(paste0(sample.name, ".barcodes.rds"))
head(cnv[,1:2])
colnames(cnv) = gsub(paste0(sample.name,"#"), "", colnames(cnv))
dim(cnv)
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dir.create(sample.name)
#1.raw correlation and CNV plot
cnv.cluster = make_correlation_heatmap(cnv, paste0(sample.name,"/", sample.name, ".1_1.correlation.png"), barcodes)
p = make_raw_cnv_plot(cnv, pos, cell.cluster = cnv.cluster, png.file = paste0(sample.name, "/", sample.name, ".1_2.CNVraw.png"))

ncol(cnv)
findsplitK(cnv.cluster, N = 150)
k.split = 11
cell.cluster.cut = cutree(cnv.cluster, k.split)
sort(table(cell.cluster.cut), decreasing = T)[1:4]
cell.cluster.cut[!cell.cluster.cut %in% c(1,2)] = 3
p = make_leftAnn_cnvplot(cnv, pos, paste0(sample.name, "/", sample.name,".1_3.CNVraw.clusterAnn.png"), cell.cluster = cnv.cluster, 
                         cell.cluster.cut)



cell.cluster.cut[cell.cluster.cut>2] = 2
#2.subtract CNV control
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 2, choice = "median")
cnv.cluster2 = make_correlation_heatmap(cnv2, paste0(sample.name, "/", sample.name, ".2_1.cnv2.correlation.png"), barcodes)
p = make_raw_cnv_plot(cnv2, pos, cell.cluster = cnv.cluster2, paste0(sample.name, "/", sample.name,".2_2.cnv2.png"))

ncol(cnv2)
findsplitK(cnv.cluster2, N = 150)
findsplitKrev(cnv.cluster2, N = 1600, n.cut = 1, k.max = 1800)

kk = 193
sort(table(cutree(cnv.cluster2, k = kk)), decreasing = T)[1:4]
cell.cluster.cut.2 = cutree(cnv.cluster2, k = kk)
cell.cluster.cut.2[cell.cluster.cut.2 >1] =2
table(cell.cluster.cut.2)
p = make_leftAnn_cnvplot(cnv2, pos, paste0(sample.name, "/", sample.name,".2_3.cnv2.clusterAnn.png"), cnv.cluster2, cell.cluster.cut.2)


#4.make final plot
p = make_final_cnv_plot(cnv2, pos, paste0(sample.name, "/", sample.name,".4.final.png"), cell.cluster.cut.2, k.malignant = 1)



candidate.malignant = names(which(cell.cluster.cut.2 == 1))
length(candidate.malignant)
write.table(candidate.malignant, paste0(sample.name, "/", sample.name,".4.CNVmalignant.txt"), quote = F, row.names = F, col.names = F)








#######################################################################################################################################


# Diaz data

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

# Diaz data
#SF11215	Male	46	EGAD00001005407	GBM
#SF11331	Male	55	EGAD00001005408	GBM
#SF11956	Male	63	EGAD00001005406	GBM IDHwt
#SF11979	Female	76	EGAD00001005418	GBM
#SF12017	Male	55	EGAD00001005405	GBM IDHmut

rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)
#SF11215 GBM
id = "SF11215"
cnv = readRDS(paste0("diaz/", id, ".CNVraw.rds"))
pos = readRDS(paste0("diaz/", id, ".CNVpos.rds"))
barcodes = readRDS(paste0("diaz/", id, ".barcodes.rds"))
head(cnv[,1:2])
colnames(cnv) = gsub(paste0(id, "#"), "", colnames(cnv))
dim(cnv)
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dim(cnv)
dir.create(paste0("diaz/", id))
#1.raw correlation and CNV plot
cnv.cluster = make_correlation_heatmap(cnv, paste0("diaz/",id, "/", id, ".1_1.correlation.png"), barcodes, plot.dend = F)
p = make_raw_cnv_plot(cnv, pos, cnv.cluster, png.file = paste0("diaz/",id, "/", id, ".1_2.CNVraw.png"))

ncol(cnv)
findsplitK(cnv.cluster, N = 50)

cell.cluster.cut = cutree(cnv.cluster, k = 3)
table(cell.cluster.cut)
cell.cluster.cut[cell.cluster.cut>2] = 3
p = make_leftAnn_cnvplot(cnv, pos, png.file = paste0("diaz/",id, "/", id, ".1_3.CNVraw.clusterAnn.png"), 
                         cnv.cluster, cell.cluster.cut )


#subtract control
cell.cluster.cut[cell.cluster.cut>2] = 2
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 2, choice = "median")
cnv.cluster2 = make_correlation_heatmap(cnv2, paste0("diaz/",id, "/", id, ".2_1.cnv2.cor.png"), barcodes)
p = make_raw_cnv_plot(cnv2, pos, cnv.cluster2, png.file = paste0("diaz/",id, "/", id, ".2_2.cnv2.png"))


findsplitKrev(cnv.cluster2, N = 109, n.cut = 1, k.max = 108)
k.split = 98
sort(table(cutree(cnv.cluster2, k = k.split)), decreasing = T)[1:4]
cell.cluster.cut.2 = cutree(cnv.cluster2, k = k.split)
cell.cluster.cut.2[cell.cluster.cut.2>2] = 3

p = make_leftAnn_cnvplot(cnv2, pos, paste0("diaz/",id, "/", id,".2_3.cnv2.clusterAnn.png"),
                         cnv.cluster2, cell.cluster.cut.2)


cell.cluster.cut.2[cell.cluster.cut.2>2] = 2
p = make_final_cnv_plot(cnv2, pos, paste0("diaz/", id, "/", id, ".4.final.png"), 
                        cell.cluster.cut.2, k.malignant = 1)


candidate.malignant = names(which(cell.cluster.cut == 1))
length(candidate.malignant)
write.table(candidate.malignant, paste0("diaz/", id, "/", id, ".4.CNVmalignant.txt"), quote = F, row.names = F, col.names = F)

#######################################################################################################################################

rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)
#SF11331 GBM
id = "SF11331"
cnv = readRDS(paste0("diaz/", id, ".CNVraw.rds"))
pos = readRDS(paste0("diaz/", id, ".CNVpos.rds"))
barcodes = readRDS(paste0("diaz/", id, ".barcodes.rds"))
head(cnv[,1:2])
colnames(cnv) = gsub(paste0(id, "#"), "", colnames(cnv))
dim(cnv)
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dir.create(paste0("diaz/", id))
#1.raw correlation and CNV plot
cnv.cluster = make_correlation_heatmap(cnv, paste0("diaz/",id, "/", id, ".1_1.correlation.png"), barcodes, plot.dend = F)
p = make_raw_cnv_plot(cnv, pos, cnv.cluster, png.file = paste0("diaz/",id, "/", id, ".1_2.CNVraw.png"))
ncol(cnv)

#only 13 few cells

#######################################################################################################################################

rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)
#SF11956 GBM
id = "SF11956"
cnv = readRDS(paste0("diaz/", id, ".CNVraw.rds"))
pos = readRDS(paste0("diaz/", id, ".CNVpos.rds"))
barcodes = readRDS(paste0("diaz/", id, ".barcodes.rds"))
head(cnv[,1:2])
colnames(cnv) = gsub(paste0(id, "#"), "", colnames(cnv))
dim(cnv)
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dir.create(paste0("diaz/", id))
#1.raw correlation and CNV plot
cnv.cluster = make_correlation_heatmap(cnv, paste0("diaz/",id, "/", id, ".1_1.correlation.png"), barcodes)
p = make_raw_cnv_plot(cnv, pos, cnv.cluster, png.file = paste0("diaz/",id, "/", id, ".1_2.CNVraw.png"))

findsplitK(cnv.cluster, N = 10)
cell.cluster.cut = cutree(cnv.cluster, k = 4)
table(cell.cluster.cut)
cell.cluster.cut[cell.cluster.cut>3] = 3
p = make_leftAnn_cnvplot(cnv, pos, png.file = paste0("diaz/",id, "/", id, ".1_3.CNVraw.clusterAnn.png"), 
                         cnv.cluster, cell.cluster.cut)

cell.cluster.cut[cell.cluster.cut>2] = 2
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 2, choice = "median")
cnv.cluster2 = make_correlation_heatmap(cnv2, paste0("diaz/",id, "/", id, ".2_1.cnv2.cor.png"), barcodes)
p = make_raw_cnv_plot(cnv2, pos, cnv.cluster2, png.file = paste0("diaz/",id, "/", id, ".2_2.cnv2.png"))

ncol(cnv2)
findsplitKrev(cnv.cluster2, N = 1400, n.cut = 1, k.max = 200)
k.split = 116
sort(table(cutree(cnv.cluster2, k = k.split)), decreasing = T)[1:4]
cell.cluster.cut.2 = cutree(cnv.cluster2, k = k.split)
cell.cluster.cut.2[cell.cluster.cut.2>2] = 3

p = make_leftAnn_cnvplot(cnv2, pos, paste0("diaz/",id, "/", id,".2_3.cnv2.clusterAnn.png"),
                         cnv.cluster2, cell.cluster.cut.2)


cell.cluster.cut.2[cell.cluster.cut.2>2] = 2
p = make_final_cnv_plot(cnv2, pos, paste0("diaz/", id, "/", id, ".4.final.png"), 
                        cell.cluster.cut.2, k.malignant = 1)

##maybe some cells with CNA were treated as non CNA cells

candidate.malignant = names(which(cell.cluster.cut.2 == 1))
length(candidate.malignant)
write.table(candidate.malignant, paste0("diaz/", id, "/", id, ".4.CNVmalignant.txt"), quote = F, row.names = F, col.names = F)

#######################################################################################################################################

rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)
#SF11979 GBM
id = "SF11979"

 message("NO cells left for this sample!!!!")
#######################################################################################################################################

rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)
#SF12017 GBM IDH mutant !!!!!!!!!!!!!!!
id = "SF12017"

#######################################################################################################################################















rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)
#SF11612 Oligodendroglioma recurrent
id = "SF11612"
cnv = readRDS(paste0("diaz/", id, ".CNVraw.rds"))
pos = readRDS(paste0("diaz/", id, ".CNVpos.rds"))
barcodes = readRDS(paste0("diaz/", id, ".barcodes.rds"))
head(cnv[,1:2])
colnames(cnv) = gsub(paste0(id, "#"), "", colnames(cnv))
dim(cnv)
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dir.create(paste0("diaz/", id))
#1.raw correlation and CNV plot
cnv.cluster = make_correlation_heatmap(cnv, paste0("diaz/",id, "/", id, ".1_1.correlation.png"), barcodes)
p = make_raw_cnv_plot(cnv, pos, cnv.cluster, png.file = paste0("diaz/",id, "/", id, ".1_2.CNVraw.png"))
findsplitK(cnv.cluster)
k.split = 46
cell.cluster.cut = cutree(cnv.cluster, k = k.split)
sort(table(cell.cluster.cut), decreasing = T)[1:4]
cell.cluster.cut[cell.cluster.cut>4] = 5
p = make_leftAnn_cnvplot(cnv, pos, png.file = paste0("diaz/",id, "/", id, ".1_3.CNVraw.clusterAnn.png"), 
                         cnv.cluster, cell.cluster.cut)

cell.cluster.cut[cell.cluster.cut>2] = 1
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 1, choice = "median")
cnv.cluster2 = make_correlation_heatmap(cnv2, paste0("diaz/",id, "/", id, ".2_1.cnv2.cor.png"), barcodes)
p = make_raw_cnv_plot(cnv2, pos, cnv.cluster2, png.file = paste0("diaz/",id, "/", id, ".2_2.cnv2.png"))

ncol(cnv2)
findsplitKrev(cnv.cluster2, N = 1800, n.cut = 1, k.max = 325)
k.split = 70
sort(table(cutree(cnv.cluster2, k = k.split)), decreasing = T)[1:4]
cell.cluster.cut.2 = cutree(cnv.cluster2, k = k.split)
cell.cluster.cut.2[cell.cluster.cut.2>2] = 3
table(cell.cluster.cut.2)
p = make_leftAnn_cnvplot(cnv2, pos, paste0("diaz/",id, "/", id,".2_3.cnv2.clusterAnn.png"),
                         cnv.cluster2, cell.cluster.cut.2)


cell.cluster.cut.2[cell.cluster.cut.2 != 2] = 1
p = make_final_cnv_plot(cnv2, pos, paste0("diaz/", id, "/", id, ".4.final.png"), 
                        cell.cluster.cut.2, k.malignant = 2)

##some CNA/normal cells not clear

candidate.malignant = names(which(cell.cluster.cut.2 == 2))
length(candidate.malignant)
write.table(candidate.malignant, paste0("diaz/", id, "/", id, ".4.CNVmalignant.txt"), quote = F, row.names = F, col.names = F)

#######################################################################################################################################





rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)
#SF11949 O IDH
id = "SF11949"
cnv = readRDS(paste0("diaz/", id, ".CNVraw.rds"))
pos = readRDS(paste0("diaz/", id, ".CNVpos.rds"))
barcodes = readRDS(paste0("diaz/", id, ".barcodes.rds"))
head(cnv[,1:2])
colnames(cnv) = gsub(paste0(id, "#"), "", colnames(cnv))
dim(cnv)
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dir.create(paste0("diaz/", id))
#1.raw correlation and CNV plot
cnv.cluster = make_correlation_heatmap(cnv, paste0("diaz/",id, "/", id, ".1_1.correlation.png"), barcodes)
p = make_raw_cnv_plot(cnv, pos, cnv.cluster, png.file = paste0("diaz/",id, "/", id, ".1_2.CNVraw.png"))

ncol(cnv)
findsplitK(cnv.cluster, N = 10, n.cut = 3, k.init = 6)
k.split = 15
cell.cluster.cut = cutree(cnv.cluster, k = k.split)
sort(table(cell.cluster.cut), decreasing = T)[1:4]
cell.cluster.cut[cell.cluster.cut>2] = 3
p = make_leftAnn_cnvplot(cnv, pos, png.file = paste0("diaz/",id, "/", id, ".1_3.CNVraw.clusterAnn.png"), 
                         cnv.cluster, cell.cluster.cut)

cell.cluster.cut[cell.cluster.cut != 1] = 2
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 2, choice = "median")
cnv.cluster2 = make_correlation_heatmap(cnv2, paste0("diaz/",id, "/", id, ".2_1.cnv2.cor.png"), barcodes)
p = make_raw_cnv_plot(cnv2, pos, cnv.cluster2, png.file = paste0("diaz/",id, "/", id, ".2_2.cnv2.png"))

ncol(cnv2)
findsplitKrev(cnv.cluster2, N = 695, n.cut = 1, k.max = 62)
k.split = 27
sort(table(cutree(cnv.cluster2, k = k.split)), decreasing = T)[1:4]
cell.cluster.cut.2 = cutree(cnv.cluster2, k = k.split)
cell.cluster.cut.2[cell.cluster.cut.2>2] = 2
table(cell.cluster.cut.2)
p = make_leftAnn_cnvplot(cnv2, pos, paste0("diaz/",id, "/", id,".2_3.cnv2.clusterAnn.png"),
                         cnv.cluster2, cell.cluster.cut.2)


p = make_final_cnv_plot(cnv2, pos, paste0("diaz/", id, "/", id, ".4.final.png"), 
                        cell.cluster.cut.2, k.malignant = 1)


candidate.malignant = names(which(cell.cluster.cut.2 == 1))
length(candidate.malignant)
write.table(candidate.malignant, paste0("diaz/", id, "/", id, ".4.CNVmalignant.txt"), quote = F, row.names = F, col.names = F)

#######################################################################################################################################

rm(cnv, pos, cnv2, UC, cnv.cluster, cnv.cluster2, cell.cluster.cut, cell.cluster.cut.2, cn.upper, umi.lower, 
   umi.upper, k.split, p, barcodes, barcodes2, candidate.malignant)
#SF11964 A IDH
id = "SF11964"
cnv = readRDS(paste0("diaz/", id, ".CNVraw.rds"))
pos = readRDS(paste0("diaz/", id, ".CNVpos.rds"))
barcodes = readRDS(paste0("diaz/", id, ".barcodes.rds"))
head(cnv[,1:2])
colnames(cnv) = gsub(paste0(id, "#"), "", colnames(cnv))
dim(cnv)
barcodes = barcodes[match(colnames(cnv), barcodes$barcode),]
barcodes$promoter.ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1); 

plot(density(log10(barcodes$passed_filters)))
plot(density(as.vector(cnv)))
cnv[cnv < -3] = -3
cnv[cnv > 3] = 3

dir.create(paste0("diaz/", id))
#1.raw correlation and CNV plot
cnv.cluster = make_correlation_heatmap(cnv, paste0("diaz/",id, "/", id, ".1_1.correlation.png"), barcodes)
p = make_raw_cnv_plot(cnv, pos, cnv.cluster, png.file = paste0("diaz/",id, "/", id, ".1_2.CNVraw.png"))
ncol(cnv)
findsplitK(cnv.cluster, N = 10)
k.split = 6
cell.cluster.cut = cutree(cnv.cluster, k = k.split)
sort(table(cell.cluster.cut), decreasing = T)[1:4]
cell.cluster.cut[cell.cluster.cut>2] = 2
p = make_leftAnn_cnvplot(cnv, pos, png.file = paste0("diaz/",id, "/", id, ".1_3.CNVraw.clusterAnn.png"), 
                         cnv.cluster, cell.cluster.cut)

table(cell.cluster.cut)
cnv2 = subtractControl(cnv, cell.cluster.cut, k.control = 2, choice = "median")
cnv.cluster2 = make_correlation_heatmap(cnv2, paste0("diaz/",id, "/", id, ".2_1.cnv2.cor.png"), barcodes)
p = make_raw_cnv_plot(cnv2, pos, cnv.cluster2, png.file = paste0("diaz/",id, "/", id, ".2_2.cnv2.png"))


ncol(cnv2)
findsplitKrev(cnv.cluster2, N = 3120, n.cut = 1, k.max = 100)
k.split = 54
sort(table(cutree(cnv.cluster2, k = k.split)), decreasing = T)[1:4]
cell.cluster.cut.2 = cutree(cnv.cluster2, k = k.split)
cell.cluster.cut.2[cell.cluster.cut.2>2] = 2
table(cell.cluster.cut.2)
p = make_leftAnn_cnvplot(cnv2, pos, paste0("diaz/",id, "/", id,".2_3.cnv2.clusterAnn.png"),
                         cnv.cluster2, cell.cluster.cut.2)


p = make_final_cnv_plot(cnv2, pos, paste0("diaz/", id, "/", id, ".4.final.png"), 
                        cell.cluster.cut.2, k.malignant = 1)


candidate.malignant = names(which(cell.cluster.cut.2 == 1))
length(candidate.malignant)
write.table(candidate.malignant, paste0("diaz/", id, "/", id, ".4.CNVmalignant.txt"), quote = F, row.names = F, col.names = F)

#######################################################################################################################################


