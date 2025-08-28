setwd("/Users/wuz6/Documents/Project/25.pHGG_spatial/01.MP_pHGG_GBM/")
###compare pHGG, GBM, IDH sampling resulted MPs

mp1 <- read.csv("00.GBM.sampleSize13.seed123.MP.csv", head =T)
colnames(mp1) <- paste0("GBM.s1:",1:ncol(mp1))
mp2 <- read.csv("00.GBM.sampleSize13.seed12345.MP.csv", head =T)
colnames(mp2) <- paste0("GBM.s2:",1:ncol(mp2))
mp3 <- read.csv("00.GBM.sampleSize13.seed1234567.MP.csv", head =T)
colnames(mp3) <- paste0("GBM.s3:",1:ncol(mp3))

mp4 <- read.csv("00.pHGG.sampleSize13.seed123.MP.csv", head =T)
colnames(mp4) <- paste0("pHGG.s1:",1:ncol(mp4))
mp5 <- read.csv("00.pHGG.sampleSize13.seed12345.MP.csv", head =T)
colnames(mp5) <- paste0("pHGG.s2:",1:ncol(mp5))
mp6 <- read.csv("00.pHGG.sampleSize13.seed1234567.MP.csv", head =T)
colnames(mp6) <- paste0("pHGG.s3:",1:ncol(mp6))

mp7 <- read.csv("01.IDH_0-25-12_14-14-4/01.MP_initial.csv", head =T)
colnames(mp7) <- paste0("IDHmut:",1:ncol(mp7))

combined.mp <- cbind(mp1, mp2, mp3,mp4,mp5,mp6,mp7)
colnames(combined.mp)

# calculate similarity between programs
nmf_intersect <- apply(combined.mp , 2, function(x) apply(combined.mp , 2, function(y) length(intersect(x,y)))) 

# hierarchical clustering of the similarity matrix 
nmf_intersect_hc_all <- hclust(as.dist(max(nmf_intersect)-nmf_intersect), method="average") 
nmf_intersect_hc_all <- reorder(as.dendrogram(nmf_intersect_hc_all), colMeans(nmf_intersect))
nmf_intersect         <- nmf_intersect[order.dendrogram(nmf_intersect_hc_all), order.dendrogram(nmf_intersect_hc_all)]


pheatmap::pheatmap(nmf_intersect)

library(ComplexHeatmap)
pair.col <- RColorBrewer::brewer.pal(12, "Paired")
RColorBrewer::display.brewer.pal(12, "Paired")

cancer.col <- pair.col[c(2,4,9)]
names(cancer.col) <- c("pHGG", "GBM", "IDHmut")
Heatmap(nmf_intersect)

cancer.type <- sapply(colnames(nmf_intersect), function(x) unlist(strsplit(x, "[\\.:]"))[1])
table(cancer.type)

row.ann = rowAnnotation(cancer = cancer.type, col = list(cancer = cancer.col))
col.ann = HeatmapAnnotation(cancer = cancer.type, col = list(cancer = cancer.col))

h = Heatmap(nmf_intersect, name = "Jaccard\nIndex", cluster_rows = F, cluster_columns = F, 
            right_annotation = row.ann, bottom_annotation = col.ann,
            show_row_names = T, show_column_names = T, col = custom_magma,
            #top_annotation = ha, 
            border_gp = gpar(col = "black", lty = 1),
            width = ncol(nmf_intersect)*unit(4, "mm"), 
            height = nrow(nmf_intersect)*unit(4, "mm"), raster_device = "jpeg",
            use_raster = TRUE, raster_resize_mat = TRUE)

pdf("00.IDH-GBM-pHGG.sampling13.MP.pdf", width = 15, height = 15, useDingbats = F)
draw(h)
dev.off()




library(biomaRt)  

mart = useMart("ensembl")
mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

genes <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name","transcript_biotype"), 
                        filters = c("transcript_biotype","chromosome_name"),
                        values = list("protein_coding",mart = mart))
####                        





####final MPs
mp.GBM <- readRDS("01.GBM_0-25-12_14-14-4/MP_list.rds")
names(mp.GBM) <- paste0("GBM-", 1:length(mp.GBM))
mp.pHGG <- readRDS("01.pHGG_0-25-12_14-14-4/MP_list.rds")
names(mp.pHGG) <- paste0("pHGG-", 1:length(mp.pHGG))

mp.both <- readRDS("MP_list.rds")
names(mp.both) <- paste0("combined-", 1:length(mp.both))

combined.mp <- do.call(cbind, c(mp.GBM, mp.pHGG, mp.both))

# calculate similarity between programs
nmf_intersect <- apply(combined.mp , 2, function(x) apply(combined.mp , 2, function(y) length(intersect(x,y)))) 

# hierarchical clustering of the similarity matrix 
nmf_intersect_hc_all <- hclust(as.dist(max(nmf_intersect)-nmf_intersect), method="average") 
nmf_intersect_hc_all <- reorder(as.dendrogram(nmf_intersect_hc_all), colMeans(nmf_intersect))
nmf_intersect         <- nmf_intersect[order.dendrogram(nmf_intersect_hc_all), order.dendrogram(nmf_intersect_hc_all)]


heatmap(nmf_intersect)

library(ComplexHeatmap)

Heatmap(nmf_intersect)

h = Heatmap(nmf_intersect, name = "Jaccard\nIndex", cluster_rows = F, cluster_columns = F,
            show_row_names = T, show_column_names = T, col = custom_magma,
            #top_annotation = ha, 
            border_gp = gpar(col = "black", lty = 1),
            width = ncol(m1)*unit(0.5, "mm"), 
            height = nrow(m1)*unit(0.5, "mm"), raster_device = "jpeg",
            use_raster = TRUE, raster_resize_mat = TRUE)

pdf("01.MP.combined.vs.GBM.vs.pHGG.pdf", width = 15, height = 15, useDingbats = F)
draw(h)
dev.off()


############
# calculate similarity between programs
combined.mp <- do.call(cbind,mp.both)
nmf_intersect <- apply(combined.mp , 2, function(x) apply(combined.mp , 2, function(y) length(intersect(x,y))))

# hierarchical clustering of the similarity matrix 
nmf_intersect_hc_all <- hclust(as.dist(max(nmf_intersect)-nmf_intersect), method="average") 
nmf_intersect_hc_all <- reorder(as.dendrogram(nmf_intersect_hc_all), colMeans(nmf_intersect))
nmf_intersect         <- nmf_intersect[order.dendrogram(nmf_intersect_hc_all), order.dendrogram(nmf_intersect_hc_all)]


heatmap(nmf_intersect)

library(ComplexHeatmap)

Heatmap(nmf_intersect)

h = Heatmap(nmf_intersect, name = "Jaccard\nIndex", cluster_rows = F, cluster_columns = F,
            show_row_names = T, show_column_names = T, col = custom_magma,
            #top_annotation = ha, 
            border_gp = gpar(col = "black", lty = 1),
            width = ncol(m1)*unit(0.5, "mm"), 
            height = nrow(m1)*unit(0.5, "mm"), raster_device = "jpeg",
            use_raster = TRUE, raster_resize_mat = TRUE)

draw(h)
#pdf("01.MP.combined.vs.GBM.vs.pHGG.pdf", width = 15, height = 15, useDingbats = F)
draw(h)
dev.off()



###

mp.both <- readRDS("MP_list.rds")
names(mp.both) <- paste0("combined-", 1:length(mp.both))
mp.cell.GBM = data.frame(readxl::read_xlsx("/Users/wuz6/Documents/Project/data/signatures/C.Greenwald2024Cell.xlsx"), stringsAsFactors = F)
mp.cell.GBM = as.vector(mp.cell.GBM)
combined.mp <- do.call(cbind,c(mp.cell.GBM, mp.both))
nmf_intersect <- apply(combined.mp , 2, function(x) apply(combined.mp , 2, function(y) length(intersect(x,y)))) 

# hierarchical clustering of the similarity matrix 
nmf_intersect_hc_all <- hclust(as.dist(max(nmf_intersect)-nmf_intersect), method="average") 
nmf_intersect_hc_all <- reorder(as.dendrogram(nmf_intersect_hc_all), colMeans(nmf_intersect))
nmf_intersect         <- nmf_intersect[order.dendrogram(nmf_intersect_hc_all), order.dendrogram(nmf_intersect_hc_all)]


heatmap(nmf_intersect)

library(ComplexHeatmap)

Heatmap(nmf_intersect)

h = Heatmap(nmf_intersect, name = "Jaccard\nIndex", cluster_rows = F, cluster_columns = F,
            show_row_names = T, show_column_names = T, col = custom_magma,
            #top_annotation = ha, 
            border_gp = gpar(col = "black", lty = 1),
            width = ncol(m1)*unit(0.5, "mm"), 
            height = nrow(m1)*unit(0.5, "mm"), raster_device = "jpeg",
            use_raster = TRUE, raster_resize_mat = TRUE)
pdf("01.MP.combined.vs.cellGBM.pdf", width = 15, height = 15, useDingbats = F)
draw(h)
dev.off()
