library(ggplot2)
library(scalop)
library(Seurat)

#source("/Users/wuz6/Documents/Project/25.pHGG_spatial/01.MP2025HE/utils/utils-hca.R")
#rlang::env_unlock(env = asNamespace('scalop'))
#rlang::env_binding_unlock(env = asNamespace('scalop'))
#assign('.hca', .hca, envir = asNamespace('scalop'))
#rlang::env_binding_lock(env = asNamespace('scalop'))
#rlang::env_lock(asNamespace('scalop'))
#remotes::install_github("jlaffy/scalop")
setwd("/Users/wuz6/Documents/Project/25.pHGG_spatial/01.MP2025HE/")
cancer.col = paletteer::paletteer_d("ggthemes::Tableau_20")[c(1,3,5)]
cancer.col <- as.character(cancer.col)
names(cancer.col) = c("GBM", "pHGG", "IDHmut")

cancer.col2 = paletteer::paletteer_d("ggthemes::Tableau_20")[c(1,3,5,6)]
cancer.col2 <- as.character(cancer.col2)
names(cancer.col2) = c("GBM", "pHGG", "IDH-A", "IDH-O")

###score MP in spatial data
MP_list = readRDS(paste0("01.all_0-30-10_10-12-12-4/MP_list.rds"))

####
meta <- (read.delim("../00.data/00.data.folder.csv", header = T, sep = ",")) #$Sample
#meta <- meta[grep("PHE", meta$Sample, value = F, invert = T),] #exclude PHE #can include here
meta$Tumor3 <- ifelse(grepl("AIDH|A IDH", meta$Tumor2), "IDH-A", 
                            ifelse(grepl("OIDH|O IDH", meta$Tumor2), "IDH-O", meta$Tumor))
table(meta$Tumor3)

samples_names = meta$Sample
samples_names
set.seed(1234567)
length(samples_names)


mp.cell.GBM <- data.frame(readxl::read_xlsx("/Users/wuz6/Documents/Project/data/signatures/C.Greenwald2024Cell.xlsx"), stringsAsFactors = F)
mp.cell.GBM <- as.vector(mp.cell.GBM)
names(mp.cell.GBM) <- paste0("Greenwald_", names(mp.cell.GBM))
names(mp.cell.GBM)
mp.cell.GBM <- lapply(mp.cell.GBM, function(x){
  x[!is.na(x) & x != ""]
})


#nFeature_Spatial.cut.off <- 800; #0
nCount_Spatial.cut.off <- 1000
data.home <- "/Users/wuz6/Documents/Project/25.pHGG_spatial/00.data/"

names(MP_list) <- c("Vasc", "Oligo", "Neuron", "Mac", "MES.Hyp", "MES", "InflammatoryMac", "AC", "Cellcycle", "unknown1",
                             "Interferon", "Perivascular", "Neuron2", "Neuron3", "Prolif.Metab", "unknown2", "NPC", "OPC.NPC", "Reactive.Ast", "EMT") #pHGG+IDH+GBM
names(MP_list)

MP_list.named <- c(MP_list, mp.cell.GBM)
names(MP_list.named)
##
... <- lapply(samples_names, function(k){
    #Check if directory exists
   if(! file.exists(paste0(data.home, k))){
     message("Missing sample: ", k) 
   }
  1
})

MP.scores = lapply(samples_names, function(k){
  #k = "UKF243"
  message(k)
  spatial_loc_file = paste0(data.home, k, "/outs/spatial/tissue_positions.csv")
  if(!file.exists(spatial_loc_file)){
    spatial_loc_file = paste0(data.home, k, "/outs/spatial/tissue_positions_list.csv.gz")
    if(!file.exists(spatial_loc_file)) spatial_loc_file = paste0(data.home, k, "/outs/spatial/tissue_positions_list.csv")
    tissue_pos <- read.csv(spatial_loc_file, head =F)
    #barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres
    colnames(tissue_pos) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
  }else{
    tissue_pos <- read.csv(spatial_loc_file, head =T)
  }
  
  
  sp.data <- Load10X_Spatial(paste0(data.home, k, "/outs/"), 
                            filename = "filtered_feature_bc_matrix.h5", assay = "Spatial")
  sp.data@meta.data$orig.ident <- k
  sp.data <- PercentageFeatureSet(sp.data, pattern = "^MT-", col.name = "percent.mt")
  sp.data <- subset(sp.data, subset = nCount_Spatial > nCount_Spatial.cut.off & percent.mt < 20)# & nFeature_Spatial > nFeature_Spatial.cut.off )
  sp.data <- NormalizeData(sp.data, scale.factor = 1e6)
  ###
  n1 = ncol(sp.data@meta.data) + 1
  #sp.data = AddMetaData(sp.data, x)
  sig.score = sigScores(data.matrix(sp.data[['Spatial']]$data), MP_list.named, expr.center = TRUE, conserved.genes = 0.5)
  ###
  #xy = sp.data@images$slice1@coordinates
  #xy <- GetTissueCoordinates(sp.data)[,1:2] #pixel coordinates
  
  arrayrow_col <- tissue_pos[match(rownames(sig.score), tissue_pos$barcode), c("array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")]
  cbind(rownames(sig.score), sp.data@meta.data, arrayrow_col, sig.score)
})

MP.scores = do.call(rbind, MP.scores)
head(MP.scores, 2)
colnames(MP.scores)[1:2] <- c("barcode", "sample")
nrow(MP.scores)


library(corrplot)
heat.col = as.character(paletteer::paletteer_c("grDevices::Red-Green", 30))
colnames(MP.scores)

MP.scores.corr = cor(MP.scores[,-c(1:9)])


MP.scores.corr.mean = lapply(unique(MP.scores$sample), function(ss){
  y = cor(MP.scores[MP.scores$sample == ss,-c(1:9)])
  y
})
MP.scores.corr.mean <- purrr::reduce(.x = MP.scores.corr.mean,.f = `+`)/length(MP.scores.corr.mean)
dim(MP.scores.corr.mean)
MP.scores.corr.mean[1:4,1:4]


#pdf(paste0("01.", study.tumor,"_0-30-10_10-12-12-4/01.rawMPscore.correlation.pdf"), width = 7, height = 7, useDingbats = F)
pdf(paste0("01.all_0-30-10_10-12-12-4/01.rawMPscore.correlation.pdf"), width = 14, height = 14, useDingbats = F)
corrplot(MP.scores.corr, method = "color", order = "hclust", title = "correlation across all spots",
         col = heat.col, tl.col = "black") #addrect = 8, 
corrplot(MP.scores.corr.mean, method = "color", order = "hclust", title = "mean correlation across samples",
         col = heat.col,tl.col = "black")
dev.off()

colnames(MP.scores)[duplicated(colnames(MP.scores))]

#non_cell_cycle_names <- names(MP_list.named)[-7]
non_cell_cycle_names <- names(MP_list)
non_cell_cycle_names
spot.max.M = apply(MP.scores[,non_cell_cycle_names], 1, function(x)  non_cell_cycle_names[which.max(x)])
#spot.second.M = apply(mp.scores, 1, function(x) colnames(mp.scores)[order(x, decreasing = T)[2]])
#spot.max.M = data.frame(MP = spot.max.M, MP2 =spot.second.M, stringsAsFactors = T)
MP.scores$MPclass = spot.max.M
MP.scores$cancer = meta$Tumor[match(MP.scores$sample, meta$Sample)]
table(MP.scores$cancer, spot.max.M)

sample.mp.freq = data.frame(table(MP.scores$sample, MP.scores$MPclass)/rowSums(table(MP.scores$sample, MP.scores$MPclass)))
head(sample.mp.freq)
x = tidyr::pivot_wider(sample.mp.freq, names_from = "Var2",values_from = "Freq")


#

manualcolors1<-c('black','forestgreen', 'red2', 'orange', 'cornflowerblue', 
                 'magenta', 'darkolivegreen4', 'indianred1', 'tan4', 'darkblue', 
                 'mediumorchid1','firebrick4',  'yellowgreen', 'lightsalmon', 'tan3',
                 "tan1",'darkgray', 'wheat4', '#DDAD4B', 'chartreuse', 
                 'seagreen1', 'moccasin', 'mediumvioletred', 'seagreen','cadetblue1',
                 "darkolivegreen1" ,"tan2" ,   "tomato3" , "#7CE3D8","gainsboro")

#tableau_colors
manualcolors <- c("#ffffff",
                  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B",
                  "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF", "#9C4D5C", "#2A7F3A",
                  "#F5B8B3", "#B3B8F0", "#D4E1F4", "#F0B2B8", "#9CBBE8", "#6B8D6F",
                  "#E1A43D", "#8C8B9C", "#D7A9B5", "#6D8DA7", "#B9C4B8", "#C5A3B6",
                  "#3C8D5B", "#E4B735", "#B37D9C", "#A7B9E5", "#F0C3B8", "#6B5F5C"
)


mp.col = manualcolors[2:30]
names(mp.col) = names(MP_list)

x <- data.frame(x)
rownames(x) <- x[,1]

hc <- hclust(dist(x[,-1]))

p = ggplot(sample.mp.freq, aes(x = Var1, y = Freq, fill = Var2))+
  geom_bar(stat = "identity")+
  geom_tile(inherit.aes = F, data = meta, aes(x = Sample, y = -0.05, fill = Tumor, height = 0.08), show.legend = T)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5) )+
  scale_fill_manual(values = c(mp.col, cancer.col))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(limits = hc$labels[hc$order])+
  labs(x = "", fill = "")
p 

pdf(paste0("01.all_0-30-10_10-12-12-4/01.sample.MPfraction.bar.pdf"), width = 12, height = 6, useDingbats = F)
p
dev.off()



if(1){
  purity1 <- readRDS("01.purity/01.inhouse.ESTIMATE.scores.rds")
  purity1 <- do.call(rbind, purity1)
  purity1$sampleName[purity1$sampleName == "Y963B11"] <- "Y963_B11" 
  purity1$sampleName[purity1$sampleName == "Y963B8"] <- "Y963_B8"
  purity2 <- readRDS("01.purity/01.Ravi2022CancerCell.ESTIMATE.scores.rds")
  purity2 <- do.call(rbind, purity2)
  purity3 <- readRDS("01.purity/01.Alissa.Greenwald.2024.GBM.Cell.ESTIMATE.scores.rds")
  purity3 <- do.call(rbind, purity3)
  purity3$sampleName <- sapply(purity3$sampleName, function(x) strsplit(x, "\\_")[[1]][2])
  purity <- rbind(purity1, purity2, purity3)
  rm(purity1, purity2, purity3)
}
purity$barcode1 <- paste0(purity$sampleName, "_", purity$sample)
unique(purity$sampleName)
unique(MP.scores$sample[!MP.scores$sample %in% purity$sampleName])

MP.scores$barcode1 <- paste0(MP.scores$sample, "_", MP.scores$barcode)
MP.scores$purity <- purity$purity[match(MP.scores$barcode1, purity$barcode1)]
plot(density(MP.scores$purity, na.rm = T))
###purity
purity.mean <- aggregate(purity ~ MPclass, data = MP.scores, FUN = mean)

p <- ggplot(MP.scores, aes(x = MPclass, y = purity, fill = MPclass)) +
  geom_violin(show.legend = F)+
  #geom_jitter(width = 0.2, size = .1, color = "gray")+
  geom_boxplot(width = 0.2, outlier.shape = NA)+
  scale_fill_manual(values = mp.col, name = "")+
  scale_x_discrete(limits = purity.mean$MPclass[order(purity.mean$purity)])+
  #scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(x = "")
p

pdf("01.MP-purity.violin.pdf", width = 8, height = 6, useDingbats = F)
p
dev.off()

if(1){
  write.csv(MP.scores, file = paste0("01.all_0-30-10_10-12-12-4/01.MPscores.csv"), quote = F, row.names = F)
}


### Plot spatial MP class for each sample
use.dir <- paste0("01.all_0-30-10_10-12-12-4/01.initial.mp.assignment.plot/")
dir.create(use.dir)

... <- lapply(1:nrow(meta), function(i){
     #i = 2
     k <- meta$Sample[i]
     if(!grepl("zh", k)){
       return(1)
     } 
     data.lab <- meta$Source[i]
     use.pt.size.factor = 3
     #use.pt.size.factor = 2
     if(data.lab == "Aldape Lab") use.pt.size.factor <- 2
     message(k)
     
     mp.metadata = MP.scores[MP.scores$sample == k, c("barcode", "MPclass")]
     rownames(mp.metadata) <- mp.metadata$barcode
     mp.metadata <- mp.metadata[,"MPclass", drop = F]
     
     sp.data <- Load10X_Spatial(paste0(data.home, k, "/outs/"), 
                                filename = "filtered_feature_bc_matrix.h5", assay = "Spatial")
     #x <- GetTissueCoordinates(sp.data)
     #message("Range x: ", paste0(range(x[,1]), collapse = "-"), "; Range y: ", paste0(range(x[,2]), collapse = "-"))
     #only use spot in mp.metadata
     sp.data <- sp.data[,rownames(mp.metadata)]
     sp.data <- AddMetaData(sp.data, mp.metadata)
     
     p <- SpatialDimPlot(sp.data, group.by = "MPclass", image.alpha = 0.4, image.scale = "lowres", pt.size.factor = 3) &
       scale_fill_manual(values = mp.col, name = "") &
       labs(title = k)
     #p
      xx <- MP.scores[MP.scores$sample ==k,]
      if(grepl("zh", k)){
        p2 = ggplot(xx, aes(x = -array_row, y = array_col, col = MPclass))+
          geom_point(size=2) + 
          labs(title= k, col = "") +
          scale_color_manual(values = mp.col, name = "")+
          theme_void()+
          theme(aspect.ratio = 1)+
          guides(colour = guide_legend(override.aes = list(size=5)))
        
      }else{
        p2 = ggplot(xx, aes(x = array_col, y = -array_row, col = MPclass))+
          geom_point(size=2) + 
          labs(title= k, col = "") +
          scale_color_manual(values = mp.col, name = "")+
          theme_void()+
          theme(aspect.ratio = 1)+
          guides(colour = guide_legend(override.aes = list(size=5)))
      }
    
    #p2
    pdf(paste0(use.dir, "/01.",k, ".GBM-MP.pdf"), width = 14, height = 7, useDingbats = F)
    print(cowplot::plot_grid(p, p2))
    dev.off()
    1;
})


##second score?
#MP.scores = read.csv(file = paste0("01.", study.tumor, "_0-30-10_10-12-12-4/01.MPscores.csv"), check.names = F)

spot.score.max <- apply(MP.scores[,non_cell_cycle_names], 1, max)
spot.score.max.name <- apply(MP.scores[,non_cell_cycle_names], 1, function(x) non_cell_cycle_names[which.max(x)])
spot.score.second <- apply(MP.scores[,non_cell_cycle_names], 1, function(x) sort(x, decreasing = T)[2] )
spot.score.second.names <- apply(MP.scores[,non_cell_cycle_names], 1, function(x) non_cell_cycle_names[order(x, decreasing = T)[2]] )

table(spot.score.max.name, MP.scores$sample)

quantile(spot.score.max)
plot(density(spot.score.max, na.rm = T), main = "MP max score")
lines(density(spot.score.second, na.rm = T), col = "blue")

plot(density(spot.score.max-spot.score.second, na.rm = T), main = "MP max score")
###plot MP score heatmap
m.score <- lapply(non_cell_cycle_names, function(mp){
  #x <- MP.scores[MP.scores$MPclass == mp, names(MP_list.named)]
  x <- MP.scores[MP.scores$MPclass == mp,]
  x <- x[order(x[[mp]], decreasing = T),]
  x
})
head(m.score[[1]])
m.score <- do.call(rbind, m.score)

reds30 <- rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 30))

min(m.score[,unique(MP.scores$MPclass)])
max(m.score[,unique(MP.scores$MPclass)])
colnames(MP.scores)
#png(paste0("01.", study.tumor, "_0-30-10_10-12-12-4/01.MP.score.heatmap.png"), width = 10, height = 5, unit = "in", res = 300)
#png(paste0("01.", study.tumor, "_0-30-10_10-12-12-4/01.GBM-MP.score.heatmap.png"), width = 10, height = 5, unit = "in", res = 300)
#Heatmap(t(m.score[,c(non_cell_cycle_names, "Cellcycle")]), name = " ",
h <- Heatmap(t(m.score[,c(non_cell_cycle_names)]), name = " ",         
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F, col = circlize::colorRamp2(c(-1, 0, 1), c("#1F78B4", "white", "red")),
        #width = ncol(m.score)*unit(2, "mm"), height = nrow(m.score)*unit(2, "mm"),
        raster_device = "jpeg", use_raster = TRUE, raster_resize_mat = TRUE, raster_quality = 10,
        row_title = " ", column_title = " ")
png("01.all_0-30-10_10-12-12-4/01.MP.score.heatmap.png", width = 10, height = 5, unit = "in", res = 300)
draw(h)
dev.off()






