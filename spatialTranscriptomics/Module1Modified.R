library(reticulate)
library(Seurat)
#BiocManager::install("jpmam1/scalop")
library(scalop)
library(leiden)
library(RColorBrewer)
library(scales)
library(reshape2)
library(scales)
library(NMF)
library(MetBrewer)
library(colorspace)
library(tibble)
library(dplyr)
library(data.table)
library(stringr)
library(readr)
library(Matrix)
library(bigmemory)
library(doMC)
library(patchwork)
library(ggplot2)
library(paletteer)
library(viridis)
library(ComplexHeatmap)
rm(list=ls())

source("/Users/wuz6/Documents/Project/25.pHGG_spatial/01.MP/utils/robust_nmf_programs.R")
mat2list <- function(x) {
  stopifnot(is.matrix(x))
  lapply(seq_len(ncol(x)), function(i) x[, i])
}

library(scalop)
source("/Users/wuz6/Documents/Project/25.pHGG_spatial/01.MP/utils/utils-hca.R")
rlang::env_unlock(env = asNamespace('scalop'))
rlang::env_binding_unlock(env = asNamespace('scalop'))
assign('.hca', .hca, envir = asNamespace('scalop'))
#assign('.hca_dist', .hca_dist, envir = asNamespace('scalop'))
rlang::env_binding_lock(env = asNamespace('scalop'))
rlang::env_lock(asNamespace('scalop'))

###
setwd("/Users/wuz6/Documents/Project/25.pHGG_spatial/01.MP_pHGG_GBM/")

# Custom color palette
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

###### Per sample Leiden clustering##############################


samples_names <- (read.delim("../00.data/00.data.folder.csv", header = T, sep = ",")) #$Sample
table(samples_names$Tumor)
samples_names = samples_names$Sample[samples_names$Tumor == "GBM"] #GBM
#samples_names = samples_names$Sample[samples_names$Tumor == "pHGG"] #pHGG
set.seed(1234567)
samples_names <- sample(samples_names, 13)

#samples_names <- samples_names$Sample[samples_names$Tumor == "IDHmut"]
samples_names
samples_names1 <- sapply(gsub("Y963_", "Y963", samples_names), function(x){stringr::str_split(x, pattern = "\\_")[[1]][[1]]})
length(samples_names)

sig_th <- .005
samples_names
# set parameters 
mp_num <- 14 ###???

#distinct16_pal<-c("#11A579","#F2B701","#66C5CC","#80BA5A","#F6CF71","#7F3C8D","#CF1C90","#3969AC","#f97b72","#E73F74","#4b4b8f","#ACA4E2","#31C53F","#B95FBB","#D4D915","#28CECA")
#dim2use_list <- c(8,20,19,9,12,14,16,10,10,14,8,10,12,16,7,5,14,16,4,12,6,6,8,12,13,10) # set by jackstraw

cluster.path = "/Users/wuz6/Documents/Project/25.pHGG_spatial/01.MP/01.banksyClusterMarkers/"
leiden_programs <- sapply(c(1:length(samples_names)), function(i){
  #i = 1
  samp_name <- stringr::str_split(samples_names[[i]], pattern = "\\_")[[1]][1]
  cat(samples_names[i]," ", samp_name, "\n")
  cluster_file = paste0(cluster.path, "/01.", samples_names[i], ".AllBanksyMarkers.csv.gz")
  #cluster.res = read.csv(paste0(cluster.path, "/01.", samples_names[i], ".AllBanksycluster.csv"), header = T, row.names = 1)
  #cluster.res = read.csv(paste0(cluster.path, "/01.", samples_names[i], ".Seuratcluster.csv"), header = T, row.names = 1)
  spatial_sample_programs = read.csv(cluster_file, header = T, row.names = 1)
  #spatial_sample_programs = read.csv(paste0(cluster.path, "/01.", samples_names[i], ".AllSeuratMarkers.csv.gz"), header = T, row.names = 1)
  spatial_sample_programs = spatial_sample_programs[spatial_sample_programs$p_val_adj< sig_th,]
  
  spatial_sample_programs$sample = samp_name
  spatial_sample_programs$res = gsub("[A-Za-z\\_]+", "", spatial_sample_programs$res)
  spatial_sample_programs$res = gsub("^\\.", "", spatial_sample_programs$res)
  spatial_sample_programs$resCluster = paste0(samp_name, "_", spatial_sample_programs$res, "_C", spatial_sample_programs$cluster)
  #table(spatial_sample_programs$resCluster)
  
  genes_num <- table(spatial_sample_programs$resCluster)
  keep.cluster = names(which(genes_num >= 60))
  spatial_sample_programs = spatial_sample_programs[spatial_sample_programs$resCluster %in% keep.cluster,]
  
  x = spatial_sample_programs[,c("gene", "resCluster", "avg_log2FC")]
  x = data.frame(tidyr::pivot_wider(x, names_from = resCluster, values_from = avg_log2FC, values_fill = 0))
  rownames(x) = x[,1]
  x = x[,-1]
  x
})
names(leiden_programs) = samples_names1 #name is required for robust_nmf_programs
#saveRDS(leiden_programs, file = "pHGG_GBM.rawLeiden.rds")


leiden_programs_sig <- lapply(leiden_programs, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:60]))
# for each sample, select robust NMF programs (i.e. observed using different ranks in the same sample), 
#remove redundancy due to multiple ranks, and apply a filter based on the similarity to programs from other samples. 
leiden_filter_all <- robust_nmf_programs(leiden_programs_sig, intra_min = 0, intra_max = 25, inter_filter=T, inter_min = 12)

length(leiden_filter_all) #

###
leiden_programs_sig <- lapply(leiden_programs_sig, function(x) x[, is.element(colnames(x), leiden_filter_all),drop=F])
leiden_programs_sig <- do.call(cbind, leiden_programs_sig)
dim(leiden_programs_sig) 
#saveRDS(leiden_programs_sig, file = "pHGG_GBM.robustLeiden-0-25-12.rds")


###### Generation of NMF + Leiden GBM spatial metaprograms#########

#

#robust_nmf_programs function
## Create list of NMF matrics where each sample is an entry
path <- "/Users/wuz6/Documents/Project/25.pHGG_spatial/01.MP/01.NMF/"
sample_ls <- list.files(path)
sample_ls
#sample_ls <- grep("UKF|GSM", sample_ls, value = T) #exlude pHGG
#sample_ls <- grep("UKF|GSM", sample_ls, value = T, invert = T) #pHGG
#sample_ls <- grep(paste0(c("Y963",samples_names1), collapse = "|"),sample_ls, value = T)
sample_ls
sample_ls <- grep(paste0(samples_names, collapse="|"), sample_ls, value = T)
length(sample_ls)
length(samples_names)
## Create list of NMF matrics where each sample is an entry
prog_genes_ls <- list()
for(i in seq_along(sample_ls)) {
  nmf_obj <- readRDS(paste(path, sample_ls[[i]], sep = "/"))
  samp_name <- stringr::str_split(gsub("Y963_", "Y963", sample_ls[[i]]), pattern = "[.]")[[1]][1]
  samp_name <- gsub("GSM\\d+\\_", "", samp_name)
  message(samp_name)
  nmf_mats <- c()
  for(j in names(nmf_obj$fit)) {
    get_basis <- basis(nmf_obj$fit[[j]])
    colnames(get_basis)  <- paste0(samp_name, ".", j, ".", 1:j)
    nmf_mats <- cbind(nmf_mats, get_basis)
  }
  prog_genes_ls[[i]] <- nmf_mats
  names(prog_genes_ls)[i] <- paste0(samp_name)
  rm(nmf_obj, nmf_mats, get_basis)
}
Genes_nmf_w_basis <- do.call(c, list(prog_genes_ls))

#saveRDS(Genes_nmf_w_basis, file = "pHGG_GBM_wBasis.rda")

background.genes = unique(unlist(lapply(prog_genes_ls, function(x) rownames(x))))
#saveRDS(background.genes, file = "background.genes.rds")
#length(background.genes) #

# Find robust NMFs
# get gene programs (top 50 genes by NMF score)
nmf_programs_sig <- lapply(prog_genes_ls, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:60]))
length(nmf_programs_sig) #
# for each sample, select robust NMF programs (i.e. observed using different ranks in the same sample), 
#remove redundancy due to multiple ranks, and apply a filter based on the similarity to programs from other samples. 

###intra_min: keep those detected multi-times
###intra_max: remove redundancy
###inter_min: keep those detected multi-samples
nmf_filter_all <- robust_nmf_programs(nmf_programs_sig, intra_min = 0, intra_max = 25, inter_filter=T, inter_min = 12)
#
nmf_programs_sig <- lapply(nmf_programs_sig, function(x) x[, is.element(colnames(x), nmf_filter_all),drop=F])
nmf_programs_sig <- do.call(cbind, nmf_programs_sig)
dim(nmf_programs_sig)
#nmf programs
#saveRDS(nmf_programs_sig, file = "GBM_all_nmf_robust.rds")
#nmf_programs_sig = readRDS("GBM_all_nmf_robust.rds")
#saveRDS(nmf_programs_sig, file = "pHGG_all_nmf_robust.rds")
#saveRDS(nmf_programs_sig, file = "pHGG_GBM_all_nmf_robust-35-10-10.rds")
#saveRDS(nmf_programs_sig, file = "pHGG_GBM_all_nmf_robust-0-25-12.rds")
#nmf_programs_sig = readRDS("pHGG_all_nmf_robust.rds")









#combine all signatures
#all_leiden_programs <- readRDS("all_leiden_programs.rds")
dim(leiden_programs_sig) #
nmf_programs_sig <- cbind(nmf_programs_sig, leiden_programs_sig)
dim(nmf_programs_sig) #

#nmf_programs_sig = readRDS("all_leiden_NMFrobust_programs.rds")
#nmf_programs_sig = readRDS("pHGGGBM_all_BanksyMulti_NMFrobust_programs.rds")


###Zhichao, exlcude programs with high number of MT- genes.
mt.cnt = apply(nmf_programs_sig, 2, function(x){
  length(grep("^MT", x))
})
pg = unique(as.character(nmf_programs_sig))
n.mt.pg = length(grep("^MT", pg))
pg = length(pg)
mt.phyper = phyper(mt.cnt -1 , n.mt.pg,  pg- n.mt.pg, 60, lower.tail = F)
sum(mt.phyper < 0.05) #
nmf_programs_sig = nmf_programs_sig[, -which(mt.phyper <0.05)]
dim(nmf_programs_sig) #

if(0){
  saveRDS(nmf_programs_sig, file = "pHGGGBM_all_BanksyMulti_NMFrobust-0-25-12_programs60.rds")
}






##########################################################################################
##########################################################################################
#nmf_programs_sig = readRDS(file = "pHGGGBM_all_BanksyMulti_NMFrobust_programs60.rds")
#leiden_programs = readRDS("pHGG_GBM.rawLeiden.rds")
#Genes_nmf_w_basis = readRDS("pHGG_GBM_wBasis.rda")


# calculate similarity between programs
nmf_intersect <- apply(nmf_programs_sig , 2, function(x) apply(nmf_programs_sig , 2, function(y) length(intersect(x,y)))) 

# hierarchical clustering of the similarity matrix 
nmf_intersect_hc_all <- hclust(as.dist(max(nmf_intersect)-nmf_intersect), method="average") 
nmf_intersect_hc_all <- reorder(as.dendrogram(nmf_intersect_hc_all), colMeans(nmf_intersect))
nmf_intersect         <- nmf_intersect[order.dendrogram(nmf_intersect_hc_all), order.dendrogram(nmf_intersect_hc_all)]

#nmf_intersect<-readRDS("MP/NMF/nmf_intersect_124.RDS")
#nmf_programs_sig<-readRDS("MP/NMF/nmf_programs_sig_124.RDS")


### use a clustering approach that updates MPs in each iteration 
nmf_intersect_KEEP    <- nmf_intersect
nmf_programs_sig_KEEP <- nmf_programs_sig




###
### Parameters (later change to function form)v1-keep!
Min_intersect_initial <- 14  # the minimal intersection cutoff for defining the Founder NMF program of a cluster
Min_intersect_cluster <- 14 # the minimal intersection cuttof for adding a new NMF to the forming cluster 
Min_group_size        <- 4   # the minimal group size to consider for defining the Founder_NMF of a MP  #This will determine the number of MPs
##Min_group_size = 4, 
##Min_group_size = 3, 

nmf_intersect = nmf_intersect_KEEP
nmf_programs_sig = nmf_programs_sig_KEEP
dim(nmf_programs_sig) #GBM 50x 317 #pHGG 50x270


#write.csv(nmf_intersect, file = "nmf_intersect.test.csv", quote = F)

##for each column, number of overlap programs (overlaped with > Min_intersect_initial )
Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)
plot(density(Sorted_intersection))
abline(v = 14, lty = 2, col = "blue");
quantile(Sorted_intersection)

Cluster_list <- list()   ### Every entry contains the NMFs of a chosec cluster
k <- 1
Curr_cluster <- c()
MP_list      <- list()

while (Sorted_intersection[1]>Min_group_size) {
  message("K", k)
  #
  Curr_cluster <- c(Curr_cluster , names(Sorted_intersection[1]))
  
  ### intersection between all remaining NMFs and Genes in MP 
  Genes_MP                   <- nmf_programs_sig[,names(Sorted_intersection[1])] # initial genes are those in the first NMF. Genes_MP always has only 50 genes consisting of the current MP
  nmf_programs_sig           <- nmf_programs_sig[,-match(names(Sorted_intersection[1]) , colnames(nmf_programs_sig))]  # remove selected NMF
  Intersection_with_Genes_MP <- sort(apply(nmf_programs_sig, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
  NMF_history                <- Genes_MP  # has all genes in all NMFs in the current cluster, for newly defining Genes_MP after adding a new NMF 
  
  #message(names(Intersection_with_Genes_MP[1:10]), Intersection_with_Genes_MP[1:10])
  ### Create gene list - composed of intersecting genes in descending order. Update Curr_cluster each time.
  
  while ( Intersection_with_Genes_MP[1] >= Min_intersect_cluster) {   ### Define current cluster 
    #message("Intersecting ", Curr_cluster, " with ", names(Intersection_with_Genes_MP)[1], ": overlap=", Intersection_with_Genes_MP[1])
    Curr_cluster  <- c(Curr_cluster , names(Intersection_with_Genes_MP)[1])
    
    Genes_MP_temp   <- sort(table(c(NMF_history , nmf_programs_sig[,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)   ## Genes_MP is newly defined each time according to all NMFs in the current cluster 
    Genes_at_border <- Genes_MP_temp[which(Genes_MP_temp == Genes_MP_temp[50])]   ### genes with overlap equal to the 50th gene
    
    if (length(Genes_at_border)>1){
      ### Sort last genes in Genes_at_border according to maximal NMF gene scores
      ### Run over all NMF programs in Curr_cluster and extract NMF scores for each gene
      Genes_curr_NMF_score <- c()
      ###prioritize NMF clusters
      Curr_cluster1 = grep("\\_C", Curr_cluster, value = T, invert = T) #these are NMF programs
      Curr_cluster2 = grep("\\_C", Curr_cluster, value = T) #some GSM NMF programs has _
      #Curr_cluster2 = grep("\\.", Curr_cluster2, value = T, invert = T)
      #print(Curr_cluster1)
      #print(Curr_cluster2)
      for (i in Curr_cluster1) {
        curr_study           <- strsplit(i, "[.]")[[1]][[1]]
        Q                    <- Genes_nmf_w_basis[[curr_study]][ match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]])))[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))]   ,i] 
        #names(Q)             <- names(Genes_at_border[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))])  ### sometimes when adding genes the names do not appear 
        Genes_curr_NMF_score <- c(Genes_curr_NMF_score,  Q )
      }
      Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score , decreasing = TRUE)
      Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]   ## take only the maximal score of each gene - which is the first entry after sorting
      ####find leiden/banksy cluster logFC rankes
      if(length(Curr_cluster1) == 0){#no NMF, so using leiden_clustering DEGs
        for (i in Curr_cluster2) {
          curr_study           <- strsplit(i, "[_]")[[1]][[1]]
          #MKG = leiden_clustering[1, match(curr_study, samples_names)][[1]]
          Q = leiden_programs[[curr_study]][,i]
          names(Q) = rownames(leiden_programs[[curr_study]])
          Genes_curr_NMF_score = c(Genes_curr_NMF_score, Q)
        }
      }
      Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score , decreasing = TRUE)
      Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]   ## take only the maximal score of each gene - which is the first entry after sorting
      Genes_MP_temp <- c(names(Genes_MP_temp[which(Genes_MP_temp > Genes_MP_temp[50])]) , names(Genes_curr_NMF_score_sort))
      
    } else {
      Genes_MP_temp <- names(Genes_MP_temp)[1:50] 
    }
    
    NMF_history   <- c(NMF_history , nmf_programs_sig[,names(Intersection_with_Genes_MP)[1]]) 
    Genes_MP <- Genes_MP_temp[1:50]
    
    nmf_programs_sig      <- nmf_programs_sig[,-match(names(Intersection_with_Genes_MP)[1] , colnames(nmf_programs_sig))]  # remove selected NMF
    
    Intersection_with_Genes_MP <- sort(apply(nmf_programs_sig, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
    
  }
  
  Cluster_list[[paste0("Cluster_",k)]] <- Curr_cluster
  MP_list[[paste0("MP_",k)]]           <- Genes_MP
  k <- k+1
  
  
  nmf_intersect             <- nmf_intersect[-match(Curr_cluster,rownames(nmf_intersect) ) , -match(Curr_cluster,colnames(nmf_intersect) ) ]  # remove current chosen cluster
  
  Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)   ### Sort intersection of remaining NMFs not included in any of the previous clusters
  
  Curr_cluster <- c()
  print(dim(nmf_intersect)[2])
}
length(MP_list)

sapply(MP_list, function(x) length(x[!is.na(x)]))

MP <-  do.call(cbind, MP_list)
head(MP)

#write.csv(MP, file = "00.GBM.sampleSize13.seed123.MP.csv", quote = F, row.names = F)
#write.csv(samples_names, file = "00.GBM.sampleSize13.seed123.samples.csv", quote = F, row.names = F)
#write.csv(MP, file = "00.GBM.sampleSize13.seed12345.MP.csv", quote = F, row.names = F)
#write.csv(samples_names, file = "00.GBM.sampleSize13.seed12345.samples.csv", quote = F, row.names = F)
write.csv(MP, file = "00.GBM.sampleSize13.seed1234567.MP.csv", quote = F, row.names = F)
write.csv(samples_names, file = "00.GBM.sampleSize13.seed1234567.samples.csv", quote = F, row.names = F)


#write.csv(MP, file = "00.pHGG.sampleSize13.seed123.MP.csv", quote = F, row.names = F)
#write.csv(samples_names, file = "00.pHGG.sampleSize13.seed123.samples.csv", quote = F, row.names = F)
#write.csv(MP, file = "00.pHGG.sampleSize13.seed12345.MP.csv", quote = F, row.names = F)
#write.csv(samples_names, file = "00.pHGG.sampleSize13.seed12345.samples.csv", quote = F, row.names = F)
#write.csv(MP, file = "00.pHGG.sampleSize13.seed1234567.MP.csv", quote = F, row.names = F)
#write.csv(samples_names, file = "00.pHGG.sampleSize13.seed1234567.samples.csv", quote = F, row.names = F)

length(unlist(Cluster_list)) # #individual programs that been incorportaed


#### *****  Sort Jaccard similarity plot according to new clusters:
inds_sorted <- c()
for (j in 1:length(Cluster_list)){
  
  inds_sorted <- c(inds_sorted , match(Cluster_list[[j]] , colnames(nmf_intersect_KEEP)))
  
}
inds_new <- c(inds_sorted   ,   which(is.na( match(1:dim(nmf_intersect_KEEP)[2],inds_sorted)))) ### combine inds in clusters with inds without clusters
inds_new = inds_sorted #Zhichao

# plot re-ordered similarity matrix heatmap     
nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_KEEP[inds_new,inds_new]) 

custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))


nmf_intersect_meltI_NEW$sample = sapply(as.character(nmf_intersect_meltI_NEW$Var1), 
                                        function(x) unlist(strsplit(x, split = "\\.|\\_"))[1] )
nmf_intersect_meltI_NEW$type = ifelse(grepl("\\_", nmf_intersect_meltI_NEW$Var1), "Leiden", "NMF")
head(nmf_intersect_meltI_NEW)

ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))


#dir.create("01.pHGG_0-25-12_14-14-4")
#dir.create("01.GBM_0-25-12_14-14-4")
dir.create("01.IDH_0-25-12_14-14-4/")
#png("01.pHGG_0-25-12_14-14-4/01.MP_initial.jaccard.png", width = 10, height = 10, units = "in", res = 300)
#png("01.GBM_0-25-12_14-14-4//01.MP_initial.jaccard.png", width = 10, height = 10, units = "in", res = 300)
#png("01.MP_initial.jaccard.png", width = 10, height = 10, units = "in", res = 300)
png("01.IDH_0-25-12_14-14-4/01.MP_initial.jaccard.png", width = 10, height = 10, units = "in", res = 300)
ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,30), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 16, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,30), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 16, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
dev.off()


#saveRDS(Cluster_list, file = "01.pHGG_0-25-12_14-14-4/Cluster_list.rds")
#saveRDS(MP_list, file = "01.pHGG_0-25-12_14-14-4/MP_list.rds")
#saveRDS(nmf_programs_sig, file = "01.pHGG_0-25-12_14-14-4/nmf_programs_sig.rds")
#write.csv(MP, file = "01.pHGG_0-25-12_14-14-4/01.MP_initial.csv", row.names = F, quote = F)

#saveRDS(Cluster_list, file = "01.GBM_0-25-12_14-14-4//Cluster_list.rds")
#saveRDS(MP_list, file = "01.GBM_0-25-12_14-14-4//MP_list.rds")
#saveRDS(nmf_programs_sig, file = "01.GBM_0-25-12_14-14-4/nmf_programs_sig.rds")
#write.csv(MP, file = "01.GBM_0-25-12_14-14-4//01.MP_initial.csv", row.names = F, quote = F)

saveRDS(Cluster_list, file = "01.IDH_0-25-12_14-14-4/Cluster_list.rds")
saveRDS(MP_list, file = "01.IDH_0-25-12_14-14-4/MP_list.rds")
saveRDS(nmf_programs_sig, file = "01.IDH_0-25-12_14-14-4/nmf_programs_sig.rds")
write.csv(MP, file = "01.IDH_0-25-12_14-14-4/01.MP_initial.csv", row.names = F, quote = F)

#saveRDS(Cluster_list, file = "Cluster_list.rds")
#saveRDS(MP_list, file = "MP_list.rds")
#saveRDS(nmf_programs_sig, file = "01.nmf_programs_sig.rds")
#write.csv(MP, file = "01.MP_initial.csv", row.names = F, quote = F)

#names(MP_list)=c("Neuron","Vasc","MES.Hyp","Mac","OPC.AC","Oligo","LQ.Chromatin.reg","MES","Prolif.Metab","MES.Ast","Reactive.Ast","NPC","Inflammatory.Mac")

inds_new1 <- c(inds_sorted   ,   which(is.na( match(1:dim(nmf_intersect_KEEP)[2],inds_sorted)))) ### combine inds in clusters with inds without clusters
#write.csv(nmf_intersect_KEEP[inds_new1,inds_new1], file = "01.pHGG_0-25-12_14-14-4/01.MP_initial.jaccardAll.csv")
#write.csv(nmf_intersect_KEEP[inds_new1,inds_new1], file = "01.GBM_0-25-12_14-14-4//01.MP_initial.jaccardAll.csv")
write.csv(nmf_intersect_KEEP[inds_new1,inds_new1], file = "01.IDH_0-25-12_14-14-4/01.MP_initial.jaccardAll.csv")
#write.csv(nmf_intersect_KEEP[inds_new1,inds_new1], file = "01.MP_initial.jaccardAll.csv")

# plot re-ordered similarity matrix heatmap     
inds_new1 = inds_sorted #Zhichao
nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_KEEP[inds_new1,inds_new1]) 

custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))


nmf_intersect_meltI_NEW$sample = sapply(as.character(nmf_intersect_meltI_NEW$Var1), 
                                        function(x) unlist(strsplit(x, split = "\\.|\\_"))[1] )
nmf_intersect_meltI_NEW$type = ifelse(grepl("\\_", nmf_intersect_meltI_NEW$Var1), "Leiden", "NMF")
head(nmf_intersect_meltI_NEW)

ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))


#png("01.pHGG_0-25-12_14-14-4/01.MP_initial.jaccardAll.png", width = 10, height = 10, units = "in", res = 300)
#png("01.GBM_0-25-12_14-14-4//01.MP_initial.jaccardAll.png", width = 10, height = 10, units = "in", res = 300)
png("01.IDH_0-25-12_14-14-4/01.MP_initial.jaccardAll.png", width = 10, height = 10, units = "in", res = 300)
#png("01.MP_initial.jaccardAll.png", width = 10, height = 10, units = "in", res = 300)
ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,30), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 16, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,30), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 16, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
dev.off()

###
m1 = nmf_intersect_KEEP[inds_new1,rev(inds_new1)]
h = Heatmap(m1, name = "Jaccard\nIndex", cluster_rows = F, cluster_columns = F,
            show_row_names = F, show_column_names = F, col = custom_magma,
            #top_annotation = ha, 
            width = ncol(m1)*unit(1, "mm"), 
            height = nrow(m1)*unit(1, "mm"), raster_device = "jpeg",
            use_raster = TRUE, raster_resize_mat = TRUE)

draw(h)
#h 

#pdf("01.pHGG_0-25-12_14-14-4/01.MP_initial.jaccardAll.pdf", width = 25, height = 25, useDingbats = F)
#pdf("01.GBM_0-25-12_14-14-4//01.MP_initial.jaccardAll.pdf", width = 25, height = 25, useDingbats = F)
pdf("01.IDH_0-25-12_14-14-4/01.MP_initial.jaccardAll.pdf", width = 25, height = 25, useDingbats = F)
#pdf("01.MP_initial.jaccardAll.pdf", width = 25, height = 25, useDingbats = F)
draw(h)
for(k in 1:length(Cluster_list) ){
  print(names(Cluster_list)[k])
  i = which(colnames(m1) %in% Cluster_list[[k]])
  print(i)
  decorate_heatmap_body(heatmap = "Jaccard\nIndex", code = {
    i = which(colnames(m1) %in% Cluster_list[[k]])
    grid.rect(x = (mean(i) - 0.5)/ncol(m1), y =(mean(i) - 0.5)/ncol(m1),
              width = length(i) * 1/ncol(m1), height = length(i) * 1/ncol(m1), 
              gp=gpar(col= "#1F78B4", fill = NA, lwd = 1))
  })
}
decorate_heatmap_body("Jaccard\nIndex", code = {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
dev.off()




###
sapply(Cluster_list, length)


####1. check initial MP correlation
samples_names <- (read.delim("../00.data/00.data.folder.csv", header = T, sep = ","))$Sample
#samples_names = samples_names[c(1:32)] #GBM
#samples_names = samples_names[-c(1:32)]#pHGG only

#nFeature_Spatial.cut.off <- 0;
nCount_Spatial.cut.off <- 1000;
data.home = "/Users/wuz6/Documents/Project/25.pHGG_spatial/00.data/"
MP_list.named = MP_list
MP_list.named <- mp.both
#names(MP_list.named) <- c("Neuron", "Oligo", "Vasc", "MES.Hyp", "Mac", "Cycling", "Interferon", "Unknown",
#                    "Neuron.2", "Vasc.Fibroblast", "MES/Ast", "Astocyte", "Neuron.3", "Neuron.4", "InflammarotyMac.") #pHGG
#names(MP_list.named) <- c("Mac", "MES.Hyp", "Vasc", "Neuron", "Oligo", "MES", "ACOPC", "Cycling",
#                          "Prolif.Metab", "Vasc2", "EMT", "Epithelial", "NPC", "Astrocyte", "MES.Ast") #GBM
names(MP_list.named) <- c("Oligo", "Vasc", "Mac", "InflammatoryMac", "Neuron", "AC", "Neuron.2", "AC/MES") #IDH
#names(MP_list.named) <- c("Mac", "Vasc", "Neuron", "Oligo", "MES.Hyp", "MES", "Cycling", "Interferon", "ACOPC", "Necrosis", "Neuron.2",
#                          "Fibroblast","Neuron.3", "InflammatoryMac", "Neuron.4", "AC/MES", "Prolif.Metab", "NPC", "Invasiveness", "Immune")

mp.cell.GBM = data.frame(readxl::read_xlsx("/Users/wuz6/Documents/Project/data/signatures/C.Greenwald2024Cell.xlsx"), stringsAsFactors = F)
mp.cell.GBM = as.vector(mp.cell.GBM)
MP_list.named <- mp.cell.GBM

names(MP_list.named)
MP.scores = lapply(samples_names, function(k){
  #k = samples_names[1]
  #k = "AH66"
  message(k)
  k = gsub("Y963", "Y963_", k)
  sp.data = Load10X_Spatial(paste0(data.home, k, "/outs/"), 
                            filename = "filtered_feature_bc_matrix.h5", assay = "Spatial")
  sp.data@meta.data$orig.ident = k
  sp.data <- PercentageFeatureSet(sp.data, pattern = "^MT-", col.name = "percent.mt")
  sp.data <- subset(sp.data, subset = nCount_Spatial > nCount_Spatial.cut.off & percent.mt < 20); #nFeature_Spatial > nFeature_Spatial.cut.off & 
  sp.data = NormalizeData(sp.data)
  ###
  n1 = ncol(sp.data@meta.data) + 1
  #sp.data = AddMetaData(sp.data, x)
  sig.score = sigScores(data.matrix(sp.data[['Spatial']]$data), MP_list.named, expr.center = TRUE, conserved.genes = 0.5)
  ###
  #xy = sp.data@images$slice1@coordinates
  xy <- GetTissueCoordinates(sp.data)[,1:2]
  #sp.data = AddModuleScore(sp.data, features = nmf_programs_sig.list, name = "P_", verbose = F)
  #print(head(sp.data[[]][,1:10]))
  ###
  ###
  cbind(rownames(sig.score), sp.data@meta.data, xy, sig.score)
})

MP.scores = do.call(rbind, MP.scores)
head(MP.scores, 2)


library(corrplot)
heat.col = as.character(paletteer::paletteer_c("grDevices::Red-Green", 30))
MP.scores.corr = cor(MP.scores[,-c(1:7)])


MP.scores.corr.mean = lapply(unique(MP.scores$orig.ident), function(ss){
  y = cor(MP.scores[MP.scores$orig.ident == ss,-c(1:7)])
  y
})
MP.scores.corr.mean <- purrr::reduce(.x = MP.scores.corr.mean,.f = `+`)/length(MP.scores.corr.mean)
dim(MP.scores.corr.mean)



#pdf("01.pHGG_0-25-12_14-14-4//01.rawMPscore.correlation.pdf", width = 7, height = 7, useDingbats = F)
#pdf("01.GBM_0-25-12_14-14-4//01.rawMPscore.correlation.pdf", width = 7, height = 7, useDingbats = F)
pdf("01.IDH_0-25-12_14-14-4/01.rawMPscore.correlation.pdf", width = 7, height = 7, useDingbats = F)
#pdf("01.rawMPscore.correlation.pdf", width = 7, height = 7, useDingbats = F)
#pdf("01.cellGBM.MPscore.correlation.pdf", width = 7, height = 7, useDingbats = F)
corrplot(MP.scores.corr, method = "color", order = "hclust", col = heat.col, addrect = 3, tl.col = "black")
corrplot(MP.scores.corr.mean, method = "color", order = "hclust", col = heat.col, tl.col = "black")
dev.off()


#colnames(MP.scores)[11:30] = names(MP_list.named)
spot.max.M = apply(MP.scores[,names(MP_list.named)], 1, function(x) names(MP_list.named)[which.max(x)])
#spot.second.M = apply(mp.scores, 1, function(x) colnames(mp.scores)[order(x, decreasing = T)[2]])
#spot.max.M = data.frame(MP = spot.max.M, MP2 =spot.second.M, stringsAsFactors = T)
MP.scores$MPclass = spot.max.M
#MP.scores$cancer = info$Tumor[match(MP.scores$orig.ident, info$Sample)]

sample.mp.freq = data.frame(table(MP.scores$orig.ident, MP.scores$MPclass)/rowSums(table(MP.scores$orig.ident, MP.scores$MPclass)))
head(sample.mp.freq)
x = tidyr::pivot_wider(sample.mp.freq, names_from = "Var2",values_from = "Freq")


#write.csv(MP.scores, file = "01.cellpaperGBM.MPscore.csv", row.names = F, quote = F)
#write.csv(MP.scores, file = "01.initialMP.MPscore.csv", row.names = F, quote = F)

#write.csv(x, file = "01.pHGG_0-25-12_14-14-4//01.rawMP.frac.csv")
#write.csv(x, file = "01.GBM_0-25-12_14-14-4//01.rawMP.frac.csv")
write.csv(x, file = "01.IDH_0-25-12_14-14-4/01.rawMP.frac.csv")
#write.csv(x, file = "01.rawMP.frac.csv")

#mp.col = c(paletteer::paletteer_d("ggthemes::Tableau_20"), "red","black")
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


mp.col = manualcolors[2:25]
names(mp.col) = names(MP_list.named)


p = ggplot(sample.mp.freq, aes(x = Var1, y = Freq, fill = Var2))+
  geom_bar(stat = "identity")+
  #scale_fill_manual(values = MP.col)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major.y = element_blank())+
  scale_fill_manual(values = mp.col)+
  labs(x = "", fill = "")
  #scale_x_discrete(limits = sample.reordered)
p  

#pdf("01.pHGG_0-25-12_14-14-4/01.sample.MP.bar.pdf", width = 10, height = 6, useDingbats = F)
#pdf("01.GBM_0-25-12_14-14-4//01.sample.MP.bar.pdf", width = 10, height = 6, useDingbats = F)
pdf("01.IDH_0-25-12_14-14-4/01.sample.MP.bar.pdf", width = 10, height = 6, useDingbats = F)
#pdf("01.sample.MP.bar.pdf", width = 10, height = 6, useDingbats = F)
#pdf("01.sample.MPfromCellPaper.bar.pdf", width = 10, height = 6, useDingbats = F)
p
dev.off()

#dir.create("01.pHGG_0-25-12_14-14-4//01.initial.mp.assignment.plot")
#dir.create("01.GBM_0-25-12_14-14-4//01.initial.mp.assignment.plot")
dir.create("01.IDH_0-25-12_14-14-4/01.initial.mp.assignment.plot")
#dir.create("01.initial.mp.assignment.plot")
#write.csv(MP.scores, "01.initial.mp.assignment.plot/01.MP.scores.csv")
#dir.create("01.cellPaper.mp.assignment.plot")
write.csv(MP.scores, "01.IDH_0-25-12_14-14-4/01.initial.mp.assignment.plot//01.MP.scores.csv")
#write.csv(MP.scores, "01.cellPaper.mp.assignment.plot/01.MP.scores.csv")

lapply(unique(MP.scores$orig.ident), function(k){
  p = ggplot(MP.scores[MP.scores$orig.ident == k,], aes(x = y, y = -x, col = MPclass))+
    geom_point(size=2) + 
    labs(title= k, col = "") +
    scale_color_manual(values = mp.col, name = "MP")+
    theme_void()
  #pdf(paste0("01.cellPaper.mp.assignment.plot/01.",k, ".MP.pdf"), width = 7, height = 7, useDingbats = F)
  #pdf(paste0("01.pHGG_0-25-12_14-14-4//01.initial.mp.assignment.plot/01.",k, ".MP.pdf"), width = 7, height = 7, useDingbats = F)
  #pdf(paste0("01.GBM_0-25-12_14-14-4///01.initial.mp.assignment.plot/01.",k, ".MP.pdf"), width = 7, height = 7, useDingbats = F)
  pdf(paste0("01.IDH_0-25-12_14-14-4/01.initial.mp.assignment.plot/01.",k, ".MP.pdf"), width = 7, height = 7, useDingbats = F)
  #pdf(paste0("01.initial.mp.assignment.plot/01.",k, ".MP.pdf"), width = 7, height = 7, useDingbats = F)
  print(p)
  dev.off()
  1
})



###
nFeature_Spatial.cut.off <- 0;
nCount_Spatial.cut.off <- 1000;
data.home = "/Users/wuz6/Documents/Project/25.pHGG_spatial/00.data/"

#Cluster_list = readRDS(file = "01.pHGGGBM_12-12-5/Cluster_list.rds")
#nmf_programs_sig = readRDS(file = "01.pHGGGBM_15-15-5/pHGGGBM_all_BanksyMulti_NMFrobust_programs60.rds")
nmf_programs_sig = nmf_programs_sig_KEEP
dim(nmf_programs_sig)
nmf_programs_sig.list = mat2list(nmf_programs_sig)
#nmf_programs_sig.list <- nmf_programs_sig.list[-c(361,365)] #GBM, some Visium genes not in CytAssist?

#nmf_programs_sig.list <- nmf_programs_sig.list[-c(555,559)] #all, seems two gene list are unique to Frozen FFPE Visium
#nmf_programs_sig.keep = nmf_programs_sig[,-c(555,559)]

program.scores = lapply(samples_names, function(k){
  #k = samples_names[35]
  #k = "AP67"
  message(k)
  k = gsub("Y963", "Y963_", k)
  sp.data = Load10X_Spatial(paste0(data.home, k, "/outs/"), 
                            filename = "filtered_feature_bc_matrix.h5", assay = "Spatial")
  k = gsub("Y963_", "Y963", k)
  k = strsplit(k, "\\_")[[1]][1]
  sp.data@meta.data$orig.ident = k
  message(k)
  #
  if(0){ #check if gene list doesn't work to score (<50% genes available in the data)
      y = lapply(nmf_programs_sig.list, function(gg){
          x = sum(gg %in% rownames(sp.data))
          y = x/length(gg)
          y
      })
      y = unlist(y)
      print(which(y<0.5))
      return(which(y<0.5))
  }
  ###
  sp.data <- PercentageFeatureSet(sp.data, pattern = "^MT-", col.name = "percent.mt")
  sp.data <- subset(sp.data, subset = nFeature_Spatial > nFeature_Spatial.cut.off & nCount_Spatial > nCount_Spatial.cut.off & percent.mt < 20);
  sp.data = NormalizeData(sp.data)
  #sp.data = AddMetaData(sp.data, x)
  sig.score = sigScores(data.matrix(sp.data[['Spatial']]$data), nmf_programs_sig.list, expr.center = TRUE, conserved.genes = 0.5)  
  #xy = sp.data@images$slice1@coordinates
  xy <- GetTissueCoordinates(sp.data)[,1:2]
  #sp.data = AddModuleScore(sp.data, features = nmf_programs_sig.list, name = "P_", verbose = F)
  #print(head(sp.data[[]][,1:10]))
  ###
  ###
  cbind(rownames(sig.score), sp.data@meta.data, xy, sig.score)
})

#unique(unlist(program.scores))

program.scores = do.call(rbind, program.scores)
colnames(program.scores)[1:2] = c("barcode", "sample")
program.scores[1:4,1:11]

#saveRDS(program.scores, file = "01.pHGG_0-25-12_14-14-4/01.robustBanksyNMF.programScores.rds")
#saveRDS(program.scores, file = "01.GBM_0-25-12_14-14-4//01.robustBanksyNMF.programScores.rds")
saveRDS(program.scores, file = "01.robustBanksyNMF.programScores.rds")
saveRDS(nmf_programs_sig.keep, file = "01.robustBanksyNMF.programGeneList.rds")
###


x <- unlist(lapply(program.scores, length))
##
stop here


#program.scores = readRDS("01.pHGGGBM_15-15-5/01.robustBanksyNMF.programScores.rds")
###average person correlation across samples
unique(program.scores$sample)

all.cor = lapply(unique(program.scores$sample),function(kk){
    x = program.scores[program.scores$sample == kk, -c(1:10)]
    correlation = cor(x, method = "pearson")
    correlation
})

all.cor.mean <- purrr::reduce(.x = all.cor,.f = `+`)/length(all.cor)
dim(all.cor.mean)

Cluster_list = readRDS("01.pHGG_0-25-12_14-14-4//Cluster_list.rds")
mp.pp = stack(Cluster_list)

sum(!mp.pp[,1] %in% rownames(all.cor.mean))



MP0 = as.character(mp.pp[,2][match(rownames(all.cor.mean), mp.pp[,1])])
MP0 =  ifelse( is.na(MP0), "None", MP0)
MP0 = gsub("^Cluster_", "", MP0)
table(MP0)

topAnn = HeatmapAnnotation(MP0 = ifelse(is.na(MP0), "None", MP0), col = list(MP0 = mp.col))
rowAnn =rowAnnotation(MP0 = ifelse(is.na(MP0), "None", MP0), col = list(MP0 = mp.col))

hc = hclust(as.dist(1 - all.cor.mean), method = "average")
h = Heatmap(all.cor.mean, name = "pearson", show_row_names = F, show_column_names = F,
            cluster_rows = hc, cluster_columns = hc)
            #top_annotation = topAnn, left_annotation = rowAnn)
#pdf("01.pHGG_0-25-12_14-14-4/01.all.robust.program.correlation.pdf", width = 10, height = 10, useDingbats = F)
#draw(h)
#dev.off()


plot(hc)
rect.hclust(hc , k = 30, border = 2:6)
clusterCut <- stats::cutree(tree = hc, k = 30)

table(clusterCut)

xx = data.frame(clusterCut)
xx$MP0 = as.character(mp.pp[match(rownames(xx), mp.pp[,1]), 2])
xx$MP0[is.na(xx$MP0)] = "None"
xx$clusterCut = paste0("c", xx$clusterCut)
xx$MP0 = gsub("Cluster_", "", xx$MP0)
head(xx)
xx = xx[xx$MP0 != "None",]

library(ggsankey)
df = ggsankey::make_long(xx, MP0, clusterCut)
# Step 2
dagg <- df%>%
  dplyr::group_by(node)%>%
  tally()
# Step 3
df2 <- merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)
head(df2)






ggplot(df,aes(x = x, next_x = next_x, node = node, next_node = next_node, 
               fill = factor(node) ) )+
  geom_sankey(flow.alpha = 0.5,  color = "gray40", show.legend = TRUE)+
  theme_bw()+
  scale_fill_manual(values = mp.col)+
  theme(axis.title = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.position = "none")


p = ggplot(df2,aes(x = x, next_x = next_x, node = node, next_node = next_node, 
              fill = factor(node), label = paste0(node," n=", n)) )+
  geom_sankey(flow.alpha = 0.5,  color = "gray40", show.legend = TRUE)+
  geom_sankey_label(size = 3, color = "white", fill= "gray40", hjust = -0.2)+
  theme_bw()+
  scale_fill_manual(values = mp.col)+
  theme(axis.title = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.position = "none")
p  

pdf("01.correlationCluster.vs.MP0.sankey.pdf", width = 10, height = 15, useDingbats = F)
p
dev.off()


###
mp.score.cellMP <- read.csv("01.cellPaper.mp.assignment.plot/01.MP.scores.csv", head =T, row.names = 1)
head(mp.score.cellMP)
mp.score.cellMP$tumor <- 
MP.scores <- read.csv("01.initial.mp.assignment.plot/01.MP.scores.csv", head =T, row.names = 1)
identical(mp.score.cellMP[,1], MP.scores[,1])

info$Sample3 <- gsub("_","", info$Sample)

mp.score.cellMP$sample <- gsub("_", "", mp.score.cellMP$orig.ident)
MP.scores$sample <- gsub("_", "", MP.scores$orig.ident)
mp.score.cellMP$tumor = info$Tumor[match(mp.score.cellMP$sample, info$Sample)]
MP.scores$tumor = info$Tumor[match(MP.scores$sample, info$Sample)]

library(networkD3)
#links = data.frame(table(xx$MP0, xx$clusterCut), stringsAsFactors = F)
links = data.frame(table(mp.score.cellMP$MPclass, paste0("Wu:", MP.scores$MPclass)), stringsAsFactors = F)
links = data.frame(table(mp.score.cellMP$MPclass[mp.score.cellMP$tumor == "pHGG"], 
                         paste0("Wu:", MP.scores$MPclass[mp.score.cellMP$tumor == "pHGG"])), stringsAsFactors = F)
links = data.frame(table(mp.score.cellMP$MPclass[mp.score.cellMP$tumor == "GBM"], 
                         paste0("Wu:", MP.scores$MPclass[ MP.scores$tumor == "GBM"])), stringsAsFactors = F)
colnames(links) = c("source", "target", "value")
nodes <- data.frame(name=c(as.character(links$source), as.character(links$target)) %>% unique() ) #not zero indexed

links$IDsource <- match(links$source, nodes$name)-1 #zero indexed
links$IDtarget <- match(links$target, nodes$name)-1 #

head(links)



mp.nodes.col = mp.col[as.character(unique(links$source))]

other.nodes = nodes$name[!nodes$name %in% names(mp.nodes.col)]
other.nodes.col = rep("#000000", length(other.nodes))
names(other.nodes.col) = other.nodes

build_d3_palette <- function(names, colours) {
  #' convert a character vector of hex colours into a d3 scaleOrdinal palette
  if (length(names) != length(colours)) {
    stop("names and colours lengths do not match")
  }
  xx = data.frame(names, colours);
  name_list    <- jsonlite::toJSON(names)
  palette_list <- jsonlite::toJSON(colours)
  palette_text <- paste0('d3.scaleOrdinal()',
                         ' .domain(', name_list,   ')',
                         ' .range(', palette_list, ");")
  return(list(palette_text, xx))
}

use.col = build_d3_palette(as.character(c(names(mp.nodes.col), other.nodes)), c(mp.nodes.col, other.nodes.col))

p = sankeyNetwork(Links = links, Nodes = nodes, Source = 'IDsource',
              Target = 'IDtarget', Value = 'value', NodeID = 'name',
              colourScale = use.col[[1]], 
              iterations = 0, #0: auto order off
              nodeWidth = 80, fontSize = 20)

p

#saveNetwork(p, "01.correlationCluster.vs.MP0.sankey.html", selfcontained = TRUE)
#saveNetwork(p, "01.allSample.cellpaperGBM.vs.Ours.MP.sankey.html", selfcontained = TRUE)
#saveNetwork(p, "01.pHGG.cellpaperGBM.vs.Ours.MP.sankey.html", selfcontained = TRUE)
saveNetwork(p, "01.GBM.cellpaperGBM.vs.Ours.MP.sankey.html", selfcontained = TRUE)

xx = xx[xx$MP0 != "None",]

chisq.test(as.factor(xx$MP0), as.factor(xx$clusterCut))





################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
