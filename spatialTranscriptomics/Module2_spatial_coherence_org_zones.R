library(patchwork)
library(parallel)

#!!! Run Functions section first
cancer.col = RColorBrewer::brewer.pal(6, "Paired")[c(4,6)]
names(cancer.col) = c("GBM", "pHGG")
# load data ---------------------------------------------------------------
setwd("/Users/wuz6/Documents/Project/25.pHGG_spatial/01.MP_pHGG_GBM/")
source("Module2_Funs.R")
#samples_names <- (read.delim("00.sample.list.txt", header = FALSE, sep = "\t"))$V1
#sample_ls = samples_names

MP_list = readRDS("01.pHGGGBM_0-25-12_15-15-5/02.finalMP_final.rds")
MP_list
names(MP_list)[names(MP_list) == "Unknown"] = "Necrosis"
MP.col.df = data.frame(readxl::read_xlsx("01.testMP_parameter.xlsx", sheet ="Sheet3"))
MP.col = MP.col.df$col.hex
names(MP.col) = MP.col.df$ann
names(MP_list)[!names(MP_list) %in% names(MP.col)]

#MP.col = manualcolors[2:25]
names(MP_list) %in% names(MP.col)
MP.col = MP.col[names(MP_list)]
MP.col


####
#spot.MP.df = read.csv("03.allSpot.MP.csv", head =T, row.names = 1)
spot.MP.df = read.csv("01.pHGGGBM_0-25-12_15-15-5/02.MP.ProgramScores.csv", head =T, row.names = 1)
colnames(spot.MP.df)
colnames(spot.MP.df)[15]
colnames(spot.MP.df)[15] = "Necrosis"

spot.MP.df = spot.MP.df[,c(1:22,42:49)]
colnames(spot.MP.df)

unique(spot.MP.df$MPclass)
spot.MP.df$MPclass[spot.MP.df$MPclass == "Unknown"] = "Necrosis"
spot.MP.df$MP = spot.MP.df$MPclass
sample_ls =  sort(unique(spot.MP.df$sample))
####
gen_clusters <- as.character(unique(spot.MP.df$MP))
gen_clusters


win_size <- 11 # run for c(5,8,11) 
rand_num <- 100

# calculate spatial coherence by window (Alt. skip and use saved results in next section) ----------------------------------------------------------
samples_num <- c(1:length(sample_ls))

#!!!!!!!!!!!!!!!!!! ---- Biowulf--------
#all_scores <- mclapply(samples_num, all_scores_fx, mc.cores = 21) #biowulf


dir.create("Spatial_coh_zones/")
dir.create("Spatial_coh_zones/sc_windows")



# define dis-organized zones ------------------------------------------------------------

all_scores <- readRDS("Spatial_coh_zones/sc_windows/spatial_win5v3.rds")
hist(unlist(all_scores), breaks = 100, main = "radius 5: spatial coherence socres", xlab = "spatial coherence score")
th1 <- as.numeric(quantile(na.omit(unlist(all_scores)),probs = (seq(0,1,0.1)))[5])

all_zones5 <- sapply(c(1:length(sample_ls)), function(i){
  set_zones <- ifelse(all_scores[[i]] <= th1, "dis", "other")
  print(table(set_zones))
  return(set_zones)
})


all_scores <- readRDS("Spatial_coh_zones/sc_windows/spatial_win8v3.rds")
hist(unlist(all_scores)[unlist(all_scores) > 0 & unlist(all_scores) < 1], breaks = 100, main = "radius 8: spatial coherence socres", xlab = "spatial coherence score")
th1 <- as.numeric(quantile(na.omit(unlist(all_scores)),probs = (seq(0,1,0.1)))[5])

all_zones10 <- sapply(c(1:length(sample_ls)), function(i){
  set_zones <- ifelse(all_scores[[i]] <= th1, "dis", "other")
  print(table(set_zones))
  return(set_zones)
})

all_scores <- readRDS("Spatial_coh_zones/sc_windows/spatial_win11v3.rds")
hist(unlist(all_scores)[unlist(all_scores) > 0 & unlist(all_scores) < 1], breaks = 100, main = "radius 11: spatial coherence socres", xlab = "spatial coherence score")
th1 <- as.numeric(quantile(na.omit(unlist(all_scores)),probs = (seq(0,1,0.1)))[5])

all_zones15 <- sapply(c(1:length(sample_ls)), function(i){
  set_zones <- ifelse(all_scores[[i]] <= th1, "dis", "other")
  print(table(set_zones))
  return(set_zones)
})



zones_intersect <- sapply(c(1:length(sample_ls)), function(x){
  inter1 <- intersect(names(all_zones5[[x]])[all_zones5[[x]] == "dis"],names(all_zones10[[x]])[all_zones10[[x]] == "dis"])
  inter2 <- intersect(inter1,names(all_zones15[[x]])[all_zones15[[x]] == "dis"])
  return(inter2)
})


all_zones <- sapply(c(1:length(sample_ls)),function(i){
  set_zones <- ifelse(names(all_scores[[i]]) %in% zones_intersect[[i]],"dis","other")
  names(set_zones) <- names(all_scores[[i]])
  return(set_zones)
})

win_size <- 4

all_smooth <- sapply(seq_len(length(sample_ls)),function(i){
  print(sample_ls[i])
  #spots_positions <- read.csv(paste("general/GBM_data/", sample_ls[i] , "/outs/spatial/tissue_positions_list.csv", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  #row.names(spots_positions) <- spots_positions$V1
  
  #spots_clusters <- readRDS(paste("MP/mp_assign_124/", sample_ls[i], ".rds", sep = ""))
  #spots_clusters <- na.omit(spots_clusters)
  #colnames(spots_clusters) <- c("barcodes", "spot_type")
  #row.names(spots_clusters)<- spots_clusters$barcodes
  
  spots_positions = spot.MP.df[spot.MP.df$sample == sample_ls[i], c("barcode", "tissue","row", "col")]
  rownames(spots_positions) = spots_positions$barcode
  colnames(spots_positions) = c("V1", "V2","V3", "V4")
  spots_clusters = spot.MP.df[spot.MP.df$sample == sample_ls[i], c("barcode", "MP")]
  colnames(spots_clusters) <- c("barcodes", "spot_type")
  row.names(spots_clusters)<- spots_clusters$barcodes  
  
  
  neighbors_table <- prox_neighbors_table_func(spots_positions,spots_clusters)
  
  all_spots <- spots_clusters$barcodes
  
  smoothing <- sapply(all_spots, function(spot){
    win_spots <- c(spot)
    sapply(c(1:win_size), function(i){
      win_spots <<- unique(c(win_spots,unique(na.omit(as.character(neighbors_table[win_spots,])))))
    })
    win_zones <- all_zones[[i]][names(all_zones[[i]]) %in% win_spots]
    if (length(win_zones[win_zones == "dis"])/length(win_zones) >= 0.5) {return("dis")}
    
    else {return("other")}
  })
  return(smoothing)
})


# define struct zones -----------------------------------------------------
all_scores <- readRDS("Spatial_coh_zones/sc_windows/spatial_win5v3.rds")
hist(unlist(all_scores), breaks = 100)
th1 <- as.numeric(quantile(na.omit(unlist(all_scores)),probs = (seq(0,1,0.1)))[5])

all_zones5 <- sapply(seq_len(length(sample_ls)), function(i){
  set_zones <- ifelse(all_scores[[i]] >= th1, "struct","other")
  print(table(set_zones))
  return(set_zones)
})

all_scores <- readRDS("Spatial_coh_zones/sc_windows/spatial_win8v3.rds")
hist(unlist(all_scores)[unlist(all_scores) > 0 & unlist(all_scores) < 1], breaks = 100)
th1 <- as.numeric(quantile(na.omit(unlist(all_scores)),probs = (seq(0,1,0.1)))[5])

all_zones10 <- sapply(seq_len(length(sample_ls)), function(i){
  set_zones <- ifelse(all_scores[[i]] >= th1, "struct", "other")
  print(table(set_zones))
  return(set_zones)
})

all_scores <- readRDS("Spatial_coh_zones/sc_windows/spatial_win11v3.rds")
hist(unlist(all_scores)[unlist(all_scores) > 0 & unlist(all_scores) < 1], breaks = 100)
th1 <- as.numeric(quantile(na.omit(unlist(all_scores)),probs = (seq(0,1,0.1)))[5])

all_zones15 <- sapply(seq_len(length(sample_ls)), function(i){
  set_zones <- ifelse(all_scores[[i]] >= th1, "struct", "other")
  print(table(set_zones))
  return(set_zones)
})



zones_intersect <- sapply(seq_len(length(sample_ls)), function(x){
  inter1 <- intersect(names(all_zones5[[x]])[all_zones5[[x]] == "struct"],names(all_zones10[[x]])[all_zones10[[x]] == "struct"])
  inter2 <- intersect(inter1,names(all_zones15[[x]])[all_zones15[[x]] == "struct"])
  return(inter2)
})


all_zones <- sapply(seq_len(length(sample_ls)),function(i){
  set_zones <- ifelse(names(all_scores[[i]]) %in% zones_intersect[[i]],"struct","other")
  names(set_zones) <- names(all_scores[[i]])
  return(set_zones)
})

win_size <- 4
all_smooth_struct <- sapply(seq_len(length(sample_ls)),function(i){
  print(sample_ls[i])
  #spots_positions <- read.csv(paste("general/GBM_data/", sample_ls[i] , "/outs/spatial/tissue_positions_list.csv", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  #row.names(spots_positions) <- spots_positions$V1
  
  #spots_clusters <- readRDS(paste("MP/mp_assign_124/", sample_ls[i], ".rds", sep = ""))
  #spots_clusters <- na.omit(spots_clusters)
  #colnames(spots_clusters) <- c("barcodes", "spot_type")
  #row.names(spots_clusters)<- spots_clusters$barcodes
  
  spots_positions = spot.MP.df[spot.MP.df$sample == sample_ls[i], c("barcode", "tissue","row", "col")]
  rownames(spots_positions) = spots_positions$barcode
  colnames(spots_positions) = c("V1", "V2","V3", "V4")
  spots_clusters = spot.MP.df[spot.MP.df$sample == sample_ls[i], c("barcode", "MP")]
  colnames(spots_clusters) <- c("barcodes", "spot_type")
  row.names(spots_clusters)<- spots_clusters$barcodes  
  
  
  neighbors_table <- prox_neighbors_table_func(spots_positions,spots_clusters)
  
  all_spots <- spots_clusters$barcodes
  
  smoothing <- sapply(all_spots, function(spot){
    win_spots <- c(spot)
    sapply(c(1:win_size), function(i){
      win_spots <<- unique(c(win_spots,unique(na.omit(as.character(neighbors_table[win_spots,])))))
    })
    win_zones <- all_zones[[i]][names(all_zones[[i]]) %in% win_spots]
    if (length(win_zones[win_zones == "struct"])/length(win_zones) >= 0.5) {return("struct")}
    else {return("other")}
  })
  return(smoothing)
})



# final zones -------------------------------------------------------------

all_zones <- sapply(seq_len(length(sample_ls)),function(i){
  sapply(names(all_smooth[[i]]), function(s){
    if (all_smooth[[i]][s] == "dis" & all_smooth_struct[[i]][s] == "struct") {
      return("Intermediate")
    } else if (all_smooth[[i]][s] == "dis") {
      return("dis")
    } else if (all_smooth_struct[[i]][s] == "struct") {
      return("struct")
    } else {
      return("Intermediate")
    }
  })
})

names(all_zones) <- sample_ls

is_small <- sapply(all_zones, function(x){
  tx <- table(x)
  
  old_clusters <- names(tx)
  add_clusters <- c("dis","Intermediate", "struct")[!c("dis","Intermediate", "struct") %in% old_clusters]
  if (length(add_clusters) > 0) {
    tx <- c(as.numeric(tx),0)
    names(tx) <- c(old_clusters, add_clusters)
  }
  tx <- tx[sort(names(tx))]
  tx10 <- sum(tx)*0.12
  print(tx)
  return(c(tx[1]<tx10, tx[2]<tx10,tx[3]<tx10))
})
colnames(is_small) <- sample_ls

dis2use <- seq_len(length(sample_ls))[!is_small[1,]]
sturct2use <- seq_len(length(sample_ls))[!is_small[3,]] 

spatial_score_zones <- lapply(seq_len(length(sample_ls)), function(i){
  # load data
  print(sample_ls[i])
  #spots_positions <- read.csv(paste("general/GBM_data/", sample_ls[i] , "/outs/spatial/tissue_positions_list.csv", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  #row.names(spots_positions) <- spots_positions$V1
  
  #spots_clusters <- readRDS(paste("MP/mp_assign_124/", sample_ls[i], ".rds", sep = ""))
  #spots_clusters <- na.omit(spots_clusters)
  #colnames(spots_clusters) <- c("barcodes", "spot_type")
  #row.names(spots_clusters)<- spots_clusters$barcodes  
  
  
  spots_positions = spot.MP.df[spot.MP.df$sample == sample_ls[i], c("barcode", "tissue","row", "col")]
  rownames(spots_positions) = spots_positions$barcode
  colnames(spots_positions) = c("V1", "V2","V3", "V4")
  spots_clusters = spot.MP.df[spot.MP.df$sample == sample_ls[i], c("barcode", "MP")]
  colnames(spots_clusters) <- c("barcodes", "spot_type")
  row.names(spots_clusters)<- spots_clusters$barcodes  
  
  
  zone_score <- lapply(c("dis", "struct"), function(z){
    
    zspots <- na.omit(names(all_zones[[i]])[all_zones[[i]] == z]) # checkk!!!!!!!!!!!
    if (length(zspots) == 0) {
      return(NA)
    } else {
      # abundance 
      programs_comp <- sample_programs_composition(spots_clusters[spots_clusters$barcodes %in% zspots,],gen_clusters)
      
      # neighbors tables   
      
      neighbors_table <- neighbors_table_func(spots_positions,spots_clusters)
      zone_neighbors_table <- neighbors_table[zspots,] 
      
      win_rand_neighbors_table <- lapply(c(1:rand_num), function(i){
        zspots_positions <- spots_positions[spots_positions$V1 %in% zspots,]
        new_pos <- sample(zspots_positions$V1, length(zspots_positions$V1), replace = FALSE)
        pos_table <- zspots_positions
        row.names(pos_table) <- new_pos 
        pos_table$V1 <- new_pos
        pos_table$V2 <- zspots_positions[new_pos, "V2"]
        
        neighbors_table <- win_prox_neighbors_table_func(pos_table,spots_clusters[spots_clusters$barcodes %in% zspots,])
        return(neighbors_table)
      }) 
      
      # spatial coherence
      
      programs_spatial_score <- sapply(sort(gen_clusters), function(cluster){
        if (!(cluster %in% spots_clusters$spot_type[spots_clusters$barcodes %in% zspots]) | 
            length(spots_clusters$barcodes[spots_clusters$spot_type == cluster & spots_clusters$barcodes %in% zspots]) < 3) {
          prog_score <- NaN
        } else {
          program_neighbors_table = neighbors_table[row.names(neighbors_table) %in% spots_clusters$barcodes[spots_clusters$spot_type == cluster & spots_clusters$barcodes %in% zspots],]
          obs <- obs_program_spatial_score(program_neighbors_table, cluster)
          one_spatial_score <- one_val(dim(program_neighbors_table)[1])
          zero_spatial_score <- zero_val(win_rand_neighbors_table, spots_clusters[spots_clusters$barcodes %in% zspots,], cluster)
          if (obs>one_spatial_score){obs <- one_spatial_score}
          if (obs<zero_spatial_score){obs <- zero_spatial_score}
          
          prog_score <- (obs - zero_spatial_score)/(one_spatial_score - zero_spatial_score)
        }
        return(prog_score)
      })    
      return(list(programs_comp,programs_spatial_score))
    }
  })
  names(zone_score) <- c("dis", "struct")
  return(zone_score)
})



# plot Relationship between spot purity and spatial coherence -------------

#purity_score_scaled <- readRDS("CNA/mal_lev.rds")

#all_dis_purity <- sapply(seq_len(length(sample_ls)), function(i){
#  dis_purity <- purity_score_scaled[[i]][all_zones[[i]] == "dis"]
#  return(dis_purity)
#})

#hist(unlist(all_dis_purity))

#all_st_purity <- sapply(seq_len(length(sample_ls)), function(i){
#  st_purity <- purity_score_scaled[[i]][all_zones[[i]] == "struct"]
#  return(st_purity)
#})

#hist(unlist(all_st_purity))


#
zones_class_fin <- sapply(seq_len(length(sample_ls)), function(i){
  print(sample_ls[[i]])
  
  #spots_clusters <- readRDS(paste("MP/mp_assign_124/", sample_ls[[i]], ".rds", sep = ""))
  #spots_clusters <- na.omit(spots_clusters)
  #colnames(spots_clusters) <- c("barcodes", "spot_type")
  #row.names(spots_clusters)<- spots_clusters$barcodes  
  
  spots_clusters = spot.MP.df[spot.MP.df$sample == sample_ls[i], c("barcode", "MP")]
  colnames(spots_clusters) <- c("barcodes", "spot_type")
  row.names(spots_clusters)<- spots_clusters$barcodes  
  
  
  v2 <- sapply(c(1:length(all_zones[[i]])), function(s){
    if(all_zones[[i]][s] == "dis"){## & purity_score_scaled[[i]][names(all_zones[[i]][s])] < 0.5) {
      return("dis_norm")
    } else if (all_zones[[i]][s] == "dis"){## & purity_score_scaled[[i]][names(all_zones[[i]][s])] >= 0.5) {
      return("dis_mal")
    } else if (all_zones[[i]][s] == "struct"){## & purity_score_scaled[[i]][names(all_zones[[i]][s])] < 0.4) {
      return("st_norm")
    } else if (all_zones[[i]][s] == "struct"){## & purity_score_scaled[[i]][names(all_zones[[i]][s])] >= 0.4) {
      return("st_mal")
    } else {
      return(all_zones[[i]][s])
    }
  })
  names(v2) <- names(all_zones[[i]])
  
  return(v2)
})
zones_class_fin = all_zones
names(zones_class_fin) <- sample_ls


zone.colors <- c(structured_mal = "#cf5e4e",structured_norm = "#a82203", disorganized_mal = "#208cc0", disorganized_norm = "#003967",intermediate ="#f1af3a")
zone.colors <- c("disorganized" = "#208cc0", "structured" = "#cf5e4e", "intermediate" ="#f1af3a")
sapply(seq_len(length(sample_ls)), function(i){
  #i =1
  print(sample_ls[[i]])
  #spots_positions <- read.csv(paste("general/GBM_data/", sample_ls[[i]] , "/outs/spatial/tissue_positions_list.csv", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  #row.names(spots_positions) <- spots_positions$V1
  
  #spots_clusters <- readRDS(paste("MP/mp_assign_124/", sample_ls[[i]], ".rds", sep = ""))
  #spots_clusters <- na.omit(spots_clusters)
  #colnames(spots_clusters) <- c("barcodes", "spot_type")
  #row.names(spots_clusters)<- spots_clusters$barcodes  
  
  spots_positions = spot.MP.df[spot.MP.df$sample == sample_ls[i], c("barcode", "tissue","row", "col")]
  rownames(spots_positions) = spots_positions$barcode
  colnames(spots_positions) = c("V1", "V2","V3", "V4")
  spots_clusters = spot.MP.df[spot.MP.df$sample == sample_ls[i], c("barcode", "MP")]
  colnames(spots_clusters) <- c("barcodes", "spot_type")
  row.names(spots_clusters)<- spots_clusters$barcodes  
  
  
  spots_filt = spots_positions[spots_positions$V1 %in% spots_clusters$barcodes,]
  #row2plot1 <- spots_clusters$barcodes[spots_clusters$barcodes %in% names(zones_class_fin[[i]])[zones_class_fin[[i]] == "st_mal"]]
  #row2plot2 <- spots_clusters$barcodes[spots_clusters$barcodes %in% names(zones_class_fin[[i]])[zones_class_fin[[i]] == "st_norm"]]
  #row2plot3 <- spots_clusters$barcodes[spots_clusters$barcodes %in% names(zones_class_fin[[i]])[zones_class_fin[[i]] == "dis_mal"]]
  #row2plot4 <- spots_clusters$barcodes[spots_clusters$barcodes %in% names(zones_class_fin[[i]])[zones_class_fin[[i]] == "dis_norm"]]
  #spots_filt$plot <- factor(ifelse(spots_filt$V1 %in% row2plot1, "structured_mal", 
  #                                 ifelse(spots_filt$V1 %in% row2plot2, "structured_norm", 
  #                                        ifelse(spots_filt$V1 %in% row2plot3, "disorganized_mal",
  #                                               ifelse(spots_filt$V1 %in% row2plot4, "disorganized_norm","intermediate")))), levels = c("disorganized_mal","disorganized_norm", "structured_mal", "structured_norm","intermediate"))
  #spots_filt$V3_ops = -(spots_filt$V3)
  
  spots_filt$zone = zones_class_fin[[i]][match(rownames(spots_filt), names(zones_class_fin[[i]]))]
  spots_filt$plot <- factor(ifelse(spots_filt$zone == "struct", "structured", 
                                   ifelse(spots_filt$zone == "dis", "disorganized","intermediate")),
                            levels = c("disorganized", "structured", "intermediate"))
  spots_filt$V3_ops = -(spots_filt$V3)
  
  gg = ggplot(spots_filt, aes(x=V4, y=V3_ops)) + 
    geom_point(aes(col=plot), size=2) + 
    labs(title=paste(sample_ls[[i]], "spatial coherence zones"), y="pos y", x="pos x") +
    scale_color_manual(values = zone.colors, name = "zone") + 
    theme_void()
  pdf(paste0("04.spatialCoH.", sample_ls[[i]], ".pdf"), width = 5,height = 5, useDingbats = F)
  print(gg)
  dev.off()
  sample_ls[[i]]
  
})


zones_cat <- c("dis_mal","dis_norm","Intermediate","st_mal","st_norm")
zones_cat <- c("dis", "Intermediate", "struct")
zone_abund <- sapply(zones_class_fin, function(z){
  spots_zone <- as.data.frame(z)
  spots_zone <- na.omit(spots_zone)
  spots_zone$barcodes <- row.names(spots_zone)
  colnames(spots_zone) <- c("spot_type","barcodes")
  
  # abundance 
  
  zones_comp <- sample_programs_composition(spots_zone,zones_cat)
  return(zones_comp)
})
write.csv(zone_abund, file = "04.zone_abund.csv", quote = F)
saveRDS(all_zones, file = "Spatial_coh_zones/all_zones.rds")

#all_zones = readRDS("Spatial_coh_zones/all_zones.rds")

zone_abund.df = data.frame(reshape2::melt(zone_abund), stringsAsFactors = F)
zone_abund.df$Var2 = gsub("ma", "", zone_abund.df$Var2)
library(ggplot2)

x = zone_abund.df[zone_abund.df$Var1 == "struct",]
x$Var2[order(x$value)]


zone_abund.df$Var1 = ifelse(zone_abund.df$Var1 == "dis", "disorganized", 
                            ifelse(zone_abund.df$Var1 == "struct", "structured", "intermediate"))


sample.reordered = x$Var2[order(x$value, decreasing = T)]
p = ggplot(zone_abund.df, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = zone.colors)+
  scale_x_discrete(limits = sample.reordered)+
  theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1, size = 7))+
  coord_cartesian(expand = F)+labs(x="", y = "freq", fill = "zone")
p


info = read.csv("../00.data/00.data.folder.csv", head =T)[,-1]
info$Sample = gsub("ma", "", info$Sample)

sample.info = readxl::read_xlsx("../00.data/00.data.folder.xlsx")

metaprog_df <- data.frame(row.names = sample.reordered)
metaprog_df$tumor <- info$Tumor[match(sample.reordered, info$Sample)]
metaprog_df$sample <- sample.reordered
metaprog_df$region <- info$Region[match(sample.reordered, info$Sample)]

p2 <- ggplot(metaprog_df, aes(x = factor(sample, levels = sample), y = 1)) + 
  geom_tile(aes(fill = tumor), color = "black") + 
  scale_fill_manual(values = cancer.col)+
  theme(legend.position = "right", axis.ticks = element_blank(),axis.title = element_text(size=24), 
        panel.background = element_rect(fill = "white"), axis.text = element_blank(), axis.line = element_blank(), axis.text.x = element_blank()) + 
  labs(x = "", y = "", fill = "") 
p2

p3 <- ggplot(metaprog_df, aes(x = factor(sample, levels = sample), y = 1)) + 
  geom_tile(aes(fill = region), color = "black") + 
  #scale_fill_manual(values = cancer.col)+
  theme(legend.position = "right", axis.ticks = element_blank(),axis.title = element_text(size=24), 
        panel.background = element_rect(fill = "white"), axis.text = element_blank(), axis.line = element_blank(), axis.text.x = element_blank()) + 
  labs(x = "", y = "", fill = "") 
p3

pdf("04.zone_abund.pdf", width = 16, height = 8, useDingbats = F)
egg::ggarrange(p, p2, p3, ncol = 1, nrow = 3, heights = c(40, 2, 2))
dev.off()


mp.zone <- table(spot.MP.df$MPclass, spot.MP.df$zone)/rowSums(table(spot.MP.df$MPclass, spot.MP.df$zone))
mp.order <- names(sort(mp.zone[,3], decreasing = T))
mp.order

mp.zone = reshape2::melt(mp.zone)
head(mp.zone)
#reorder by % of structured


p = ggplot(mp.zone, aes(x = Var1,, y = value, fill = Var2))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = zone.colors)+
  scale_x_discrete(limits = mp.order)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  coord_cartesian(expand = F)+
  labs(x = "", y = "freq", fill = "")
pdf("04.zone_abund-program.pdf", width = 7, height = 4, useDingbats = F)
p
dev.off()

###
head(spot.MP.df)
all_zones2 = lapply(seq_len(length(all_zones)), function(x) return(cbind(names(all_zones)[x], names(all_zones[[x]]), all_zones[[x]])))

all_zones.df = do.call(rbind, all_zones2)
rownames(all_zones.df) = paste0(all_zones.df[,1], ".", all_zones.df[,2])

rownames(spot.MP.df) = gsub("ma", "", rownames(spot.MP.df))
rownames(spot.MP.df) = gsub("Y963_", "Y963", rownames(spot.MP.df))

for( i in grep("GSM", sample.info$Sample) ){
   x1 = sample.info$Sample[i]
   x2 = sample.info$SampleName[i]
   rownames(spot.MP.df) = gsub(x1, x2, rownames(spot.MP.df))
}

sum(!rownames(all_zones.df) %in% rownames(spot.MP.df))
sum(!rownames(spot.MP.df) %in% rownames(all_zones.df))

rownames(all_zones.df)[!rownames(all_zones.df) %in% rownames(spot.MP.df)][1:4]
rownames(spot.MP.df)[!rownames(spot.MP.df) %in% rownames(all_zones.df)][1:4]


spot.MP.df$zone = all_zones.df[,3][match(rownames(spot.MP.df), rownames(all_zones.df))]
spot.MP.df$zone = ifelse(spot.MP.df$zone == "dis", "disorganized", 
                         ifelse(spot.MP.df$zone == "struct", "structured", "intermediate"))

###
spot.MP.df.pHGG = spot.MP.df[spot.MP.df$cancer == "pHGG",]
spot.MP.df.GBM = spot.MP.df[spot.MP.df$cancer == "GBM",]

#x = data.frame(table(spot.MP.df$MP, spot.MP.df$zone)/rowSums(table(spot.MP.df$MP, spot.MP.df$zone)))
x1 = data.frame(table(spot.MP.df.pHGG$MP, spot.MP.df.pHGG$zone)/rowSums(table(spot.MP.df.pHGG$MP, spot.MP.df.pHGG$zone)))
x2 = data.frame(table(spot.MP.df.GBM$MP, spot.MP.df.GBM$zone)/rowSums(table(spot.MP.df.GBM$MP, spot.MP.df.GBM$zone)))

p1 = ggplot(x1, aes(x = Var1, y = Freq, fill = Var2))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = zone.colors)+
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "bottom",
                     axis.text.x = element_text(angle = 30, hjust = 1))+
  #scale_x_discrete(limits = x1$Var1[order(x1$Freq, decreasing = T)])+
  labs(x = "", fill = "", title = "pHGG")
p1

p2 = ggplot(x2, aes(x = Var1, y = Freq, fill = Var2))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = zone.colors)+
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "bottom",
                     axis.text.x = element_text(angle = 30, hjust = 1))+
  #scale_x_discrete(limits = x1$Var1[order(x1$Freq, decreasing = T)])+
  labs(x = "", fill = "", title = "GBM")
p2
cowplot::plot_grid(p1, p2)


coh_scores.win1 <- readRDS("Spatial_coh_zones/sc_windows/spatial_win1v3.rds")
coh_scores.win5 <- readRDS("Spatial_coh_zones/sc_windows/spatial_win5v3.rds")
coh_scores.win8 <- readRDS("Spatial_coh_zones/sc_windows/spatial_win8v3.rds")
coh_scores.win11 <- readRDS("Spatial_coh_zones/sc_windows/spatial_win11v3.rds")
#coh_scores <- readRDS("Spatial_coh_zones/sc_windows/GBMspatial_win1v3.rds")

#coh_scores <- readRDS("Spatial_coh_zones/sc_windows/GBMspatial_win1v3.rds")
#coh_scores <- readRDS("Spatial_coh_zones/sc_windows/Cellpaper_GBMspatial_win1v3.rds")
#coh_scores = lapply(seq_len(length(coh_scores)), function(k){y = coh_scores[[k]]; names(y) = paste0(sample_ls[k], ".", names(coh_scores[[k]])); y})
rename.coh_scores <- function(win.result){
  y <- lapply(1:length(win.result), function(k){
      y1 = win.result[[k]];
      names(y1) = paste0(names(win.result)[k], ".", names(y1));
      y1
  })
  y = unlist(y)
  y
}
coh_scores.win1 = rename.coh_scores(coh_scores.win1)
coh_scores.win5 = rename.coh_scores(coh_scores.win5)
coh_scores.win8 = rename.coh_scores(coh_scores.win8)
coh_scores.win11 = rename.coh_scores(coh_scores.win11)
coh_scores.win1[1:4]
coh_scores.win5[1:4]
coh_scores.win8[1:4]
coh_scores.win11[1:4]

spot.MP.df$coh.score.win1 = coh_scores.win1[match(rownames(spot.MP.df), names(coh_scores.win1))]
spot.MP.df$coh.score.win5 = coh_scores.win5[match(rownames(spot.MP.df), names(coh_scores.win5))]
spot.MP.df$coh.score.win8 = coh_scores.win8[match(rownames(spot.MP.df), names(coh_scores.win8))]
spot.MP.df$coh.score.win11 = coh_scores.win11[match(rownames(spot.MP.df), names(coh_scores.win11))]
sum(is.na(spot.MP.df$coh.score.win1))
sum(!is.na(spot.MP.df$coh.score.win1))
write.csv(spot.MP.df, file = "04.spatialCoHZ.MP.csv", quote = F)

#plot(density(spot.MP.df$coh.score, na.rm = T))

x = spot.MP.df[!is.na(spot.MP.df$coh.score.win5),] #[!is.na(spot.MP.df$coh.score) & spot.MP.df$sample %in% grep("GSM|UKF", unique(spot.MP.df$sample), value = T),] # spot.MP.df$cancer == "GBM" & ,]
unique(x$sample)

unique(x$coh.score.win1)

ggplot(x, aes(x = MPclass, y = coh.score.win5, fill = cancer)) +
  geom_boxplot()+
  theme_bw()

ggplot(x, aes(x = MPclass, y = coh.score.win5, fill = cancer)) +
  geom_boxplot()+
  theme_bw()


p.list = lapply(unique(x$sample), function(k){
    #k = unique(x$sample)[1]
    message(k)
    x1 = x[x$sample == k,]
    x1$GBMMPclass = gsub("Greenward.", "", x1$GBMMPclass)
    p1 = ggplot(x1, aes(x = row, y = col, col = coh.score.win5)) +
      geom_point()+
      theme_bw() + theme(panel.grid = element_blank())+
      labs(title = k)+ coord_equal()
    p2 = ggplot(x1, aes(x = row, y = col, col = MPclass)) +
      geom_point()+
      theme_bw() + theme(panel.grid = element_blank())+
      labs(title = k)+ coord_equal()+
      scale_color_manual(values = MP.col)
    pdf(paste0("Spatial_coh_zones/sc_windows/03.pHGGGBM.win5.", k, ".coh.MP.pdf"), width = 14, height = 8, useDingbats = F)
    #pdf(paste0("Spatial_coh_zones/sc_windows/04.GBM.win1.", k, ".coh.MP.pdf"), width = 14, height = 7, useDingbats = F)
    print(cowplot::plot_grid(p1,p2))
    dev.off()
    list(p1, p2)
})

#cowplot::plot_grid(plotlist = p.list[[2]])

table(x$cancer)
coh_score.mean = data.frame(aggregate(x$coh.score.win5, list(x$MPclass), mean))
coh_score.sd = data.frame(aggregate(x$coh.score.win5, list(x$MPclass), sd))
coh_score.meansd = cbind(coh_score.mean, coh_score.sd[,2])
colnames(coh_score.meansd) = c("MP", "Mean", "SD")

head(coh_score.meansd)

p = ggplot(coh_score.meansd, aes(x = reorder(MP, -Mean), y = Mean))+
  geom_point()+
  geom_errorbar(aes(x=MP, ymin=Mean-SD, ymax=Mean+SD), width=0.4, position = position_dodge(.9))+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))+
  labs(x = "metaprogram", y = "mean spatial coherence")
p

x1 = x[x$cancer == "pHGG",]
coh_score.mean.pHGG = data.frame(aggregate(x1$coh.score, list(x1$MPclass), mean))
coh_score.sd.pHGG = data.frame(aggregate(x1$coh.score, list(x1$MPclass), sd))
coh_score.meansd.pHGG = cbind(coh_score.mean.pHGG, coh_score.sd.pHGG[,2])
colnames(coh_score.meansd.pHGG) = c("MP", "Mean", "SD")
coh_score.meansd.pHGG$cancer = "pHGG"

head(coh_score.meansd.pHGG)
p1 = ggplot(coh_score.meansd.pHGG, aes(x = reorder(MP, -Mean), y = Mean))+
  geom_point()+
  geom_errorbar(aes(x=MP, ymin=Mean-SD, ymax=Mean+SD), width=0.4, position = position_dodge(.9))+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))+
  labs(x = "metaprogram", y = "mean spatial coherence", title = "pHGG")
p1


x2 = x[x$cancer == "GBM",]
coh_score.mean.GBM = data.frame(aggregate(x2$coh.score, list(x2$MPclass), mean))
coh_score.sd.GBM = data.frame(aggregate(x2$coh.score, list(x2$MPclass), sd))
coh_score.meansd.GBM = cbind(coh_score.mean.GBM, coh_score.sd.GBM[,2])
colnames(coh_score.meansd.GBM) = c("MP", "Mean", "SD")
coh_score.meansd.GBM$cancer = "GBM"

head(coh_score.meansd.GBM)
p2 = ggplot(coh_score.meansd.GBM, aes(x = reorder(MP, -Mean), y = Mean))+
  geom_point()+
  geom_errorbar(aes(x=MP, ymin=Mean-SD, ymax=Mean+SD), width=0.4)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))+
  labs(x = "metaprogram", y = "mean spatial coherence", title = "GBM")

p2

##
coh_score.meansd.GBM.pHGG = rbind(coh_score.meansd.GBM, coh_score.meansd.pHGG)

p3 = ggplot(coh_score.meansd.GBM.pHGG, aes(x = reorder(MP, -Mean), y = Mean, color = cancer, group = cancer))+
  geom_point( position = position_dodge(0.9))+
  #geom_bar( stat = "identity", alpha = 0.5, position = position_dodge()  ) +
  geom_errorbar(aes(x=MP, ymin=Mean-SD, ymax=Mean+SD), width=0.4, position = position_dodge(0.9))+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))+
  labs(x = "metaprogram", y = "mean spatial coherence")+
  scale_color_manual(values = cancer.col)
  #facet_grid(~cancer, space = "free_x")
p3

cowplot::plot_grid(p1, p2)
#pdf("01.Cellpaper.signature.win5.coh.pdf", width = 5, height = 5, useDingbats = F)
pdf("04.pHGGGBM.win5.coh.pdf", width = 9, height = 5, useDingbats = F)  
p1
p2
p3
dev.off()



plot(density(as.vector(table(spot.MP.df$MPclass, spot.MP.df$sample))))
quantile(as.vector(table(spot.MP.df$MPclass, spot.MP.df$sample)))




xx = lapply(unique(spot.MP.df$sample), function(ss){
  #ss = unique(spot.MP.df$sample)[1]
  tmp = spot.MP.df[spot.MP.df$sample ==ss, ]
  x1 = table(tmp$MP)
  mean.df = data.frame(aggregate(tmp$coh.score.win5, list(tmp$MP), function(x){mean(x, na.rm = T)}))
  colnames(mean.df) = c("MP", "coh.score")
  #sd.df = data.frame(aggregate(tmp$coh.score, list(tmp$MP), function(x){sd(x, na.rm = T)}))
  #mean.df$sd = sd.df[,2]
  mean.df$sample = ss
  mean.df$MPcnt = x1[mean.df$MP]
  ###keep 5%
  #cnt.cutoff = 0.05*sum(mean.df$MPcnt)
  #mean.df = mean.df[mean.df$MPcnt>cnt.cutoff,]
  mean.df
})
xx = do.call(rbind, xx)
xx$cancer = sample.info$Cancer[match(xx$sample, sample.info$SampleName)]
table(xx$cancer)
table(xx$MP)
head(xx)

xx.mean.sd = lapply(unique(xx$sample), function(ss){
  tmp = xx[xx$sample ==ss,]
  c(ss, mean(tmp$coh.score, na.rm = T), sd(tmp$coh.score, na.rm = T))
})
xx.mean.sd = data.frame(do.call(rbind, xx.mean.sd))
colnames(xx.mean.sd) = c("sample", "mean", "sd")
xx.mean.sd$mean = as.numeric(xx.mean.sd$mean)
xx.mean.sd$sd = as.numeric(xx.mean.sd$sd)
xx.mean.sd$cancer = sample.info$Cancer[match(xx.mean.sd$sample, sample.info$SampleName)]
head(xx.mean.sd)


###top abundant MPs
sample.mp.freq = table(spot.MP.df$sample, spot.MP.df$MPclass)/rowSums(table(spot.MP.df$sample, spot.MP.df$MPclass))

sample.mp.freq.top = lapply(1:nrow(sample.mp.freq), function(x){
  y = colnames(sample.mp.freq)[which(sample.mp.freq[x,] > 0.1)]
  print(y)
  y
})
names(sample.mp.freq.top) = rownames(sample.mp.freq)
xx2 = lapply(1:length(sample.mp.freq.top), function(k){
    ss = names(sample.mp.freq.top)[k]
    #message(ss)
    #print(sample.mp.freq.top4[k,])
    mean.df = xx[xx$sample== ss,]
    #print(head(mean.df))
    mean.df[mean.df$MP %in% sample.mp.freq.top[[k]],]
})
xx2 = do.call(rbind, xx2)

p = ggplot(xx.mean.sd, aes(x=reorder(sample, -mean), y = mean)) +
  geom_point()+
  geom_errorbar(aes(x=sample, ymin=mean-sd, ymax=mean+sd), width=0.4, position = position_dodge(.9))+
  geom_point(inherit.aes = F, data = xx2, aes(x = reorder(sample, -coh.score), y = coh.score, col = MP), size = 3, alpha = 0.8)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))+
  scale_color_manual(values = MP.col)+
  labs(x = "sample", y = "spatial coherence")+
  facet_grid(~cancer, scales = "free_x", space='free_x')
p  


pdf("04.spatial.coh.sample.pdf", width = 10, height = 5, useDingbats = F)  
p
dev.off()

###
spot.MP.df = read.csv("04.spatialCoHZ.MP.csv", head =T, row.names = 1, check.names = F)
head(spot.MP.df)

x = spot.MP.df[spot.MP.df$MPclass == "Necrosis",]
nrow(x)

library(ggplot2)
ggplot(x, aes(x =sample, y = coh.score.win5))+
  geom_boxplot()+
  theme_bw()


p = ggplot(x, aes(x = coh.score.win5))+
  geom_density(fill = "gray")+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid = element_blank())+
  coord_cartesian(expand = F)+
  labs(title = "Necorsis", x = "coherence score")

pdf("06.Necrosis.score.distribution.pdf", width = 5, height = 4, useDingbats = F)
p
dev.off()

#sum(x$coh.score.win5>0.8, na.rm = T)
spot1 <- x[which(x$coh.score.win5>0.6 & x$coh.score.win5 < 0.8),]
spot2 <- x[which(x$coh.score.win5>0.4 & x$coh.score.win5 < 0.6),]

table(spot1$cancer)
table(spot2$cancer)


###Biowulf: find genes associated with Necrosis stage
library(Seurat)
info <- read.csv("../01.preproess/00.data.folder.csv")
source("/data/wuz6/project/25.pedHGG_spatial/01.preproess/01.MP/Module2.func.R")
load("/data/wuz6/project/25.pedHGG_spatial/01.preproess/01.cluster/01.all.sample.seuratList.rda")

###
spot.MP.df = read.csv("04.spatialCoHZ.MP.csv", head =T, row.names = 1, check.names = F)
x = spot.MP.df[spot.MP.df$MPclass == "Necrosis",]
spot1 <- x[which(x$coh.score.win5>=0.4 & x$coh.score.win5 < 0.6),]
spot2 <- x[which(x$coh.score.win5>0.2 & x$coh.score.win5 < 0.4),]

#find neighbor spots
more.necrotic <- lapply(unique(spot1$sample), function(s_name){
    necrosis.spot <- spot1[spot1$sample == s_name,]
    #s_name = unique(spot1$sample)[1]
    spots_positions = spot.MP.df[spot.MP.df$sample == s_name, c("barcode", "tissue","row", "col")]
    rownames(spots_positions) = spots_positions$barcode
    colnames(spots_positions) = c("V1", "V2","V3", "V4")
    spots_clusters = spot.MP.df[spot.MP.df$sample == s_name, c("barcode", "MP")]
    colnames(spots_clusters) <- c("barcodes", "spot_type")
    row.names(spots_clusters)<- spots_clusters$barcodes
    all_win_neighbors_table <- neighbors_table_funcV2(spots_positions,spots_clusters)
    necrosis.spot.neighbors <- all_win_neighbors_table[rownames(all_win_neighbors_table) %in% necrosis.spot$barcode,]
    necrosis.spot.neighbors <- as.vector(necrosis.spot.neighbors)
    y <- spot.MP.df[spot.MP.df$sample == s_name & spot.MP.df$barcode %in% necrosis.spot.neighbors,]
    y
})
#names(more.necrotic) <- unique(spot1$sample)
more.necrotic <- do.call(rbind, more.necrotic)
more.necrotic$necrotic <- "more"

###
less.necrotic <- lapply(unique(spot2$sample), function(s_name){
  necrosis.spot <- spot2[spot2$sample == s_name,]
  #s_name = unique(spot1$sample)[1]
  spots_positions = spot.MP.df[spot.MP.df$sample == s_name, c("barcode", "tissue","row", "col")]
  rownames(spots_positions) = spots_positions$barcode
  colnames(spots_positions) = c("V1", "V2","V3", "V4")
  
  #spots_clusters <- readRDS(paste("MP/mp_assign_124/", sample_ls[i], ".rds", sep = ""))
  #spots_clusters <- na.omit(spots_clusters)
  #colnames(spots_clusters) <- c("barcodes", "spot_type")
  ###
  #spots_clusters = spot.MP.df[spot.MP.df$sample == sample_ls[i], c("barcode", "MPclass")]
  spots_clusters = spot.MP.df[spot.MP.df$sample == s_name, c("barcode", "MP")]
  colnames(spots_clusters) <- c("barcodes", "spot_type")
  row.names(spots_clusters)<- spots_clusters$barcodes
  
  
  all_win_neighbors_table <- neighbors_table_funcV2(spots_positions,spots_clusters)
  #neighbors_table <- neighbors_table_func(spots_positions,spots_clusters)
  necrosis.spot.neighbors <- all_win_neighbors_table[rownames(all_win_neighbors_table) %in% necrosis.spot$barcode,]
  necrosis.spot.neighbors <- as.vector(necrosis.spot.neighbors)
  y <- spot.MP.df[spot.MP.df$sample == s_name & spot.MP.df$barcode %in% necrosis.spot.neighbors,]
  y
})
#names(less.necrotic) <- unique(spot2$sample)
less.necrotic <- do.call(rbind, less.necrotic)
less.necrotic$necrotic = "less"
##

#all.seurat
necrotic.neighobrs <- rbind(more.necrotic, less.necrotic)


all.seurat.names <- unlist(lapply(all.seurat, function(sp.data){
  as.character(sp.data[[]]$orig.ident[1])
}))
all.seurat.names2 <- info$Sample2[match(all.seurat.names, info$Sample)]

subset.sp <- lapply(1:length(all.seurat), function(k){
    #k = 32
    s.name1 <- all.seurat.names[k]
    s.name2 <- all.seurat.names2[k]
    message(s.name1, "==", s.name2)
    use.spots <- necrotic.neighobrs[necrotic.neighobrs$sample == s.name2, c("sample","barcode","MPclass", "zone", "necrotic")]
    dup.barcode <- use.spots$barcode[duplicated(use.spots$barcode)]
    use.spots <- use.spots[!use.spots$barcode %in% dup.barcode,]
    rownames(use.spots) <- use.spots$barcode
    #sum(colnames(sp.data) %in% x$barcode)
    if(nrow(use.spots) > 0 & sum(use.spots$barcode %in% colnames(sp.data)) > 0){ ##Need to check the data
        #print(head(use.spots))
        print("has spots.")
        sp.data <- all.seurat[[k]]
        DefaultAssay(sp.data) <- "Spatial"
        y <- sp.data[,colnames(sp.data) %in% use.spots$barcode]
        y <- AddMetaData(y, use.spots)
        #print(head(y[[]]))
        return(y)
    }else{
      print("No spot.")
      return(NULL)
    }
})

names(subset.sp) <- all.seurat.names
subset.sp <- subset.sp[lengths(subset.sp)>0] #remove empty

sp.combined <- merge(subset.sp[[1]], subset.sp[-1], add.cell.ids = names(subset.sp))
sp.combined <- NormalizeData(sp.combined)
Idents(sp.combined) <- sp.combined@meta.data$necrotic
sp.combined <- JoinLayers(sp.combined)
deg <- FindMarkers(sp.combined, ident.1 = "less", ident.2 = "more", logfc.threshold = -Inf)
deg$gene <- row.names(deg)
avg.exp <- AverageExpression(sp.combined)
deg <- cbind(deg, avg.exp[[1]][deg$gene,])
write.csv(deg, file = "01.necrotic.more.less.DEG.csv", row.names=F)

setwd("/Users/wuz6/Documents/Project/25.pHGG_spatial/05.necorsis/")
