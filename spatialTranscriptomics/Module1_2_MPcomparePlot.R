#shared MPs across three cancer types
setwd("/Users/wuz6/Documents/Project/25.pHGG_spatial/01.MP2025HE")
#rm(list=ls())
cancer.col = paletteer::paletteer_d("ggthemes::Tableau_20")[c(1,3,5)]
cancer.col <- as.character(cancer.col)
names(cancer.col) = c("GBM", "pHGG", "IDHmut")

gbm <- read.csv("01.GBM_0-30-10_10-12-12-4/01.MP_initial.csv", head =T)
gbm[1:4,1:4]
colnames(gbm) <- c("Vasc", "MES.Hyp", "Neuron", "Oligo", "Mac", "MES", "AC.OPC", "Interferon", "Prolif.Metab", "Cellcycle", "InflammatoryMac")
colnames(gbm) <- paste0("GBM:", colnames(gbm))

idh <- read.csv("01.IDHmut_0-30-10_10-12-12-4/01.MP_initial.csv", head =T)
colnames(idh) <- c("Oligo", "Vasc", "Mac", "Neuron", "InflammatoryMac", "MES.Hyp", "AC", "Neuron2", "Interferon")
colnames(idh) <- paste0("IDH:", colnames(idh))

phgg <- read.csv("01.pHGG_0-30-10_10-12-12-4/01.MP_initial.csv", head =T)
colnames(phgg) <- c("Vasc", "Neuron", "Oligo", "Mac", "MES.Hyp", "unknown", "Interferon", 
                    "Cellcycle", "InflammatoryMac", "Reactive.Ast", "Neuron2", "MES")
colnames(phgg) <- paste0("pHGG:", colnames(phgg))

all.glioma <- read.csv("01.all_0-30-10_10-12-12-4/01.MP_initial.csv", head =T)
colnames(all.glioma) <- c("Vasc", "Oligo", "Neuron", "Mac", "MES.Hyp", "MES", "InflammatoryMac", "AC", "Cellcycle", "unknown1",
                        "Interferon", "Perivascular", "Neuron2", "Neuron3", "Prolif.Metab", "unknown2", "NPC", "OPC.NPC", "Reactive.Ast", "EMT") 
colnames(all.glioma) <- paste0("Combined:", colnames(all.glioma))

###
signature.path = "/Users/wuz6/Documents/Project/data/signatures/"
C.Greenwald2024GBM = data.frame(readxl::read_xlsx(file.path(signature.path, "C.Greenwald2024Cell.xlsx"), sheet = 1))
C.Greenwald2024IDH = data.frame(readxl::read_xlsx(file.path(signature.path, "C.Greenwald2024Cell.xlsx"), sheet = 2))
head(C.Greenwald2024GBM)
head(C.Greenwald2024IDH)
sum(colnames(C.Greenwald2024IDH) %in% colnames(C.Greenwald2024GBM))
colnames(C.Greenwald2024GBM) <- paste0("Greenwald-GBM:", colnames(C.Greenwald2024GBM))
colnames(C.Greenwald2024IDH) <- paste0("Greenwald-IDH:", colnames(C.Greenwald2024IDH))
###
Ravi2023 = data.frame(readxl::read_xlsx(file.path(signature.path, "Ravi2023CancerCell.xlsx"), sheet = 1, skip = 1))
head(Ravi2023)
Ravi.list <- lapply(unique(Ravi2023$Regional.Programm), function(x) Ravi2023$Gene[Ravi2023$Regional.Programm == x])
names(Ravi.list) <- unique(Ravi2023$Regional.Programm)
Ravi2023.df <- data.frame(do.call(cbind, Ravi.list))
colnames(Ravi2023.df) <- paste0("Ravi:", colnames(Ravi2023.df))

###
combined <- cbind(gbm, idh, phgg, all.glioma, C.Greenwald2024GBM, C.Greenwald2024IDH)
combined[1:2,]

signature.intersect <- apply(combined, 2, function(x){apply(combined, 2, function(y){ length(intersect(x, y))}) })
signature.intersect2 <- apply(combined, 2, function(x){apply(Ravi2023.df, 2, function(y){ length(intersect(x, y)) }) })
signature.intersect3 <- apply(Ravi2023.df, 2, function(x){apply(combined, 2, function(y){ length(intersect(x, y)) }) })
signature.intersect4 <- apply(Ravi2023.df, 2, function(x){apply(Ravi2023.df, 2, function(y){ length(intersect(x, y)) }) })
signature.intersect[1:5, 1:5]
signature.intersect2[1:4,1:4]
signature.intersect3[1:4,1:4]
signature.intersect4[1:4,1:4]

signature.intersect <- rbind(signature.intersect, signature.intersect2)
dim(signature.intersect)
dim(signature.intersect3)
dim(signature.intersect4)

signature.intersect <- cbind(signature.intersect, rbind(signature.intersect3, signature.intersect4))

hc <- hclust(as.dist(max(signature.intersect)-signature.intersect), method="average") 
hc$labels[hc$order]

hc1 <- reorder(as.dendrogram(hc), colMeans(signature.intersect))

signature.intersect         <- signature.intersect[order.dendrogram(hc1), order.dendrogram(hc1)]
signature.intersect[1:4,1:4]


# plot re-ordered similarity matrix heatmap     
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

df <- reshape2::melt(signature.intersect) 
head(df)
df$sample = sapply(as.character(df$Var1), function(x) unlist(strsplit(x, split = ":"))[1] )

n.p <- length(unique(df$Var1))
n.p

pair.col <- RColorBrewer::brewer.pal(12, "Paired")
library(ggplot2)

y.mp <- data.frame(mp = hc$labels[hc$order])
y.mp$cancer <- ifelse(grepl("^IDH|Greenwald-IDH", y.mp$mp), "IDHmut",
                      ifelse(grepl("^pHGG", y.mp$mp), "pHGG", 
                             ifelse(grepl("^Combined", y.mp$mp), "Combined","GBM")))
#View(y.mp)
y.mp$color <- cancer.col[match(y.mp$cancer, names(cancer.col))]
y.mp$color[y.mp$cancer == "Combined"] <- "black"
p <- ggplot(data = df, aes(x=Var1, y=Var2, fill=100*value/(100-value))) + 
  geom_tile() + 
  scale_x_discrete(limits=hc$labels[hc$order], expand=c(0,0)) +
  scale_y_discrete(limits=hc$labels[hc$order], expand=c(0,0), position = "right")+
  #scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(0,50), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 25, oob=squish, name="Similarity\n(Jaccard index)")  +
  #geom_hline(yintercept = 1:n.p, color = rep(pair.col[1:2], round(n.p/2))[1:n.p], size = 0.1) +
  geom_hline(yintercept = 1:nrow(y.mp), color = y.mp$color, size = 0.1) +
  #geom_vline(xintercept = 1:nrow(y.mp), color = y.mp$color, size = 0.1) +
  #theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  #theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))+
  theme_bw()+
  theme(panel.grid = element_blank(), aspect.ratio = 1,
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.text = element_text(size = 5),
        axis.text.y = element_text(color = y.mp$color, size = 5))+
  #coord_equal()+
  labs(x = "", y = "")
p

ggsave(p, file = "01.MP_shared_across_cancer_types.pdf", width = 10, height = 8, units = "in", device = "pdf")  



