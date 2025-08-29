GBM_path_subtypes <- TCGAbiolinks::TCGAquery_subtype(tumor = "gbm")
LGG_path_subtypes <- TCGAbiolinks::TCGAquery_subtype(tumor = "lgg")

identical(colnames(GBM_path_subtypes), colnames(LGG_path_subtypes))

subtypes <- rbind(GBM_path_subtypes, LGG_path_subtypes)

table(subtypes$IDH.status, subtypes$IDH.codel.subtype)

data <- read.table("../../00.data/TCGA.GBMLGG.sampleMap_HiSeqV2.gz", row.names = 1, header = 1, check.names = F)
data[1:4,1:4]
data <- data[rownames(data) %in% rownames(seurat), ]

sum(colnames(data) %in% subtypes$patient)

table(substr(colnames(data), 14,16))

exclude.idx <- which(substr(colnames(data), 14,16) == "11")
data <- data[,-exclude.idx]

data.patient <- substr(colnames(data), 1,12)
sum(duplicated(data.patient)) #20, can ignore

###
subtypes$patient <- as.character(subtypes$patient)
aidh.idx <- which(data.patient %in% subtypes$patient[subtypes$IDH.status == "Mutant" & 
                                                       subtypes$IDH.codel.subtype == "IDHmut-non-codel"])

oidh.idx <- which(data.patient %in% subtypes$patient[subtypes$IDH.status == "Mutant" & 
                                                       subtypes$IDH.codel.subtype == "IDHmut-codel"])

gbm.idx <- which(data.patient %in% subtypes$patient[subtypes$IDH.status == "WT" & 
                                                      subtypes$Grade == "G4"])

length(aidh.idx) #268
length(oidh.idx) #172
length(gbm.idx) #148



res1 <- res[rownames(res) %in% rownames(data) & res$avg_log2FC > 0.5, ]
sub.data <- data[rownames(res1), c(gbm.idx, aidh.idx, oidh.idx)]

col.ann <- data.frame(tumor = c(rep("GBM", length(gbm.idx)), rep("A IDH", length(aidh.idx)), 
                                rep("O IDH", length(oidh.idx))))
rownames(col.ann) <- colnames(sub.data)
gene.ann <- res1[,6, drop = F]
head(gene.ann)

dev.off()
pdf("all.glioma/01.spatialTumorTypeDEG.BlukTCGA.clusterHeatmap.pdf", width = 10, height = 15, useDingbats = F)
pheatmap(sub.data,
         #cluster_rows = gene.cluster, cluster_cols = roi.cluster,
         #scale = "row",
         annotation_col = col.ann,
         annotation_row = gene.ann,
         annotation_colors = list(tumor = cancer.col),
         show_rownames = F, show_colnames = F, 
         color = red_blue_20)

dev.off()




###

gbm <- lapply(1:nrow(data), function(k){
  t <- t.test(as.numeric(data[k,gbm.idx]), as.numeric(data[k, c(oidh.idx, aidh.idx)]))
  c(t$p.value, t$estimate)
  
})
gbm <- do.call(rbind, gbm)
gbm <- data.frame(apply(gbm, 2, as.numeric))
colnames(gbm) <- c("P.value", "GBM", "IDH")
rownames(gbm) <- rownames(data)
gbm$FDR <- p.adjust(gbm$P.value, method = "fdr")
gbm$FC <- gbm$GBM/gbm$IDH
gbm <- gbm[!is.na(gbm$P.value), ]
head(gbm)
sum(gbm$FDR<0.05)
sum(gbm$FDR<0.05 & (gbm$FC > 2| gbm$FC < 0.5))

###
idh.ao <- lapply(1:nrow(data), function(k){
  t <- t.test(as.numeric(data[k,aidh.idx]), as.numeric(data[k, oidh.idx]))
  c(t$p.value, t$estimate)
  
})
idh.ao <- do.call(rbind, idh.ao)
idh.ao <- data.frame(apply(idh.ao, 2, as.numeric))
colnames(idh.ao) <- c("P.value", "AIDH", "OIDH")
rownames(idh.ao) <- rownames(data)
idh.ao$FDR <- p.adjust(idh.ao$P.value, method = "fdr")
idh.ao$FC <- idh.ao$AIDH/idh.ao$OIDH
idh.ao <- idh.ao[!is.na(idh.ao$P.value), ]
head(idh.ao)
sum(idh.ao$FDR<0.05)
sum(idh.ao$FDR<0.05 & (idh.ao$FC > 2| idh.ao$FC < 0.5))

###
idh.a.over <- rownames(idh.ao)[idh.ao$FDR<0.05 & idh.ao$FC > 2]
idh.a.under <- rownames(idh.ao)[idh.ao$FDR<0.05 & idh.ao$FC < 0.5]
gbm.over <- rownames(gbm)[gbm$FDR<0.05 & gbm$FC > 2]
gbm.under <- rownames(gbm)[gbm$FDR<0.05 & gbm$FC < 0.5]


###
length(gbm.idh.genes)
length(idh.ao.genes)
length(intersect(idh.ao.genes, gbm.idh.genes))

write.csv(gbm, file = "GBM.vs.IDH.csv", quote = F)
write.csv(idh.ao, file = "IDHA.vs.IDHO.csv", quote = F)

save(idh.a.over, idh.a.under, gbm.over, gbm.under, file = "00.GBM.IDH.A.O.DEG.rda")


sub.data <- data[rownames(gene.ann), c(gbm.idx, aidh.idx, oidh.idx)]

col.ann <- data.frame(tumor = c(rep("GBM", length(gbm.idx)), rep("A IDH", length(aidh.idx)), 
                                                        rep("O IDH", length(oidh.idx))))
rownames(col.ann) <- colnames(sub.data)

pdf("all.glioma/01.tcgaDEG.BlukTCGA.clusterHeatmap.pdf", width = 10, height = 15, useDingbats = F)
pheatmap(sub.data,
         #cluster_rows = gene.cluster, cluster_cols = roi.cluster,
         #scale = "row",
         annotation_col = col.ann,
         annotation_row = gene.ann,
         annotation_colors = list(tumor = cancer.col),
         show_rownames = F, show_colnames = F, 
         color = red_blue_20)

dev.off()

