###Impute ligand and receptor expression
library(glmnet)
setwd("/Users/wuz6/Documents/Project/08.spatialTranscriptome/01.RNA/04.LR")

red_blue_20 = c("#124984","#1f63a8","#2f79b5","#3f8ec0","#5fa5cd","#87beda","#a7d0e4",
                "#c5dfec","#dbeaf2","#edf2f5","#f9f0eb","#fbe3d4","#fbd0b9","#f7b799",
                "#f09c7b","#e17860","#d25849","#c13639","#ae172a","#8a0b25");

LR <- read.table("https://fantom.gsc.riken.jp/5/suppl/Ramilowski_et_al_2015/data/PairsLigRec.txt", head =T, sep = "\t")
table(LR$Pair.Source)
LR <- LR[LR$Pair.Source == "known",] #n = 849
#LR$Ligand.ApprovedSymbol
#LR$Receptor.ApprovedSymbol
load("../01.QCfirst/CTA.allGlioma.seuratObj.rda")
seurat <- subset(seurat, location == "TC")
unique(seurat@meta.data$location)

sum(rownames(seurat) %in% LR$Ligand.ApprovedSymbol)
sum(rownames(seurat) %in% LR$Receptor.ApprovedSymbol)
sum(LR$Ligand.ApprovedSymbol %in% rownames(seurat) & 
      LR$Receptor.ApprovedSymbol %in% rownames(seurat)) #152


idx <- which(LR$Ligand.ApprovedSymbol %in% rownames(seurat) & 
      LR$Receptor.ApprovedSymbol %in% rownames(seurat)) #152

LR[idx[1],]
tumor.class <- unique(seurat@meta.data$tumor)[1:3]

LR.cor <- lapply(tumor.class, function(j){
    sub.seurat <- subset(seurat, tumor == j)
    message(j)
    yy <- lapply(idx, function(k){
      ligand <- LR$Ligand.ApprovedSymbol[k]
      receptor <- LR$Receptor.ApprovedSymbol[k]
      exp.ligand <- mean(sub.seurat@assays$RNA@data[ligand,])
      exp.receptor <- mean(sub.seurat@assays$RNA@data[receptor,])
      y <- cor.test(sub.seurat@assays$RNA@data[ligand,], sub.seurat@assays$RNA@data[receptor,], method = "spearman")
      c(j, ligand, receptor, exp.ligand, exp.receptor, y$p.value, y$estimate)
    })
    yy <- do.call(rbind, yy)
    yy
})
LR.cor <- data.frame(do.call(rbind, LR.cor))
head(LR.cor)

LR.cor[,4:7] <- apply(LR.cor[,4:7], 2, as.numeric)
colnames(LR.cor) <- c("Cancer", "Ligand", "Receptor", "LigandExpression","ReceptorExpression","P", "Rho")
head(LR.cor)


ggplot(LR.cor, aes(x = LigandExpression, y = ReceptorExpression, col = Rho))+
  geom_point()+
  theme_bw()+
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)


ggplot(LR.cor, aes(x = LigandExpression, y = ReceptorExpression, col = Rho))+
  geom_point()+
  theme_bw()+
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)

####
cta <- read.csv("../01.edgeR.QCfirst/00.IDHastrocytoma.CTA.cpm.filtered.csv", head =T, row.names = 1, check.names = F)
cta[1:4,1:4]
cta.keep.genes <- rownames(cta)


##CGGA data
rna <- read.table("../../00.data/TCGALGG/CGGA.mRNAseq_693.RSEM-genes.20200506.txt", head =T, row.names = 1)
clin <- read.delim("../../00.data/TCGALGG/CGGA.mRNAseq_693_clinical.20200506.txt", head =T, stringsAsFactors = F, sep = "\t")
clin <- clin[clin$IDH_mutation_status =="Mutant" & clin$X1p19q_codeletion_status == "Non-codel",]
rna <- rna[, colnames(rna) %in% clin$CGGA_ID]
rna <- DGEList(counts=rna)
rna <- cpm(rna, log = T)
rna[1:4,1:4]
sum(!d1[,1] %in% rownames(rna))
d1[,1][!d1[,1] %in% rownames(rna)]

grep("TRGC", rownames(rna), value = T)


overlap.genes <- intersect(cta.keep.genes, rownames(rna))
length(overlap.genes)


LR <- LR[which(LR$Ligand.ApprovedSymbol %in% rownames(rna) & LR$Receptor.ApprovedSymbol %in% rownames(rna)),]
nrow(LR) #674
LR.genes <- unique(c(LR$Ligand.ApprovedSymbol, LR$Receptor.ApprovedSymbol))
length(LR.genes) # 546

x <- t(rna[overlap.genes,])
x[1:4,1:4]

k.fold[[1]]

#1. test on CTA genes
set.seed(123)
k.fold = caret::createFolds(1:ncol(rna), k = 10, returnTrain = T);

test.cta.res <- lapply(1:length(overlap.genes), function(k){
#test.cta.res <- lapply(1:5, function(k){
    
      x1 <- t(rna[overlap.genes[-k],])
    y1 <- rna[overlap.genes[k],]
    #
    cv.correlation <- lapply(k.fold, function(idx){
      #idx <- k.fold[[1]]
      cvfit <- cv.glmnet(x1[idx, ], y1[idx])
      y.p <- predict(cvfit, newx = x1[-idx,], s = "lambda.min")
      correlation <- cor(y.p[,1], y1[-idx], method = "pearson", use = "pairwise.complete.obs")
      correlation
    })
    cv.correlation.median <- median(unlist(cv.correlation))
    cvfit <- cv.glmnet(x1, y1)
    p <- predict(cvfit, newx = t(cta[overlap.genes[-k],]), s = "lambda.min")
    correlation <- cor(as.vector(p[,1], "numeric"), t(cta[overlap.genes[k],]), method = "pearson", use = "pairwise.complete.obs")
    return(c(cv.correlation.median, correlation))
})

test.cta.correlation<- data.frame(do.call(rbind, test.cta.res))
colnames(test.cta.correlation) <- c("CV", "CTA.test")
rownames(test.cta.correlation) <- overlap.genes
head(test.cta.correlation)

write.csv(test.cta.correlation, file = "01.CTAgenes.test.pearson.csv", quote = F)
test.cta.correlation <- read.csv("01.CTAgenes.test.pearson.csv", head =T, row.names = 1)
#test.cta.correlation[,2][test.cta.correlation[,2] < 0 ] = 0
p <- ggplot(test.cta.correlation, aes(x = CV, y = CTA.test)) + 
  geom_density_2d_filled(na.rm = T)+
  geom_point(na.rm = T, size = 0.1, alpha = 0.5)+ 
  #geom_bin2d(bins = 70)+
  geom_vline(xintercept = 0.8, linetype = "dashed")+
  #scale_fill_gradient2(high = red_blue_20[20], low = red_blue_20[1], mid = red_blue_20[10], midpoint = 6)+
  theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank())+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  labs(x = "Cross Validataion on CGGA-693", y = "Test on CTA")
p

pdf("01.CTAgenes.test.pearson.pdf", width = 5, height = 4, useDingbats = F)
p
dev.off()


#2. test on CTA non-present LR genes
set.seed(123)
k.fold = caret::createFolds(1:nrow(x), k = 10, returnTrain = T);

all.fit <- lapply(LR.genes, function(gene.symbol){
  #gene.symbol <- LR.genes[1]
  y <- t(rna[gene.symbol,])
  cv.correlation <- lapply(k.fold, function(idx){
    #idx <- k.fold[[1]]
    cvfit <- cv.glmnet(x[idx,], y[idx])
    y.p <- predict(cvfit, newx = x[-idx,], s = "lambda.min")
    correlation <- cor(y.p[,1], y[-idx], method = "pearson")
    correlation
    })
  cv.correlation.median <- median(unlist(cv.correlation))
  cvfit <- cv.glmnet(x, y)
  return(list(gene = gene.symbol, median.correlation= cv.correlation.median, fit = cvfit))
})
save(all.fit, file = "CTApanel.LRfit.rda")

load("CTApanel.LRfit.rda")
all.correlation<- data.frame(Gene = unlist(lapply(all.fit, function(x) x[[1]])),
                             Pearson = unlist(lapply(all.fit, function(x) x[[2]])) )
head(all.correlation)

n.good <- sum(all.correlation$Pearson > 0.8, na.rm = T)

p <- ggplot(all.correlation, aes(x = 'x', y = Pearson)) +
  geom_boxplot(width = 0.5, fill= "grey", na.rm = T, outlier.size = .5)+
  geom_hline(yintercept = 0.8, linetype = "dashed")+
  geom_text(x = 0.6, y = 0.85, label = paste0("N > 0.8: ", n.good))+
  theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank(),
                     axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  #scale_x_continuous(expand = c(0, 0))+
  #scale_y_continuous(expand = c(0, 0))
  labs(x = "")
p

pdf("01.LRgenes.train.pearson.pdf", width = 5, height = 4, useDingbats = F)  
p
dev.off()

