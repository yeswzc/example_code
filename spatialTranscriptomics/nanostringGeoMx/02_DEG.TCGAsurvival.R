setwd("/Users/wuz6/Documents/Project/08.spatialTranscriptome/01.RNA/02.DEG.samplewise/")
library(TCGAbiolinks)
library("survival")
library("survminer")

lgg_clin <- TCGAquery_subtype(tumor = "lgg")
lgg_clin <- lgg_clin[lgg_clin$IDH.status == "Mutant" & lgg_clin$X1p.19q.codeletion == "non-codel",]
nrow(lgg_clin)

rna <- read.table("../../00.data/TCGALGG/TCGA.LGG.sampleMap_HiSeqV2.gz", head =T, row.names = 1, check.names = F)
rna[1:4,1:4]

lgg_clin[1:4,1:4]

colnames(rna) <- sapply(colnames(rna), function(x) substr(x, 1,12))
colnames(rna)
overlap <- intersect(colnames(rna), lgg_clin$patient)

lgg_clin <- lgg_clin[lgg_clin$patient %in% overlap,]
rna <- rna[,colnames(rna) %in% overlap]

dim(rna)
dim(lgg_clin)

survival.data <- subset(lgg_clin, select = c("Survival..months.", "Vital.status..1.dead."))
head(survival.data)

###

spatial_deg <- read.csv("04.all.DEGs.csv", head =T, row.names = 1)
head(spatial_deg)

tc.up <- rownames(spatial_deg)[spatial_deg$TC == "up"]
tc.up <- tc.up[tc.up %in% rownames(rna)]
pheatmap(rna[tc.up,], show_colnames = F, scale = "row")

tc.down <- rownames(spatial_deg)[spatial_deg$TC == "down"]
tc.down <- tc.down[tc.down %in% rownames(rna)]
pheatmap(rna[tc.down,], show_colnames = F, scale = "row")

knn1 <- kmeans(t(rna[tc.down,]), 2)
knn1$cluster
lgg_clin$km <- knn1$cluster[match(lgg_clin$patient, names(knn1$cluster))]
fit <- survfit(Surv(Survival..months., Vital.status..1.dead.) ~ km, data = lgg_clin)
print(fit)
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


###
ie.up <- rownames(spatial_deg)[spatial_deg$IE == "up"]
ie.up <- ie.up[ie.up %in% rownames(rna)]
pheatmap(rna[ie.up,], show_colnames = F, scale = "row")

ie.down <- rownames(spatial_deg)[spatial_deg$IE == "down"]
ie.down <- ie.down[ie.down %in% rownames(rna)]
pheatmap(rna[ie.down,], show_colnames = F, scale = "row")

knn1 <- kmeans(t(rna[ie.down,]), 2)
knn1$cluster
lgg_clin$km <- knn1$cluster[match(lgg_clin$patient, names(knn1$cluster))]
fit <- survfit(Surv(Survival..months., Vital.status..1.dead.) ~ km, data = lgg_clin)
print(fit)
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


###
tc.ie.up <- rownames(spatial_deg)[spatial_deg$TCvsIE == "up"]
tc.ie.up <- tc.ie.up[tc.ie.up %in% rownames(rna)]
pheatmap(rna[tc.ie.up,], show_colnames = F, scale = "row")

tc.ie.down <- rownames(spatial_deg)[spatial_deg$TCvsIE == "down"]
tc.ie.down <- tc.ie.down[tc.ie.down %in% rownames(rna)]
pheatmap(rna[tc.ie.down,], show_colnames = F, scale = "row")

knn1 <- kmeans(t(rna[tc.ie.down,]), 2)
knn1$cluster
lgg_clin$km <- knn1$cluster[match(lgg_clin$patient, names(knn1$cluster))]
fit <- survfit(Surv(Survival..months., Vital.status..1.dead.) ~ km, data = lgg_clin)
print(fit)
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


###
head(spatial_deg)
ie.tc.both.up <- spatial_deg$gene[spatial_deg$IE=="up" & spatial_deg$TC == "up"]
ie.tc.both.up <- ie.tc.both.up[ie.tc.both.up%in% rownames(rna)]
pheatmap(rna[ie.tc.both.up,], show_colnames = F, scale = "row")
