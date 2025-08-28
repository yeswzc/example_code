#R/3.6 biowulf

library(cisTopic)
library(Rtsne)
library(densityClust)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

xsp_file = "GBM.LSI.inhouse.PeakCluster.xsp_final.rda"
output_prefix = "GBM.LSI.inhouse.malignant.cisTopic"

if(0){
load(xsp_file)
#sink(paste0(output_prefix,".cisTopic.log"))
print(x.sp);
cat("data loaded\n\n")
#select malignent cells
malignant.clusters = c(0, 2,3,5,7:13, 15, 18:22)
idx = which(x.sp@cluster %in% malignant.clusters)
if(length(idx) > 0){x.sp = x.sp[idx,]; rm(idx);} else{stop("No cells selected...Quit...\n");}
colnames(x.sp@pmat) = as.character(x.sp@peak);
rownames(x.sp@pmat) = x.sp@barcode;
colnames(x.sp@pmat) = gsub("[b']+", "", colnames(x.sp@pmat))
#
cisTopicObject = createcisTopicObject(t(x.sp@pmat), project.name = "GBM")
#
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(2, 5, 10:25, 30, 35, 40), tmp = paste0(output_prefix,"Tmp"),
                                   seed=123, nCores=20, iterations = 500, addModels=FALSE)

save(cisTopicObject, file = paste0(output_prefix, ".cisTopicObject.rda"))
# For WarpLDA
cisTopicObject@calc.params[['runWarpLDAModels']]$seed <- seed
cisTopicObject@calc.params[['runWarpLDAModels']]$iterations <- iterations
cisTopicObject@calc.params[['runWarpLDAModels']]$alpha <- alpha
cisTopicObject@calc.params[['runWarpLDAModels']]$alphaByTopic <- alphaByTopic
cisTopicObject@calc.params[['runWarpLDAModels']]$beta <- beta

png(paste0(output_prefix, ".modelSelection.png"), width = 12, height = 12, units = "in", res = 300)
par(mfrow=c(3,3))
cisTopicObject <- selectModel(cisTopicObject, type='maximum')
cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
cisTopicObject <- selectModel(cisTopicObject, type='derivative')
dev.off())
}
#cisTopicObject@selected.model$topics
#row as topics, columns are peak contribution to the topic


#A. Identification of cell states using the cell-cisTopic distributions
cisTopicObject <- runtSNE(cisTopicObject, target='cell', seed=123, pca=FALSE, method='Probability')
#cisTopicObject <- runUmap(cisTopicObject, target='cell', seed=123, method='Probability')
#cisTopicObject = runPCA(cisTopicObject, target = 'cell', method = "Z-score", seed = 123)

cellassign <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
set.seed(123)
cat("Density cluster based on TSNE\n...")
DR <- Rtsne(t(cellassign), pca=F)
DRdist <- dist(DR$Y)
dclust <- densityClust(DRdist,gaussian=T)
## Distance cutoff calculated to ?
dclust <- findClusters(dclust, rho = 50, delta = 2.5)
# Check thresholds
output_prefix = "test"
png(paste0(output_prefix, ".densityClust.cutoff.png"), width = 7, height = 4, units = "in", res= 300)
options(repr.plot.width=6, repr.plot.height=6)
plot(dclust$rho,dclust$delta,pch=20,cex=0.6,xlab='rho', ylab='delta')
points(dclust$rho[dclust$peaks],dclust$delta[dclust$peaks],col="red",pch=20,cex=0.8)
text(dclust$rho[dclust$peaks]-2,dclust$delta[dclust$peaks]+1.5,labels=dclust$clusters[dclust$peaks])
abline(v=50)
abline(h=2.5)
dev.off()


#Add sample name
sample.name = sapply(cisTopicObject@cell.names, function(x) as.character(unlist(strsplit(x, "\\#"))[1]) )
# Add cluster information
densityClust <- dclust$clusters
meta.data = data.frame(densityClust = densityClust, sample = sample.name)
rownames(meta.data) <- cisTopicObject@cell.names
meta.data[,1] <- as.factor(meta.data[,1])
meta.data$sample = sample(1:10, size = nrow(meta.data), replace = T)
meta.data[,2] <- as.factor(meta.data[,2])

cisTopicObject <- addCellMetadata(cisTopicObject, meta.data);
rm(meta.data,densityClust, dclust, DRdist, DR, cellassign)



#
png(paste0(output_prefix, ".plotFeatures.png"), width = 10, height = 10, units = "in", res = 300)
par(mfrow=c(2,2))
plotFeatures(cisTopicObject, method='tSNE', target='cell', topic_contr=NULL, 
             colorBy=c('nCounts', 'nAcc','densityClust', 'graphBasedClusters_CRA'), 
             cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low="#0a3b70", 
             col.mid='floralwhite', col.high="#760521", intervals=10)

dev.off()

png(paste0(output_prefix, ".topicHeatmap.png"), width = 15, height = 8, units = "in", res= 300)
cellTopicHeatmap(cisTopicObject, method='Z-score', colorBy=c('densityClust', 'sample'), 
                 col.low = "#0a3b70", col.mid = "floralwhite", col.high = "#760521");
dev.off()

##color by topic score
png(paste0(output_prefix, ".topic.tSNE.png"), width = 15, height = 40, units = "in", res= 300)
par(mfrow=c(5,5))
plotFeatures(cisTopicObject, method='tSNE', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
dev.off()




#Enrichment of epigenomic signatures in the cells
#Not applicatable/lack of signature:




### B. Analysis of the regulatory topics
#### Defining topics
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
# These files can be uploaded in IGV or UCSC for visualisation.
dir.create(paste0(output_prefix, '.output/'))
getBigwigFiles(cisTopicObject, path= paste0(output_prefix, '.output/cisTopics_asBW'), seqlengths=seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene))

# the most contributing regions in a topic
pdf(paste0(output_prefix, ".topicsScore.pdf"))
par(mfrow=c(5,4))
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.99, plot=TRUE)
dev.off()
#save  regions sets selected and distributions for each cisTopic 
getBedFiles(cisTopicObject, path= paste0(output_prefix, '.output/cisTopics_asBed'))

#### Topic visualization
cisTopicObject <- runtSNE(cisTopicObject, target='region', perplexity=200, check_duplicates=FALSE)
#plotFeatures(cisTopicObject, method='tSNE', target='region', topic_contr=NULL, colorBy=c('nCells'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)
#par(mfrow=c(2,5))
#plotFeatures(cisTopicObject, method='tSNE', target='region', topic_contr='Z-score', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1')


#### Enrichment of epigenomic signatures
#Required signatures....


#### Annotation of topics to genes and GO terms
cisTopicObject <- annotateRegions(cisTopicObject, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb='org.Hs.eg.db')
#plot region type annotations (promoter, downstream...)
png(paste0(output_prefix, ".topicLocation.png"), width = 10, height = 7, units = "in",res = 300)
signaturesHeatmap(cisTopicObject, selected.signatures = 'annotation')
dev.off()
png(paste0(output_prefix, ".topicLocationtSNE.png"), width = 5, height = 4, units = "in",res = 300)
plotFeatures(cisTopicObject, method='tSNE', target='region', topic_contr=NULL, colorBy=c('annotation'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)
dev.off()

#Identify enriched GO terms per topic,
date()
cisTopicObject <- GREAT(cisTopicObject, genome='hg38', fold_enrichment=2, geneHits=1, sign=0.05, request_interval=10)

#export GREAT 
dir.create(paste0(output_prefix, '.output/rGREAT/'))
.. = lapply(names(cisTopicObject@binarized.rGREAT), function(topic){
    GO.bp = cisTopicObject@binarized.rGREAT[[topic]][['GO Biological Process']]
    write.csv(GO.bp, file = paste0(output_prefix, '.output/rGREAT/',topic,".BP.csv"), row.names = F)
    0;
})
#head(cisTopicObject@binarized.rGREAT$Topic1[['GO Biological Process']])

date()
#topics=c(6,10),
png(paste0(output_prefix,".GREAT.png"), width = 15, height = 7, units = "in", res = 300)
ontologyDotPlot(cisTopicObject, top=5,  var.y='name', order.by='Binom_Adjp_BH')
dev.off()

save(cisTopicObject, file = paste0(output_prefix,".final.cisTopicObject.rda"))

#(Transcription factor) motif enrichment
#hg19 only

#Gene accessibility scores
#already has cicero