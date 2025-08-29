library(glmnet)
library(GenomicRanges)
#library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(parallel)


###
promoter.info = "/data/wuz6/project/MEP2020/bin/data/hg19_promoter.txt.gz"


rna = read.table("TCGA.GBM.HiSeqV2.gz", head =T, row.names = 1, check.names = F)
sum(is.na(rna))


epic.ann = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
methy.beta = read.table("TCGA.GBM.HumanMethylation450.gz", head = T, row.names=1, check.names = F)
#
#sort(colSums(is.na(methy.beta)), decreasing=T)[1:20]
# unique(colSums(is.na(methy.beta)))
#remove probes that has NA 
methy.beta = methy.beta[rowSums(is.na(methy.beta))==0,]
#remove probes that removed in EPIC array
methy.beta = methy.beta[rownames(methy.beta) %in% epic.ann$Name, ]
#log2(betai/(1 - betai)).
methy.M = log2(methy.beta/(1-methy.beta))

max(methy.M)
min(methy.M)

###prepare 
use.samples = intersect(colnames(methy.M), colnames(rna))
length(use.samples) #65?
#keep duplicated samples (some patient can have multi-samples) for now
rna = rna[,use.samples]
rna.unique.cnt = apply(rna, 1, function(x) length(unique(x)))
rna = rna[rna.unique.cnt > 20,]

methy.M = methy.M[,use.samples]





###get methylation locations
epic.ann = epic.ann[epic.ann$Name %in% rownames(methy.M),]
probe.range <- GRanges(seqnames = epic.ann$chr, 
                       ranges = IRanges(epic.ann$pos, width = 1), 
                       strand = epic.ann$strand, 
                       name = epic.ann$Name)

###
promoter.info.df = read.table(promoter.info, head =T, sep = "\t")

sum(promoter.info.df$gene.name %in% rownames(rna)) #18115
promoter.info.df = promoter.info.df[promoter.info.df$gene.name %in% rownames(rna), ]

###start from here
set.seed(123)
n.folds = 20
my.folds = sample(rep(seq(n.folds), length = length(use.samples)))
my.fit = vector(mode='list', length= nrow(promoter.info.df))
use.range = 10*1000*1000; #10Mb #in fact it will be 20Mb in total

names(my.fit) = promoter.info.df$gene.name

glmnet.res = mclapply(1:nrow(promoter.info.df), function(k){
  #
  my.gene = promoter.info.df$gene.name[k]
  chr = promoter.info.df$chr[k]
  start = promoter.info.df$start[k]
  end = promoter.info.df$end[k]
  strand =promoter.info.df$strand #this doesn't matter
  #
  gene.range = GRanges(seqnames = Rle(chr),
                  ranges = IRanges(start = start-use.range , start + use.range)) #not using start-, end+
  ovl = findOverlaps(gene.range, probe.range)
  use.probes = probe.range$name[subjectHits(ovl)]
  #
  if(length(use.probes) == 0){return(NULL)}
  
  y = t(rna[my.gene,,drop = F])
  x = t(methy.M[use.probes,])
  cv.fit = cv.glmnet(x, y, keep=TRUE, alpha=0.5, foldid= my.folds)
  ###
  k.lambda.min = which(cv.fit$lambda == cv.fit$lambda.min)
  if(var(cv.fit$fit.preval[,k.lambda.min]) != 0 ) {
    R = cor(cv.fit$fit.preval[, k.lambda.min],y)
    R2 = ifelse(R>0, R^2, 0)
  } else {
    R2 = 0
  }
  ###
  return(list( c(my.gene, R2), cv.fit))
}, mc.cores = 10)
#
names(glmnet.res) = promoter.info.df$gene.name
save(glmnet.res, file = "GBM.cvfit.rda")
###

cv.r2 = lapply(glmnet.res, function(x) x[[1]])
cv.r2 = data.frame(do.call(rbind, cv.r2))
cv.r2[,2] = as.numeric(as.character(cv.r2[,2]))
sum(cv.r2[,2] > 0.4)
#[1] 166