#!/public/agis/ruanjue_group/wuzhichao/software/R/R-3.1.3/bin/R
args<-commandArgs()
if (length(args) < 8) {
	cat("\n\tUsage: Rscript ldHeatmap.R <plink.matrix.ld> <.map or .bim> <SNP> <SNP> <SNP ...>
	Calculate the R2 using: plink --file ped --r2 square --out plink.matrix.ld
	<PeakSNP> 
	.raw.tiff Heatmap highlighted SNPs region with r2>=0.6 against PeakSNP
	.rmDev.tiff, remove edege SNPs with distance >20 Kb with its neighbor
	The 'LDheatmap' package is required in /public/agis/ruanjue_group/wuzhichao/software/R/R-3.1.3/bin/Rscript\n
")
    q(save="no", status=1)
}

library(LDheatmap)

ldResult <- args[6] #plink --r2 square, matrix
map <- args[7]	#plink map file or bim file

SNP <- NULL
for (i in 8:length(args) ){
	SNP <- c(SNP,args[i])
}

ldData <- read.table(ldResult,head=F)
mk <- read.table(map,head=F)
dist <- as.vector(mk[,4])
mk <- as.vector(mk[,2])
colnames(ldData) <- mk #required for SNP.name and etc
row.names(ldData) <- mk

SNP.POS <- which(mk==SNP)

pic <- paste(ldResult,".LDheatmap.naming.tiff",sep="")
tiff(filename=pic,width=12,height=12,units="in",compression="lzw",res=300)
par(omi=c(0.5,0.5,0,0))

MyHeatmap <- LDheatmap(as.matrix(ldData),genetic.distances=dist,flip=TRUE,color=heat.colors(20),SNP.name=c(SNP))

dev.off()
