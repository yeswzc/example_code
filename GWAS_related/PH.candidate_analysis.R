##PH
#install.packages("ggsignif")
setwd("/Users/wuzhichao/Project/01.Xu300Resequencing/07.GWAS.refR498/PH/chr12.24M")
source("/Users/wuzhichao/Desktop/multiplot.R")
source("/Users/wuzhichao/Project/01.Xu300Resequencing/07.GWAS.refR498/haplotype_analysis.R")
library(genemodel)
source("/Users/wuzhichao/Desktop/genemodel.R")


#####chr12.24M===========================================================================================================
#========================================================================================================================
d1<-read.table("OsR498G1222051500.OsIAA3.hap.stat.txt",head=T)
pdf("OsR498G1222051500.OsIAA3.hap.stat.pie.pdf",width = 10,height = 5)
make_haplotype_pie(d1)
dev.off()

d1<-read.table("OsR498G1222051500.OsIAA3.hap.PX2017.txt",head=T)
#pdf("OsR498G1222051500.OsIAA3.hap.PX2017.pdf",width = 5,height = 5)
make_aov_box(d1,"PX2017","PH")
dev.off()

par(mfrow=c(2,1))
#OsR498G1222051500.01.T01
OsR498G1222051500.01.T01 = data.frame( type=c("5' utr","coding_region","coding_region","coding_region","3' utr","exon","exon","exon","intron","intron"),
                                       coordinates=c("1-130","131-282","381-743","1040-1118","1119-2018","1-282","381-743","1040-2018","283-380","744-1039") )

#OsR498G1222051500.01.T02
OsR498G1222051500.01.T02 = data.frame( type=c("5' utr","coding_region","coding_region","coding_region","exon","exon","exon","intron","intron"),
                                       coordinates=c("1-130","131-282","339-743","1040-1118","1-282","339-743","1040-1118","283-338","744-1039") )
pdf("OsR498G1222051500.OsIAA3.genemodel.pdf",width = 10,height = 5)
genemodel.plot(model=OsR498G1222051500.01.T01,start=24271170, bpstop=24273187, orientation="forward",xaxis=T)
mutation.plot(24271523,24271523,text = "A/G")
genemodel.plot(model=OsR498G1222051500.01.T02,start=24271170, bpstop=24273187, orientation="forward",xaxis=F)
mutation.plot(24271523,24271523,text = "A/G")


#candidate 2
d1<-read.table("OsR498G1222052600.bZIP.hap.stat.txt",head=T)
pdf("OsR498G1222052600.bZIP.hap.stat.pie.pdf",width = 10,height = 5)
make_haplotype_pie(d1)
dev.off()

d1<-read.table("OsR498G1222052600.bZIP.hap.PX2017.txt",head=T)
d1<-d1[d1$HET!="HET",]
pdf("OsR498G1222052600.bZIP.hap.PX2017.pdf",width = 5,height = 5)
make_aov_box(d1,"PX2017","PH")
dev.off()

d1<-read.table("OsR498G1222052600.bZIP.hap.HN2017.txt",head=T)
d1<-d1[d1$HET!="HET",]
pdf("OsR498G1222052600.bZIP.hap.PX2017.pdf",width = 5,height = 5)
make_aov_box(d1,"HN2017","PH")
dev.off()


#####chr12.21M===========================================================================================================
#========================================================================================================================
setwd("/Users/wuzhichao/Project/01.Xu300Resequencing/07.GWAS.refR498/PH/chr12.21M/")
d1<-read.table("OsR498G1221868900.Fbox.hap.stat.txt",head=T)
pdf("OsR498G1221868900.Fbox.hap.stat.pie.pdf",width = 10,height = 5)
make_haplotype_pie(d1)
dev.off()
d1<-read.table("OsR498G1221868900.Fbox.hap.HN2018.txt",head=T)
make_aov_box(d1,"HN2018")

#####chr12.16M===========================================================================================================
#========================================================================================================================
setwd("/Users/wuzhichao/Project/01.Xu300Resequencing/07.GWAS.refR498/PH/chr12.16M/")
d1<-read.table("OsR498G1221568200.DEC.hap.stat.txt",head=T)
pdf("OsR498G1221568200.DEC.hap.stat.pdf",width = 10,height = 5)
make_haplotype_pie(d1)
dev.off()
d1<-read.table("OsR498G1221568200.DEC.hap.HN2017.txt",head=T)
d1 <- d1[d1$HET!="HET",]
pdf("OsR498G1221568200.DEC.hap.HN2017.nonHET.pdf",width = 5,height = 5)
make_aov_box(d1,"HN2017")
dev.off()

#OsR498G1221568200.01.T01
OsR498G1221568200.01.T01 = data.frame( type=c("5' utr","coding_region","coding_region","3' utr","exon","exon","intron"),
                                       coordinates=c("4771-5562","415-1027","4106-4770","1-414","1-1027","4106-5562","1028-4105") )
pdf("OsR498G1221568200.DEC.genemodel.pdf",width = 10,height = 3)
genemodel.plot(model=OsR498G1221568200.01.T01,start=16273747, bpstop=16279308, orientation="reverse",xaxis=T)
mutation.plot(16274608,16274608,text = "C/T")
dev.off()


###chr11.22M=========================================================================================================================
###==================================================================================================================================
setwd("/Users/wuzhichao/Project/01.Xu300Resequencing/07.GWAS.refR498/PH/chr11.22M/")
d1<-read.table("OsR498G1120343900.Pish.hap.stat.txt",head=T)
pdf("OsR498G1120343900.Pish.hap.stat.pie.pdf",width = 10,height = 5)
make_haplotype_pie(d1)
dev.off()

d1<-read.table("OsR498G1120343900.Pish.hap.HN2018.txt",head=T)
pdf("OsR498G1120343900.Pish.hap.HN2018.pdf",width = 5,height = 5)
make_aov_box(d1,"HN2018")
dev.off()
#####
d1<-read.table("OsR498G1120341600.OsWAK.hap.stat.txt",head=T)
pdf("OsR498G1120341600.OsWAK.hap.stat.pdf",width = 10,height = 5)
make_haplotype_pie(d1)
dev.off()
d1<-read.table("OsR498G1120341600.OsWAK.hap.HN2018.txt",head=T)
pdf("OsR498G1120341600.OsWAK.hap.HN2018.pdf",width = 5,height = 5)
make_aov_box(d1,"HN2018")
dev.off()
#####
d1<-read.table("OsR498G1120297500.c1.hap.stat.txt",head=T)
make_haplotype_pie(d1) #too many haplotypes
d1<-read.table("OsR498G1120297500.c1.hap.HN2018.txt",head=T)
d1<-d1[d1$HET!="HET",]
make_aov_box(d1,"HN2018") #no significance,even removed heterozygous accessions
#####
d1<-read.table("OsR498G1120330600.c2.hap.stat.txt",header = T)
pdf("OsR498G1120330600.c2.hap.stat.pdf",width = 10,height = 5)
make_haplotype_pie(d1)
dev.off()
d1<-read.table("OsR498G1120330600.c2.hap.PX2017.txt",head=T)
pdf("OsR498G1120330600.c2.hap.PX2017.pdf",width = 5,height = 5)
make_aov_box(d1,"PX2017")
dev.off()
#####
d1<-read.table("OsR498G1120331000.c3_maf.hap.stat.txt",head=T)
pdf("OsR498G1120331000.c3_maf.hap.stat.pie.pdf",width = 10,height = 5)
make_haplotype_pie(d1)
dev.off()
d1<-read.table("OsR498G1120331000.c3_maf.hap.PX2017.txt",head=T)
pdf("OsR498G1120331000.c3_maf.hap.PX2017.pdf",width = 5,height = 5)
make_aov_box(d1,"PX2017")
dev.off()

###chr11.2.3M========================================================================================================================
###==================================================================================================================================
setwd("/Users/wuzhichao/Project/01.Xu300Resequencing/07.GWAS.refR498/PH/chr11_2M//")
d1<-read.table("OsR498G111942790.ThrProt.hap.stat.txt",head=T)
pdf("OsR498G111942790.ThrProt.hap.stat.pdf",width = 10,height = 5)
make_haplotype_pie(d1)
dev.off()
d1<-read.table("OsR498G111942790.ThrProt.hap.HN2018.txt",head=T)
pdf("OsR498G111942790.ThrProt.hap.HN2018.txt",height = 5,width = 5)
make_aov_box(d1,"HN2018")
dev.off()


###chr1.26M==========================================================================================================================
###==================================================================================================================================
setwd("/Users/wuzhichao/Project/01.Xu300Resequencing/07.GWAS.refR498/PH/chr1.26M/")
d1<-read.table("OsR498G0101599100.CBS.hap.stat.txt",head=T)
pdf("OsR498G0101599100.CBS.hap.stat.pdf",width = 10,height = 5)
make_haplotype_pie(d1)
dev.off()
d1<-read.table("OsR498G0101599100.CBS.hap.PX2017.txt",head=T)
pdf("OsR498G0101599100.CBS.hap.PX2017.pdf",width = 5,height = 5)
make_aov_box(d1,"PX2017")
dev.off()
#####
d1<-read.table("OsR498G0101597100.hap.stat.txt",head=T)
pdf("OsR498G0101597100.hap.stat.pdf",width = 10,height = 5)
make_haplotype_pie(d1)
dev.off()
d1<-read.table("OsR498G0101597100.hap.PX2017.txt",head=T)
pdf("OsR498G0101597100.hap.PX2017.pdf",width = 5,height = 5)
make_aov_box(d1,"PX2017")
dev.off()



###chr1.26M==========================================================================================================================
###==================================================================================================================================
setwd("/Users/wuzhichao/Project/01.Xu300Resequencing/07.GWAS.refR498/PH/chr1.13M/")
d1<-read.table("OsR498G0100861300.H.hap.stat.txt",head=T)
pdf("OsR498G0100861300.H.hap.stat.pie.pdf",width = 10,height = 5)
make_haplotype_pie(d1)
dev.off()
d1<-read.table("OsR498G0100861300.H.hap.HN2018.txt",head=T)
make_aov_box(d1,"HN2018")
####
d1<-read.table("OsR498G0100849000.GX2.hap.stat.txt",head=T)
pdf("OsR498G0100849000.GX2.hap.stat.pie.pdf",width = 10,height = 5)
make_haplotype_pie(d1)
dev.off()
d1<-read.table("OsR498G0100849000.GX2.hap.HN2018.txt",head=T)
pdf("OsR498G0100849000.GX2.hap.HN2018.pdf",width = 5,height = 5)
make_aov_box(d1,"HN2018")
dev.off()
###########
d1<-read.table("OsR498G0100882700.c18ABA.hap.stat.txt",head=T)
pdf("OsR498G0100882700.c18ABA.hap.stat.pie.pdf",width = 10,height = 5)
make_haplotype_pie(d1)
dev.off()
d1<-read.table("OsR498G0100882700.c18ABA.hap.HN2018.txt",head=T)
pdf("OsR498G0100882700.c18ABA.hap.HN2018.pdf",width = 5,height = 5)
make_aov_box(d1,"HN2018")
dev.off()
###########
d1<-read.table("OsR498G0100882300.c17AMP.hap.stat.txt",head=T)
pdf("OsR498G0100882300.c17AMP.hap.stat.pie.pdf",width = 10,height = 5)
make_haplotype_pie(d1)
dev.off()
d1<-read.table("OsR498G0100882300.c17AMP.hap.HN2018.txt",head=T)
make_aov_box(d1,"HN2018")
