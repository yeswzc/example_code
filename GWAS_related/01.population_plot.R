setwd("/Users/wuzhichao/Documents/Project/01.Xu300Resequencing/04.populationStructure/04.LD")
d1<-read.table("01.IRRI.PLINK.ld.decay",head=T)
d2<-read.table("02.Landrace.PLINK.ld.decay",head=T)
d3<-read.table("03.Primary.PLINK.ld.decay",head=T)
d4<-read.table("04.Modern.PLINK.ld.decay",head=T)
d5<-read.table("05.whole.PLINK.ld.decay",head=T)
d6<-read.table("06.nonAdmix.PLINK.ld.decay",head=T)


d1$Population="IRRI"
d2$Population="Landrace"
d3$Population="Primary"
d4$Population="Modern"
d5$Population="Whole"
d6$Population="nonAdmix"

d<-rbind(d1,d2,d3,d4,d5,d6)
d$Population<-as.factor(d$Population)
library(ggplot2)
#pdf("LDdecay.pdf",width = 6,height = 5)
ggplot(data=d,aes(x=Distance.kb.,y=LDdecay))+geom_line(aes(colour=Population))+
   scale_color_manual(values = c("Landrace"="#f1956c","IRRI"="#6177ab","Primary"="#cf037c","Modern"="#a483c4",
                                 "Whole"="black","nonAdmix"="#40ad26"))+
  theme_bw()+theme(legend.position="left",legend.title = element_blank(),panel.border =element_rect(size=0.1),
                   axis.text = element_text(size=20),axis.title = element_text(size=25),legend.text = element_text(size=20))+
  xlab("Distance (Kb)")+labs( y=expression( r^2) ) 
#dev.off()



##PCA
library(ggplot2)
setwd("/Users/wuzhichao/Documents/Project/01.Xu300Resequencing/04.populationStructure/01.PCA")
d<-read.csv("01.rice298.allSNP.PCA.evec.csv",head=F,comment.char = "#")
head(d)
pdf("01.rice298.allSNP.PCA.evec.PH.revise.pdf",width = 6,height = 5)
ggplot(d)+geom_point(aes(x=V2,y=V3,color=V12))+
  scale_color_manual(values =c("admix"="#757575","Landrace"="#f1956c","IRRI"="#6177ab","Primary"="#cf037c","Modern"="#a483c4"))+
  theme_bw()+xlab("PC 1")+ylab("PC 2")+
  theme(legend.position = "none",axis.title = element_text(size=25),panel.grid = element_blank())

ggplot(d)+geom_point(aes(x=V4,y=V5,color=V12))+
  scale_color_manual(values =c("admix"="#757575","Landrace"="#f1956c","IRRI"="#6177ab","Primary"="#cf037c","Modern"="#a483c4"))+
  theme_bw()+xlab("PC 3")+ylab("PC 4")+
  theme(legend.position = "none",axis.title = element_text(size=25),panel.grid = element_blank())

ggplot(d)+geom_point(aes(x=V5,y=V6,color=V12))+
  scale_color_manual(values =c("admix"="#757575","Landrace"="#f1956c","IRRI"="#6177ab","Primary"="#cf037c","Modern"="#a483c4"))+
  theme_bw()+xlab("PC 4")+ylab("PC 5")+
  theme(legend.position = "none",axis.title = element_text(size=25),panel.grid = element_blank())
dev.off()


#tSNE test
setwd("/Users/wuzhichao/Documents/Project/01.Xu300Resequencing/04.populationStructure/01.PCA/")
library(Rtsne) # Load package
d<-read.csv("01.rice298.allSNP.PCA.evec.csv",head=T)
set.seed(42) # Sets seed for reproducibility
tsne_out <- Rtsne(as.matrix(d[,2:4]),pca = F,theta = 0,max_iter = 2500) # Run TSNE
plot(tsne_out$Y,col=d$Pop) # Plot the result






##Tree and admimxture were draw using iTol
setwd("/Users/wuzhichao/Documents/Project/01.Xu300Resequencing/04.populationStructure/02.admixture/")
d<-read.csv("01.seed.1234567.CV.error",head=F,sep=" ")
head(d)
d$size<-16
d$size[9]=15 #"#6177ab" #K8 lowest
d$size[5]=17 #"#f1956c" #K4
d$size<-as.numeric(d$size)
pdf("01.seed.1234567.CV.error.pdf",width = 7,height = 5)
ggplot(data=d,aes(x=V1,y=V2))+geom_line()+geom_point(shape=d$size,size=3)+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12),expand = c(0,0.5))+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid = element_blank(),
                   axis.line = element_line(colour = "black"))+
  xlab("K")+ylab("CV Error")
#+geom_vline(xintercept = 4,colour="grey",linetype="dashed")
#+geom_rect(aes(ymin=0.475,ymax=0.495,xmin=3.8,xmax=4.2),fill="grey")
  
dev.off()  

#landrace	#f1956c
#IRRI	#6177ab
#Primary	#5a6727
#Modern	#a483c4
#mid-term #40ad26
#restor #0094bd
#other1 #450044

rm(d)
d<-read.csv("01.rice.2986.ADMIXTURE.run20.clumpp.csv",head=T)
pdf("01.rice.2986.ADMIXTURE.run20.clumpp.pdf",width=10,height = 4)
par(mfrow=c(6,1),mar=c(0,1,0,0),mgp=c(0,1,0))
barplot(t(d[,3:4]),col =c("#6177ab","#f1956c") ,border = NA,space = 0,yaxt="n",xaxt="n",cex.lab =1,xaxs="i",ylab="K2")
barplot(t(d[,5:7]),col =c("#6177ab","#f1956c","#a483c4") ,border = NA,space = 0,yaxt="n",xaxt="n",cex.lab =1,xaxs="i",ylab="K3")
barplot(t(d[,8:11]),col =c("#f1956c","#5a6727","#a483c4","#6177ab") ,border = NA,space = 0,yaxt="n",xaxt="n",cex.lab =1,xaxs="i",ylab="K4")
barplot(t(d[,12:16]),col =c("#40ad26","#a483c4","#5a6727","#f1956c","#6177ab") ,border = NA,space = 0,yaxt="n",xaxt="n",cex.lab =1,xaxs="i",ylab="K5")
barplot(t(d[,17:22]),col =c("#0094bd","#5a6727","#a483c4","#40ad26","#f1956c","#6177ab") ,border = NA,space = 0,yaxt="n",xaxt="n",cex.lab =1,xaxs="i",ylab="K6")
barplot(t(d[,23:29]),col =c("#6177ab","#5a6727","#f1956c","#0094bd","#a483c4","#450044","#40ad26"),
        border = NA,space = 0,yaxt="n",xaxt="n",cex.lab =1,xaxs="i",ylab="K7")
dev.off()

#====================Tree====================


library(ape)

colD<-read.csv(
  "/Users/wuzhichao/Documents/Project/01.Xu300Resequencing/04.populationStructure/02.admixture/02.0.ADMIXTURE.csv",head=T)

edgeColor<-function(x){
  #tree<-tre
  tree<-x
  edge<-tree$edge
  edge<-cbind(edge,"grey",0.2,0.2)
  
  for (i in 1:length(colD[,1])) {
    #i=1
    ind<-colD[i,1]
    ind<-sub("IRIS_313-","IRIS",ind)
    j<-which.edge(tree,ind)
    edge[j,3]<-as.character(colD$Colour[i])
    cc<-paste(ind,edge[j,3],sep=" ")
    #print(cc)
    edge[j,4]=2
    edge[j,5]=2
  }
  return(edge)
}


setwd("/Users/wuzhichao/Documents/Project/01.Xu300Resequencing/04.populationStructure/")
tre<-read.tree("01.randomSNP.fastTree.nwk")
edgeNJ<-edgeColor(tre)


plot(tre,type="p",show.tip.label=T,use.edge.length=T,edge.color=edgeNJ[,3],edge.width=1) #tip.color=tip[,2],

tre<-read.tree("01.randomSNP.fastTree.GTR.nwk")
edgeNJ<-edgeColor(tre)
plot(tre,type="p",show.tip.label=T,use.edge.length=T,edge.color=edgeNJ[,3],edge.width=1) #tip.color=tip[,2],


