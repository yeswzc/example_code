library(ggplot2)
library(reshape2)
source("/Users/wuzhichao/Desktop/multiplot.R")

mytheme <- theme_bw()+theme(legend.position="top",
                            panel.border=element_blank(),
                            panel.grid.major=element_line(linetype="dashed"),
                            panel.grid.minor=element_blank(),
                            plot.title=element_text(size=32,colour="#000000"), 
                            legend.text=element_text(size=32,colour="#000000"),
                            legend.key=element_blank(),
                            legend.key.size = unit(2.5,'cm'),
                            legend.key.height = unit(0.5,'cm'),
                            legend.title = element_blank(),
                            axis.text=element_text(size=32,colour="#000000"),
                            strip.text=element_text(size=32,colour = c("#000000")),
                            strip.background=element_blank()
)
mytheme2 <- theme_bw()+theme(legend.position="top",
                            #panel.border=element_blank(),
                            #panel.grid.major=element_line(linetype="dashed"),
                            panel.grid.minor=element_blank(),
                            plot.title=element_text(size=32,colour="#000000",hjust = 0.5), 
                            legend.text=element_text(size=32,colour="#000000"),
                            legend.key=element_blank(),
                            legend.key.size = unit(2.5,'cm'),
                            legend.key.height = unit(0.5,'cm'),
                            legend.title = element_blank(),
                            axis.text=element_text(size=32,colour="#000000"),
                            axis.title = element_text(size=32,colour="#000000"),
                            #strip.text=element_text(size=16,colour = c("#000000")),
                            #strip.background=element_blank(),
                            
)

mytheme2<-theme_bw()+
  theme(axis.title = element_text(size=32),axis.text = element_text(size=32),legend.position = "none")

pie_theme = mytheme+theme(axis.text=element_blank(),
                        axis.ticks=element_blank(),
                        axis.title=element_blank(),
                        panel.grid.major=element_blank())

colour7 <- c("#cf037c","#a483c4","#f1956c","#40ad26","#6177ab","#0094bd","#de1d2a")
mycolour_7 <- scale_fill_manual(values=colour7)
mycolour_k4 <- scale_fill_manual(values = c("Landrace"="#f1956c","IRRI"="#6177ab","Primary"="#cf037c","Modern"="#a483c4","admix"="#b4b4b4"))
mycolour_k4 <- scale_fill_manual(values = c("#6177ab","#f1956c","#cf037c","#a483c4","#b4b4b4"))


make_haplotype_pie<-function(d1){
  library(grid)
  d1 = d1[d1$N>2,]
  d1$HapID<-as.factor(as.character(d1$HapID))
  d1$POP<-factor(d1$POP,levels=c("IRRI","Landrace","Primary","Modern","admix"))
  d1$Ratio2=NA
  d1$POP2=NA
  j=0
  for(i in levels(d1$POP)){
    j=j+1
    i_index<-which(d1$POP==i)
    hapN<-sum(d1$N[i_index])
    d1$Ratio2[i_index] = d1$N[i_index]/hapN
    #hapN=ceiling(hapN/2) #This is error,
    d1$POP2[i_index]=paste(0,j,d1$POP[i_index],"\nn=",hapN,sep="")
  }
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2,1)))
  vplayout = function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  p1<-ggplot(d1)+geom_bar(aes(x="",y=Ratio,fill=HapID),stat = "identity",width = 3)+ coord_polar(theta="y")+mycolour_7+
    facet_wrap(~POP,nrow=1)+xlab("")+ylab("")+pie_theme+theme(legend.position="bottom")
  print(p1,vp=vplayout(1,1))
  p2<-ggplot(d1)+geom_bar(aes(x="",y=Ratio2,fill=HapID),stat = "identity",width = 3)+ coord_polar(theta="y")+mycolour_7+
    facet_wrap(~POP2,nrow =1)+xlab("")+ylab("")+theme_bw()+pie_theme+theme(legend.position="bottom")
  print(p2,vp=vplayout(2,1))
}

autoY<-function(pheno,my_compa){
  maxY <- max(pheno)
  minY = min(pheno,na.rm = T)
  one<-(maxY-minY)/10
  annY=NULL
  for (i in 1:length(my_compa)){
    j=maxY+one*i
    annY<-c(annY,j)
  }
  return(annY)
}





################################################################################################################
library(ggsignif)

make_box_plot<-function(d1){
  d1<-d1[!is.na(d1$GLWR),]
  d1 = d1[d1$N>2,]
  d1$HAPID<-as.factor(paste("Hap",as.character(d1$HAPID),sep=""))
  d1$POP<-factor(d1$POP,levels=c("IRRI","Landrace","Primary","Modern","admix"))
#ggplot(d1,aes(x=HAPID,y=GLWR,fill=HAPID))+geom_boxplot()+mytheme2+mycolour_7+xlab("")+ylab("GLWR (2018 HN)")

  my_compa<-combn(levels(d1$HAPID),2,simplify = F) ##
  annY<-autoY(d1$GLWR,my_compa)
  p<-ggplot(d1,aes(x=HAPID,y=GLWR,fill=HAPID))+geom_boxplot()+mytheme2+mycolour_7+xlab("")+ylab("GLWR (2018 HN)")+
    geom_signif(comparisons =my_compa,test = "wilcox.test",y_position = annY)+
    scale_x_discrete(position = "top")+
    theme(axis.text.y = element_text(hjust = 0))#+coord_flip()
  #geom_violin(draw_quantiles=c(0.5))
  return(p)
}




make_aov_box<-function(d1,Title,phenotypeID){
  #phenotypeID="PH"
  d1 =d1[!is.na(d1[[phenotypeID]]),]
  d1 = d1[d1$N>2,]
  d1$HAPID<-factor(paste("Hap",as.character(d1$HAPID),sep=""))
  #d1$HAP.MEAN=NA
  HapMean<-rep(NA,length(levels(d1$HAPID)))
  HapLevel=levels(d1$HAPID)
  for(i in 1:length(HapLevel)){
    i_index = which(d1$HAPID == HapLevel[i] )
    meanValue = mean(d1[[phenotypeID]][i_index])
    #d1$HAP.MEAN[i_index]=meanValue
    HapMean[i] = meanValue
  }
  HapMeanOrder = order(HapMean,decreasing = T)
  HapMeanOrder
  HapLevel[HapMeanOrder]
  d1$HAPID = factor(as.character(d1$HAPID),levels = HapLevel[HapMeanOrder] ) #if not coord_flip, use this

  #BR.t<-bartlett.test(GLWR~HAPID,data=d1) #if p-value !< 0.05 or 0.01, 不符合齐性检验
  #KR.t<-kruskal.test(GLWR~HAPID,data=d1)
  AOV <- aov(d1[[phenotypeID]]~d1$HAPID)
  TK.t<-TukeyHSD(AOV) #P adj < 0.05 => significant
  TK.t=data.frame( TK.t$`d1$HAPID`)
  
  Signif.level<-rep(NA,length(HapLevel)-1)
  
  for(i in 2:length(HapMeanOrder)){
    a<-HapMeanOrder[i-1]
    b<-HapMeanOrder[i]
    hap1 = HapLevel[a]
    hap2 = HapLevel[b]
    comp = paste(hap1,"-",hap2,sep="")
    testID = which(rownames(TK.t) == comp)
    if(length(testID)==0){
      comp = paste(hap2,"-",hap1,sep="")
      testID = which(rownames(TK.t) == comp)
    }
    Pvalue <- TK.t$p.adj[testID]
    if(Pvalue < 0.01){
      Signif.level[i-1]="**"
      cat(comp,Pvalue,"**\n")
    }else if(Pvalue < 0.05){
      Signif.level[i-1]="*"
      cat(comp,Pvalue,"*\n")
    }else{
      Signif.level[i-1]="n.s."
      cat(comp,Pvalue,"n.s.\n")
    }
  }
  #mycolour_hap <- scale_fill_manual(values=colour7[HapMeanOrder]) #if not coord_flip
  d1$HAPID = factor(as.character(d1$HAPID),levels = rev(HapLevel[HapMeanOrder]) ) #reverse HapID for good coord_flip()
  mycolour_hap <- scale_fill_manual(values=rev(colour7[HapMeanOrder]))
  #HapLevel=HapLevel[HapMeanOrder]
  Signif.data = data.frame(HAPID=HapLevel[HapMeanOrder][-length(HapLevel)],Signif=Signif.level)
  textY<-max(d1[[phenotypeID]])+mean(d1[[phenotypeID]])/20
  p<-ggplot(d1)+geom_boxplot(aes_string(x='HAPID',y=phenotypeID,fill='HAPID'))+mycolour_hap+labs(title=Title,x="")+
    geom_text(data=Signif.data,aes(x=HAPID,y=textY,label=Signif),size=6)+mytheme2+theme(legend.position = "none")+
    coord_flip()
    
  return(p)
}








