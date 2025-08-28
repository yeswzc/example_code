
#
setwd("/Users/wuzhichao/Project/xoo/GWAS.LD/")
gwa = data.frame(readxl::read_xlsx("../Supplemental Data Set 5_snp.xlsx", skip = 3, col_names =F))
head(gwa)
str(gwa)

xian.specific = unique(gwa[,4][which(gwa[,6] != "-" & gwa[,7] == "-")])
geng.specific = unique(gwa[,4][which(gwa[,6] == "-" & gwa[,7] != "-")])

length(xian.specific)
length(geng.specific)

#sample size:
#701 total
#419 indca
#216 jap
#244 overseas
#136 china_landrace
#319 china_improve
#455 china
f = read.table("../TableS5.snps.P.list.txt.snp.uniq.AF.xls", head = T)
f = f[,c(1:3,6:9)]
head(f)
max(f[,5])
max(f[,7])

f[,4] = ifelse(f[,4]>0.5, 1-f[,4], f[,4])
f[,6] = ifelse(f[,6]>0.5, 1-f[,6], f[,6])
#sapply(list, function)

xian.specific = gsub("Chr", "rs", xian.specific)
geng.specific = gsub("^Chr", "rs", geng.specific)
n.xian = 419
n.geng = 216

res1= lapply(1:length(xian.specific), function(k){
  print(k)
  #k = 238
  i = which(f$SNP == xian.specific[k])
  n10 = round(n.xian*(1-f[i,5]),0)
  n11 = round(n10*f[i,4], 0);
  n12 = n10-n11;
  n20 = round(n.geng*(1-f[i,7]), 0)
  n21 = round(n20*f[i,6], 0);
  n22 = n20-n21
  #if(n22<0) n22 = 0
  
  M = matrix(c(n11, n12, n21, n22), nrow = 2);
  
  rs.fisher = fisher.test(M, alternative = "greater");
  return(c(rs.fisher$p.value, rs.fisher$estimate))
});

res1 = data.frame(do.call(rbind, res1))
colnames(res1)[1] = "p"
res1$snp = xian.specific
head(res1)

max(res1$p)
nrow(res1) #3411
sum(res1$p <0.05); #3189 #94.5%
###在3411 个Xian specifc GWA SNP里面，3189 (94.%) 在籼稻群体里面显著高(fisher.exact test, p<0.05)。

res2= lapply(1:length(geng.specific), function(k){
  print(k)
  #k = 1
  i = which(f$SNP == geng.specific[k])
  n10 = round(n.xian*(1-f[i,5]),0)
  n11 = round(n10*f[i,4], 0);
  n12 = n10-n11;
  n20 = round(n.geng*(1-f[i,7]), 0)
  n21 = round(n20*f[i,6], 0);
  n22 = n20-n21
  #if(n22<0) n22 = 0
  
  M = matrix(c(n11, n12, n21, n22), nrow = 2);
  
  rs.fisher = fisher.test(M, alternative = "less");
  return(c(rs.fisher$p.value, rs.fisher$estimate))
});

res2 = data.frame(do.call(rbind, res2))
res2$snp = geng.specific
head(res2)
colnames(res2)[1] = "p"

nrow(res2) #245
sum(res2$p <0.05) #118 #48.2%
#在297个Geng specific GWA SNP里面 147个(49.5%)在Geng稻群体频率显著高(fisher.exact test, p<0.05)



idx1 = which(f$SNP %in% xian.specific)
idx2 = which(f$SNP %in% geng.specific)



library(ggplot2)

data1 = f[idx1,]
identical(data1$SNP, res1$snp)
data1 = cbind(data1, res1)
data1$sig = as.character(ifelse(data1$p <0.05, "p <0.05", "p >= 0.05"))

p1 = ggplot(data1, aes(x = indca.Allele1.freq, y = jap.Allele1.freq , color = sig)) +
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c("black", "grey"))+
  theme(aspect.ratio = 1)+
  labs(x = "MAF (Xian)", y = "MAF (Geng)")


###
data2 = f[idx2,]
identical(data2$SNP, res2$snp)
data2 = cbind(data2, res2)
data2$sig = as.character(ifelse(data2$p <0.05, "p <0.05", "p >= 0.05"))
p2 = ggplot(data2, aes(x = indca.Allele1.freq, y = jap.Allele1.freq , color = sig)) +
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c("black", "grey"))+
  theme(aspect.ratio = 1)+
  labs(x = "MAF (Xian)", y = "MAF (Geng)")

#install.packages("cowplot")
pdf("XG.specific.GWA.SNP.MAF.test.pdf", width = 8, height = 4, useDingbats = F)
cowplot::plot_grid(p1, p2, nrow = 1)
dev.off()


ggplot(data2, aes(x = indca.missing, y = jap.missing , color = sig)) +
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c("black", "grey"))+
  theme(aspect.ratio = 1)+
  labs(x = "MAF (Xian)", y = "MAF (Geng)")

data1$pop = "Xian"
data2$pop = "Geng"
data = rbind(data1, data2)

#write.table(data[,-11], file = "XG.specific.GWA.SNP.fisher.xls", quote = F, row.names = F)


###
#write.table(unique(c(xian.specific, geng.specific)), file = "XG.specific.gwaSNP.txt", row.names = F, col.names = F, quote = F)
#awk '{print $0"\tXian"}' //analysis/disk5/wuzhichao/01.Project/06.Xoo/03.indca/01.emmax/00.pheno/00.pheno.txt>XG.phenotype.txt
#awk '{if(NR>1) print $0"\tGeng"}' /analysis/disk5/wuzhichao/01.Project/06.Xoo/04.jap/01.emmax/00.pheno/00.pheno.txt >> XG.phenotype.txt
#plink --bfile /analysis/disk5/wuzhichao/01.Project/06.Xoo/01.Geno/public/rice701.resequencing --keep XG.phenotype.txt --extract XG.specific.gwaSNP.txt --recode 12 --out XG.specific
#on server
gwa = data.frame(readxl::read_xlsx("Supplemental Data Set 5_snp.xlsx", skip = 3, col_names =F))
#XG.specific.ped
#XG.specific.fam
#B001 B001 0 0 0 -9 2 2
phenotype = read.table("XG.phenotype.txt", head = T, row.names = 1)
genotype = read.table("XG.specific.ped", head = F, stringsAsFactors = F, row.names = 1)[,-c(1:5)];
x1 = ncol(genotype);x1;
geno1 = genotype[,seq(1, x1, by = 2)]
geno2 = genotype[,seq(2, x1, by = 2)]

genotype = geno1 + geno2; rm(geno1, geno2)
genotype.snp = read.table("XG.specific.map", head =F , stringsAsFactors = F)
colnames(genotype) = genotype.snp[,2]; rm(genotype.snp);
phenotype = phenotype[rownames(genotype),];


head(gwa)
gwa[,1] = as.character(gwa[,1]);
gwa[,4] = gsub("^Chr", "rs", gwa[,4]);



res = lapply(colnames(genotype), function(snp){
  #print(snp)
    #snp = xian.specific[k];
    #snp = colnames(genotype)[66]
  #snp = "rs7_183195"
  #snp = "rs2_18431061"
  #snp = "rs2_9385406"
    index = which(gwa[,4] == snp);
    res = lapply(index, function(i){
        xoo = gwa[i,1];
        data = data.frame(pheno = phenotype[[xoo]], geno = genotype[[snp]], pop = phenotype[["pop"]]);
        data = data[data[,2]!=0 & data[,2] != 3,]
        data$geno = factor(data$geno)
        tt = table(data$pop, data$geno);
        if(length(unique(data$geno)) == 1) return(c(xoo, snp, "fix1"));# fixed genotype at two pop
        if(sum(tt[,1]) < 6 | sum(tt[,2]<6) ) return(c(xoo, snp, "fix2")); #low homo genotyping rate
        if(min(tt) == 0) return(c(xoo, snp, "fix3")); #fixed at at least one population
        test = summary(aov(pheno~geno*pop,data=data));
        xx = grep("geno:pop", rownames(test[[1]]))
        p = test[[1]][["Pr(>F)"]][xx]
        if(length(xx) == 0) p = "NA"
        return(c(xoo, snp, p));
    })
    res = do.call(rbind, res);
    res;
});



res = data.frame(do.call(rbind, res));
colnames(res) = c("xoo", "snp", "pvalue")
write.table(res, file = "XG.specific.noHet.2way.avo.xls", sep = "\t", quote=F, row.names=F)




rm(list=ls())
##fisher test if xian/geng specific GWA SNP was in higher frequency
gwa = data.frame(readxl::read_xlsx("Supplemental Data Set 5_snp.xlsx", skip = 3, col_names =F))
gwa[,1] = as.character(gwa[,1]);
gwa[,4] = gsub("^Chr", "rs", gwa[,4]);
####
xian.specific = unique(gwa[,4][which(gwa[,6] != "-" & gwa[,7] == "-")])
geng.specific = unique(gwa[,4][which(gwa[,6] == "-" & gwa[,7] != "-")])


#XG.specific.ped
#XG.specific.fam
#B001 B001 0 0 0 -9 2 2
phenotype = read.table("XG.phenotype.txt", head = T, row.names = 1)
genotype = read.table("XG.specific.ped", head = F, stringsAsFactors = F, row.names = 1)[,-c(1:5)];


x1 = ncol(genotype);x1;
geno1 = genotype[,seq(1, x1, by = 2)]
geno2 = genotype[,seq(2, x1, by = 2)]
#genotype = rbind(geno1, geno2);

genotype.snp = read.table("XG.specific.map", head =F , stringsAsFactors = F)
colnames(geno1) = genotype.snp[,2]; 
colnames(geno2) = genotype.snp[,2]; 
rm(genotype.snp, genotype);
geno.xg = as.character(phenotype$pop[match(rownames(geno1), rownames(phenotype))])

res1 = lapply(xian.specific, function(snp){
  print(snp)
    #snp = "rs1_10323754"
    index1 = which(colnames(geno1) == snp);
    data = matrix(c(geno1[,index1], geno2[,index1], geno.xg), ncol = 3);
    data = data[data[,1] == data[,2] & data[,1] != 0,];
    data = matrix(c(data[,1], data[,2], data[,3], data[,3]), ncol = 2);
    data = data[data[,1] != 0,]
    M = table(data[,1], data[,2]);
    M = M[, c("Geng", "Xian"), drop = F];
    #print(M)
    if(nrow(M)!= 2) return(c(snp, NA, NA, 0, sum(M[,1]), 0, sum(M[,2]) ) )
    M = apply(M, 2, sort);
    #if(ncol(M) != 2) return(c(snp, NA, M[1,2], sum(M[,1]), M[1,2], sum(M[,2])))
    

    rs.fisher = fisher.test(M, alternative = "less");
    return(c(snp, rs.fisher$p.value, rs.fisher$estimate, M[1,2], sum(M[,1]), M[1,2], sum(M[,2])));
})
res1 = data.frame(do.call(rbind, res1))
res1$specific = "xian"

res2 = lapply(geng.specific, function(snp){
  print(snp)
  #snp ="rs1_19157735" 
  index1 = which(colnames(geno1) == snp);
  data = matrix(c(geno1[,index1], geno2[,index1], geno.xg), ncol = 3);
  data = data[data[,1] == data[,2] & data[,1] != 0,];
  data = matrix(c(data[,1], data[,2], data[,3], data[,3]), ncol = 2);
  data = data[data[,1] != 0,]
  M = table(data[,1], data[,2]);
  M = M[, c("Geng", "Xian"), drop = F];
  #print(M)
  if(nrow(M)!= 2) return(c(snp, NA, NA, 0, sum(M[,1]),0, sum(M[,2]) ) )
  M = apply(M, 2, sort);
  #if(ncol(M) != 2) return(c(snp, NA, M[1,2], sum(M[,1]), M[1,2], sum(M[,2])))
  
  
  rs.fisher = fisher.test(M, alternative = "greater");
  return(c(snp, rs.fisher$p.value, rs.fisher$estimate, M[1,2], sum(M[,1]), M[1,2], sum(M[,2])));
})
res2 = data.frame(do.call(rbind, res2))
res2$specific = "geng"

res = rbind(res1, res2);
colnames(res) = c("SNP", "Pvalue", "OddsRatio", "Geng.minor.cnt", "Geng.major.cnt", "Xian.minor.cnt", "Xian.major.cnt", "XGspecific")
write.table(res, file = "XG.specific.GWA.SNP.fisher.noHet.xls", quote = F, row.names = F)

