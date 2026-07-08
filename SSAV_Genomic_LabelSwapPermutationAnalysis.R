library(ggplot2)#I use this
library(plotrix) #for like confidence intervals
library(boot)#for bootstrapping
library(lme4)#I use this

library(ggpubr)#??
library(tidyverse)
library(rapport)
library(vcd)
library(grid)##???
library(rstatix)# I think I do use this

library(scales)
library(readr)#reading excel files
#library(tidymodels)###???



#--------IDENTIFY GENE PROPERTIES OF SAMPLE PERMUTATION SWAP------#


#permutation_groups=combn(c(1,2,3,4,5,6),3,simplify=FALSE)

###---G130----####
##IMPORT DATA####
#Import the file with gene features
SBGE.SSS.total.data <- read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SBGEandSSSdataForMBE.csv", header = TRUE)
SBGE.table<-SBGE.SSS.total.data[,-1] #SBGE.table is just abbreivated from the original file name

#Important file output from bedtools intersect that maps all SNPs present at generation G130 (coverage 40 to 500) to corresponding genes
#this contains multiple entries per SNP (entries are 'per bedtools feature')
backgroundFeatures=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130a_cov40-500v2_Background_genes_short.txt", sep=",", header = TRUE)
backgroundFeatures<-backgroundFeatures[-c(1)]

##add genes and information present in backgroundFeatures to SBGE.table 
overlap.geneIDX=SBGE.table[[1]]%in%backgroundFeatures[["geneID"]] 
#obtain a list where chromosome information is missing for a gene
missing.chrm.IDX=is.na(SBGE.table[,"chrm"]) 
helpIDX=overlap.geneIDX&missing.chrm.IDX 
help.gene.list=SBGE.table[helpIDX,"geneID"] 
#pull out the gene&chrm pairs for genes that appear in BOTH SBGE.table and backgroundFeatures, but where chromosome is missing
chrm.list=unique(subset(backgroundFeatures[,c("geneID","chrm")],backgroundFeatures[,"geneID"]%in%help.gene.list))
SBGE.subset=SBGE.table[helpIDX,c("geneID","chrm")]
merged.chrms=merge(chrm.list,SBGE.subset,by="geneID",all.x=FALSE)
SBGE.table[helpIDX,"chrm"]=merged.chrms[[2]]# 

#####Import the permutation file (NOTE ONLY INCLUDE SNPS ON CHR 2)####
##notes 37,084 SNPs (FDR < 0.05), corresponding to 2127 genes AFTER FILTERING for MAF
##--RAW----> 19094 SNPs occur in genes/map to genes summed across all
#Perm 1= swap 126   "sig" SNPs= 3034/175 906       |
#Perm 2= swap 134  "sig" SNPs=  3051/175 906       
#Perm 3= swap 135  "sig" SNPs= 5219/175 906        |

#Perm 4= swap 126   "sig" SNPs= 5900/175 906       |1113 genes
#Perm 5= swap 134  "sig" SNPs=  5199/175 906       |1038 genes
#Perm 6= swap 135  "sig" SNPs= 2832/175 906        |698 genes
#Perm 7= swap 136. "sig" SNPs= 2760/175906        |667 genes
#Perm 8= swap 145# "sig" SNPs= 2790/175906        |691 genes
#Perm 9= swap 146# "sig" SNPs= 2743/175906        |653 genes
#Perm 10=swap 156# "sig" SNPs= 4487/175906        |963 genes

##--Segregating in 4 populations
#Perm 4= swap 126   "sig" SNPs= 990       | genes
#Perm 5= swap 134  "sig" SNPs=  633      | genes
#Perm 6= swap 135  "sig" SNPs=  54       | genes
#Perm 7= swap 136. "sig" SNPs=  44       | genes
#Perm 8= swap 145# "sig" SNPs=  54       | genes
#Perm 9= swap 146# "sig" SNPs=  55       | genes
#Perm 10=swap 156# "sig" SNPs=  637       | genes

#--Seg4 AND MAF5 filter-->only 1470 across ALL permutations map to SNPs
#Perm 4= swap 126   "sig" SNPs= 834      |312 genes
#Perm 5= swap 134  "sig" SNPs=  503      |242 genes
#Perm 6= swap 135  "sig" SNPs=  33       |23 genes
#Perm 7= swap 136. "sig" SNPs=  33       |21 genes
#Perm 8= swap 145# "sig" SNPs=  37       |18 genes
#Perm 9= swap 146# "sig" SNPs=  32       |16 genes
#Perm 10=swap 156# "sig" SNPs=  507      |209 genes

LableSwap.SNPS=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130_permutated_sample1_3newCMHchr2.csv.",sep=",", header = TRUE)

# 
# #note SNPs are derrived from label swap on "SSAV_G130_MAF_cleaned_allFreq.csv"
# LableSwap.SNPS_part1=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130_permutated_sample4_6newCMHchr2.csv.",sep=",", header = TRUE)
# #remove first empty row and column X (not sure what it is)
# LableSwap.SNPS_part1=LableSwap.SNPS_part1[-1,-1]
# 
# LableSwap.SNPS_part2=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130_permutated_sample7_10newCMHchr2.csv.",sep=",", header = TRUE)
# #remove extra column and row
# LableSwap.SNPS_part2=LableSwap.SNPS_part2[-1,-1]
# #p.adjust should be FDR
# LableSwap.SNPS=rbind(LableSwap.SNPS,LableSwap.SNPS_part1,LableSwap.SNPS_part2)
# #setting it this way to save some memory
# 
# LableSwap.SNPS_part1=0

###
##calc exact # sig SNPs per Perm Swap (without ANY FILTERING)
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 1","p.adjust"]<0.05)#
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 2","p.adjust"]<0.05)#
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 3","p.adjust"]<0.05)# 
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 4","p.adjust"]<0.05)#
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 5","p.adjust"]<0.05)#
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 6","p.adjust"]<0.05)# 
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 7","p.adjust"]<0.05)#
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 8","p.adjust"]<0.05)#
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 9","p.adjust"]<0.05)# 
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 10","p.adjust"]<0.05)#
#length(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 7","p.adjust"]<0.05)#175906

# ##Examining overlap of SNPs between permutations sets
# P45s=(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 4","p.adjust"]<0.05)&(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 5","p.adjust"]<0.05)
# #overlap SNPs=3042->out of 5900 and 5199
# P78s=(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 7","p.adjust"]<0.05)&(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 8","p.adjust"]<0.05)
# #overlap SNPs=2058->out of 2760 and 2790
# P910s=(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 9","p.adjust"]<0.05)&(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 10","p.adjust"]<0.05)
# #overlap SNPs=2134

# ##--determining that there are a specific set of ~1800 SNPs that appear in all dataset.
# sum(P45s&P78s)#1837
# sum(P45s&P78s&P910s)#1829

##------MAF and PopSegregating Filtering----####

#imports table containg allele frequencies for each SNP in each sample type for each population
G130.MAF=read.table("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130_MAF_cleaned_MAF1_Freq.csv",header=TRUE)
G130.chr2=subset(G130.MAF,chrm=="2L"|chrm=="2R")
#We filter out SNPs that segregate in less than half the populations
#this calculates how many populations a SNP segrates in
G130.MAF$Seg.numR=rep(0,length(G130.MAF$ID))#number of Red populations where a variant is segregating
for(col in seq(5,15,2)){#this iterates over the frequency columns with Red Sample types
  #we only look at red sample types because a SNP can't be segregating in N without also segregating in R
  seg.IDX=G130.MAF[,col]!=1
  G130.MAF[seg.IDX,"Seg.numR"]=G130.MAF[seg.IDX,"Seg.numR"]+1}

LableSwap.SNPS$PopSeg.num4=FALSE
LableSwap.SNPS[LableSwap.SNPS$ID%in%G130.MAF[G130.MAF$Seg.numR>3,"ID"],"PopSeg.num4"]=TRUE

#Report just the SegPop # filter
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 4"&LableSwap.SNPS$PopSeg.num4,"p.adjust"]<0.05)#990
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 5"&LableSwap.SNPS$PopSeg.num4,"p.adjust"]<0.05)#633
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 6"&LableSwap.SNPS$PopSeg.num4,"p.adjust"]<0.05)#54
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 7"&LableSwap.SNPS$PopSeg.num4,"p.adjust"]<0.05)#44
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 8"&LableSwap.SNPS$PopSeg.num4,"p.adjust"]<0.05)#54
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 9"&LableSwap.SNPS$PopSeg.num4,"p.adjust"]<0.05)#55
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 10"&LableSwap.SNPS$PopSeg.num4,"p.adjust"]<0.05)#637

###MAF Filter: Use NR.avg
LableSwap.SNPS$MAF5=FALSE
LableSwap.SNPS[LableSwap.SNPS$ID%in%G130.MAF[G130.MAF$MAF.NRavg>0.05,"ID"],"MAF5"]=TRUE

sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 4"&LableSwap.SNPS$PopSeg.num4&LableSwap.SNPS$MAF5,"p.adjust"]<0.05)#834
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 5"&LableSwap.SNPS$PopSeg.num4&LableSwap.SNPS$MAF5,"p.adjust"]<0.05)#834
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 6"&LableSwap.SNPS$PopSeg.num4&LableSwap.SNPS$MAF5,"p.adjust"]<0.05)#834
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 7"&LableSwap.SNPS$PopSeg.num4&LableSwap.SNPS$MAF5,"p.adjust"]<0.05)#834
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 8"&LableSwap.SNPS$PopSeg.num4&LableSwap.SNPS$MAF5,"p.adjust"]<0.05)#834
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 9"&LableSwap.SNPS$PopSeg.num4&LableSwap.SNPS$MAF5,"p.adjust"]<0.05)#834
sum(LableSwap.SNPS[LableSwap.SNPS$PermutationSet=="Perm 10"&LableSwap.SNPS$PopSeg.num4&LableSwap.SNPS$MAF5,"p.adjust"]<0.05)#834


####======GENE LEVEL ANALYSIS======####

#this ignores different feature types for simplicity
LabelSwap.genes=merge(LableSwap.SNPS,unique(backgroundFeatures[,c("ID","geneID")]),by="ID")
LabelSwap.genes$Candidate=LabelSwap.genes$p.adjust<0.05&LabelSwap.genes$MAF5&LabelSwap.genes$PopSeg.num4#19094 (sum all)-->with MAF=4226--->with MAF AND Seg4 1470

###-----Put the Perm Candidates in TABLE----
SBGE.table$Perm4.Can=SBGE.table$geneID%in%LabelSwap.genes[LabelSwap.genes$PermutationSet=="Perm 4"&LabelSwap.genes$Candidate==TRUE,"geneID"]
SBGE.table$Perm5.Can=SBGE.table$geneID%in%LabelSwap.genes[LabelSwap.genes$PermutationSet=="Perm 5"&LabelSwap.genes$Candidate==TRUE,"geneID"]
SBGE.table$Perm6.Can=SBGE.table$geneID%in%LabelSwap.genes[LabelSwap.genes$PermutationSet=="Perm 6"&LabelSwap.genes$Candidate==TRUE,"geneID"]
SBGE.table$Perm7.Can=SBGE.table$geneID%in%LabelSwap.genes[LabelSwap.genes$PermutationSet=="Perm 7"&LabelSwap.genes$Candidate==TRUE,"geneID"]
SBGE.table$Perm8.Can=SBGE.table$geneID%in%LabelSwap.genes[LabelSwap.genes$PermutationSet=="Perm 8"&LabelSwap.genes$Candidate==TRUE,"geneID"]
SBGE.table$Perm9.Can=SBGE.table$geneID%in%LabelSwap.genes[LabelSwap.genes$PermutationSet=="Perm 9"&LabelSwap.genes$Candidate==TRUE,"geneID"]
SBGE.table$Perm10.Can=SBGE.table$geneID%in%LabelSwap.genes[LabelSwap.genes$PermutationSet=="Perm 10"&LabelSwap.genes$Candidate==TRUE,"geneID"]
#SBGE.table$Perm7.Can=SBGE.table$geneID%in%LabelSwap.genes[LabelSwap.genes$PermutationSet=="Perm 7","geneID"]

###calculate the # sig genes per Perm Swap (without ANY FILTERING)
 sum(SBGE.table$Perm4.Can,na.rm=TRUE)#Raw=1113->MAF5&Seg4 filter=312
 sum(SBGE.table$Perm5.Can)#Raw=1038->MAF5&Seg4 filter=242
 sum(SBGE.table$Perm6.Can)#Raw=1698->MAF5&Seg4 filter=23
 sum(SBGE.table$Perm7.Can)#Raw=1667->MAF5&Seg4 filter=21
 sum(SBGE.table$Perm8.Can)#Raw=1691->MAF5&Seg4 filter=18
 sum(SBGE.table$Perm9.Can)#Raw=1653->MAF5&Seg4 filter=16
 sum(SBGE.table$Perm10.Can)#Raw=1953->MAF5&Seg4 filter=209
 
 P45=SBGE.table$Perm5.Can&SBGE.table$Perm4.Can#Raw=761 genes->MAF5&Seg4 filter 59
 P78=SBGE.table$Perm7.Can&SBGE.table$Perm8.Can#Raw=552 genes->MAF5&Seg4 filter 2
 P910=SBGE.table$Perm9.Can&SBGE.table$Perm10.Can#Raw=552 ->MAF5&Seg4 filter 4
 sum(P45&SBGE.table$Perm6.Can)#Raw=573 ->MAF5&Seg4 filter 2
 sum(SBGE.table$Perm6.Can&P78)#Raw=542 ->MAF5&Seg4 filter 0
 sum(P45&P78)#Raw=522->MAF5&Seg4 filter 1
 sum(P45&P78&P910&SBGE.table$Perm6.Can)#Raw=520->MAF5&Seg4 filter 0
 
 
 
##--REMOVE GENES THAT DO NOT APPEAR IN OUR SEQUENCING by assign NA
SBGE.table[!SBGE.table$geneID%in%LabelSwap.genes[LabelSwap.genes$PermutationSet=="Perm 4","geneID"],"Perm4.Can"]=NA
SBGE.table[!SBGE.table$geneID%in%LabelSwap.genes[LabelSwap.genes$PermutationSet=="Perm 5","geneID"],"Perm5.Can"]=NA
SBGE.table[!SBGE.table$geneID%in%LabelSwap.genes[LabelSwap.genes$PermutationSet=="Perm 6","geneID"],"Perm6.Can"]=NA
SBGE.table[!SBGE.table$geneID%in%LabelSwap.genes[LabelSwap.genes$PermutationSet=="Perm 7","geneID"],"Perm7.Can"]=NA
SBGE.table[!SBGE.table$geneID%in%LabelSwap.genes[LabelSwap.genes$PermutationSet=="Perm 8","geneID"],"Perm8.Can"]=NA
SBGE.table[!SBGE.table$geneID%in%LabelSwap.genes[LabelSwap.genes$PermutationSet=="Perm 9","geneID"],"Perm9.Can"]=NA
SBGE.table[!SBGE.table$geneID%in%LabelSwap.genes[LabelSwap.genes$PermutationSet=="Perm 10","geneID"],"Perm10.Can"]=NA
#SBGE.table$Perm7.Can=SBGE.table$geneID%in%LabelSwap.genes[LabelSwap.genes$PermutationSet=="Perm 7","geneID"]


#####----Graph SBGE Distribution----##

#quick way
graph.diff.density<-function(data,IDX,col,legends,leg_pos){
  
  x <-data[IDX,] #SBGE.chr2$totalExonLength
  x<-x[!is.na(x[[col]]),col]
  y <-data[!IDX,]
  y<-y[!is.na(y[[col]]),col]
  #generate main LINE density
  u<-range(c(x, y))
  dSig.chr2<-density(x,from=u[1], to=u[2]) #works
  dNSig.chr2<-density(y,from=u[1], to=u[2]) #works
  ##dd_xy<-dx$y -dy$y
  diff.SigNon <-dSig.chr2$y - dNSig.chr2$y
  
  #Generate Simple Kernel Density Estimate uing the default R function
  dSig.fit2 <- replicate(1000,{
    #Sample with replacement (for bootstrap from original dataset). Save the resample to x
    samp <- sample(x, replace=TRUE)                    
    #Generate the density from the resampled dataset, and extract y coordinates to generate variablity bands
    density(samp, from=min(dSig.chr2$x), to=max(dSig.chr2$x))$y}) 
  #note!!!.................................mac(X value data))
  #.........................................................^ end DENSITY
  #.........................................................fit stores y ONLY
  #density(x, from=u[1], to=u[2])}) 
  #Apply the quantile function to the y coordinates to get the bounds of the polygon to be drawn on the y axis?
  dSig.fit3 <- apply(dSig.fit2, 1, quantile, c(0.025,0.975) )
  
  dNSig.fit2 <- replicate(1000,{
    samp <- sample(y, replace=TRUE)  
    #NOTE, this is NOT the harmonized range!
    density(samp, from=min(dNSig.chr2$x), to=max(dNSig.chr2$x))$y}) 
  dNSig.fit3 <- apply(dNSig.fit2, 1, quantile, c(0.025,0.975) )
  
  #generate CI for differences
  test=dSig.fit2 - dNSig.fit2
  diff.CI=apply(dSig.fit2-dNSig.fit2, 1, quantile, c(0.025,0.975) )
  diff.CI.1=apply(test, 1, quantile, c(0.025,0.975) )
  lower_graph=min(diff.SigNon)*1.2
  max=range(c(dSig.fit2,dNSig.fit2))
  upper_graph=max[2]*0.4
  
  plot(dSig.chr2, col=2,ylab="Density",xlim=c(u[1],u[2]),xlab=col,
       main=paste('Difference in',col, "Distribution" ),ylim=c(lower_graph,upper_graph))#,xlim=c(-10,15))
  #plot(dNSig.chr2 , col=1, ylim=c(0, .20), main='Difference in SBGE Distribution', xlab='',ylab="Density",xlim=c(-10,15))
  polygon( c(dNSig.chr2$x, rev(dNSig.chr2$x)), c(dNSig.fit3[1,], rev(dNSig.fit3[2,])),
           col='lightgrey', density = -0.5, border=F)
  #polygon( c(dSig.chr2$x, rev(dSig.chr2$x)), c(dSig.fit3[1,], rev(dSig.fit3[2,])),
  #         col='lightblue', density = -0.5, border=F)
  #same x as normal plots but different y
  polygon( c(dSig.chr2$x, rev(dSig.chr2$x)), c(diff.CI[1,], rev(diff.CI[2,])),
           col="pink", density = -0.5, border=F)
  abline(h=0, lty=3, col=8) # lty=2, lwd=2,
  lines(dNSig.chr2, col=1,lwd=2)
  lines(dSig.chr2, col="blue",lwd=2)
  lines(dSig.chr2$x, diff.SigNon , col="pink", lty=1, lwd=2)  ## <---------------- difference#
  #mtext(sprintf('N(x) = %s  Bandwidth(x) = %s', dSig.chr2$n, signif(dSig.chr2$bw, 3)), 1, 2)
  #mtext(sprintf('N(y) = %s  Bandwidth(y) = %s', dNSig.chr2$n, signif(dNSig.chr2$bw, 3)), 1, 3)
  legend(leg_pos, legend=legends, col=c(2,1,3), 
         lty=c(1, 1, 2), lwd=c(1, 1, 2))
}#Difference, False, True

#Perm 7 graph
SBGE.chr2=subset(SBGE.table,chrm=="2L"|chrm=="2R")

graph.diff.density(SBGE.table[!is.na(SBGE.table$Perm4.Can),],SBGE.table[!is.na(SBGE.table$Perm4.Can),"Perm4.Can"],"Whole.SBGE.Osada",c("Diff","background Perm4"," Cand Perm4"),"topright")
graph.diff.density(SBGE.table[!is.na(SBGE.table$Perm5.Can),],SBGE.table[!is.na(SBGE.table$Perm5.Can),"Perm5.Can"],"Whole.SBGE.Osada",c("Diff","background Perm5","Cand Perm5"),"topright")

graph.diff.density(SBGE.table[!is.na(SBGE.table$Perm6.Can),],SBGE.table[!is.na(SBGE.table$Perm6.Can),"Perm6.Can"],"Whole.SBGE.Osada",c("Diff","background Perm6","Cand Perm6"),"topright")
graph.diff.density(SBGE.table[!is.na(SBGE.table$Perm7.Can),],SBGE.table[!is.na(SBGE.table$Perm7.Can),"Perm7.Can"],"Whole.SBGE.Osada",c("Diff","background Perm7"," Cand Perm 7"),"topright")
graph.diff.density(SBGE.table[!is.na(SBGE.table$Perm8.Can),],SBGE.table[!is.na(SBGE.table$Perm8.Can),"Perm8.Can"],"Whole.SBGE.Osada",c("Diff","background Perm8","Cand Perm 8"),"topright")
graph.diff.density(SBGE.table[!is.na(SBGE.table$Perm9.Can),],SBGE.table[!is.na(SBGE.table$Perm9.Can),"Perm9.Can"],"Whole.SBGE.Osada",c("Diff","background Perm9","Cand Perm 9"),"topright")

graph.diff.density(SBGE.table[!is.na(SBGE.table$Perm10.Can),],SBGE.table[!is.na(SBGE.table$Perm10.Can),"Perm10.Can"],"Whole.SBGE.Osada",c("Diff","background Perm10","Cand Perm 10"),"topright")


#preliminary look at rmf



####proper graphing####
col="Whole.SBGE.Osada"
null=all.shifted.genes
#
Real.Cand <-SBGE.chr2[SBGE.chr2$Candidates,] 
Real.Cand<-Real.Cand[!is.na(Real.Cand[[col]]),col]
#This gets the united range for ALL genes on chromosome 2 so range will be equal and comparable
u=range(SBGE.chr2[!is.na(SBGE.chr2[[col]]),col])
#the real SBGE density is obtained
Real.Cand.dens<-density(Real.Cand,from=u[1], to=u[2])
#the real distribution is bootstrapped to obtain 95% confidence intervals
Real.Cand.dens.fit2 <- replicate(1000,{
  #Sample with replacement (for bootstrap from original dataset). Save the resample to x
  samp <- sample(Real.Cand, replace=TRUE)                    
  #Generate the density from the resampled dataset, and extract y coordinates to generate variablity bands
  density(samp, from=min(Real.Cand.dens$x), to=max(Real.Cand.dens$x))$y}) 
Real.Cand.dens.fit3 <- apply(Real.Cand.dens.fit2, 1, quantile, c(0.025,0.975) )

#this creates the fake SBGE density distribution
Fake.density=null_density_dist(SBGE.chr2,null,col,u)
#the confidence intervals for the null distribution are obtained--> for each x value 
#(SBGE bin as determined by density function) the 95% CI are taken from the y (density) values
Fake.density.fit3<-apply(Fake.density, 1, quantile, c(0.025,0.975) )
Fake.density.fit2<-apply(Fake.density, 1, quantile, c(0.5 ))
#adjust Fake.Candidates.fit2 so that it integrates to 1 and is equal to the real distribution
Cand.x <- Real.Cand.dens$x  ## 512 evenly spaced points on [min(x) - 3 * d$bw, max(x) + 3 * d$bw]
bin.x.width <- Cand.x [2L] - Cand.x[1L]  ## spacing / bin size
all.y <- Real.Cand.dens$y  ## 512 density values for SBGE
curve.area.R <- sum(all.y) * bin.x.width 
curve.area.M<-sum(Fake.density.fit2)*bin.x.width
#the corrected median of the null distribution
Adjust.Fake.density.fit2<-Fake.density.fit2*(curve.area.R/curve.area.M)

#Taking the difference between the null distribution and the candiates
single.diff=Real.Cand.dens$y-Fake.density.fit2
single.diff.CI.1=Real.Cand.dens$y-Fake.density.fit3[1,]
single.diff.CI.2=Real.Cand.dens$y-Fake.density.fit3[2,]

#data ranges for graphing purposes
max=range(c(Real.Cand.dens.fit2,Fake.density))
upper_graph=max[2]
lower_graph=min(single.diff.CI.2=Real.Cand.dens$y-Fake.density.fit3[2,])
# 
# png(
#   "Genomic_Fig1v2.png",
#   width     = 2.75,
#   height    = 2.75,
#   units     = "in",
#   res       = 1200,
#   pointsize = 4)
plot(Real.Cand.dens, col="white",ylab="Density",xlab="Sex-Biased Gene Expression (SBGE)",xlim=c(-10,15),#xlim=c(u[1],u[2]),
     main=paste('Difference in Whole Body SBGE Distribution' ),ylim=c(-0.025,upper_graph*0.9),
     cex.axis=1.2, cex.lab=1.4,cex.main=1.6)
polygon( c(Real.Cand.dens$x, rev(Real.Cand.dens$x)), c(Fake.density.fit3[1,], rev(Fake.density.fit3[2,])),
         col='grey', density = -0.5, border=F)
#polygon( c(Real.Cand.dens$x, rev(Real.Cand.dens$x)),  c(single.diff.CI.1, rev(single.diff.CI.2)),
#         col='pink', density = -0.5, border=F)
abline(h=0, lty=3, col=8) 
lines(Real.Cand.dens, col="blue",lwd=2)
lines(Real.Cand.dens$x,Adjust.Fake.density.fit2, col="black",lwd=2)
#lines(Real.Cand.dens$x,single.diff, col="darkred",lwd=2,lty=2)
#legend("topright", legend=c("Candidates","Null Distribution","Candidates-Null"), col=c("blue","grey","pink"), 
legend("topright", legend=c("Candidates","Null Distribution"), col=c("blue","grey"), 
       lty=c(1, 1, 2), lwd=c(1, 1, 2),cex=1.4)
# dev.off()



######-----Permuting Generation 7 Data-----#####
#####----------create G7 Permutation Data------
##Assumes from scratch
#note R.avg and N.avg are averages per sample per site!
##other if minimum coverage is set to 2 then have ~500k sites becaus
#stars 1.7 mill
G7.all<-read.table("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7a_ALL_cov2_rc",header=TRUE,sep="",
                   col.names=c("chrm","pos","rc","allele_count","allele_states	",
                   "deletion_sum","snp_type","major_alleles(maa)","minor_alleles(mia)",
                   "maa_1	","maa_2","maa_3","maa_4","maa_5","maa_6","maa_7",
                   "maa_8	","maa_9","maa_10","maa_11","maa_12","mia_1",
                   "mia_2","mia_3","mia_4","mia_5","mia_6","mia_7","mia_8",
                   "mia_9","mia_10","mia_11","mia_12"))
#maa_1 -> major allle, complementary to mia_1 -> minor allele
#goes N1, R1, N2, R2
#keep only the minor allele frequencies and remove "major allele" since this is inverse
G7.all = G7.all[,c(1:10,22:33)]


###------FILTERING STEp------
#due to server issues, filter this here

#renaming col in case deleted
colnames(G7.all)[1]="chrm"
G7.all=subset(G7.all,chrm=="2L"|chrm=="2R")#600k
colnames(G7.all)[2]="pos"
colnames(G7.all)[4]="allele_count"
G7.all=subset(G7.all,allele_count==2)
#only 560k sites remain


#I tell if the allele is segregating in at least 3 populations, minaor allele COUNT
#must be greater than 3
#Do this using "Red" frequency only since this is already half NonRed
G7.all$Total.CountA= as.numeric( str_split_fixed(G7.all[,"mia_2"],"/",n=2)[,1])  +   as.numeric( str_split_fixed(G7.all[,"mia_4"],"/",n=2)[,1])+
 as.numeric( str_split_fixed(G7.all[,"mia_6"],"/",n=2)[,1])  +   as.numeric( str_split_fixed(G7.all[,"mia_8"],"/",n=2)[,1])+
as.numeric( str_split_fixed(G7.all[,"mia_10"],"/",n=2)[,1])  +   as.numeric( str_split_fixed(G7.all[,"mia_12"],"/",n=2)[,1])

#this is a proxy for a seg3 filter
G7.all=G7.all[G7.all$Total.CountA>3,]#213k sites remain

G7.all$Total.Coverage=as.numeric( str_split_fixed(G7.all[,"mia_2"],"/",n=2)[,2])  +   as.numeric( str_split_fixed(G7.all[,"mia_4"],"/",n=2)[,2])+
  as.numeric( str_split_fixed(G7.all[,"mia_6"],"/",n=2)[,2])  +   as.numeric( str_split_fixed(G7.all[,"mia_8"],"/",n=2)[,2])+
  as.numeric( str_split_fixed(G7.all[,"mia_10"],"/",n=2)[,2])  +   as.numeric( str_split_fixed(G7.all[,"mia_12"],"/",n=2)[,2])

#MAF filter x>0.05
G7.all=G7.all[G7.all$Total.CountA/G7.all$Total.Coverage>0.05,]#197k sites remain

#Add ID
G7.all$ID=paste(G7.all$chrm,G7.all$pos)

#------unstacking process
#split treatments
#25 should be ID column
G7.Non=G7.all[,c(1,2,10,12,14,16,18,20,25)]
G7.Non$Treatment="NonRed"
colnames(G7.Non)[3]="Pop1"
colnames(G7.Non)[4]="Pop2"
colnames(G7.Non)[5]="Pop3"
colnames(G7.Non)[6]="Pop4"
colnames(G7.Non)[7]="Pop5"
colnames(G7.Non)[8]="Pop6"
#25 should be ID column
G7.Red=G7.all[,c(1,2,11,13,15,17,19,21,25)]
G7.Red$Treatment="Red"
colnames(G7.Red)[3]="Pop1"
colnames(G7.Red)[4]="Pop2"
colnames(G7.Red)[5]="Pop3"
colnames(G7.Red)[6]="Pop4"
colnames(G7.Red)[7]="Pop5"
colnames(G7.Red)[8]="Pop6"

G7.data=G7.all[1,c("chrm","pos","ID")]
G7.data$Treatment="Wack"
G7.data$Population="Wack"
G7.data$CountA=0
G7.data$CountB=0
##unstack populations----this will take a minute
for(pop in c("Pop1","Pop2","Pop3","Pop4","Pop5","Pop6")){
  #pulling out pop information
  temp.Red=G7.Red[,c("chrm","pos","ID","Treatment")]
  temp.Red$Population=pop
  temp.Red$CountA=as.numeric(str_split_fixed(G7.Red[,pop],"/",n=2)[,1])
  temp.Red$CountB=as.numeric(str_split_fixed(G7.Red[,pop],"/",n=2)[,2])
  #now NonRed
  temp.Non=G7.Non[,c("chrm","pos","ID","Treatment")]
  temp.Non$Population=pop
  temp.Non$CountA=as.numeric(str_split_fixed(G7.Non[,pop],"/",n=2)[,1])
  temp.Non$CountB=as.numeric(str_split_fixed(G7.Non[,pop],"/",n=2)[,2])
  
  G7.data=rbind(G7.data,rbind(temp.Red,temp.Non))}
#remove fake row. G7.data now has 2.4 mill rows (but still only 197k sites)
G7.data=G7.data[-1,]

pval_per_SNP_G7<-function(data.MAF,SNP.IDX){
  #make this a bit smaller so in linear model not indexing as large a dataset each time
  smaller.data=data.MAF[data.MAF$ID%in%SNP.IDX,]  
  #save every site
  G7.pval=smaller.data[,c("chrm","pos","ID")]
  #this gets 1 hit per site
  G7.pval=unique(G7.pval)
  G7.pval$pval=1
  #data for a single SNP
  for(SNP in SNP.IDX){
    #print(SNP)
    #note, even though I am only 
    temp.model=glm(cbind(CountA,CountB)~Treatment+Population,data=smaller.data[smaller.data$ID==SNP,],family="quasibinomial")
    #this will return the pvalue for every SNP from all treatment and population (whether it's different from 1)
    temp.ps=coef(summary(temp.model))[,4]
    #this returns second value on the list which is whether there is a different bt Red/Non
    temp.p=temp.ps[2]
    G7.pval[G7.pval$ID==SNP,"pval"]=temp.p}
  
  
  return(G7.pval)}

#the SNP ID should be from G7.all (197k IDs) because it's the IDs to include. 
#test=pval_per_SNP_G7(G7.data,G7.all[1:5000,"ID"])
#5000=21 secs. 10k=66 seconds
#test2=pval_per_SNP_G7(G7.data,G7.all[1:25000,"ID"])
#25k=5 minutes 44 seconds

# G7.P0.p1=pval_per_SNP_G7(G7.data,G7.all[1:25000,"ID"])
# G7.P0.p2=pval_per_SNP_G7(G7.data,G7.all[25001:50000,"ID"])
# G7.P0.p3=pval_per_SNP_G7(G7.data,G7.all[50001:75000,"ID"])
# G7.P0.p4=pval_per_SNP_G7(G7.data,G7.all[75001:100000,"ID"])
# G7.P0.p5=pval_per_SNP_G7(G7.data,G7.all[100001:125000,"ID"])
# G7.P0.p6=pval_per_SNP_G7(G7.data,G7.all[125001:150000,"ID"])
# G7.P0.p7=pval_per_SNP_G7(G7.data,G7.all[150001:175000,"ID"])
# G7.P0.p8=pval_per_SNP_G7(G7.data,G7.all[175001:length(G7.all$ID),"ID"])
# G7.P0.sigInfo=bind_rows(G7.P0.p1,G7.P0.p2,G7.P0.p3,G7.P0.p4,G7.P0.p5,G7.P0.p6,G7.P0.p7,G7.P0.p8)
# G7.P0.sigInfo$p.adjust=p.adjust(G7.P0.sigInfo$pval,method="BH")
# #write.csv(G7.P0.sigInfo,"C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Real_pval.csv")
# sum(G7.P0.sigInfo$p.adjust<0.05)#6734 SNPs
# 

#---------Label Swap Permutation---------
#Permutate all the populations first and then run the SNP significance (just for ease of record keeping)
#eventually will move
# #Perm 1 = population 123 are swapped
# G7.Perm1=G7.data
# G7.Perm1[G7.Perm1$Population%in%c("Pop1","Pop2","Pop3")&G7.Perm1$Treatment=="NonRed","Treatment"]="Non1"
# G7.Perm1[G7.Perm1$Population%in%c("Pop1","Pop2","Pop3")&G7.Perm1$Treatment=="Red","Treatment"]="NonRed"
# G7.Perm1[G7.Perm1$Population%in%c("Pop1","Pop2","Pop3")&G7.Perm1$Treatment=="Non1","Treatment"]="Red"
# #Perm 2 = population 124 are swapped
# G7.Perm2=G7.data
# G7.Perm2[G7.Perm2$Population%in%c("Pop1","Pop2","Pop4")&G7.Perm2$Treatment=="NonRed","Treatment"]="Non1"
# G7.Perm2[G7.Perm2$Population%in%c("Pop1","Pop2","Pop4")&G7.Perm2$Treatment=="Red","Treatment"]="NonRed"
# G7.Perm2[G7.Perm2$Population%in%c("Pop1","Pop2","Pop4")&G7.Perm2$Treatment=="Non1","Treatment"]="Red"
# #Perm 3 = population 125 are swapped
# G7.Perm3=G7.data
# G7.Perm3[G7.Perm3$Population%in%c("Pop1","Pop2","Pop5")&G7.Perm3$Treatment=="NonRed","Treatment"]="Non1"
# G7.Perm3[G7.Perm3$Population%in%c("Pop1","Pop2","Pop5")&G7.Perm3$Treatment=="Red","Treatment"]="NonRed"
# G7.Perm3[G7.Perm3$Population%in%c("Pop1","Pop2","Pop5")&G7.Perm3$Treatment=="Non1","Treatment"]="Red"
# #Perm 4 = population 126 are swapped
# G7.Perm4=G7.data
# G7.Perm4[G7.Perm4$Population%in%c("Pop1","Pop2","Pop6")&G7.Perm4$Treatment=="NonRed","Treatment"]="Non1"
# G7.Perm4[G7.Perm4$Population%in%c("Pop1","Pop2","Pop6")&G7.Perm4$Treatment=="Red","Treatment"]="NonRed"
# G7.Perm4[G7.Perm4$Population%in%c("Pop1","Pop2","Pop6")&G7.Perm4$Treatment=="Non1","Treatment"]="Red"
# #Perm 5 = population 134 are swapped
# G7.Perm5=G7.data
# G7.Perm5[G7.Perm5$Population%in%c("Pop1","Pop3","Pop4")&G7.Perm5$Treatment=="NonRed","Treatment"]="Non1"
# G7.Perm5[G7.Perm5$Population%in%c("Pop1","Pop3","Pop4")&G7.Perm5$Treatment=="Red","Treatment"]="NonRed"
# G7.Perm5[G7.Perm5$Population%in%c("Pop1","Pop3","Pop4")&G7.Perm5$Treatment=="Non1","Treatment"]="Red"
# 
# #Perm 6 = population 135 are swapped
# G7.Perm6=G7.data
# G7.Perm6[G7.Perm6$Population%in%c("Pop1","Pop3","Pop5")&G7.Perm6$Treatment=="NonRed","Treatment"]="Non1"
# G7.Perm6[G7.Perm6$Population%in%c("Pop1","Pop3","Pop5")&G7.Perm6$Treatment=="Red","Treatment"]="NonRed"
# G7.Perm6[G7.Perm6$Population%in%c("Pop1","Pop3","Pop5")&G7.Perm6$Treatment=="Non1","Treatment"]="Red"
# # #Perm 7 = population 136 are swapped
# G7.Perm7=G7.data
# G7.Perm7[G7.Perm7$Population%in%c("Pop1","Pop3","Pop6")&G7.Perm7$Treatment=="NonRed","Treatment"]="Non1"
# G7.Perm7[G7.Perm7$Population%in%c("Pop1","Pop3","Pop6")&G7.Perm7$Treatment=="Red","Treatment"]="NonRed"
# G7.Perm7[G7.Perm7$Population%in%c("Pop1","Pop3","Pop6")&G7.Perm7$Treatment=="Non1","Treatment"]="Red"
#Perm 8 = population 145 are swapped
# G7.Perm8=G7.data
# G7.Perm8[G7.Perm8$Population%in%c("Pop1","Pop4","Pop5")&G7.Perm8$Treatment=="NonRed","Treatment"]="Non1"
# G7.Perm8[G7.Perm8$Population%in%c("Pop1","Pop4","Pop5")&G7.Perm8$Treatment=="Red","Treatment"]="NonRed"
# G7.Perm8[G7.Perm8$Population%in%c("Pop1","Pop4","Pop5")&G7.Perm8$Treatment=="Non1","Treatment"]="Red"
#Perm 9 = population 146 are swapped
# G7.Perm9=G7.data
# G7.Perm9[G7.Perm9$Population%in%c("Pop1","Pop4","Pop6")&G7.Perm9$Treatment=="NonRed","Treatment"]="Non1"
# G7.Perm9[G7.Perm9$Population%in%c("Pop1","Pop4","Pop6")&G7.Perm9$Treatment=="Red","Treatment"]="NonRed"
# G7.Perm9[G7.Perm9$Population%in%c("Pop1","Pop4","Pop6")&G7.Perm9$Treatment=="Non1","Treatment"]="Red"
# #Perm 10 = population 156 are swapped
# G7.Perm10=G7.data
# G7.Perm10[G7.Perm10$Population%in%c("Pop1","Pop5","Pop6")&G7.Perm10$Treatment=="NonRed","Treatment"]="Non1"
# G7.Perm10[G7.Perm10$Population%in%c("Pop1","Pop5","Pop6")&G7.Perm10$Treatment=="Red","Treatment"]="NonRed"
# G7.Perm10[G7.Perm10$Population%in%c("Pop1","Pop5","Pop6")&G7.Perm10$Treatment=="Non1","Treatment"]="Red"


# G7.P1.p1=pval_per_SNP_G7(G7.Perm1,G7.all[1:25000,"ID"])
# G7.P1.p2=pval_per_SNP_G7(G7.Perm1,G7.all[25001:50000,"ID"])
# G7.P1.p3=pval_per_SNP_G7(G7.Perm1,G7.all[50001:75000,"ID"])
# G7.P1.p4=pval_per_SNP_G7(G7.Perm1,G7.all[75001:100000,"ID"])
# G7.P1.p5=pval_per_SNP_G7(G7.Perm1,G7.all[100001:125000,"ID"])
# G7.P1.p6=pval_per_SNP_G7(G7.Perm1,G7.all[125001:150000,"ID"])
# G7.P1.p7=pval_per_SNP_G7(G7.Perm1,G7.all[150001:175000,"ID"])
# G7.P1.p8=pval_per_SNP_G7(G7.Perm1,G7.all[175001:length(G7.all$ID),"ID"])
# G7.P1.sigInfo=bind_rows(G7.P1.p1,G7.P1.p2,G7.P1.p3,G7.P1.p4,G7.P1.p5,G7.P1.p6,G7.P1.p7,G7.P1.p8)
# G7.P1.sigInfo$p.adjust=p.adjust(G7.P1.sigInfo$pval,method="BH")
# write.csv(G7.P1.sigInfo,"C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm1_pval.csv")
# sum(G7.P1.sigInfo$p.adjust<0.05)

# 
# G7.P2.p1=pval_per_SNP_G7(G7.Perm2,G7.all[1:25000,"ID"])
# G7.P2.p2=pval_per_SNP_G7(G7.Perm2,G7.all[25001:50000,"ID"])
# G7.P2.p3=pval_per_SNP_G7(G7.Perm2,G7.all[50001:75000,"ID"])
# G7.P2.p4=pval_per_SNP_G7(G7.Perm2,G7.all[75001:100000,"ID"])
# G7.P2.p5=pval_per_SNP_G7(G7.Perm2,G7.all[100001:125000,"ID"])
# G7.P2.p6=pval_per_SNP_G7(G7.Perm2,G7.all[125001:150000,"ID"])
# G7.P2.p7=pval_per_SNP_G7(G7.Perm2,G7.all[150001:175000,"ID"])
# G7.P2.p8=pval_per_SNP_G7(G7.Perm2,G7.all[175001:length(G7.all$ID),"ID"])
# G7.P2.sigInfo=bind_rows(G7.P2.p1,G7.P2.p2,G7.P2.p3,G7.P2.p4,G7.P2.p5,G7.P2.p6,G7.P2.p7,G7.P2.p8)
# G7.P2.sigInfo$p.adjust=p.adjust(G7.P2.sigInfo$pval,method="BH")
# write.csv(G7.P2.sigInfo,"C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm2_pval.csv")
# sum(G7.P2.sigInfo$p.adjust<0.05)
# 
# G7.P3.p1=pval_per_SNP_G7(G7.Perm3,G7.all[1:25000,"ID"])
# G7.P3.p2=pval_per_SNP_G7(G7.Perm3,G7.all[25001:50000,"ID"])
# G7.P3.p3=pval_per_SNP_G7(G7.Perm3,G7.all[50001:75000,"ID"])
# G7.P3.p4=pval_per_SNP_G7(G7.Perm3,G7.all[75001:100000,"ID"])
# G7.P3.p5=pval_per_SNP_G7(G7.Perm3,G7.all[100001:125000,"ID"])
# G7.P3.p6=pval_per_SNP_G7(G7.Perm3,G7.all[125001:150000,"ID"])
# G7.P3.p7=pval_per_SNP_G7(G7.Perm3,G7.all[150001:175000,"ID"])
# G7.P3.p8=pval_per_SNP_G7(G7.Perm3,G7.all[175001:length(G7.all$ID),"ID"])
# G7.P3.sigInfo=bind_rows(G7.P3.p1,G7.P3.p2,G7.P3.p3,G7.P3.p4,G7.P3.p5,G7.P3.p6,G7.P3.p7,G7.P3.p8)
# G7.P3.sigInfo$p.adjust=p.adjust(G7.P3.sigInfo$pval,method="BH")
# write.csv(G7.P3.sigInfo,"C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm3_pval.csv")
# sum(G7.P3.sigInfo$p.adjust<0.05)

# G7.P4.p1=pval_per_SNP_G7(G7.Perm4,G7.all[1:25000,"ID"])
# G7.P4.p2=pval_per_SNP_G7(G7.Perm4,G7.all[25001:50000,"ID"])
# G7.P4.p3=pval_per_SNP_G7(G7.Perm4,G7.all[50001:75000,"ID"])
# G7.P4.p4=pval_per_SNP_G7(G7.Perm4,G7.all[75001:100000,"ID"])
# G7.P4.p5=pval_per_SNP_G7(G7.Perm4,G7.all[100001:125000,"ID"])
# G7.P4.p6=pval_per_SNP_G7(G7.Perm4,G7.all[125001:150000,"ID"])
# G7.P4.p7=pval_per_SNP_G7(G7.Perm4,G7.all[150001:175000,"ID"])
# G7.P4.p8=pval_per_SNP_G7(G7.Perm4,G7.all[175001:length(G7.all$ID),"ID"])
# G7.P4.sigInfo=bind_rows(G7.P4.p1,G7.P4.p2,G7.P4.p3,G7.P4.p4,G7.P4.p5,G7.P4.p6,G7.P4.p7,G7.P4.p8)
# G7.P4.sigInfo$p.adjust=p.adjust(G7.P4.sigInfo$pval,method="BH")
# write.csv(G7.P4.sigInfo,"C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm4_pval.csv")
# sum(G7.P4.sigInfo$p.adjust<0.05)

# G7.P5.p1=pval_per_SNP_G7(G7.Perm5,G7.all[1:25000,"ID"])
# G7.P5.p2=pval_per_SNP_G7(G7.Perm5,G7.all[25001:50000,"ID"])
# G7.P5.p3=pval_per_SNP_G7(G7.Perm5,G7.all[50001:75000,"ID"])
# G7.P5.p4=pval_per_SNP_G7(G7.Perm5,G7.all[75001:100000,"ID"])
# G7.P5.p5=pval_per_SNP_G7(G7.Perm5,G7.all[100001:125000,"ID"])
# G7.P5.p6=pval_per_SNP_G7(G7.Perm5,G7.all[125001:150000,"ID"])
# G7.P5.p7=pval_per_SNP_G7(G7.Perm5,G7.all[150001:175000,"ID"])
# G7.P5.p8=pval_per_SNP_G7(G7.Perm5,G7.all[175001:length(G7.all$ID),"ID"])
# G7.P5.sigInfo=bind_rows(G7.P5.p1,G7.P5.p2,G7.P5.p3,G7.P5.p4,G7.P5.p5,G7.P5.p6,G7.P5.p7,G7.P5.p8)
# G7.P5.sigInfo$p.adjust=p.adjust(G7.P5.sigInfo$pval,method="BH")
# write.csv(G7.P5.sigInfo,"C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm5_pval.csv")
# sum(G7.P5.sigInfo$p.adjust<0.05)
# 
# G7.P6.p1=pval_per_SNP_G7(G7.Perm6,G7.all[1:25000,"ID"])
# G7.P6.p2=pval_per_SNP_G7(G7.Perm6,G7.all[25001:50000,"ID"])
# G7.P6.p3=pval_per_SNP_G7(G7.Perm6,G7.all[50001:75000,"ID"])
# G7.P6.p4=pval_per_SNP_G7(G7.Perm6,G7.all[75001:100000,"ID"])
# G7.P6.p5=pval_per_SNP_G7(G7.Perm6,G7.all[100001:125000,"ID"])
# G7.P6.p6=pval_per_SNP_G7(G7.Perm6,G7.all[125001:150000,"ID"])
# G7.P6.p7=pval_per_SNP_G7(G7.Perm6,G7.all[150001:175000,"ID"])
# G7.P6.p8=pval_per_SNP_G7(G7.Perm6,G7.all[175001:length(G7.all$ID),"ID"])
# G7.P6.sigInfo=bind_rows(G7.P6.p1,G7.P6.p2,G7.P6.p3,G7.P6.p4,G7.P6.p5,G7.P6.p6,G7.P6.p7,G7.P6.p8)
# G7.P6.sigInfo$p.adjust=p.adjust(G7.P6.sigInfo$pval,method="BH")
# write.csv(G7.P6.sigInfo,"C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm6_pval.csv")
# sum(G7.P6.sigInfo$p.adjust<0.05)

# G7.P7.p1=pval_per_SNP_G7(G7.Perm7,G7.all[1:25000,"ID"])
# G7.P7.p2=pval_per_SNP_G7(G7.Perm7,G7.all[25001:50000,"ID"])
# G7.P7.p3=pval_per_SNP_G7(G7.Perm7,G7.all[50001:75000,"ID"])
# G7.P7.p4=pval_per_SNP_G7(G7.Perm7,G7.all[75001:100000,"ID"])
# G7.P7.p5=pval_per_SNP_G7(G7.Perm7,G7.all[100001:125000,"ID"])
# G7.P7.p6=pval_per_SNP_G7(G7.Perm7,G7.all[125001:150000,"ID"])
# G7.P7.p7=pval_per_SNP_G7(G7.Perm7,G7.all[150001:175000,"ID"])
# G7.P7.p8=pval_per_SNP_G7(G7.Perm7,G7.all[175001:length(G7.all$ID),"ID"])
# G7.P7.sigInfo=bind_rows(G7.P7.p1,G7.P7.p2,G7.P7.p3,G7.P7.p4,G7.P7.p5,G7.P7.p6,G7.P7.p7,G7.P7.p8)
# G7.P7.sigInfo$p.adjust=p.adjust(G7.P7.sigInfo$pval,method="BH")
# write.csv(G7.P7.sigInfo,"C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm7_pval.csv")
# # sum(G7.P7.sigInfo$p.adjust<0.05)
# 
# G7.P8.p1=pval_per_SNP_G7(G7.Perm8,G7.all[1:25000,"ID"])
# G7.P8.p2=pval_per_SNP_G7(G7.Perm8,G7.all[25001:50000,"ID"])
# G7.P8.p3=pval_per_SNP_G7(G7.Perm8,G7.all[50001:75000,"ID"])
# G7.P8.p4=pval_per_SNP_G7(G7.Perm8,G7.all[75001:100000,"ID"])
# G7.P8.p5=pval_per_SNP_G7(G7.Perm8,G7.all[100001:125000,"ID"])
# G7.P8.p6=pval_per_SNP_G7(G7.Perm8,G7.all[125001:150000,"ID"])
# G7.P8.p7=pval_per_SNP_G7(G7.Perm8,G7.all[150001:175000,"ID"])
# G7.P8.p8=pval_per_SNP_G7(G7.Perm8,G7.all[175001:length(G7.all$ID),"ID"])
# G7.P8.sigInfo=bind_rows(G7.P8.p1,G7.P8.p2,G7.P8.p3,G7.P8.p4,G7.P8.p5,G7.P8.p6,G7.P8.p7,G7.P8.p8)
# G7.P8.sigInfo$p.adjust=p.adjust(G7.P8.sigInfo$pval,method="BH")
# write.csv(G7.P8.sigInfo,"C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm8_pval.csv")
# sum(G7.P8.sigInfo$p.adjust<0.05)

# G7.P9.p1=pval_per_SNP_G7(G7.Perm9,G7.all[1:25000,"ID"])
# G7.P9.p2=pval_per_SNP_G7(G7.Perm9,G7.all[25001:50000,"ID"])
# G7.P9.p3=pval_per_SNP_G7(G7.Perm9,G7.all[50001:75000,"ID"])
# G7.P9.p4=pval_per_SNP_G7(G7.Perm9,G7.all[75001:100000,"ID"])
# G7.P9.p5=pval_per_SNP_G7(G7.Perm9,G7.all[100001:125000,"ID"])
# G7.P9.p6=pval_per_SNP_G7(G7.Perm9,G7.all[125001:150000,"ID"])
# G7.P9.p7=pval_per_SNP_G7(G7.Perm9,G7.all[150001:175000,"ID"])
# G7.P9.p8=pval_per_SNP_G7(G7.Perm9,G7.all[175001:length(G7.all$ID),"ID"])
# G7.P9.sigInfo=bind_rows(G7.P9.p1,G7.P9.p2,G7.P9.p3,G7.P9.p4,G7.P9.p5,G7.P9.p6,G7.P9.p7,G7.P9.p8)
# G7.P9.sigInfo$p.adjust=p.adjust(G7.P9.sigInfo$pval,method="BH")
# write.csv(G7.P9.sigInfo,"C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm9_pval.csv")
# sum(G7.P9.sigInfo$p.adjust<0.05)
# 
# G7.P10.p1=pval_per_SNP_G7(G7.Perm10,G7.all[1:25000,"ID"])
# G7.P10.p2=pval_per_SNP_G7(G7.Perm10,G7.all[25001:50000,"ID"])
# G7.P10.p3=pval_per_SNP_G7(G7.Perm10,G7.all[50001:75000,"ID"])
# G7.P10.p4=pval_per_SNP_G7(G7.Perm10,G7.all[75001:100000,"ID"])
# G7.P10.p5=pval_per_SNP_G7(G7.Perm10,G7.all[100001:125000,"ID"])
# G7.P10.p6=pval_per_SNP_G7(G7.Perm10,G7.all[125001:150000,"ID"])
# G7.P10.p7=pval_per_SNP_G7(G7.Perm10,G7.all[150001:175000,"ID"])
# G7.P10.p8=pval_per_SNP_G7(G7.Perm10,G7.all[175001:length(G7.all$ID),"ID"])
# G7.P10.sigInfo=bind_rows(G7.P10.p1,G7.P10.p2,G7.P10.p3,G7.P10.p4,G7.P10.p5,G7.P10.p6,G7.P10.p7,G7.P10.p8)
# G7.P10.sigInfo$p.adjust=p.adjust(G7.P10.sigInfo$pval,method="BH")
# write.csv(G7.P10.sigInfo,"C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm10_pval.csv")
# sum(G7.P10.sigInfo$p.adjust<0.05)



##---combine into one file and FULLY MAF filter (NR.MAF avg as with G130)
# 
# G7.real=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Real_pval.csv")
# G7.P1=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm1_pval.csv")
# G7.P2=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm2_pval.csv")
# G7.P3=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm3_pval.csv")
# G7.P4=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm4_pval.csv")
# G7.P5=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm5_pval.csv")
# G7.P6=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm6_pval.csv")
# G7.P7=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm7_pval.csv")
# G7.P8=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm8_pval.csv")
# G7.P9=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm9_pval.csv")
# G7.P10=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_Perm10_pval.csv")
# 
# #all of the p.adjust should be in the same order so there doesn't need to be any merging
# G7.real$Perm1=G7.P1$p.adjust
# G7.real$Perm2=G7.P2$p.adjust
# G7.real$Perm3=G7.P3$p.adjust
# G7.real$Perm4=G7.P4$p.adjust
# G7.real$Perm5=G7.P5$p.adjust
# G7.real$Perm6=G7.P6$p.adjust
# G7.real$Perm7=G7.P7$p.adjust
# G7.real$Perm8=G7.P8$p.adjust
# G7.real$Perm9=G7.P9$p.adjust
# G7.real$Perm10=G7.P10$p.adjust

# write.csv(G7.real,"C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_allPerms_pval.csv")
# 



# G7.all$ID=paste(G7.all$chrm,G7.all$pos)
# G7.all$Sig=G7.all$ID%in%G7.real[G7.real$Sig,"ID"]
# 
# 
# test=G7.all[G7.all$Sig,]
#####-----actually analyze and filter Permutation data----####

#the significance data. Because This was generated on my computer, columns were removed to reduce file sizes
#note this was run on cov2 filter
G7.data=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_allPerms_pval.csv",header=TRUE)
#remove misc column
G7.data=G7.data[,-1]
G7.data=G7.data[,-1]
#import other files with more information.

# ##this is the file I used for G7 MAF binning for D.simulans earlier investigation. Note this is run cov5 filter
# #contains the N.avg and R.avg ferquency
# G7.MAF<-read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_totalCount_cov5_w_avg.txt")
# G7.MAF$N.avg=1-G7.MAF$N.avg
# G7.MAF$R.avg=1-G7.MAF$R.avg
# G7.MAF$NR.avg=apply(G7.MAF[,c("N.avg","R.avg")],1,mean)
# G7.MAF$Total.Major.Count= as.numeric( str_split_fixed(G7.MAF[,"Red"],"/",n=2)[,1])  +   as.numeric( str_split_fixed(G7.MAF[,"NonRed"],"/",n=2)[,1])
# G7.MAF$Total.Coverage= as.numeric( str_split_fixed(G7.MAF[,"Red"],"/",n=2)[,2])  +   as.numeric( str_split_fixed(G7.MAF[,"NonRed"],"/",n=2)[,2])
# G7.MAF$Minor.Count=G7.MAF$Total.Coverage-G7.MAF$Total.Major.Count
# #now only select those on chrm2
# G7.MAF=subset(G7.MAF,chrm=="2L"|chrm=="2R")#99 468
# 
# #sum(G7.MAF$Minor.Count<3)
# #Do an additional MAF5 filteirng pass. The intial MAF5 filtering pass was done on combined counts, not averaged frequencies
# #which means it does not remove sites for which one Sample Type (Red or NonRed) was low/basically 0
# G7.double.MAF=G7.data[G7.data$ID%in%G7.MAF[G7.MAF$NR.avg>0.05,"ID"],]

#sum(G7.data$ID%in%G7.MAF[,"ID"])#only 67 508 sites 
#this is likely because
#apply a segregating 3 filter

G7.all<-read.table("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7a_ALL_cov2_rc",header=TRUE,sep="",
                   col.names=c("chrm","pos","rc","allele_count","allele_states	",
                               "deletion_sum","snp_type","major_alleles(maa)","minor_alleles(mia)",
                               "maa_1	","maa_2","maa_3","maa_4","maa_5","maa_6","maa_7",
                               "maa_8	","maa_9","maa_10","maa_11","maa_12","mia_1",
                               "mia_2","mia_3","mia_4","mia_5","mia_6","mia_7","mia_8",
                               "mia_9","mia_10","mia_11","mia_12"))
#maa_1 -> major allle, complementary to mia_1 -> minor allele
#goes N1, R1, N2, R2
#keep only the minor allele frequencies and remove "major allele" since this is inverse
G7.all = G7.all[,c(1:10,22:33)]
G7.all=subset(G7.all,chrm=="2L"|chrm=="2R")#600k
G7.all=subset(G7.all,allele_count==2)
#only 560k sites remain

#Do actual seg3 filtering and actual coverage (per pop) filtering
#minor alle counts start a column 11 (col 10 is major allele)
G7.all$Pop.SegR=0
G7.all$Pop.SegN=0
G7.all$cov5=TRUE
for(col in seq(from=11,to=22,by=2)){
  G7.all[ as.numeric(str_split_fixed(G7.all[,col],"/",n=2)[,1])>0,"Pop.SegN"]=  G7.all[ as.numeric(str_split_fixed(G7.all[,col],"/",n=2)[,1])>0,"Pop.SegN"]+1
  G7.all[ as.numeric(str_split_fixed(G7.all[,col+1],"/",n=2)[,1])>0,"Pop.SegR"]=  G7.all[ as.numeric(str_split_fixed(G7.all[,col+1],"/",n=2)[,1])>0,"Pop.SegR"]+1
  G7.all[as.numeric( str_split_fixed(G7.all[,col+1],"/",n=2)[,2])<5|
           as.numeric( str_split_fixed(G7.all[,col+1],"/",n=2)[,2])<5,"cov5"]=FALSE
}

#check that G7.MAF and G7.all are showing the same data
sum(!G7.all$cov5)
G7.all$ID=paste(G7.all$chrm,G7.all$pos)
# sum(G7.all[G7.all$cov5,"ID"]%in%G7.MAF[,"ID"])#99k pass filter. where are the 

#Seg3 filter proper. Note this is kind of an average so not exact 
#ideally R and N would be added per population (N1 and R1 count as 1) but Pop.Seg.avg isn't used as the filter
G7.all$Pop.Seg.avg=(G7.all$Pop.SegR+G7.all$Pop.SegN)/2
#how many fail seg3 filter?
sum(G7.all$Pop.SegR>3)
#sum(G7.all$Pop.SegN<3)

#Do this using "Red" frequency only since this is already half NonRed
G7.all$Red.MinorCount=as.numeric( str_split_fixed(G7.all[,"mia_2"],"/",n=2)[,1])  +   as.numeric( str_split_fixed(G7.all[,"mia_4"],"/",n=2)[,1])+
  as.numeric( str_split_fixed(G7.all[,"mia_6"],"/",n=2)[,1])  +   as.numeric( str_split_fixed(G7.all[,"mia_8"],"/",n=2)[,1])+
  as.numeric( str_split_fixed(G7.all[,"mia_10"],"/",n=2)[,1])  +   as.numeric( str_split_fixed(G7.all[,"mia_12"],"/",n=2)[,1])

G7.all$Red.Coverage=as.numeric( str_split_fixed(G7.all[,"mia_2"],"/",n=2)[,2])  +   
  +as.numeric( str_split_fixed(G7.all[,"mia_4"],"/",n=2)[,2])+
  as.numeric( str_split_fixed(G7.all[,"mia_6"],"/",n=2)[,2])  +   as.numeric( str_split_fixed(G7.all[,"mia_8"],"/",n=2)[,2])+
  as.numeric( str_split_fixed(G7.all[,"mia_10"],"/",n=2)[,2])  +   as.numeric( str_split_fixed(G7.all[,"mia_12"],"/",n=2)[,2])

G7.all$NonRed.MinorCount=as.numeric( str_split_fixed(G7.all[,"mia_1"],"/",n=2)[,1])  +  
  as.numeric( str_split_fixed(G7.all[,"mia_3"],"/",n=2)[,1])+  as.numeric( str_split_fixed(G7.all[,"mia_5"],"/",n=2)[,1])  + 
  as.numeric( str_split_fixed(G7.all[,"mia_7"],"/",n=2)[,1])+ as.numeric( str_split_fixed(G7.all[,"mia_9"],"/",n=2)[,1])+ 
  as.numeric( str_split_fixed(G7.all[,"mia_11"],"/",n=2)[,1])  
G7.all$NonRed.Coverage=as.numeric( str_split_fixed(G7.all[,"mia_1"],"/",n=2)[,2])  +   
  as.numeric( str_split_fixed(G7.all[,"mia_3"],"/",n=2)[,2])+
  as.numeric( str_split_fixed(G7.all[,"mia_5"],"/",n=2)[,2])  +   as.numeric( str_split_fixed(G7.all[,"mia_7"],"/",n=2)[,2])+
  as.numeric( str_split_fixed(G7.all[,"mia_9"],"/",n=2)[,2])  +   as.numeric( str_split_fixed(G7.all[,"mia_11"],"/",n=2)[,2])




  
#add relevant information to G7.data
G7.data=merge(G7.data,G7.all[,c("ID","Pop.SegR","Pop.SegN","cov5","NonRed.MinorCount","Red.MinorCount","NonRed.Coverage","Red.Coverage")],by="ID",all.x=TRUE)
G7.data$R.avg=G7.data$Red.MinorCount/G7.data$Red.Coverage
G7.data$N.avg=G7.data$NonRed.MinorCount/G7.data$NonRed.Coverage

#write.csv(G7.data,"C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_allPerms_pval_info.csv")
G7.data=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_allPerms_pval_info.csv",header=TRUE)
#after filtering
G7.data$NR.avg=apply(G7.data[,c("R.avg","N.avg")],1,mean)
sum(G7.data$NR.avg>0.05)

#filter without cov5
G7.filtered=G7.data[G7.data$Pop.SegR>3&G7.data$NR.avg>0.05,]

sum(G7.filtered$p.adjust<0.05)
sum(G7.filtered$Perm1<0.05)
sum(G7.filtered$Perm2<0.05)
sum(G7.filtered$Perm3<0.05)
sum(G7.filtered$Perm4<0.05)
sum(G7.filtered$Perm5<0.05)
sum(G7.filtered$Perm6<0.05)
sum(G7.filtered$Perm7<0.05)
sum(G7.filtered$Perm8<0.05)
sum(G7.filtered$Perm9<0.05)
sum(G7.filtered$Perm10<0.05)
# 
# G7.filtered=G7.filtered[G7.filtered$cov5,]
# sum(G7.filtered$p.adjust<0.05)
# sum(G7.filtered$Perm1<0.05)
# sum(G7.filtered$Perm2<0.05)
# sum(G7.filtered$Perm3<0.05)
# sum(G7.filtered$Perm4<0.05)
# sum(G7.filtered$Perm5<0.05)
# sum(G7.filtered$Perm6<0.05)
# sum(G7.filtered$Perm7<0.05)
# sum(G7.filtered$Perm8<0.05)
# sum(G7.filtered$Perm9<0.05)
# sum(G7.filtered$Perm10<0.05)
# 

####-------- G7---Now examine gene level differences---###

SBGE.SSS.total.data <- read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SBGEandSSSdataForMBE.csv", header = TRUE)
SBGE.table<-SBGE.SSS.total.data[,-1] #SBGE.table is just abbreivated from the original file name

#Important file output from bedtools intersect that maps all SNPs present at generation G130 (coverage 40 to 500) to corresponding genes
#this contains multiple entries per SNP (entries are 'per bedtools feature')
backgroundFeatures=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130a_cov40-500v2_Background_genes_short.txt", sep=",", header = TRUE)
backgroundFeatures<-backgroundFeatures[-c(1)]

##add genes and information present in backgroundFeatures to SBGE.table 
overlap.geneIDX=SBGE.table[[1]]%in%backgroundFeatures[["geneID"]] 
#obtain a list where chromosome information is missing for a gene
missing.chrm.IDX=is.na(SBGE.table[,"chrm"]) 
helpIDX=overlap.geneIDX&missing.chrm.IDX 
help.gene.list=SBGE.table[helpIDX,"geneID"] 
#pull out the gene&chrm pairs for genes that appear in BOTH SBGE.table and backgroundFeatures, but where chromosome is missing
chrm.list=unique(subset(backgroundFeatures[,c("geneID","chrm")],backgroundFeatures[,"geneID"]%in%help.gene.list))
SBGE.subset=SBGE.table[helpIDX,c("geneID","chrm")]
merged.chrms=merge(chrm.list,SBGE.subset,by="geneID",all.x=FALSE)
SBGE.table[helpIDX,"chrm"]=merged.chrms[[2]]# 

#note there will be some duplicate IDs if a site overlaps 2 genes
test.genes=merge(G7.filtered,unique(backgroundFeatures[,c("ID","geneID")]),by="ID",all.x=TRUE)
length(unique(test.genes[,"geneID"]))
       
length(unique(test.genes[test.genes$p.adjust<0.05,"geneID"]))
length(unique(test.genes[test.genes$Perm1<0.05,"geneID"]))
length(unique(test.genes[test.genes$Perm2<0.05,"geneID"]))
length(unique(test.genes[test.genes$Perm3<0.05,"geneID"]))
length(unique(test.genes[test.genes$Perm4<0.05,"geneID"]))
length(unique(test.genes[test.genes$Perm5<0.05,"geneID"]))
length(unique(test.genes[test.genes$Perm6<0.05,"geneID"]))
length(unique(test.genes[test.genes$Perm7<0.05,"geneID"]))
length(unique(test.genes[test.genes$Perm8<0.05,"geneID"]))
length(unique(test.genes[test.genes$Perm9<0.05,"geneID"]))
length(unique(test.genes[test.genes$Perm10<0.05,"geneID"]))


(unique(test.genes[test.genes$p.adjust<0.05,"geneID"]))
(unique(test.genes[test.genes$Perm1<0.05,"geneID"]))
(unique(test.genes[test.genes$Perm2<0.05,"geneID"]))
(unique(test.genes[test.genes$Perm3<0.05,"geneID"]))
(unique(test.genes[test.genes$Perm4<0.05,"geneID"]))
(unique(test.genes[test.genes$Perm5<0.05,"geneID"]))
(unique(test.genes[test.genes$Perm6<0.05,"geneID"]))
(unique(test.genes[test.genes$Perm7<0.05,"geneID"]))
(unique(test.genes[test.genes$Perm8<0.05,"geneID"]))
(unique(test.genes[test.genes$Perm9<0.05,"geneID"]))
(unique(test.genes[test.genes$Perm10<0.05,"geneID"]))

#G7.real.genes=unique(backgroundFeatures[backgroundFeatures$ID%in%G7.P0.sigInfo[G7.P0.sigInfo$p.adjust<0.05,"ID"],"geneID"])

#populations combined
#####
#G7.MAF<-read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_totalCount_cov5_w_avg.txt")


#####







