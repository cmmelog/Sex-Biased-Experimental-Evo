library(tidyverse)
library(plotrix) #for like confidence intervals
library(sm) #for density plot comparsons
library(boot)#for bootstrapping
library(rstatix)
library(ggplot2)

#############==============================########################
######========= STEP ONE: Importing DATA  =======#####

##===Import a dataframe containing glm tested SNPs which identifies candidate SNPs
newtest=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130_newCMHchr2_full_noMAF.csv")
newtest=newtest[,-1]#remove "X" row
newtest$FDR5=newtest$p.adjust<0.05

#===important allele frequency data from generation 130===
G130.MAF=read.table("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130_MAF_cleaned_MAF1_Freq.csv",sep="",header=TRUE)

newtest$MAF5.G130=newtest$ID%in%G130.MAF[G130.MAF$MAF.Navg>0.05|G130.MAF$MAF.Ravg>0.05,"ID"]
newtest$MAF1.G130=newtest$ID%in%G130.MAF[G130.MAF$MAF.Navg>0.01|G130.MAF$MAF.Ravg>0.01,"ID"]

##===G7 MAF data===
G7.MAF<-read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_totalCount_cov5_w_avg.txt")
G7.MAF$N.avg=1-G7.MAF$N.avg
G7.MAF$R.avg=1-G7.MAF$R.avg
newtest$MAF5.G7=newtest$ID%in%G7.MAF[G7.MAF$Navg>0.05|G7.MAF$Ravg>0.05,"ID"]
newtest$MAF1.G7=newtest$ID%in%G7.MAF[G7.MAF$Navg>0.01|G7.MAF$Ravg>0.01,"ID"]

#import Drosophila simulans ( signor et al 2018 Call80, DP10, N=182) reference that has been converted to Drosophila melanogaster coordinates using liftOver
sim.frq3=read.table("C:\\Users\\user\\Documents\\Documents\\Genomics\\Dsim_dmel6_freqRI_Call80DP10.csv",header=TRUE)
Dsim.chr2=subset(sim.frq3,Chrom.r6=="2L"|Chrom.r6=="2R")
sample.Dsim=Dsim.chr2[,c("Chrom.r6","Pos.r6","MAF_sim","A1","A2")]
colnames(sample.Dsim)[1]="chrm"
colnames(sample.Dsim)[2]="pos"


#####=========================================
#####=============STEP TWO: BINNING=========####
#Add candidate information (passed FDR5 threshold) to each allele frequency data frame
G7.MAF$FDR5=G7.MAF$ID%in%newtest[newtest$FDR5,"ID"]
G130.MAF$FDR5=G130.MAF$ID%in%newtest[newtest$FDR5,"ID"]
G7.MAF$in.G130=G7.MAF$ID%in%G130.MAF$ID
#determine the average allele frequency
G7.MAF$NR.avg=(G7.MAF$R.avg+G7.MAF$N.avg)/2
G7.chr2=subset(G7.MAF,chrm=="2L"|chrm=="2R")
G130.chr2=subset(G130.MAF,chrm=="2L"|chrm=="2R")

#this function will downsample non-candidate SNPs so that the allele frequency distribution
#is the same as candidate SNPs
bin_by_MAF<-function(data,col,bin_size){
  #this creates a list of IDs where the portion of sampling will match FDR5
  down.sampled.nonC=c()
for(MAF in seq(0,0.5,bin_size)){
  window.subset=data[data[col]>MAF&data[col]<(MAF+bin_size),]
  cand.num=sum(window.subset$FDR5)
  non.csample=sample(window.subset[!window.subset$FDR5,"ID"],cand.num)
  down.sampled.nonC=append(down.sampled.nonC,non.csample)
  #  print(paste("MAF window=",MAF,"\n",length(window.cand$ID), "candidates",length(window.cand[window.cand$FDR5,"ID"])))
  print(paste("MAF=",MAF,"cand.num",cand.num))
  }
  return(down.sampled.nonC)
}


#####================================================
#####========= STEP THREE: Dsim FREQ =============#####
#create a SNP ID
sample.Dsim$ID=paste(sample.Dsim$chrm,sample.Dsim$pos)

#determine the trans-pecific status of each SNP-->see among our SNPs if they are also present
#in Dsimulans
G7.chr2$in.Dsim=G7.chr2$ID%in%sample.Dsim$ID
G130.chr2$in.Dsim=G130.chr2$ID%in%sample.Dsim$ID
#Make sure Dsim SNPs also pass the MAF cut-off
G7.chr2$in.Dsim.filt=G7.chr2$ID%in%sample.Dsim[sample.Dsim$MAF_sim>0.05,"ID"]
G130.chr2$in.Dsim.filt=G130.chr2$ID%in%sample.Dsim[sample.Dsim$MAF_sim>0.05,"ID"]

#downsample to meet MAF. Potential filering that identies if the same allele exists 
#at the same SNP position between Dsim and generation 7 or 130 data
G7.chr2$same.allele=FALSE
for(overlap.SNP in G7.chr2[G7.chr2$in.Dsim,"ID"]){
  if(grepl(sample.Dsim[sample.Dsim$ID==overlap.SNP,"A1"],
           G7.chr2[G7.chr2$ID==overlap.SNP,"allele_states"])){
    if(grepl(sample.Dsim[sample.Dsim$ID==overlap.SNP,"A2"],
             G7.chr2[G7.chr2$ID==overlap.SNP,"allele_states"])){
      #only assign as TRUE (the same) if BOTH A1 and A2 are present in allele_states
      G7.chr2[G7.chr2$ID==overlap.SNP,"same.allele"]=TRUE
    }}}

G130.chr2$same.allele=FALSE
for(overlap.SNP in G130.chr2[G130.chr2$in.Dsim,"ID"]){
  if(grepl(sample.Dsim[sample.Dsim$ID%in%overlap.SNP,"A1"],
           G130.chr2[G130.chr2$ID%in%overlap.SNP,"allele_states"])){
    if(grepl(sample.Dsim[sample.Dsim$ID%in%overlap.SNP,"A2"],
             G130.chr2[G130.chr2$ID%in%overlap.SNP,"allele_states"])){
      #only assign as TRUE (the same) if BOTH A1 and A2 are present in allele_states
      G130.chr2[G130.chr2$ID%in%overlap.SNP,"same.allele"]=TRUE
    } }} 

G7.chr2$same.Dsim=G7.chr2$in.Dsim&G7.chr2$same.allele #63/364
G130.chr2$same.Dsim=G130.chr2$in.Dsim&G130.chr2$same.allele# 679/3750
#if filtering for only sites where the same allele is present
G7.chr2$in.Dsim.sfilt=G7.chr2$same.allele& G7.chr2$ID%in%sample.Dsim[sample.Dsim$MAF_sim>0.05,"ID"]
G130.chr2$in.Dsim.sfilt=G130.chr2$same.allele&G130.chr2$ID%in%sample.Dsim[sample.Dsim$MAF_sim>0.05,"ID"]

test=bin_by_MAF(G7.chr2,"NR.avg",0.02) #0.4 cut-off
#downsample background SNPs for Generation 7
G7.downsample=G7.chr2[G7.chr2$ID%in%test|G7.chr2$FDR5,]
#downsample background SNPs for generation 130
G130.test=bin_by_MAF(G130.chr2,"MAF.NRavg",0.02) #0.4 cut-off
G130.downsample=G130.chr2[G130.chr2$ID%in%G130.test|G130.chr2$FDR5,]

#test the enrichment of whether trans-specific SNPs are more likely to occur in Candidates or not
#each allele frequency data frame (G7 or G130) has already been determined to be in Dsimulans or not
test=table(G130.downsample[G130.downsample$MAF.NRavg>0.05,c("FDR5","in.Dsim")])
test
chisq.test(test)
test[3]/sum(test[1,])*100
test[4]/sum(test[2,])*100


