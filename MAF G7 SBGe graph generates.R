#MAF G7 SBGe graph generates

library(tidyverse)
library(plotrix) #for like confidence intervals
library(sm) #for density plot comparsons
library(boot)#for bootstrapping
library(rstatix)
library(ggplot2)

G7.MAF<-read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_totalCount_cov5_w_avg.txt")
G7.MAF$N.avg=1-G7.MAF$N.avg
G7.MAF$R.avg=1-G7.MAF$R.avg
G7.MAF$MAF.NRavg=(G7.MAF$N.avg+G7.MAF$R.avg)/2
G7.MAF$ID=paste(G7.MAF$chrm,G7.MAF$pos)

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


#now I have to determine candidates as from G130

####==========================Candidate FDR5===================#######
#this imports all SNPs on chromosome 2 that have undergone significance testing undering script SSAV_frequencyDifferences
#since the calculation of significane is independent of other SNPs, filtering afterward will not affect raw p.value
newtest=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130_newCMHchr2_full_noMAF.csv")
newtest=newtest[,-1]#remove "X" row
#old p.adjust
newtest$FDR5=newtest$p.adjust<0.05


#any gene with 1 significant SNP becomes a Candidate gene: this pulls out unique diverged SNPs
genes.1SNP=unique(backgroundFeatures[backgroundFeatures$ID%in%newtest$ID,"geneID"])
#this pulls out unique genes that pass the FDR5 cut-off (pre-filtering)
FDR5.genes.1SNP=unique(backgroundFeatures[backgroundFeatures$ID%in%newtest[newtest$FDR5,"ID"],"geneID"])

#####=========================Minor Allele Filtering at Generation 130=#####
#imports table containg allele frequencies for each SNP in each sample type for each population
G130.MAF=read.table("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130_MAF_cleaned_MAF1_Freq.csv",header=TRUE)

#We filter out SNPs that segregate in less than half the populations
#this calculates how many populations a SNP segrates in
G130.MAF$Seg.numR=rep(0,length(G130.MAF$ID))#number of Red populations where a variant is segregating
for(col in seq(5,15,2)){#this iterates over the frequency columns with Red Sample types
  #we only look at red sample types because a SNP can't be segregating in N without also segregating in R
  seg.IDX=G130.MAF[,col]!=1
  G130.MAF[seg.IDX,"Seg.numR"]=G130.MAF[seg.IDX,"Seg.numR"]+1}
#specifying the cut off of how many populations a SNP must segrate in
newtest$seg3.threshold=newtest$ID%in%G130.MAF[G130.MAF$Seg.numR>3,"ID"]
seg3threshold.genes.1SNP=unique(backgroundFeatures[backgroundFeatures$ID%in%newtest[newtest$seg3.threshold,"ID"],"geneID"])#

#Filter of the Minor Allele Frequency determined only using Red genotypes
newtest$R.MAF5=newtest$ID%in%G130.MAF[G130.MAF$MAF.Ravg>0.05,"ID"]#
R.MAF5.genes.1SNP=unique(backgroundFeatures[backgroundFeatures$ID%in%newtest[newtest$R.MAF5,"ID"],"geneID"])#
#Filter of the Minor Allele Frequency determined only using both genotypes
newtest$N.MAF5=newtest$ID%in%G130.MAF[G130.MAF$MAF.NRavg>0.05,"ID"]#
N.MAF5.genes.1SNP=unique(backgroundFeatures[backgroundFeatures$ID%in%newtest[newtest$N.MAF5,"ID"],"geneID"])#

G130.MAF$FDR5=G130.MAF$ID%in%newtest[newtest$FDR5,"ID"]

#specfication of all genes containing significant SNPs after both filtering steps and having a p.adjust<0.05
Candidates.gene.1SNP=unique(backgroundFeatures[backgroundFeatures$ID%in%newtest[newtest$N.MAF5&newtest$seg3.threshold&newtest$FDR5,"ID"],"geneID"])#
#actual candidates can only be on chromosome 2 due to the experimental selection protocol
SBGE.chr2=subset(SBGE.table,chrm=="2L"|chrm=="2R")
SBGE.chr2$Candidates=SBGE.chr2$geneID%in%Candidates.gene.1SNP
#######----ACTUAL MAF BINNING-----

#create the little bin again?

#graph
for(freq in seq(0.05,0.45,by=0.1)){
  #an IDX for SNPs with a minor allele frequency (calculated as average between Red and NonRed sample types) within specified window
  windowMAF.IDX=G7.MAF[,"MAF.NRavg"]>freq&G7.MAF[,"MAF.NRavg"]<freq+0.01
  #generate an index to convert SNP ID from G7 MAF dataset to "backgroundFeatures" dataset which links SNPs to genes they occur in
  windowgene1.IDX=backgroundFeatures$ID%in%G7.MAF[windowMAF.IDX,"ID"]
  windowgene.list=unique(backgroundFeatures[windowgene1.IDX,"geneID"])
  #subset the data frame SBGE chr2 to only include genes which had a SNP occuring in the minor allele window
  window.SBGE=SBGE.chr2[SBGE.chr2$geneID%in%windowgene.list,]
  #print(paste(sum(window.SBGE$Candidates),freq))
  
  col="Whole.SBGE.Osada"
  #generate bootstrapped density distribution of SBGE for candidate genes
  #generate a data.frame of SBGE.chr2 with only candidates
  Real.Cand <-window.SBGE[window.SBGE$Candidates,] 
  Real.Cand<-Real.Cand[!is.na(Real.Cand[[col]]),col]
  #This gets the united range for ALL genes on chromosome 2 so range will be equal and comparable
  u=range(window.SBGE[!is.na(window.SBGE[[col]]),col])
  #the real candidate SBGE density is obtained
  Real.Cand.dens<-density(Real.Cand,from=u[1], to=u[2])
  #the candidate distribution is bootstrapped to obtain 95% confidence intervals
  Real.Cand.dens.fit2 <- replicate(1000,{
    #Sample with replacement (for bootstrap from original dataset). Save the resample to x
    samp <- sample(Real.Cand, replace=TRUE)                    
    #Generate the density from the resampled dataset, and extract y coordinates to generate variablity bands
    density(samp, from=min(Real.Cand.dens$x), to=max(Real.Cand.dens$x))$y}) 
  #this represents the confidence intervals for the candidate distribution
  Real.Cand.dens.fit3 <- apply(Real.Cand.dens.fit2, 1, quantile, c(0.025,0.975) )
  
  #repeat the same step as above using Non-Candidate/Background genes
  Fake.Cand <-window.SBGE[!window.SBGE$Candidates,] 
  Fake.Cand<-Fake.Cand[!is.na(Fake.Cand[[col]]),col]
  #the real SBGE density is obtained
  Fake.Cand.dens<-density(Fake.Cand,from=u[1], to=u[2])
  #the real distribution is bootstrapped to obtain 95% confidence intervals
  Fake.Cand.dens.fit2 <- replicate(1000,{
    #Sample with replacement (for bootstrap from original dataset). Save the resample to x
    samp <- sample(Fake.Cand, replace=TRUE)                    
    #Generate the density from the resampled dataset, and extract y coordinates to generate variablity bands
    density(samp, from=min(Fake.Cand.dens$x), to=max(Fake.Cand.dens$x))$y}) 
  Fake.density.fit3 <- apply(Fake.Cand.dens.fit2, 1, quantile, c(0.025,0.975) )
  
  #Taking the difference between the background distribution and the candiates distribution
  single.diff=Real.Cand.dens$y-Fake.Cand.dens$y
  single.diff.CI.1=Real.Cand.dens$y-Fake.density.fit3[1,]
  single.diff.CI.2=Real.Cand.dens$y-Fake.density.fit3[2,]
  
  #data ranges for graphing purposes
  max=range(c(Real.Cand.dens.fit2,Fake.Cand.dens.fit2))
  upper_graph=max[2]
  lower_graph=min(single.diff.CI.2=Real.Cand.dens$y-Fake.density.fit3[2,])
  
  plot(Real.Cand.dens, col="white",ylab="Density",xlab="Sex-Biased Gene Expression(SBGE)",xlim=c(-10,15),#xlim=c(u[1],u[2]),
       main=paste('Difference in SBGE Distribution',freq,"to",(freq+0.1)),ylim=c(lower_graph*0.95,upper_graph*0.95))
  #polygon( c(Real.Cand.dens$x, rev(Real.Cand.dens$x)), c(Fake.density.fit3[1,], rev(Fake.density.fit3[2,])),
  #         col='grey', density = -0.5, border=F)
  polygon( c(Real.Cand.dens$x, rev(Real.Cand.dens$x)),  c(single.diff.CI.1, rev(single.diff.CI.2)),
           col='pink', density = -0.5, border=F)
  abline(h=0, lty=3, col=8) 
  lines(Real.Cand.dens, col="blue",lwd=2)
  #the Background distribution
  lines(Real.Cand.dens$x,Fake.Cand.dens$y, col="black",lwd=2)
  legend("topright", legend=c("Candidates","Background","Difference"), col=c("blue","black","pink"), 
         lty=c(1, 1, 2), lwd=c(1, 1, 2))
  
}
