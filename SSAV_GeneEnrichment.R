library(ggplot2)#
library(plotrix) #
library(boot)#
library(lme4)#
library(ggpubr)#
library(tidyverse)
library(rapport)
library(grid)##
library(rstatix)#
library(readr)#

####################
#-----FUNCTIONS----#
####################
bootmean.func <- function(x,i){mean(x[i])}
bootmean.func1<-function(x,i){
  test=boot(x[i],bootmean.func,R=10000)
  return(test[[1]])
}
boot.err.func1<-function(x,i){
  test=boot(x[i],bootmean.func,R=10000)
  #test1=boot(SBGE.cat.table[SBGE.cat.table$SBGEcat.body.Osada=="MB","rmf"],bootmean.func,R=10000)
  meanCI=t.test(test$t, conf.level = 0.95)$conf.int #allegedly this gives 95% CI of mean
  bootCI=boot.ci(boot.out = test, type = c( "perc"))#this gives percentile based CI
  CI=c(bootCI$percent[4],bootCI$percent[5]) #roughly encapsulated 2*standard error diff from mean
  #CI=quantile(test$t,c(0.025,0.975))
  #return(apply(test$t,2,sd)[1])
  return(CI)
}
#I do use this
boot_ftable<-function(data1,col1,col2){
  data_table=table(data1[,c(col1,col2)])
  fraction.table=data.frame(fraction=c(0),category=c(0))
  names=rownames(data_table)
  for (i in 0:dim(data_table)[1]){
    #print(paste(i,test[i,2],sum(test[i,2])))
    fraction.table[i,"fraction"]=sum(data_table[i,2])/sum(data_table[i,])
    fraction.table[i,"category"]=names[i]
  }
  return(fraction.table)}
#pretty sure I don't call these. Delete later.
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
  lower_graph=min(diff.SigNon)*1.1
  max=range(c(dSig.fit2,dNSig.fit2))
  upper_graph=max[2]
  
  plot(dSig.chr2, col=2,ylab="Density",xlim=c(u[1],u[2]),xlab=col,
       main=paste('Difference in',col, "Distribution" ),ylim=c(lower_graph,upper_graph))#,xlim=c(-10,15))
  #plot(dNSig.chr2 , col=1, ylim=c(0, .20), main='Difference in SBGE Distribution', xlab='',ylab="Density",xlim=c(-10,15))
  polygon( c(dNSig.chr2$x, rev(dNSig.chr2$x)), c(dNSig.fit3[1,], rev(dNSig.fit3[2,])),
           col='lightgrey', density = -0.5, border=F)
  polygon( c(dSig.chr2$x, rev(dSig.chr2$x)), c(dSig.fit3[1,], rev(dSig.fit3[2,])),
           col='lightgreen', density = -0.5, border=F)
  #same x as normal plots but different y
  polygon( c(dSig.chr2$x, rev(dSig.chr2$x)), c(diff.CI[1,], rev(diff.CI[2,])),
           col=2, density = -0.5, border=F)
  abline(h=0, lty=3, col=8) # lty=2, lwd=2,
  lines(dNSig.chr2, col=1,lwd=2)
  lines(dSig.chr2, col=3,lwd=2)
  lines(dSig.chr2$x, diff.SigNon , col="darkred", lty=1, lwd=2)  ## <---------------- difference#
  #mtext(sprintf('N(x) = %s  Bandwidth(x) = %s', dSig.chr2$n, signif(dSig.chr2$bw, 3)), 1, 2)
  #mtext(sprintf('N(y) = %s  Bandwidth(y) = %s', dNSig.chr2$n, signif(dNSig.chr2$bw, 3)), 1, 3)
  legend(leg_pos, legend=legends, col=c(2,1,3), 
         lty=c(1, 1, 2), lwd=c(1, 1, 2))
}#Difference, False, True
graph.point.distribution<-function(point,distribution,col,title2,legend){
  plotting_values=data.frame(means=distribution)
  x_title=paste("mean",col)
  
  ggplot(plotting_values,aes(x=means))+ggtitle(title2)+theme_bw() +
    #xlim(c(1.265,1.335))+
    theme(text = element_text(size = 16),legend.position="top")+
    theme(plot.title = element_text(size = 22)) + xlab(x_title)+
    #geom_density()+
    geom_histogram(bins=25,colour="black",fill="lightgrey",show.legend=TRUE)+
    geom_vline(aes(xintercept=point), color="blue", linetype="dashed", size=1,show.legend=TRUE)}

#This converts per feature information to per gene
#takes a bedtools intersect output txt file which lists SNPs per "Feature"
#function requires some format editing from output file including creating a column "geneID" containing the Flybase ID name
unstack.geneTable<-function(table){ 
  #pulls out the genes
  genes=unique(table[["geneID"]])
  #this creates a blank table where columns respond to quetions of interst about the SNPS
  #is.exon would be: are the SNPs here exons
  unstack.table<-data.frame(geneID=character(),chrm=character(),is.CDS=logical(),is.exon=logical(),is.5UTR=logical(),is.3UTR=logical(),is.pseudogene=logical(),
                            is.stop_codon=logical(),is.start_codon=logical(),is.RNA=logical(),num.positions=double())
  #for each gene, the SNPs/features corresponding to that gene are pulled out and information about the them is recorded
  for (i.gene in genes){ 
    current.gene=subset(table,geneID==i.gene)
    num.positions=length(unique(current.gene[['position']]))
    featureType=unique(current.gene[['featureType']])
    is.CDS=FALSE
    if ("CDS"%in%featureType) {is.CDS=TRUE}
    is.exon=FALSE
    if ("exon"%in%featureType) {is.exon=TRUE}
    is.5UTR=FALSE
    if ("5UTR"%in%featureType) {is.5UTR=TRUE}
    is.3UTR=FALSE
    if ("3UTR"%in%featureType) {is.3UTR=TRUE}
    is.pseudogene=FALSE
    if ("pseudogene="%in%featureType) {is.pseudogene==TRUE}
    is.stop_codon=FALSE
    if ("stop_codon"%in%featureType) {is.stop_codon=TRUE}
    is.start_codon=FALSE
    if ("start_codon"%in%featureType) { is.start_codon=TRUE}
    is.RNA=FALSE #setting defeault value
    any.RNA=str_detect(featureType,"RNA")#this detects which elements and where are RNA
    if( TRUE%in%any.RNA) {#this will be true if ANY TRUE in above 
      is.RNA=TRUE }
    new.entry=c(i.gene,current.gene[1,1],is.CDS,is.exon,is.5UTR,is.3UTR,is.pseudogene,is.stop_codon,is.start_codon,is.RNA,num.positions)
    #f = rbind(df, output)}
    unstack.table=rbind(unstack.table,new.entry)
  }
  colnames(unstack.table)<-c("geneID","chrm","is.CDS","is.exon","is.5UTR",'is.3UTR',"is.pseudogene","is.stop_codon","is.start_codon","is.RNA","num.positions")
  return(unstack.table)
}

##############################-
####----Import Files------####
##############################-
#import Amerdeep Singh's chart SBGEandSSSdataForMBE constructed from Singh and Agrawal(2023)
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



#=============================================================#
####==========================Candidate FDR5===================#######
#this imports all SNPs on chromosome 2 that have undergone significance testing undering script SSAV_frequencyDifferences
#since the calculation of significane is independent of other SNPs, filtering afterward will not affect raw p.value
newtest=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130_newCMHchr2_full_noMAF.csv")
newtest=newtest[,-1]#remove "X" row
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

#######-----Creation of Null Distributions#####
#Significant Filtered SNPs are shifted using ShiftPosition.R script to create 1000 permutations of the same number of filtered SNPs
#these SNPs are then converted to genes using bedtools intersect and put into a table using PermToGeneTable.R
all.shifted.genes=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130_permshiftCHR2_NR.MAF5_filt_THENSHIFT_geneTable.csv")
all.shifted.genes=all.shifted.genes[,-1]#removes row names
#The first permutation represents the actual candidate diverged SNPs and is removed
all.shifted.genes$Perm..0.5<-NULL

#XXXXXX
null_dist.PermNum<-function(data,null,col,title){
  #Cand.mean=mean(data[cand.IDX,col],na.rm=TRUE)
  #pulls out each permutation number (minus the original data at 0)
  permutations=unique(null[null$PermNum>0,"PermNum"])
  null.means=c()
  count=0
  for(perm in permutations){
    null.subset=subset(null,PermNum==perm)
    data.subset=data[data$geneID%in%null.subset$geneID,col]
    subset.mean=mean(data.subset,na.rm=TRUE)
    count=count+1
    null.means[count]=subset.mean
  }
  
  #graph.point.distribution(Cand.mean,null.means,col,paste(title,"compared to null distribution"),
  #                         c("Candidate mean"))
  return=null.means
}

#this returns the null distribution of a given metric (here specified by col). One value of col
#is calculated per permutated
#null refers to the permutated gene list
null_dist<-function(data,null,col){
  #pulls out each permutation number (minus the original data at 0)
  permutations=colnames(null)
  null.means=c()
  count=0
  for(perm in permutations){
    subset.IDX=null[[perm]]
    data.subset=data[data$geneID%in%null[subset.IDX,"geneID"],col]
    subset.mean=mean(data.subset,na.rm=TRUE)
    count=count+1
    null.means[count]=subset.mean
  }
  
  #graph.point.distribution(Cand.mean,null.means,col,paste(title,"compared to null distribution"),
  #                         c("Candidate mean"))
  #there's a random null and it says the length is 101 instead of 100
  return=null.means[-1]
}
#the null density distribution: returns one density distribution of the specified column per permutation
#since the x values are the same, this ONLY returns the y values of density
null_density_dist<-function(data,null,col,u){
  #u is the united range
  Fake.density=c()
  #this ensures it doesn't use geneID 
  total.perm.Num=colnames(null[,-1])
  
  
  for(perm in total.perm.Num){
    print(perm)
    #null[[perm]] is the gene IDX and the first col geneID has the genes
    temp.genes=null[null[[perm]],"geneID"]
    #print(paste("sum",sum(null[[perm]])))
    #temp.genes=null[[perm]]
    temp.Fake.Cand <-data[data$geneID%in%temp.genes,] #SBGE.chr2$totalExonLength
    temp.Fake.Cand<-temp.Fake.Cand[!is.na(temp.Fake.Cand[[col]]),col]
    #density for temp sample
    temp.density=density(temp.Fake.Cand,from=u[1], to=u[2])$y
    #add the y values to a growing list
    #is this rbind or cbind
    #print(temp.density[1])
    Fake.density=cbind(Fake.density,temp.density)
  }
  
  return(Fake.density)
}

######---------------------------------------------------------------------------------#
#######-----Comparison of Null to Diverged Genes: Sex-Bias Gene Distribution---------#####

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

#this creates the fake SBGE density distribution (from the null)
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

plot(Real.Cand.dens, col="white",ylab="Density",xlab="Sex-Biased Gene Expression(SBGE)",xlim=c(-10,15),#xlim=c(u[1],u[2]),
     main=paste('Difference in',col,' SBGE Distribution' ),ylim=c(lower_graph*0.95,upper_graph*0.95))
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
              
              lty=c(1, 1, 2), lwd=c(1, 1, 2))


####-----------------SBGE Enrichment per SBGE category---####
#examining diverged candidate encirhment over discrete sex biased gene categories 
#for our data set get the fraction of significant candidate diverged genes in each category
read.cand=(table(SBGE.chr2[SBGE.chr2$Candidates,"SBGEcat.body.Osada"])
           /sum(table(SBGE.chr2[SBGE.chr2$Candidates,"SBGEcat.body.Osada"])))*100
#convert to data.frame for ease of plotting
plot.cand=as.data.frame((read.cand))
colnames(plot.cand)=c("SBGE.cat","means")
plot.cand$SBGE.cat<-factor(plot.cand$SBGE.cat,levels=c("extFB","sFB","FB","UB","MB","sMB","extMB"))

#for the null distribution of candidate lists, find the null fraction for each SBGE category
null=all.shifted.genes
#record a list of all the permutations (the column names)
permutations=colnames(all.shifted.genes)
#remove the "geneID" column
permutations=permutations[2:1000]
SBGE.category=c("extFB","sFB","FB","UB","MB","sMB","extMB")
null.means=data.frame()
count=0
#calling this variable data for ease of use
data=SBGE.chr2
#ordering the variables
data$SBGEcat.body.Osada<-factor(data$SBGEcat.body.Osada,levels=c("extFB","sFB","FB","UB","MB","sMB","extMB"))
#for each permutation in all.shifted.genes this calculates the fraction of genes in each SBGE category
for(perm in permutations){
  #take only the genes that appear in a given permutation
  subset.IDX=null[[perm]]
  data.subset=data[data$geneID%in%null[subset.IDX,"geneID"],]
  #calculate the fraction genes within each SBGE for that fraction
  temp.percent=(table(data.subset$SBGEcat.body.Osada)/
                  sum(table(data.subset$SBGEcat.body.Osada)))*100
  #because this is 7 values in an order, we want it to be in a list not a df
  null.means=rbind(null.means,temp.percent)
  count=count+1}
#rename the column names of null means so they correspond to the categories rather than numbers
colnames(null.means)=SBGE.category

#unstack the null means for ease of graphing later
unstack.SBGE=data.frame(means=c(),SBGE.cat=c())
for(cat in 1:length(null.means)){
  temp.distr=data.frame(means=null.means[[cat]],SBGE.cat=rep(SBGE.category[cat],length(null.means[[cat]])))
  unstack.SBGE=rbind(unstack.SBGE,temp.distr)  }
#ordering the categories
unstack.SBGE$SBGE.cat<-factor(unstack.SBGE$SBGE.cat,levels=c("extFB","sFB","FB","UB","MB","sMB","extMB"))

#graphing the distribution of SBGE fractions compared to the real diverged gene mean
ggplot(unstack.SBGE,aes(x=SBGE.cat,y=means,fill=SBGE.cat))+
  scale_fill_manual(values=c("firebrick4","red2","orange","lightgoldenrod","palegreen1","steelblue2","royalblue4"))+
  geom_violin(trim=FALSE,position=position_dodge(),alpha=0.5,scale="width",linewidth=0.6,colour="black")  +
  ylab("mean % of genes")+xlab("SBGE categories")+scale_y_continuous(breaks=seq(0,27,5))+
  geom_point(inherit.aes=FALSE,data=plot.cand, 
             aes(x=SBGE.cat,y=abs(means),fill=SBGE.cat),
             stat="identity",position=position_nudge(0),size=3,pch=21,fill="black")+
  theme_classic() +theme(text = element_text(size = 16),legend.position="none")+
  theme(plot.title = element_text(size = 20)) +
  ggtitle("% of genes: null vs candidates") 

######-----------SBGE Distribution by MAF window---------

#For each minor allele freq bin, generate the same graph was done for the entire chr 2 dataset
for(freq in seq(0.05,0.45,by=0.05)){
  #an IDX for SNPs with a minor allele frequency (calculated as average between Red and NonRed sample types) within specified window
  windowMAF.IDX=G130.MAF[,"MAF.NRavg"]>freq&G130.MAF[,"MAF.NRavg"]<freq+0.05
  #generate an index to convert SNP ID from G130 MAF dataset to "backgroundFeatures" dataset which links SNPs to genes they occur in
  windowgene1.IDX=backgroundFeatures$ID%in%G130.MAF[windowMAF.IDX,"ID"]
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
       main=paste('Difference in SBGE Distribution',freq,"to",(freq+0.05)),ylim=c(lower_graph*0.95,upper_graph*0.95))
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


#######------Comparison of Null to Diverged Genes:rmf (general)-------------------------####
#examining whether diverged genes have higher rmf
#get the null ditribution of rmf values
rmf.null=null_dist(SBGE.chr2[,],all.shifted.genes,"rmf")
#this specifies the p<0.05 cut-off for the null distribution
one.tailed.cutoff=quantile(rmf.null,0.95) 
two.tailed.cutoff=quantile(rmf.null,0.975)
actual.mean=mean(SBGE.chr2[SBGE.chr2$Candidates,"rmf"],na.rm=TRUE)
#reformats the data fr plotting with ggplot
rmf.null.for.graph=data.frame(rmf.means=rmf.null,perm=seq(1:length(rmf.null)),toy=rep(TRUE,length(rmf.null)))

#PLOTTING
ggplot(rmf.null.for.graph,aes(x=toy,y=rmf.means,fill=toy),fill="grey")+
  geom_violin(trim=TRUE,position=position_dodge())  +
  ylab("mean rmf")+ylim(0.3,0.5)+
  geom_point(inherit.aes=FALSE,data=SBGE.chr2[SBGE.chr2$Candidates,], 
             aes(x=Candidates,y=actual.mean),
             stat="identity",position=position_nudge(0),size=5,pch=21,fill="blue")+
  theme_bw() +theme(text = element_text(size = 16),legend.position="none")+
  theme(plot.title = element_text(size = 20)) +
  ggtitle("Mean rmf:\n Null vs Diverged genes") 

#####----------------rmf per SBGE catgory==============
#Examining whether candidates have elevated when broken down in discrete sex-biased expression categories
#generate rmf, null, for SBGE
#this represents the diverged candidates mean rmf per category
rmf.by.SBGE=data.frame(cat=c(),rmf.mean=c())
#will will be the null distribution of each category, a list of distributions
rmf.by.SBGE.null=c()
count=1
col="rmf"
#for each of the seven SBGE categories in the data, calculate the null distribution of mean rmf
for(SBGE.cat in c("extFB","sFB","FB","UB","MB","sMB","extMB")){
 
  rmf.by.SBGE[count,"cat"]=SBGE.cat
  #take the genes within the SBGE category
  SBGE.sub=subset(SBGE.chr2[,],SBGEcat.body.Osada==SBGE.cat)
  #calculate the actual mean rmf of diverged candidate genes in those categories
  rmf.by.SBGE[count,"rmf.mean"]=mean(SBGE.sub[SBGE.sub$Candidates,col],na.rm=TRUE)
  #for the specific SBGE category generate the null distribution of rmf
  temp.distr=null_dist(SBGE.sub,all.shifted.genes,col)
  #add the null distribution to the list of null distributions
  rmf.by.SBGE.null[count]=data.frame(distr=temp.distr)
  count=count+1
}

#unstack the values of mean rmf fore each SBGE category so that they can be graphed
#unstack rmf.by.SBGE2
SBGE.category=c("extFB","sFB","FB","UB","MB","sMB","extMB")
unstack.null.rmfBySBGE=data.frame(rmf.means=c(),SBGE.cat=c())
for(cat in 1:length(rmf.by.SBGE.null)){
  temp.distr=data.frame(rmf.means=rmf.by.SBGE.null[[cat]],SBGE.cat=rep(SBGE.category[cat],length(rmf.by.SBGE.null[[cat]])))
  unstack.null.rmfBySBGE=rbind(unstack.null.rmfBySBGE,temp.distr)  }
#ordering the categories
unstack.null.rmfBySBGE$SBGE.cat<-factor(unstack.null.rmfBySBGE$SBGE.cat,levels=c("extFB","sFB","FB","UB","MB","sMB","extMB"))

#generating the upper cutoff for rmf
rmf.by.SBGE.cutoffs=lapply(rmf.by.SBGE.null,quantile,(0.975))
rmf.by.SBGE.cutoffs.lower=lapply(rmf.by.SBGE.null,quantile,(0.025))

#PLOTTING:graphing a comparison of the null distribution of rmf per each SBGE category
ggplot(unstack.null.rmfBySBGE,aes(x=SBGE.cat,y=rmf.means,fill=SBGE.cat))+
  scale_fill_brewer(palette="OrRd",direction=-1)+
  #scale_fill_brewer(palette="RdYlBu")+
  geom_violin(trim=FALSE,position=position_dodge())  +
  ylab("mean rmf")+xlab("SBGE categories")+

  theme_bw() +theme(text = element_text(size = 16),legend.position="none")+
  theme(plot.title = element_text(size = 20)) +
  ggtitle("rmf by SBGE categories: null vs candidates") 

#Creating a combined figure of general rmf and SBGE
rmf.null.for.graph$SBGE.cat="ALL"
unstack.null.rmfBySBGE=rbind(unstack.null.rmfBySBGE,rmf.null.for.graph[,c("rmf.means","SBGE.cat")])
unstack.null.rmfBySBGE$SBGE.cat<-factor(unstack.null.rmfBySBGE$SBGE.cat,levels=c("ALL","extFB","sFB","FB","UB","MB","sMB","extMB"))

#Add the general mean rmf for all candidatess
rmf.by.SBGE[8,"rmf.mean"]=as.numeric(mean(SBGE.chr2[SBGE.chr2$Candidates,"rmf"],na.rm=TRUE))
rmf.by.SBGE[8,"cat"]="ALL"

#PLOTTING
ggplot(unstack.null.rmfBySBGE,aes(x=SBGE.cat,y=rmf.means,fill=SBGE.cat))+
  scale_fill_manual(values=c("darkgrey","firebrick4","red2","orange","lightgoldenrod","lightskyblue1","steelblue3","royalblue4"))+
  geom_violin(trim=FALSE,position=position_dodge(),alpha=0.8,linewidth=0.55)  +
  ylab("mean rmf")+xlab("SBGE categories")+
  geom_point(inherit.aes=FALSE,data=rmf.by.SBGE, 
             aes(x=cat,y=abs(rmf.mean),fill=cat),
             stat="identity",position=position_nudge(0),size=3,pch=21,fill="black")+
  theme_bw() +theme(text = element_text(size = 16),legend.position="none")+
  theme(plot.title = element_text(size = 20)) +geom_vline(xintercept=1.6,colour="black")+
  ggtitle("rmf by SBGE categories: null vs candidates") 


####------Test overlap with Grieshop et al., 2025-----###
#examining overlap in transcriptomic and genomic candidates derrived by the same populations
M.Liu.data<-read.table("C:\\Users\\user\\Documents\\Documents\\Genomics\\All.geno_candidates.tsv")
#reformate column names to be consistent
colnames(M.Liu.data)=M.Liu.data[1,]
M.Liu.data=M.Liu.data[-1,]

Liu.data=M.Liu.data[,c("FlyBaseID","A.m.Sig","A.f.Sig","Sig")]
colnames(Liu.data)[1]="geneID"
colnames(Liu.data)[4]="DE.genes"

SBGE.chr2=merge(SBGE.chr2,Liu.data,by.x="geneID",all.x=TRUE)
#actual test of enrichment for diverged Candidates in DE genes (Grieshop et al., 2025 Candidates)
chisq.test(table(SBGE.chr2[,c("Candidates","DE.genes")]))


#determining the bootstraped fraction of overlap between DE genes and diverged Candidates
#Fraction of DE/Not that are genomic candidates
DE_boot_fraction=replicate(1000,{
  #Sample with replacement (for bootstrap from original dataset).
  resampID <- sample(SBGE.chr2[,"geneID"], replace=TRUE)
  resamp=SBGE.chr2[SBGE.chr2$geneID%in%resampID,]
  #Generate the distribution from the resampled dataset, and extract y coordinates to generate variablity bands
  boot_ftable(resamp,"Sig.transcript","Candidates")$fraction})
row.names(DE_boot_fraction)=c("not.DE","DE.Genes")

#fraction of Genomic Candidates/Background that are DE
G_boot_fraction=replicate(1000,{
  #Sample with replacement (for bootstrap from original dataset).
  resampID <- sample(SBGE.chr2[,"geneID"], replace=TRUE)
  resamp=SBGE.chr2[SBGE.chr2$geneID%in%resampID,]
  #generate the fraction table for Candidates
  boot_ftable(resamp,"Candidates","Sig.transcript")$fraction})
row.names(G_boot_fraction)=c("Background","Candidates")

####--------------------------------SIDU------------------####
#examining enrichment for genes with known sex specific splicing
#a generate chi-squared test between Candidate/Non-Candidates genes 
SDIU.test= table(SBGE.chr2[,c("Candidates","SDIU.body.sig")])
chisq.test(SDIU.test)

#calculates the fraction across categories
boot_ftable<-function(data1,col1,col2){
  data_table=table(data1[,c(col1,col2)])
  fraction.table=data.frame(fraction=c(0),category=c(0))
  names=rownames(data_table)
  #get the fraction of category per each one
  for (i in 0:dim(data_table)[1]){
    fraction.table[i,"fraction"]=sum(data_table[i,2])/sum(data_table[i,])
    fraction.table[i,"category"]=names[i]
  }
  return(fraction.table)}

#this will generate row1=FALSE (not sig) row2=TRUE (in Candidate) % of genes overlap with Candidates
boot_fraction=replicate(1000,{
  resampID <- sample(SBGE.chr2[,"geneID"], replace=TRUE)
  resamp=SBGE.chr2[SBGE.chr2$geneID%in%resampID,]
  #Generate the density from the resampled dataset, and extract y coordinates to generate variablity bands
  boot_ftable(resamp,"Candidates","SDIU.body.sig")$fraction})

SDIU.enrich=data.frame(category=c("Background","Candidate"),CI_up=apply(boot_fraction,1,quantile,probs=0.975),
                       CI_low=apply(boot_fraction,1,quantile,probs=0.025),mean_SDIU=boot_ftable(SBGE.chr2,"Candidates","SDIU.body.sig")$fraction)

#measures general likelihood of being SDIU
chance=table(SBGE.chr2[,"SDIU.body.sig"])["TRUE"]/sum(table(SBGE.chr2[,"SDIU.body.sig"]))
boot_chance=replicate(1000,{
  resampID <- sample(SBGE.chr2[,"geneID"], replace=TRUE)
  resamp=SBGE.chr2[SBGE.chr2$geneID%in%resampID,]
  #Generate the density from the resampled dataset, and extract y coordinates to generate variablity bands
  table(resamp[,"SDIU.body.sig"])["TRUE"]/sum(table(resamp[,"SDIU.body.sig"]))})
#####same as above but slightly different???  
boot_fraction=replicate(1000,{
  resampID <- sample(SBGE.chr2[,"geneID"], replace=TRUE)
  resamp=SBGE.chr2[SBGE.chr2$geneID%in%resampID,]
  #Generate the density from the resampled dataset, and extract y coordinates to generate variablity bands
  boot_ftable(resamp,"SDIU.body.sig","Candidates")$fraction})

SDIU.enrich2=data.frame(category=c("Not SDIU","SDIU"),CI_up=apply(boot_fraction,1,quantile,probs=0.975),
                        CI_low=apply(boot_fraction,1,quantile,probs=0.025),mean_SDIU=boot_ftable(SBGE.chr2,"SDIU.body.sig","is.Sig")$fraction)


enrichment_plotting=SDIU.enrich2
enrichment_plotting=rbind(enrichment_plotting,c("Chance",quantile(boot_chance,0.975),
                                                quantile(boot_chance,0.025),chance))

ggplot(enrichment_plotting,
       aes(x=category,fill=category,y=as.numeric(mean_SDIU)),ylim=c(0,0.6))+
  #geom_point(stat="identity",position=position_dodge(0.9),size=6,pch=21)+
  geom_bar(stat="identity",position=position_dodge(0.9),size=1,
           colour="black")+
  scale_y_continuous("fraction that are SDIU",limits = c(0, 0.75)) +
  scale_colour_manual(values=c("black","black","black","black","black","black","black","black"))+
  scale_fill_brewer(palette="Set1")+ 
  ggtitle("(chr2 only) w/95% CI\nSDIU splice enrichment")+theme_bw() + 
  theme(text = element_text(size = 16),legend.position="bottom")+
  theme(plot.title = element_text(size = 22)) + 
  geom_errorbar(aes(ymin=as.numeric(CI_low), ymax=as.numeric(CI_up)),
                width=.4,position=position_dodge(.9),show.legend = FALSE)


##################################################################-
