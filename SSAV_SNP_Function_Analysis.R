library(tidyverse)
library(boot)

library(rapport)
library(lme4)#requires package Matrix
library(vcd)#requires package grid
library(car)
library(rstatix)
library(ggplot2)

##===STEP ONE===IMPORT DATA
#
SBGE.table <- read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SBGEandSSSdataForMBE.csv", header = TRUE)
SBGE.table <-SBGE.table[,-1]


backgroundFeatures=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130a_cov40-500v2_Background_genes_short.txt", sep="", header = FALSE)
colnames(backgroundFeatures)<-c('chrm','position','featureType','gene.check','geneID')
backgroundFeatures<-subset(unique(backgroundFeatures), gene.check=="gene_id"&featureType!="mRNA"&chrm!="rDNA"&chrm!="Unmapped_Scaffold_8_D1580_D1567") #any use is in looking for mRNA and gene DON'T overlap?
backgroundFeatures<-backgroundFeatures[-c(4)]
backgroundFeatures$ID=paste(backgroundFeatures$chrm,backgroundFeatures$position)

#there are some extra filtering steps in the old version of this file but I don't used them
G130.MAF<- read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130_MAF_cleaned_allFreq.csv",header=TRUE,sep=",")
G7.data<- read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G7_totalCount_cov5.txt", header = TRUE,sep="")
colnames(G7.data)[1]= "chrm"

##Add EXTRA INFO TO SBGE
#ADD chrm info (from backgroundFeatures) to SBGE.table IF OVERLAP with geneID
#DOES NOT update Is.A or Is.X however!
overlap.geneIDX=SBGE.table[[1]]%in%backgroundFeatures[["geneID"]] #13875
missing.chrm.IDX=is.na(SBGE.table[,"chrm"]) #the length is 13875 sum is 4030 ->now only 1847 are na
helpIDX=overlap.geneIDX&missing.chrm.IDX #sum 4030
help.gene.list=SBGE.table[helpIDX,"geneID"]
chrm.list=unique(subset(backgroundFeatures[,c("geneID","chrm")],backgroundFeatures[,"geneID"]%in%help.gene.list))#4030
SBGE.subset=SBGE.table[helpIDX,c("geneID","chrm")]
merged.chrms=merge(chrm.list,SBGE.subset,by="geneID",all.x=FALSE)
SBGE.table[helpIDX,"chrm"]=merged.chrms[[2]]# 


#this one is referenced
newtest=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130_newCMHchr2_full_noMAF.csv")#702 876
newtest=newtest[,-1]#remove "X" row
newtest$FDR5=newtest$p.adjust<0.05


#SNPs IDs from G130 (position on chromosome) were given to the variant effector predictor https://useast.ensembl.org/Tools/VEP
#this was done after triallelic sites were filtered out but before any subsequent filtering

##===STEP TWO===Obtain VEP files
#taken from G130 data. Split into parts due to file size
VEP_2LA=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\VEP\\VEP_SSAV_G130_2LA.txt",header=TRUE,sep="")
VEP_2LB=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\VEP\\VEP_SSAV_G130_2LB.txt",header=TRUE,sep="")
VEP_2RA=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\VEP\\VEP_SSAV_G130_2RA.txt",header=TRUE,sep="")
VEP_2RB=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\VEP\\VEP_SSAV_G130_2RB.txt",header=TRUE,sep="")

#subset necessary columns so they can be combined, 
VEP_2LA=VEP_2LA[,c("Location","Allele","Consequence","Gene","BIOTYPE","DISTANCE")]
VEP_2LB=VEP_2LB[,c("Location","Allele","Consequence","Gene","BIOTYPE","DISTANCE")]
VEP_2RA=VEP_2RA[,c("Location","Allele","Consequence","Gene","BIOTYPE","DISTANCE")]
VEP_2RB=VEP_2RB[,c("Location","Allele","Consequence","Gene","BIOTYPE","DISTANCE")]

#create an ID that can be compared across all files
create_ID<-function(location_info){
  change1=str_split_i(location_info,"-",1)
  return(gsub(":"," ",change1))}
VEP_2LA$ID=create_ID(VEP_2LA$Location)
VEP_2LB$ID=create_ID(VEP_2LB$Location)
VEP_2RA$ID=create_ID(VEP_2RA$Location)
VEP_2RB$ID=create_ID(VEP_2RB$Location)

##filter from G130
##G130 MAF has following filtering steps:
#A)No tri-alleic states
#B)Variants present in ALL populations
VEP_2LA=VEP_2LA[VEP_2LA$ID%in%G130.MAF$ID,]#
VEP_2LB=VEP_2LB[VEP_2LB$ID%in%G130.MAF$ID,]#
VEP_2RA=VEP_2RA[VEP_2RA$ID%in%G130.MAF$ID,]#
VEP_2RB=VEP_2RB[VEP_2RB$ID%in%G130.MAF$ID,]#
#combine files into one
VEP_chr2=rbind(VEP_2LA,VEP_2LB,VEP_2RA,VEP_2RB)

#write.csv(VEP_chr2,"C:\\Users\\user\\Documents\\Documents\\Genomics\\VEP_SSAV_G130_2ALL_subcol_noMAF_filt.txt")

####===STEP THREE====FILTERING VEP RESULTS========####
VEP_chr2=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\VEP_SSAV_G130_2ALL_subcol_noMAF_filt.txt",header=TRUE)
VEP_chr2=VEP_chr2[,c("ID","Allele","Consequence","BIOTYPE","Gene")]

#designate if a SNP is a candidate diverged SNP
VEP_chr2$is.Sig=FALSE
VEP_chr2[VEP_chr2$ID%in%newtest[newtest$FDR5,"ID"],"is.Sig"]=TRUE

#remove duplicate rows
VEP_chr2.filt=VEP_chr2[!duplicated(VEP_chr2),]
VEP_chr2.filt=VEP_chr2.filt[VEP_chr2.filt$ID%in%G130.MAF$ID,]#610 844

##remove transposable elements "FBti"
VEP_chr2.filt=VEP_chr2.filt[!grepl("FBti",VEP_chr2.filt$Gene),]#450 114-->5.5 mill
Multi.Cons.IDX=grepl(",",VEP_chr2.filt$Consequence)
VEP_chr2.filt[Multi.Cons.IDX,"Consequence"]=str_split_i(VEP_chr2.filt[Multi.Cons.IDX,"Consequence"],",",1)

####### TOTAL VEP HIERARCHY IS AS FOLLOWS
#transcript_ablation  > splice_acceptor_variant  > splice_donor_variant >  stop_gained  > frameshift_variant
#stop_lost start_lost > transcript_amplification > feature_elongation   > feature_truncation > inframe_insertion
#inframe_deletion     > missense_variant         >protein_altering_variant > splice_donor_5th_base_variant
#splice_region_variant > splice_donor_region_variant > splice_polypyrimidine_tract_variant > incomplete_terminal_codon_variant
#start_retained_variant > stop_retained_variant  > synonymous_variant > coding_sequence_variant > mature_miRNA_variant
#5_prime_UTR_variant  > 3_prime_UTR_variant  > non_coding_transcript_exon_variant > intron_variant
# NMD_transcript_variant > non_coding_transcript_variant > coding_transcript_variant > upstream_gene_variant
#downstream_gene_variant > TFBS_ablation > TFBS_amplification > TF_binding_site_variant > regulatory_region_ablation
#regulatory_region_amplification > regulatory_region_variant > intergenic_variant > sequence_variant


####----This section concatenates named groups and removes Consequences where the numbers are too low

#splice>stop>frameshift>missense>protein alt variant>start>syn>coding seq>5UTR>3UTR>intron>upstream>downstream>intergenic
#collapse some of the above categories
VEP_chr2.filt[grepl("splice",VEP_chr2.filt$Consequence),"Consequence"]="splice"
VEP_chr2.filt[grepl("start",VEP_chr2.filt$Consequence),"Consequence"]="start"
VEP_chr2.filt[grepl("stop",VEP_chr2.filt$Consequence),"Consequence"]="stop"


VEP_chr2.filt[grepl("protein_altering_variant",VEP_chr2.filt$Consequence),"Consequence"]="missense_variant" #"FBgn0031288" "FBgn0031288" "FBgn0031288" "FBgn0031585" "FBgn0085424" "FBgn0085472"
VEP_chr2.filt[grepl("coding_sequence_variant",VEP_chr2.filt$Consequence),"Consequence"]="missense_variant" #"FBgn0031288" "FBgn0031288" "FBgn0031288" "FBgn0031585" "FBgn0085424" "FBgn0085472"
VEP_chr2.filt[grepl("inframe",VEP_chr2.filt$Consequence),"Consequence"]="inframe_indel"
#Remove the confusing " non_coding_transcript_exon_variant A sequence variant that changes non-coding exon sequence in a non-coding transcript
#basically if something changes RNA or whatever
VEP_chr2.filt=VEP_chr2.filt[VEP_chr2.filt$Consequence!="non_coding_transcript_exon_variant",]
#write.csv(VEP_chr2.filt,"C:\\Users\\user\\Documents\\Documents\\Genomics\\VEP_filter_intermediate_1.txt",row.names=FALSE)


#######===============================================================================####
#######------------------------------Final Filtering and Graphing---------------######
#######=============================================================================######


#this is the latest VEP filtered intermediate file that does all the above filtering steps so I don't have to run it each time
VEP_chr2.filt=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\VEP_filter_intermediate_1.txt")
VEP_chr2.filt[grepl("protein_altering_variant",VEP_chr2.filt$Consequence),"Consequence"]="missense_variant" #"FBgn0031288" "FBgn0031288" "FBgn0031288" "FBgn0031585" "FBgn0085424" "FBgn0085472"
VEP_chr2.filt[grepl("coding_sequence_variant",VEP_chr2.filt$Consequence),"Consequence"]="missense_variant" #"FBgn0031288" "FBgn0031288" "FBgn0031288" "FBgn0031585" "FBgn0085424" "FBgn0085472"

#removes duplicates again
VEP_chr2.filt=VEP_chr2.filt[!duplicated(VEP_chr2.filt[,c("ID","Consequence","Gene","BIOTYPE")]),]#removes 359 262

#collapsing removes about14 SNPs
VEP_chr2.filt=VEP_chr2.filt[!duplicated(VEP_chr2.filt[,c("ID","Consequence","Gene")]),]#biotype doesn't change
#VEP_chr2.filt=VEP_chr2.filt[!duplicated(VEP_chr2.filt[,c("ID","Consequence")]),]#doesn't really matter if per gene
#there is massive overlap so put these in one category for now. May change later
#VEP_chr2.filt[grepl("downstream_gene_variant",VEP_chr2.filt$Consequence),"Consequence"]="updownstream_gene_variant"
#VEP_chr2.filt[grepl("upstream_gene_variant",VEP_chr2.filt$Consequence),"Consequence"]="updownstream_gene_variant"
VEP_chr2.filta=VEP_chr2.filt

VEP_chr2.filt=VEP_chr2.filt[!duplicated(VEP_chr2.filt[,c("ID","Consequence")]),]
#359 247->276 754 SNPs
ordered.Cons=c("splice","stop","frameshift_variant","infram_indel","missense_variant","start","synonymous_variant",
               "5_prime_UTR_variant","3_prime_UTR_variant","non_coding_transcript_exon_variant","intron_variant",
               "upstream_gene_variant","downstream_gene_variant","intergenic_variant")


#Remove duplicate Consequences per ID
VEP_chr2.filt=VEP_chr2.filt[!duplicated(VEP_chr2.filt[,c("ID","Consequence")]),]#182 883 (one for every ID)


#this is where the distinction candidate/non candidates made again
newtest=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130_newCMHchr2_full_noMAF.csv")#
newtest=newtest[,-1]#remove "X" row
newtest$FDR5=newtest$p.adjust<0.05
VEP_chr2.filt$is.Sig=VEP_chr2.filt$ID%in%newtest[newtest$FDR5,"ID"]
#This would be where the MAF comes from (but is not 100% MAF filtered due to differenced between using N, R or NR for MAF filtering)
G130.MAF=read.table("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130_MAF_cleaned_MAF1_Freq.csv",header=TRUE)
G130.MAF$Seg.numR=rep(0,length(G130.MAF$ID))#number of Red populations where a variant is segregating

#this makes sure that SNPs are segregating appropriately in all populations
for(col in seq(5,15,2)){#doing cols 4-15 by 1 would give seg in EVERY population
  #only look at red because can't be segregating in N without R (sine the Red sample type is Red/NonRed heterozygotes)
  seg.IDX=G130.MAF[,col]!=1
  G130.MAF[seg.IDX,"Seg.numR"]=G130.MAF[seg.IDX,"Seg.numR"]+1}


#VEP_chr2.filta=VEP_chr2.filt[!duplicated(VEP_chr2.filt[,c("ID","Consequence","Gene")]),]#biotype doesn't change
#turn upsteam and downsteam of a gene into one category. A SNP occuring upstream of one gene and downstream of another only gets counted once
VEP_chr2.filt[grepl("stream",VEP_chr2.filt$Consequence),"Consequence"]="updown_stream"
VEP_chr2.filt=VEP_chr2.filt[!duplicated(VEP_chr2.filt[,c("ID","Consequence")]),]

#this function bootstraps the production of a table (how many counts in each category) used to commute a fraction
boot_ftable<-function(data1,col1,col2){
  #the columns are is.Sig and Consequence
  data_table=table(data1[,c(col1,col2)])
  fraction.table=data.frame(fraction=c(0),category=c(0))
  names=rownames(data_table)
  for (i in 0:dim(data_table)[1]){
    #print(paste(i,test[i,2],sum(test[i,2])))
    fraction.table[i,"fraction"]=sum(data_table[i,2])/sum(data_table[i,])
    fraction.table[i,"category"]=names[i]
  }
  return(fraction.table)}


#filter out consequences with low SNP counts, 
VEP_chr2.filt_sub1=subset(VEP_chr2.filt[,c("is.Sig","Consequence","ID","BIOTYPE")],Consequence!="stop"&Consequence!="inframe_indel"&Consequence!="start"&Consequence!="frameshift_variant")
VEP_chr2.filt_sub=VEP_chr2.filt_sub1[VEP_chr2.filt_sub1$ID%in%G130.MAF[,"ID"],]#30,000 not in G130.MAF
VEP_chr2.filt_sub2=merge(VEP_chr2.filt_sub,G130.MAF[,c("MAF.NRavg","ID","MAF.Ravg")],by.x="ID",by.y="ID",all.x=TRUE)
VEP_chr2.filt_sub2=VEP_chr2.filt_sub2[!duplicated(VEP_chr2.filt_sub2[,c("ID","Consequence","BIOTYPE")]),]

##########FILTERING MAF AND FOR SEGREGATING IN HALF THE POPULATIONS
VEP_chr2.filt_sub3=VEP_chr2.filt_sub2[VEP_chr2.filt_sub2$ID%in%G130.MAF[G130.MAF$Seg.numR>3&G130.MAF$MAF.NRavg>0.05,"ID"],]#30,000 not in G130.MAF


#####GRAPH GENIC VS INTER_GENIC 

#Checking the break down of all categories across biotypes
#table(VEP_chr2.filt_sub3[,c("Consequence","BIOTYPE")])[,1]
VEP_chr2.filt_sub3$genic="genic"
VEP_chr2.filt_sub3[VEP_chr2.filt_sub3$Consequence=="intergenic_variant","genic"]="intergenic"
VEP_chr2.filt_sub3[VEP_chr2.filt_sub3$Consequence=="non_coding_transcript_exon_variant","genic"]="non-coding exon variant"

#get confidence intervals around the mean for genic vs non-genic %candidates
#get the mean: note this applies boot fraction per ro
VEP_inter_vs_genic=replicate(1000,{
  #Sample with replacement (for bootstrap from original dataset).
  resampID <- sample(VEP_chr2.filt_sub3[,"ID"], replace=TRUE)
  resamp=VEP_chr2.filt_sub3[VEP_chr2.filt_sub3$ID%in%resampID,]
  
  #Generate the fraction from the resampled dataset, and extract y coordinates to generate fraction
  boot_ftable(resamp,"genic","is.Sig")$fraction})
#order is "genic (not protein coding include tRNA""intergentic"  "non-coding transcription variant"


#take the CI around the mean fraction of candidates in each category
VEP_inter_vs_genic_CI=data.frame(category=c("Genic","Intergenic","\nNonCoding Transcription Exon"),CI_up=apply(VEP_inter_vs_genic,1,quantile,probs=0.975),
                            CI_low=apply(VEP_inter_vs_genic,1,quantile,probs=0.025))
#the mean (using bootftable)
VEP_inter_vs_genic_CI$mean=boot_ftable(VEP_chr2.filt_sub3,"genic","is.Sig")$fraction

#only take the first 2 categories, genic and intergenic. NonCoding Transcription Exon excluded e
x=ggplot(VEP_inter_vs_genic_CI[1:2,],
         aes(x=category,fill=category,y=mean,label=category))+
  #geom_point(stat="identity",position=position_dodge(0.9),size=6,pch=21)+
  geom_bar(stat="identity",position=position_dodge(0.9),lwd=1,colour="black")+
  scale_y_continuous("fraction of Candidates",limits = c(0, max(VEP_inter_vs_genic_CI$CI_up)*1.08)) +
  scale_colour_manual(values=c("black","black","black","black","black","black","black","black"))+
  scale_fill_manual(values=c("lightgoldenrod","burlywood4"))+
  xlab("Variant Effects") +  
  ggtitle(paste("chr2","NRMAF(seg4)"))+theme_bw() + 
  theme(text = element_text(size = 16),legend.position="bottom")+
  theme(plot.title = element_text(size = 22)) + 
  geom_errorbar(aes(ymin=as.numeric(CI_low), ymax=as.numeric(CI_up)),lwd=.5,
                width=.4,position=position_dodge(.9),show.legend = FALSE)

print(x)


#Only take SNPs that affect protein coding regions (including updown stream of protein coding regions)
MAF.title="(int1)NRMAF(seg4): protein coding"
#only take SNPs that occur in (or around in the case of updownstream) protein coding genes
VEP_chr2.filt_sub=VEP_chr2.filt_sub3[VEP_chr2.filt_sub3$BIOTYPE=="protein_coding",]#VEP_chr2.filt_sub3$BIOTYPE=="protein_coding"

#Take the bootstraped fraction of candidates in each functional category
VEP_boot_fraction=replicate(1000,{
  #Sample with replacement (for bootstrap from original dataset).
  resampID <- sample(VEP_chr2.filt_sub[,"ID"], replace=TRUE)
  resamp=VEP_chr2.filt_sub[VEP_chr2.filt_sub$ID%in%resampID,]
  
  #Generate the fraction from the resampled dataset, and extract y coordinates to generate fraction
  boot_ftable(resamp,"Consequence","is.Sig")$fraction})

names=row.names(table(VEP_chr2.filt_sub[,c("Consequence","is.Sig")]))
#renaming this for brevity for labelling. ORDER IS FIXED
names=c("3'UTR","5'UTR","intron","missense","splice",
        #"non_coding_exon","\nsynonymous","\nupdown_stream","intergenic")
         "synonymous","updown_stream")

#get the confidence interval for the mean fraction of candidates in each category
VEP_chr2.filt_CI=data.frame(category=names,CI_up=apply(VEP_boot_fraction,1,quantile,probs=0.975),
                            CI_low=apply(VEP_boot_fraction,1,quantile,probs=0.025))

VEP_chr2.filt_CI$mean=boot_ftable(VEP_chr2.filt_sub,"Consequence","is.Sig")$fraction

#fix the order of factors
VEP_chr2.filt_CI_p<-VEP_chr2.filt_CI  %>%
  mutate(category=fct_relevel(category,
   #"\nupdown_stream","5'UTR","splice","intron","missense","\nsynonymous","3'UTR"
   "missense","synonymous","splice","intron","5'UTR","3'UTR","updown_stream",                          
  #"splice","missense","\nsynonymous","5'UTR","3'UTR","intron","\nupdown_stream","intergenic","\nnon_coding_exon"
  ))
VEP_chr2.filt_CI=VEP_chr2.filt_CI_p
#generate a line which is the UNWEIGHTED average of fraction of Candidates within each category (GEN fraction is 0.1091133)
unweight_AvgFractionCandidate=mean(VEP_chr2.filt_CI$mean)+0.0005#0.1066425-->differs by 0.0025or 0.25%

VEP_plot=ggplot(VEP_chr2.filt_CI,aes(x=category,fill=category,y=mean,label=category))+
geom_bar(stat="identity",position=position_dodge(0.9),lwd=1,colour="black")+
  scale_y_continuous("fraction of Candidates",limits = c(0, max(VEP_chr2.filt_CI$CI_up)*1.02)) +
  scale_colour_manual(values=c("black","black","black","black","black","black","black","black"))+
  scale_fill_manual(values=c("darkseagreen","palegreen","skyblue1","lightblue3","plum2","mediumpurple3","lightyellow3"))+
  xlab("Variant Effects") +theme_classic() + 
  theme(text = element_text(size = 16),legend.position="right")+ 
  geom_errorbar(aes(ymin=as.numeric(CI_low), ymax=as.numeric(CI_up)),lwd=.5,
                width=.4,position=position_dodge(.9),show.legend = FALSE)+
  geom_abline(slope=0,lwd=.5,colour="black",intercept=unweight_AvgFractionCandidate)

print(VEP_plot)


