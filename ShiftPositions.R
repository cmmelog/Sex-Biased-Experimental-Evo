
## This file is for shifting SNP positions.
## It assumes their are 3 segments we care about
## Chrm 2L  (from first variable SNP to last)
## Chrm 2R - segment 1
## Chrm 2R - segment 2
### The 2 "segments" on 2R are separated by the "DsRed region", which we should exclude
#cov40.500.SNPS.chr2=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130_newCMHchr2_full.csv")
#cov40.500.SNPS.chr2=read.csv("/plas1/c.melo.gavin/SSAV/SSAV_G130_newCMHchr2_full.csv")
#test=subset(cov40.500.SNPS.chr2,chrm=="2R"&pos>1180500)

## The position boundaries for each of these segments.
## You can just use the ==first and last usable SNPs== within each region as the boundaries
## They are currently set with dummy values for testing
#this.start2L <- 11  ## change to true value
#this.end2L <- 100 ## change to true value
#this.start2R.seg1 <- 1001 ## change to true value
#this.end2R.seg1 <- 1100 ## change to true value
#this.start2R.seg2 <- 2001 ## change to true value
#this.end2R.seg2 <- 2100 ## change to true value

this.start2L <- 99  ## 2L 99
this.end2L <- 23513704 ## 2L 23513704
this.start2R.seg1 <- 1001 ## 2R 67
this.end2R.seg1 <- 10697465 ## 2R 10697465 last before "DSRED cut-off 10,700,000"->then goes to 2R 10700002
this.start2R.seg2 <- 1189744## 2R 1189744  (jump from 1184223 to 1189744 to 1190252)
this.end2R.seg2 <- 25286461## 2R 25286461

## total combined tract length is: 
totalTractLength = (this.end2L - this.start2L + 1) + (this.end2R.seg1 - this.start2R.seg1 + 1) + (this.end2R.seg2 - this.start2R.seg2 + 1)
#the shift  must be <= totalLength


ConvertToNewPosition<-function(chrm, pos, shift){
  ## shift must be between 0 and total length of all chromosome segments
  start2L <- this.start2L
  end2L <- this.end2L
  start2R.seg1 <- this.start2R.seg1
  end2R.seg1 <- this.end2R.seg1
  start2R.seg2 <- this.start2R.seg2
  end2R.seg2 <- this.end2R.seg2
  
  ## map to continuous space
  last2LPosInContSpace = end2L - start2L + 1
  first2R.seg1PosInContSpace = last2LPosInContSpace + 1
  last2R.seg1PosInContSpace = first2R.seg1PosInContSpace + (end2R.seg1 - start2R.seg1)
  first2R.seg2PosInContSpace = last2R.seg1PosInContSpace + 1
  last2R.seg2PosInContSpace = first2R.seg2PosInContSpace + (end2R.seg2 - start2R.seg2)
  
  #y is new position
  y=0
  if(chrm == "2L") y = pos - start2L + 1 + shift
  if(chrm == "2R" & pos <= end2R.seg1) y = first2R.seg1PosInContSpace + (pos - start2R.seg1) + shift
  if(chrm == "2R" & pos > end2R.seg1) y = first2R.seg2PosInContSpace + (pos - start2R.seg2) + shift
  
  ## now we effectively make the y space circular
  ## this fails if y is > last2R.seg2PosInContSpace, which it should never be if "shift" is not larger
  ## than total length of track
  if(y > last2R.seg2PosInContSpace) y = y - last2R.seg2PosInContSpace 
  
  # map back to chromosomal structure
  new.chrm = NA
  new.pos = NA
  if(y <= last2LPosInContSpace) {
    new.chrm = "2L"
    new.pos = start2L + y - 1
  }
  if (y <= last2LPosInContSpace) {
    new.chrm = "2L"
    new.pos = start2L + y - 1
  } else if (y <= last2R.seg1PosInContSpace) {
    new.chrm = "2R"
    new.pos = start2R.seg1 + (y - first2R.seg1PosInContSpace) 
  } else if (y <= last2R.seg2PosInContSpace) {
    new.chrm = "2R"
    new.pos = start2R.seg2 + (y - first2R.seg2PosInContSpace)
  }

  return(c(new.chrm, new.pos))
}


## Takes data frame with col 1 for chromosome and col 2 for position
## Returns same type of data frame with shifted positions
## It is ***important*** that the data frame you pass to this function
## only contains SNPs within the designated segments 

MakeDataFrameOfShiftedPositions<-function(df, shift){
  df.shifted = as.data.frame(t(sapply(1:nrow(df), function(x) ConvertToNewPosition(df[x, 1], df[x, 2], shift))))
  #the capital is called in the other function so can't change it
  names(df.shifted) = c("Chrm", "Pos")
  df.shifted$Pos = as.numeric(df.shifted$Pos)
  return(df.shifted)
}


{
## testing it
#  test.shift= 100
#  fake.sites.test = data.frame(chrm = c(rep("2L", 3), rep("2R", 6)), pos = c(11, 22, 100, 1001, 1090, 1100, 2001, 2080, 2100))
#  fake.sites.test.shifted = MakeDataFrameOfShiftedPositions(fake.sites.test, test.shift)
#  cbind(fake.sites.test, fake.sites.test.shifted) 

## Each permutation of the data should shift the data by drawing a "shift" value using
#test.shifts=sample(1:totalTractLength, 100)
#fake.sites.test = data.frame(chrm = c(rep("2L", 3), rep("2R", 6)), pos = c(11, 22, 100, 1001, 1090, 1100, 2001, 2080, 2100))
#perm.test=fake.sites.test
#for(shift in test.shifts){
#  test.shift= shift
#  fake.sites.test.shifted = MakeDataFrameOfShiftedPositions(fake.sites.test, test.shift)
#  perm.test=cbind(perm.test, fake.sites.test.shifted) 
#}
}

test.shifts=sample(1:totalTractLength, 4)
sites.test = data.frame(chrm = cov40.500.SNPS.chr2$chrm, pos = cov40.500.SNPS.chr2$pos)
sites.test=sites.test[,c("chrm","pos")]
perm.test=sites.test
for(shift in test.shifts){
  test.shift= shift
  sites.test.shifted = MakeDataFrameOfShiftedPositions(sites.test, test.shift)
  perm.test=cbind(perm.test, sites.test.shifted) 
}

#write.csv(perm.test,"/plas1/c.melo.gavin/SSAV/SSAV_G130_permshiftCHR2_part1.csv",row.names = FALSE)

test.perm.shift=read.csv("C:\\Users\\user\\Documents\\Documents\\Genomics\\SSAV_G130_permshiftCHR2_part3.csv")

perm.shift4BED=test.perm.shift[,c(1,2)]
perm.shift4BED$end=perm.shift4BED$pos+1
perm.shift4BED$ShiftPerm=0 #setting permutation 0->this is the original set
for(col in  seq.int(1,length(test.perm.shift),2)){
  print(col)
  if(col!=1){
    #SNPs=length(test.perm.shift$chrm)
    temp.perm=test.perm.shift[,c(col,col+1),]
    colnames(temp.perm)[1]="chrm"
    colnames(temp.perm)[2]="pos"
    temp.perm$end=temp.perm$pos+1
    temp.perm$ShiftPerm=col-2
    perm.shift4BED=rbind(perm.shift4BED,temp.perm)
  }
  
}
