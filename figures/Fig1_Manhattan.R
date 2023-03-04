library(data.table)
library(ggplot2)

#2) manhattan plot part
X <- NULL
clumpDat <- NULL

for(chr in 1:23){
    print(chr)
    X <- rbind(X,fread(paste0("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/Menopause_LinearEvaluated/Menopause_df2_atMax_45_52_",chr,".res")))
}


clumpDat <- NULL
for(chr in 1:23){
    dat <- fread(paste0("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/Menopause_atMax_",chr,"_clumping_p1_df2.clumped"))
    splitsNeeded <- max(dat$TOTAL)+1

    datSplit <- str_split_fixed(dat$SP2, pattern=",",n=splitsNeeded)
    datSplit[datSplit == "" | datSplit == "NONE"] <- NA
    datSplit <- as.data.table(datSplit)
    datSplit[,indexSNP:=dat$SNP]
    
    datLong <- melt(datSplit,id.vars="indexSNP")
    datLong[,value:=gsub("(1)","",value,fixed=T)]
    datLong[,variable:=NULL]

    indexSNPData <- copy(datLong)
    indexSNPData[,value:=indexSNP]
    indexSNPData <- unique(indexSNPData)

    #Add the index SNPs to the main file
    datLong <- rbind(datLong,indexSNPData)
    datLong <- datLong[!is.na(value),]

    datLong <- merge(datLong,dat[,list(SNP,CHR,BP,P)],by.x="indexSNP",by.y="SNP",all.x=T)

    clumpDat <- rbind(clumpDat,datLong)
}


#Read the novel/previously found results + clumps
replDat <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/ReplicatedNovelAssociationsatMax_45_52.csv")
prevDat <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/PreviousAssociationsatMax_45_52.csv")

replDat <- replDat[,list(SNP)]
prevDat <- prevDat[,list(Variant)]
names(prevDat) <- "SNP"
replDat[,specType:="Novel discovery"]
prevDat[,specType:="Previous discovery"]

mergDat <- rbind(replDat,prevDat)
mergDat[specType == "Previous discovery",specType:=NA]

clumpDat2 <- merge(clumpDat,mergDat,by.x="indexSNP",by.y="SNP",all.x=T)
clumpDat3 <- clumpDat2[!is.na(specType),list(value,specType)]

X <- merge(X,clumpDat3,by.x="SNP",by.y="value",all.x=T)

X_all_sub <- X[P<0.005,]
X_all_sub <- X_all_sub[order(CHR,BP),]

fwrite(X_all_sub,file="/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/MarieDat.txt")
X_all_sub <- fread("MarieDat.txt")
X_all_sub[,max_bp:=max(BP),by="CHR"]
toAdd <- unique(X_all_sub[,list(CHR,max_bp)])
toAdd[,toAdd:=lag(max_bp)]
toAdd[is.na(toAdd),toAdd:=0]
toAdd[,cumToAdd:=cumsum(as.numeric(toAdd))]

X_all_sub <- merge(X_all_sub,toAdd[,list(CHR,cumToAdd)],by="CHR",all.x=T)
X_all_sub[,bpAdj:=BP+cumToAdd]

axis_set <- X_all_sub %>% 
  group_by(CHR) %>% 
  summarize(center = mean(bpAdj))

X_all_sub[,CHR2:=" "]
X_all_sub[(CHR)%%2 == 0,CHR2:="  "]

#X_all_sub[!is.na(specType),CHR2:=specType]
