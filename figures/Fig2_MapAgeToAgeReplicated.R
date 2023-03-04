library(data.table)

if(F){
dat <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/replication/PerAgeSlopeMaxTResults")
bimEst <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/replication/PerAgeSlopeMaxTResults.bim")

dat <- merge(dat,bimEst[,list(V2,V5)],by.x="SNP",by.y="V2",all.x=T)
setnames(dat,"V5","effAll")

t_grid <- seq(41,55,2)
replicatedResults <- NULL
for(t in t_grid){
        #Take SNPs to test
        snpUse <- NULL
        for(chr in 1:23){
            fileName <- paste0("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/Menopause_p1_df2_",t,"_",chr,"_COJO.jma.cojo")
            if(file.exists(fileName)){
                tmp <- fread(fileName)[,list(Chr,SNP,refA,b,p)]
                snpUse <- rbind(snpUse,tmp)
            }
        }
        datUse <- dat[SNP %in% snpUse$SNP,]
        datUse[,tUse:=t]
        datUse[,se1:= abs(b1/qnorm(0.5 * p1))]
        datUse[,se0:= abs(b0/qnorm(0.5 * p0))]

        datUse[,bUse:= b0 + b1 * (tUse - 49)]
        datUse[,seUse:=sqrt(se0**2 + (tUse-49)**2 * se1**2 + 2*(tUse-49)*cov_b0_b1)]
        datUse[,chisq:=(bUse/seUse)**2]
        datUse[,pUse:=1-pchisq(chisq,df=2)]
        datUse <- datUse[,list(SNP,tUse,bUse,pUse,effAll)]
        datUse <- merge(datUse,snpUse,by="SNP")
        datUse[,bUse2:=bUse]
        datUse[effAll != refA,bUse2:=bUse * (-1)]
        datUse2 <- datUse[bUse2 * b > 0 & pUse < 0.05,]
        #fwrite(datUse,paste0("/data/MenopauseAnalysis/MenopausePerAge/PerAgeReplication",t))
        replicatedResults <- rbind(replicatedResults,datUse2)

}



res_list <- list()
t_grid <- seq(41,53,2)  #age2 =t+2

for(i in 1:length(t_grid)){
    useAge1 <- t_grid[i]
    useAge2 <- useAge1 + 2
    replicated1 <- replicatedResults[tUse == useAge1,]
    replicated2 <- replicatedResults[tUse == useAge2,]

    allResults1 <- NULL
    for(chr in 1:23){
        fileName <- paste0("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/Menopause_p1_df2_",useAge1,"_",chr,"_COJO.jma.cojo")
        if(file.exists(fileName)){
            dat <- fread(fileName)[,list(Chr,SNP)]
            dat <- dat[SNP %in% replicated1$SNP,]   #Use only replicated loci
            allResults1 <- rbind(allResults1,dat)
        }
    }

    allResults1[,isAge1:=T]

    allResults2 <- NULL
    for(chr in 1:23){
        fileName <- paste0("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/Menopause_p1_df2_",useAge2,"_",chr,"_COJO.jma.cojo")
        if(file.exists(fileName)){
            dat <- fread(fileName)[,list(Chr,SNP)]
            dat <- dat[SNP %in% replicated2$SNP,]   #Use only replicated loci
            allResults2 <- rbind(allResults2,dat)
        }
    }

    #Read clumps for age2
    clumpDat2 <- NULL
    for(chr in 1:23){
        fileName <- paste0("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/Menopause_",useAge2,"_",chr,"_clumping_p1_df2.clumped")
        if(file.exists(fileName)){
            dat <- fread(fileName)
            splitsNeeded <- max(dat$TOTAL)+1

            datSplit <- stringr::str_split_fixed(dat$SP2, pattern=",",n=splitsNeeded)
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

            datLong <- merge(datLong,dat[,list(SNP,CHR,P)],by.x="indexSNP",by.y="SNP",all.x=T)

            clumpDat2 <- rbind(clumpDat2,datLong)
        }
    }
    #Keep only the significant ones (after COJO)
    clumpDat2 <- clumpDat2[indexSNP %in% allResults2$SNP,]

    #Check if the index SNP in age1 is the clump in age2
    clumpDat2 <- merge(clumpDat2,allResults1,all.x=T,by.x=c("value","CHR"),by.y=c("SNP","Chr"))
    allResults2_usedAge1 <- clumpDat2[isAge1 ==T,]

    allResults2 <- merge(allResults2,allResults2_usedAge1[,list(indexSNP,value)],by.x="SNP",by.y="indexSNP",all.x=T)
    res_list[[i]] <- allResults2
}
#This section is for smart data set merging


#Let's map one by one
baseResult <- NULL
    for(chr in 1:23){
        fileName <- paste0("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/Menopause_p1_df2_",41,"_",chr,"_COJO.jma.cojo")
        if(file.exists(fileName)){
            dat <- fread(fileName)[,list(SNP)]
            dat <- dat[SNP %in% replicatedResults[tUse == 41,]$SNP,]   #Use only replicated loci
            baseResult <- rbind(baseResult,dat)
        }
}
names(baseResult) <- c("SNP_41")
#names(baseResult) <- c("SNP_43","Chr","SNP_41")

#But for now let's do it manually
i <- 1
datUse <- res_list[[i]][,list(SNP,value)]
nameVec <- c(paste0("SNP_",t_grid[i]+2), paste0("SNP_",t_grid[i]))
names(datUse) <- nameVec
datUse[is.na(SNP_41),SNP_41:=paste0("missing41_",1:.N)]
baseResult <- merge(baseResult,datUse,all.x=T,all.y=T,by=c(paste0("SNP_",t_grid[i])))
baseResult <- unique(baseResult)

i <- 2
datUse <- res_list[[i]][,list(SNP,value)]
nameVec <- c(paste0("SNP_",t_grid[i]+2), paste0("SNP_",t_grid[i]))
names(datUse) <- nameVec
datUse[is.na(SNP_43),SNP_43:=paste0("missing43_",1:.N)]
baseResult <- merge(baseResult,datUse,all.x=T,all.y=T,by=c(paste0("SNP_",t_grid[i])))

i <- 3
datUse <- res_list[[i]][,list(SNP,value)]
nameVec <- c(paste0("SNP_",t_grid[i]+2), paste0("SNP_",t_grid[i]))
names(datUse) <- nameVec
datUse[is.na(SNP_45),SNP_45:=paste0("missing45_",1:.N)]
baseResult <- merge(baseResult,datUse,all.x=T,all.y=T,by=c(paste0("SNP_",t_grid[i])))

i <- 4
datUse <- res_list[[i]][,list(SNP,value)]
nameVec <- c(paste0("SNP_",t_grid[i]+2), paste0("SNP_",t_grid[i]))
names(datUse) <- nameVec
datUse[is.na(SNP_47),SNP_47:=paste0("missing_47",1:.N)]
baseResult <- merge(baseResult,datUse,all.x=T,all.y=T,by=c(paste0("SNP_",t_grid[i])))

i <- 5
datUse <- res_list[[i]][,list(SNP,value)]
nameVec <- c(paste0("SNP_",t_grid[i]+2), paste0("SNP_",t_grid[i]))
names(datUse) <- nameVec
datUse[is.na(SNP_49),SNP_49:=paste0("missing_49",1:.N)]
baseResult <- merge(baseResult,datUse,all.x=T,all.y=T,by=c(paste0("SNP_",t_grid[i])))

i <- 6
datUse <- res_list[[i]][,list(SNP,value)]
nameVec <- c(paste0("SNP_",t_grid[i]+2), paste0("SNP_",t_grid[i]))
names(datUse) <- nameVec
datUse[is.na(SNP_51),SNP_51:=paste0("missing_51",1:.N)]
baseResult <- merge(baseResult,datUse,all.x=T,all.y=T,by=c(paste0("SNP_",t_grid[i])))

i <- 7
datUse <- res_list[[i]][,list(SNP,value)]
nameVec <- c(paste0("SNP_",t_grid[i]+2), paste0("SNP_",t_grid[i]))
names(datUse) <- nameVec
datUse[is.na(SNP_53),SNP_53:=paste0("missing_53",1:.N)]
baseResult <- merge(baseResult,datUse,all.x=T,all.y=T,by=c(paste0("SNP_",t_grid[i])))

#code the results that are missing as NAs
baseResult[grepl("missing",SNP_41),SNP_41:=NA]
baseResult[grepl("missing",SNP_43),SNP_43:=NA]
baseResult[grepl("missing",SNP_45),SNP_45:=NA]
baseResult[grepl("missing",SNP_47),SNP_47:=NA]
baseResult[grepl("missing",SNP_49),SNP_49:=NA]
baseResult[grepl("missing",SNP_51),SNP_51:=NA]
baseResult[grepl("missing",SNP_53),SNP_53:=NA]

#
bRes <- baseResult[,list(SNP_41,SNP_43,SNP_45,SNP_47,SNP_49,SNP_51,SNP_53,SNP_55)]

datCombine <- NULL
#How many new significant at 41
s41 <- length(unique(bRes$SNP_41)) -1 #-1 for NA
datCombine <- rbind(datCombine,c(s41,41,41))

#remaining significant at 43
s43_41 <- length(unique(bRes[!is.na(SNP_43) & !is.na(SNP_41),]$SNP_43))
datCombine <- rbind(datCombine,c(s43_41,41,43))
#new significanr at 43
s43 <- length(unique(bRes[!is.na(SNP_43) & is.na(SNP_41),]$SNP_43))
datCombine <- rbind(datCombine,c(s43,43,43))

length(unique(bRes$SNP_41)) -1

#remain significant at 45
bRes_use <- bRes[!is.na(SNP_45),list(SNP_41,SNP_43,SNP_45)]
bRes_use[,naCount:=apply(bRes_use,1,function(x) sum(is.na(x)) )]
bRes_use <- bRes_use[order(naCount),]
bRes_use <- bRes_use[!duplicated(SNP_45),]

s45_43_41 <- length(unique(bRes_use[ !is.na(SNP_43) & !is.na(SNP_41),]$SNP_45))
s45_43 <- length(unique(bRes_use[!is.na(SNP_43) & is.na(SNP_41),]$SNP_45))
s45 <- length(unique(bRes_use[is.na(SNP_43) & is.na(SNP_41),]$SNP_45))

datCombine <- rbind(datCombine,c(s45_43_41,41,45))
datCombine <- rbind(datCombine,c(s45_43,43,45))
datCombine <- rbind(datCombine,c(s45,45,45))


#remain significant at 47
bRes_use <- bRes[!is.na(SNP_47),list(SNP_41,SNP_43,SNP_45,SNP_47)]
bRes_use[,naCount:=apply(bRes_use,1,function(x) sum(is.na(x)) )]
bRes_use <- bRes_use[order(naCount),]
bRes_use <- bRes_use[!duplicated(SNP_47),]

s47_45_43_41 <- length(unique(bRes_use[!is.na(SNP_47) & !is.na(SNP_45) & !is.na(SNP_43) & !is.na(SNP_41),]$SNP_47))
s47_45_43 <- length(unique(bRes_use[!is.na(SNP_47) & !is.na(SNP_45) & !is.na(SNP_43) & is.na(SNP_41),]$SNP_47))
s47_45 <- length(unique(bRes_use[!is.na(SNP_47) & !is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_47))
s47 <- length(unique(bRes_use[!is.na(SNP_47) & is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_47))
datCombine <- rbind(datCombine,c(s47_45_43_41,41,47))
datCombine <- rbind(datCombine,c(s47_45_43,43,47))
datCombine <- rbind(datCombine,c(s47_45,45,47))
datCombine <- rbind(datCombine,c(s47,47,47))

#at 49
bRes_use <- bRes[!is.na(SNP_49),list(SNP_41,SNP_43,SNP_45,SNP_47,SNP_49)]
bRes_use[,naCount:=apply(bRes_use,1,function(x) sum(is.na(x)) )]
bRes_use <- bRes_use[order(naCount),]
bRes_use <- bRes_use[!duplicated(SNP_49),]

s49_47_45_43_41 <- length(unique(bRes_use[!is.na(SNP_49) & !is.na(SNP_47) & !is.na(SNP_45) & !is.na(SNP_43) & !is.na(SNP_41),]$SNP_49))
s49_47_45_43 <- length(unique(bRes_use[!is.na(SNP_49) & !is.na(SNP_47) & !is.na(SNP_45) & !is.na(SNP_43) & is.na(SNP_41),]$SNP_49))
s49_47_45 <- length(unique(bRes_use[!is.na(SNP_49) & !is.na(SNP_47) & !is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_49))
s49_47 <- length(unique(bRes_use[!is.na(SNP_49) & !is.na(SNP_47) & is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_49))
s49 <- length(unique(bRes_use[!is.na(SNP_49) & is.na(SNP_47) & is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_49))
datCombine <- rbind(datCombine,c(s49_47_45_43_41,41,49))
datCombine <- rbind(datCombine,c(s49_47_45_43,43,49))
datCombine <- rbind(datCombine,c(s49_47_45,45,49))
datCombine <- rbind(datCombine,c(s49_47,47,49))
datCombine <- rbind(datCombine,c(s49,49,49))

#at 51
bRes_use <- bRes[!is.na(SNP_51),list(SNP_41,SNP_43,SNP_45,SNP_47,SNP_49,SNP_51)]
bRes_use[,naCount:=apply(bRes_use,1,function(x) sum(is.na(x)) )]
bRes_use <- bRes_use[order(naCount),]
bRes_use <- bRes_use[!duplicated(SNP_51),]

s51_49_47_45_43_41 <- length(unique(bRes_use[!is.na(SNP_51) & !is.na(SNP_49) & !is.na(SNP_47) & !is.na(SNP_45) & !is.na(SNP_43) & !is.na(SNP_41),]$SNP_51))
s51_49_47_45_43 <- length(unique(bRes_use[!is.na(SNP_51) & !is.na(SNP_49) & !is.na(SNP_47) & !is.na(SNP_45) & !is.na(SNP_43) & is.na(SNP_41),]$SNP_51))
s51_49_47_45 <- length(unique(bRes_use[!is.na(SNP_51) & !is.na(SNP_49) & !is.na(SNP_47) & !is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_51))
s51_49_47 <- length(unique(bRes_use[!is.na(SNP_51) & !is.na(SNP_49) & !is.na(SNP_47) & is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_51))
s51_49 <- length(unique(bRes_use[!is.na(SNP_51) & !is.na(SNP_49) & is.na(SNP_47) & is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_51))
s51 <- length(unique(bRes_use[!is.na(SNP_51) & is.na(SNP_49) & is.na(SNP_47) & is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_51))
datCombine <- rbind(datCombine,c(s51_49_47_45_43_41,41,51))
datCombine <- rbind(datCombine,c(s51_49_47_45_43,43,51))
datCombine <- rbind(datCombine,c(s51_49_47_45,45,51))
datCombine <- rbind(datCombine,c(s51_49_47,47,51))
datCombine <- rbind(datCombine,c(s51_49,49,51))
datCombine <- rbind(datCombine,c(s51,51,51))

#At 53
bRes_use <- bRes[!is.na(SNP_53),list(SNP_41,SNP_43,SNP_45,SNP_47,SNP_49,SNP_51,SNP_53)]
bRes_use[,naCount:=apply(bRes_use,1,function(x) sum(is.na(x)) )]
bRes_use <- bRes_use[order(naCount),]
bRes_use <- bRes_use[!duplicated(SNP_53),]

s53_51_49_47_45_43_41 <- length(unique(bRes_use[!is.na(SNP_53) & !is.na(SNP_51) & !is.na(SNP_49) & !is.na(SNP_47) & !is.na(SNP_45) & !is.na(SNP_43) & !is.na(SNP_41),]$SNP_53))
s53_51_49_47_45_43 <- length(unique(bRes_use[!is.na(SNP_53) & !is.na(SNP_51) & !is.na(SNP_49) & !is.na(SNP_47) & !is.na(SNP_45) & !is.na(SNP_43) & is.na(SNP_41),]$SNP_53))
s53_51_49_47_45 <- length(unique(bRes_use[!is.na(SNP_53) & !is.na(SNP_51) & !is.na(SNP_49) & !is.na(SNP_47) & !is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_53))
s53_51_49_47 <- length(unique(bRes_use[!is.na(SNP_53) & !is.na(SNP_51) & !is.na(SNP_49) & !is.na(SNP_47) & is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_53))
s53_51_49 <- length(unique(bRes_use[!is.na(SNP_53) & !is.na(SNP_51) & !is.na(SNP_49) & is.na(SNP_47) & is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_53))
s53_51 <- length(unique(bRes_use[!is.na(SNP_53) & !is.na(SNP_51) & is.na(SNP_49) & is.na(SNP_47) & is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_53))
s53 <- length(unique(bRes_use[!is.na(SNP_53) & is.na(SNP_51) & is.na(SNP_49) & is.na(SNP_47) & is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_53))
datCombine <- rbind(datCombine,c(s53_51_49_47_45_43_41,41,53))
datCombine <- rbind(datCombine,c(s53_51_49_47_45_43,43,53))
datCombine <- rbind(datCombine,c(s53_51_49_47_45,45,53))
datCombine <- rbind(datCombine,c(s53_51_49_47,47,53))
datCombine <- rbind(datCombine,c(s53_51_49,49,53))
datCombine <- rbind(datCombine,c(s53_51,51,53))
datCombine <- rbind(datCombine,c(s53,53,53))

#at 55
bRes_use <- bRes[!is.na(SNP_55),list(SNP_41,SNP_43,SNP_45,SNP_47,SNP_49,SNP_51,SNP_53,SNP_55)]
bRes_use[,naCount:=apply(bRes_use,1,function(x) sum(is.na(x)) )]
bRes_use <- bRes_use[order(naCount),]
bRes_use <- bRes_use[!duplicated(SNP_55),]

s55_53_51_49_47_45_43_41 <- length(unique(bRes_use[!is.na(SNP_55) & !is.na(SNP_53) & !is.na(SNP_51) & !is.na(SNP_49) & !is.na(SNP_47) & !is.na(SNP_45) & !is.na(SNP_43) & !is.na(SNP_41),]$SNP_55))
s55_53_51_49_47_45_43 <- length(unique(bRes_use[!is.na(SNP_55) & !is.na(SNP_53) & !is.na(SNP_51) & !is.na(SNP_49) & !is.na(SNP_47) & !is.na(SNP_45) & !is.na(SNP_43) & is.na(SNP_41),]$SNP_55))
s55_53_51_49_47_45 <- length(unique(bRes_use[!is.na(SNP_55) & !is.na(SNP_53) & !is.na(SNP_51) & !is.na(SNP_49) & !is.na(SNP_47) & !is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_55))
s55_53_51_49_47 <- length(unique(bRes_use[!is.na(SNP_55) & !is.na(SNP_53) & !is.na(SNP_51) & !is.na(SNP_49) & !is.na(SNP_47) & is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_55))
s55_53_51_49 <- length(unique(bRes_use[!is.na(SNP_55) & !is.na(SNP_53) & !is.na(SNP_51) & !is.na(SNP_49) & is.na(SNP_47) & is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_55))
s55_53_51 <- length(unique(bRes_use[!is.na(SNP_55) & !is.na(SNP_53) & !is.na(SNP_51) & is.na(SNP_49) & is.na(SNP_47) & is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_55))
s55_53 <- length(unique(bRes_use[!is.na(SNP_55) & !is.na(SNP_53) & is.na(SNP_51) & is.na(SNP_49) & is.na(SNP_47) & is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_55))
s55 <- length(unique(bRes_use[!is.na(SNP_55) & is.na(SNP_53) & is.na(SNP_51) & is.na(SNP_49) & is.na(SNP_47) & is.na(SNP_45) & is.na(SNP_43) & is.na(SNP_41),]$SNP_55))
datCombine <- rbind(datCombine,c(s55_53_51_49_47_45_43_41,41,55))
datCombine <- rbind(datCombine,c(s55_53_51_49_47_45_43,43,55))
datCombine <- rbind(datCombine,c(s55_53_51_49_47_45,45,55))
datCombine <- rbind(datCombine,c(s55_53_51_49_47,47,55))
datCombine <- rbind(datCombine,c(s55_53_51_49,49,55))
datCombine <- rbind(datCombine,c(s55_53_51,51,55))
datCombine <- rbind(datCombine,c(s55_53,53,55))
datCombine <- rbind(datCombine,c(s55,55,55))


datCombine <- as.data.table(datCombine)
names(datCombine) <- c("countSig","Significant since","t")
datCombine[,`Significant since`:=as.factor(`Significant since`)]
datCombine[,`Significant since`:=factor(`Significant since`,levels=as.character(seq(55,41,-2)))]


fwrite(datCombine,"/nfs/scistore13/robingrp/sojavee/figData/Menopause_df2_mapAgeToAge")

} # end of first if
datCombine <- fread("/nfs/scistore13/robingrp/sojavee/figData/Menopause_df2_mapAgeToAge")
datCombine[,`Significant since`:=as.factor(`Significant since`)]
datCombine[,`Significant since`:=factor(`Significant since`,levels=as.character(seq(55,41,-2)))]
library(ggpubr)
library(ggsci)

#legend.title = element_blank(),
p_mapping <- ggplot(data=datCombine,aes(x=as.factor(t),y=countSig,fill=`Significant since`)) + geom_bar(position="stack", stat="identity") + theme_bw()  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
    theme( legend.text = element_text(size=13), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
    theme(title = element_text(color="gray20", size=12)) +
    theme(axis.title = element_text(color="gray20", size=14)) +
    theme(axis.text = element_text(color="gray40",size=13)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill="white")) +
    theme(strip.text=element_text(size = 8, color = "gray10")) +
    xlab("Age") + ylab("Number of significant regions") + scale_y_continuous(breaks = seq(0,150,30))

p_mapping <- p_mapping + scale_fill_nejm()
ggsave(plot=p_mapping,filename="/nfs/scistore13/robingrp/sojavee/figures/MenopauseChangeInEffect_ReplicatedEstonia_df2.pdf",width=7,height=5)


