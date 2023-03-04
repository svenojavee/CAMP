library(data.table)
library(xtable)


dat <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/ReplicatedNovelAssociationsatMax_45_52.csv")
datFull <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/Menopause_Linear_FullCov.res")

datUse <- dat[,list(SNP,Chr,bp,refA,freq,b,p,whichMaxT)]

datUse <- merge(datUse,datFull[,list(SNP,A1,A2,B0,B1)],by="SNP",all.x=T)
datUse[,`alleles`:=paste0(A1,"/",A2)]
datUse[,whichMaxT:=as.character(format(round(whichMaxT,1),nsmall=1))]

geneDat <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/FUMA/FUMA_novelSNPs/snps.txt")
geneDat <- geneDat[,list(rsID,nearestGene,dist,func)]
datUse <- merge(datUse,geneDat,by.x="SNP",by.y="rsID",all.x=T)

datUse[,printGene:=nearestGene]
datUse[dist != "0",printGene:=paste0(nearestGene,"*")]

datUse2 <- datUse[,list(SNP,Chr,bp,alleles,freq,printGene,b,whichMaxT,p,B0,B1)]
datUse2 <- datUse2[order(Chr,bp),]

datUse2[,bp:=as.character(format(bp, nsmall=0, big.mark=","))]
datUse2[,freq:=as.character(format(round(freq,3),nsmall=3))]
datUse2[,b:=as.character(format(round(b,3),nsmall=3))]
datUse2[,B0:=as.character(format(round(B0,3),nsmall=3))]
datUse2[,B1:=as.character(format(round(B1,4),nsmall=4))]

names(datUse2) <- c("SNP","Chr","Position","Eff/Oth","MAF","Nearest gene","Maximum effect","Age-at-maximum","p-value","Effect at 49","Yearly effect size change")

datUse2[,Chr:=as.character(Chr)]

print(xtable(datUse2,digits=-2),include.rownames=F,math.style.exponents=T)

