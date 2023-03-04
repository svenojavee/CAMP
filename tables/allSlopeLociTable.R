library(data.table)
library(xtable)


dat <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/UnreplicatedNovelAssociationsSlope.csv")

datFull <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/Menopause_Linear_FullCov.res")


bim <- NULL
chrUse <- 1:23


datUse <- dat[,list(SNP,Chr,bp,refA,freq,b,p)]

datUse <- merge(datUse,datFull[,list(SNP,A1,A2,B0,SE0,B1,SE1,COV_B0_B1)],by="SNP",all.x=T)
datUse[,maxT:=49 + (B0 * COV_B0_B1 - B1 * SE0**2)/(B1 * COV_B0_B1 - B0 * SE1**2)]
datUse[,chiMaxT:=(B0 + B1 * (maxT-49))**2 / (SE0**2 + (maxT-49)**2 * SE1**2 + 2 * (maxT-49) * COV_B0_B1)]

datUse[,chi40:=(B0 + B1 * (40-49))**2 / (SE0**2 + (40-49)**2 * SE1**2 + 2 * (40-49) * COV_B0_B1)]
datUse[,chi55:=(B0 + B1 * (55-49))**2 / (SE0**2 + (55-49)**2 * SE1**2 + 2 * (55-49) * COV_B0_B1)]
datUse[,useMaxT:=maxT]
datUse[maxT < 40 | maxT > 55 & chi40 > chi55,useMaxT:=40]
datUse[maxT < 40 | maxT > 55 & chi55 > chi40,useMaxT:=55]
datUse[,chiMaxUseT:=(B0 + B1 * (useMaxT-49))**2 / (SE0**2 + (useMaxT-49)**2 * SE1**2 + 2 * (useMaxT-49) * COV_B0_B1)]
datUse[,BMaxUseT:=(B0 + B1 * (useMaxT-49)) ]
datUse[,SEMaxUseT:=sqrt(SE0**2 + (useMaxT-49)**2 * SE1**2 + 2 * (useMaxT-49) * COV_B0_B1) ]
datUse[,PMaxUseT:=pchisq(chiMaxUseT,df=2,lower.tail=F)] #Let's use df=2 instead


datUse[,`alleles`:=paste0(A1,"/",A2)]


geneDat <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/FUMA/FUMA_unreplicatedSlopeSNPs/snps.txt")
geneDat <- geneDat[,list(rsID,nearestGene,dist,func)]
datUse <- merge(datUse,geneDat,by.x="SNP",by.y="rsID",all.x=T)


datUse[SNP == "rs186475193",nearestGene:="PDE4B"]
datUse[SNP == "rs186475193",dist:="0"]
datUse[SNP == "rs12559662",dist:="0"]
datUse[SNP == "rs12559662",nearestGene:="EFHC2"]
datUse[SNP == "rs5951636",dist:="16121"]
datUse[SNP == "rs5951636",nearestGene:="MBTPS2"]

datUse[nearestGene == "STON1-GTF2A1L:GTF2A1L",nearestGene:="GTF2A1L"]
datUse[nearestGene == "AFAP1-AS1:AFAP1",nearestGene:="AFAP1"]
datUse[nearestGene == "ENTPD1-AS1:CCNJ",nearestGene:="CCNJ"]


datUse[,printGene:=nearestGene]
datUse[dist != "0",printGene:=paste0(nearestGene,"*")]

datUse2 <- datUse[,list(SNP,Chr,bp,alleles,freq,printGene,b,p,PMaxUseT)]
datUse2 <- datUse2[order(Chr,bp),]

datUse2[,bp:=as.character(bp)]
datUse2[,freq:=as.character(round(freq,3))]
datUse2[,b:=as.character(round(b,4))]

names(datUse2) <- c("SNP","Chr","Position","Eff/Oth","MAF","Nearest gene","Yearly effect change","Slope p-value","Strongest p-value")

datUse2[,Chr:=as.character(Chr)]

print(xtable(datUse2,digits=-2),include.rownames=F,math.style.exponents=T)

