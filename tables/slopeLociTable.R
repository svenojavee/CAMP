library(data.table)
library(xtable)


dat <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/ReplicatedNovelAssociationsSlope.csv")

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


geneDat <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/FUMA/FUMA_slopeSNPs/snps.txt")
geneDat <- geneDat[,list(rsID,nearestGene,dist,func)]
datUse <- merge(datUse,geneDat,by.x="SNP",by.y="rsID",all.x=T)

datUse[nearestGene == "RP11-166B2.8:GSPT1",nearestGene:="GSPT1"]
datUse[nearestGene == "RP5-864K19.4:MYCBP:GJA9",nearestGene:="MYCBP"]
datUse[SNP == "rs67596711",nearestGene:="ZNF275"]
datUse[SNP == "rs67596711",dist:="13176"]
datUse[SNP == "rs6631137",dist:="5893"]
datUse[SNP == "rs6631137",nearestGene:="GK"]
datUse[dist == "0:0:0" | dist == "0:0",dist:="0"]
datUse[,printGene:=nearestGene]
datUse[dist != "0",printGene:=paste0(nearestGene,"(",dist,")")]

datUse2 <- datUse[,list(SNP,Chr,bp,alleles,freq,printGene,b,p,PMaxUseT)]
datUse2 <- datUse2[order(Chr,bp),]

datUse2[,bp:=as.character(bp)]
datUse2[,freq:=as.character(round(freq,3))]
datUse2[,b:=as.character(round(b,4))]

names(datUse2) <- c("SNP","Chr","Position","Eff/Oth","MAF","Nearest gene","Yearly effect change","Slope p-value","Strongest p-value")

datUse2[,Chr:=as.character(Chr)]

print(xtable(datUse2,digits=-2),include.rownames=F,math.style.exponents=T)

