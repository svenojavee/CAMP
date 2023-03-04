library(survival)
library(data.table)
library(BEDMatrix)

#Specifies the hazard value at the intercept
#Only relevant for the b0 parameter interpretation
t0 <- 49
#Specify the location of intermediate results
wd <- "/location/for/intermediate/results"

args = commandArgs(trailingOnly=TRUE)
chr <- as.numeric(args[1])
snpLow <- as.numeric(args[2])
snpHigh <- as.numeric(args[3])

#Specify the location of .fam and .bim files (chr-specific)
fam <- fread(paste0("/location/of/fam/file/fileName_c",chr,".fam"))
bim <- fread(paste0("/location/of/bim/file/fileName_c",chr,".bim"))


#1) Read the covariate and phenotypic data
#cov file has 3 columns, IID, summarised covariate value, and the BayesW predictor (different for each chr)
#phen file has 3 columns, IID, FID (== IID), and the phenotype value
#fail file has 1 column of 0s and 1s (0=censoring, 1=observed value)
covDat <- fread(paste0("/location/of/cov/file/MenopauseMarginal_",chr,".cov"))
phen <- fread(paste0("/location/of/phen/file/MenopauseMarginal.phen"))
fail <- fread(paste0("location/of/fail/file/MenopauseMarginal.fail"))

phen[,fail:=fail$V1]
phen[,V2:=NULL]
phen <- merge(phen,covDat,by.x="V1",by.y="IID",all.x=T)

#Check which individuals we choose
selectInd <- fam$V1 %in% phen$V1
ids <- fam$V1[selectInd]

combineDat <- data.table(V1=ids)
combineDat[,ord:=1:.N]

combineDat <- merge(combineDat,phen,by="V1")
combineDat <- combineDat[order(ord),]
combineDat[,ord:=NULL]

combineDat[,covCombine:=scale(covCombine)]
combineDat[,genVal:=scale(genVal)]

rm(phen,fam,covDat,fail)
gc()

#2) Read the genetic data (chr-specific)
genDat <- BEDMatrix(paste0("/location/of/bed/file_c",chr))


resDat <- NULL
for(i in snpLow:snpHigh){
    print(i)
    genUse <- genDat[selectInd,i]
    combineDat[,SNP:=scale(genUse)]
    combineDatTmp <- combineDat[!is.na(SNP),]
    combineDatTmp2 <- copy(combineDatTmp)
    #We round to 0.5 years, the number of other values is very small but having them will increase the
    #run time significantly
    combineDatTmp2[fail == 1,V3:=(round(V3/0.5)*0.5)]

    m_linear <- coxph(Surv( V3, fail) ~ covCombine + genVal + SNP + tt(SNP), data=combineDatTmp2, tt=function(x, t, ...) {x*(t - t0)})
    summ <- summary(m_linear)
    coefs <- summ$coefficients
    cov_b0_b1 <- m_linear$var[3,4]

    b0 <- coefs[3,1]
    b1 <- coefs[4,1]

    z0 <- coefs[3,4]
    z1 <- coefs[4,4]

    resDat <- rbind(resDat,c(bim$V2[i],b0,b1,z0,z1,cov_b0_b1))
    t1 <- Sys.time()
    gc()
    print(Sys.time()-t1)
}


resDat2 <- as.data.table(resDat)
names(resDat2) <- c("SNP","b0","b1","t0","t1","cov_b0_b1")

fwrite(resDat2,file=paste0(wd,"/mp_",chr,"_L_",snpLow,"_H_",snpHigh),sep=" ")

print("CalculationDone")

