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
fail <- fread(paste0("/location/of/fail/file/MenopauseMarginal.fail"))

#hrt file has columns with IID, age at HRT start (NA if no HRT)
hrt <- read.table("/location/of/HRT/file/ANM_HRT.txt", header=T)
phen$order <- 1:nrow(phen)
names(phen) <- c("FID","IID","ANM","order")
dat <- merge(phen, hrt, by="IID")
dat <- dat[order(dat$order),]
HRT_dat <- dat[,c(1,5)]

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
    genUse <- genDat[selectInd,bim$V2 == useSNPsChr$SNP[i]]
    combineDat[,SNP:=scale(genUse)]
    combineDatTmp <- combineDat[!is.na(SNP),]
    combineDatTmp2 <- copy(combineDatTmp)
    combineDatTmp2[fail == 1,V3:=(round(V3/0.5)*0.5)] #rounding helps with execution speed

    #This part creates the time-varying covariate for HRT (=1, if HRT is on)
    dat1 <- copy(combineDatTmp2)
    dat2 <- copy(combineDatTmp2)
  
    dat1[,t1:=0]
    dat1[!is.na(HRTage),t2:=pmin(V3,HRTage)]
    dat1[is.na(HRTage),t2:=V3]

    dat1[HRTage == t2,useMoreData:=T]
    dat1[,isHRT:=0]
    dat1[is.na(useMoreData),eventObserved:=fail] 
    dat1[!is.na(useMoreData),eventObserved:=0] #code as censoring if that id continued in the next cohort (dat2)
  
    dat2 <- dat2[V1 %in% dat1[useMoreData == T,]$V1,] #dat2 contains people who have received HRT prior to their menopause
    dat2[,t1:=HRTage]
    dat2[,t2:=V3]
    dat2[,isHRT:=1]
    dat2[,eventObserved:=fail]
  
    finDat <- rbind(dat1[,list(V1,t1,t2,covCombine,isHRT,genVal,SNP,eventObserved)],dat2[,list(V1,t1,t2,covCombine,isHRT,genVal,SNP,eventObserved)])

    m_linear <- coxph(Surv(t1,t2,eventObserved) ~ covCombine + genVal + SNP + tt(SNP) + isHRT, cluster = V1, data=finDat, tt=function(x, t, ...) {x*(t - t0)})
   
    summ <- summary(m_linear)
    coefs <- summ$coefficients
    cov_b0_b1 <- m_linear$var[3,4]

    b0 <- coefs[3,1]
    b1 <- coefs[4,1]

    z0 <- coefs[3,4]
    z1 <- coefs[4,4]
    t1 <- Sys.time()
    gc()
    print(Sys.time()-t1)
}


resDat2 <- as.data.table(resDat)
names(resDat2) <- c("SNP","b0","b1","t0","t1","cov_b0_b1")

fwrite(resDat2,file=paste0(wd,"/mp_",chr,"_L_",snpLow,"_H_",snpHigh),sep=" ")

print("CalculationDone")

