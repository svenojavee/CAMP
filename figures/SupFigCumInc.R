library(data.table)
library(survival)
library(cmprsk)

#Let's create a cumulative incidence curve for age-at-menopause, taking into account death as competing risk

phenMenopause <- fread("/nfs/scistore13/robingrp/human_data/pheno/timeVaryingPhenotypes/MenopauseMarginal/MenopauseMarginal.phen")
failMenopause <- fread("/nfs/scistore13/robingrp/human_data/pheno/timeVaryingPhenotypes/MenopauseMarginal/MenopauseMarginal.fail")
phenMenopause[,censInd:=failMenopause$V1]

fam <- fread("/nfs/scistore13/robingrp/human_data/pheno/timeVaryingPhenotypes/MenopauseMarginal/ukb22828_c23_UKB_EST_v4.fam")
phenMenopause <- phenMenopause[V1 %in% fam$V1,]

deathDat <- fread("/nfs/scistore13/robingrp/human_data/pheno/ukb41236.csv",select=c("eid","40007-0.0"))
names(deathDat) <- c("V1","deathTime")

phenMenopause <- merge(phenMenopause,deathDat,by="V1",all.x=T)
phenMenopause[!is.na(deathTime) & censInd == 0,]
phenMenopause <- phenMenopause[order(deathTime),]

phenMenopause[censInd == 0 & !is.na(deathTime),V3:=deathTime]
phenMenopause[censInd == 0 & !is.na(deathTime),censInd:=2]

m <- cuminc(phenMenopause$V3, phenMenopause$censInd)
save(m,file="/nfs/scistore13/robingrp/sojavee/figData/MP_UK_comprisk.RData")


