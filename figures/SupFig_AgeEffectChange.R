library(data.table)
setDTthreads(5)

replDat <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/ReplicatedNovelAssociationsatMax_45_52.csv")
prevDat <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/PreviousAssociationsatMax_45_52.csv")
setnames(prevDat,"Variant","SNP")



X <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/Menopause_Linear_FullCov.res")

allDat <- rbind(replDat[,list(SNP)],prevDat[,list(SNP)])
allDat <- merge(allDat,X[,list(SNP,B0,B1)],all.x=T,by.x="SNP",by.y="SNP")

allDat[,`41`:=B0+B1*(41-49)]
allDat[,`43`:=B0+B1*(43-49)]
allDat[,`45`:=B0+B1*(45-49)]
allDat[,`47`:=B0+B1*(47-49)]
allDat[,`49`:=B0+B1*(49-49)]
allDat[,`51`:=B0+B1*(51-49)]
allDat[,`53`:=B0+B1*(53-49)]
allDat[,`55`:=B0+B1*(55-49)]

allDat[abs(`41`) < abs(`55`) ,`Effect change`:="Increasing effect" ]
allDat[abs(`41`) >= abs(`55`) ,`Effect change`:="Decreasing effect" ]

allDatOutlier <- allDat[,list(SNP,`Effect change`)]
allDat[,c("B0","B1","Effect change"):=NULL]

allDatLong <- melt(allDat,id.vars="SNP")
allDatLong[,variable:=as.numeric(as.character(variable))]

allDatLong <- merge(allDatLong,allDatOutlier,by="SNP",all.x=T)
allDatLong[`Effect change` == "Decreasing effect",lineSize:=0.1]
allDatLong[`Effect change` == "Increasing effect",lineSize:=0.5]

fwrite(allDatLong,file="/nfs/scistore13/robingrp/sojavee/figureData/SignificantEffectPlotting.txt")


#Let's make a comparison with nominally insignificant SNPs (p>0.05)
dat <- copy(X)

#Calculate at which age the chi-squared statistic is the largest
dat[,maxT:=49 + (B0 * COV_B0_B1 - B1 * SE0**2)/(B1 * COV_B0_B1 - B0 * SE1**2)]
quantile(dat$maxT,probs=c(0.05,0.1,0.9,0.95))
dat[,chiMaxT:=(B0 + B1 * (maxT-49))**2 / (SE0**2 + (maxT-49)**2 * SE1**2 + 2 * (maxT-49) * COV_B0_B1)]

#For the max values that are outside the interval [45,52], let's check the points 45 and 52
dat[,chi45:=(B0 + B1 * (45-49))**2 / (SE0**2 + (45-49)**2 * SE1**2 + 2 * (45-49) * COV_B0_B1)]
dat[,chi52:=(B0 + B1 * (52-49))**2 / (SE0**2 + (52-49)**2 * SE1**2 + 2 * (52-49) * COV_B0_B1)]

dat[,useMaxT:=maxT]
dat[(maxT < 45 | maxT > 52) & chi45 > chi52,useMaxT:=45]
dat[(maxT < 45 | maxT > 52) & chi52 > chi45,useMaxT:=52]
dat[,chiMaxUseT:=(B0 + B1 * (useMaxT-49))**2 / (SE0**2 + (useMaxT-49)**2 * SE1**2 + 2 * (useMaxT-49) * COV_B0_B1)]
dat[,BMaxUseT:=(B0 + B1 * (useMaxT-49)) ]
dat[,SEMaxUseT:=sqrt(SE0**2 + (useMaxT-49)**2 * SE1**2 + 2 * (useMaxT-49) * COV_B0_B1) ]

dat[,PMaxUseT:=pchisq(chiMaxUseT,df=2,lower.tail=F)] #Let's use df=2 instead
#Keep only the nominally insignificant results
dat2 <- dat[PMaxUseT >0.05,list(SNP,B0,B1)]
rm(X,dat)
gc()
gridVals <- seq(41,52,0.5)
fullCJ <- CJ(gridVals,dat2$SNP)
names(fullCJ) <- c("age","SNP")

fullCJ <- merge(fullCJ,dat2,by="SNP",all.x=T)
fullCJ[,b:=B0+B1*(age-49)]
fullCJ[,B0:=NULL]
fullCJ[,B1:=NULL]
fullCJ[,SNP:=NULL]

fullCJ[,q25:=round(quantile(b,probs=0.25),5),by=age]
fullCJ[,q75:=round(quantile(b,probs=0.75),5),by=age]
gc()
fullCJ[,q025:=round(quantile(b,probs=0.025),5),by=age]
fullCJ[,q975:=round(quantile(b,probs=0.975),5),by=age]
gc()
fullCJ[,q995:=round(quantile(b,probs=0.995),5),by=age]
fullCJ[,q005:=round(quantile(b,probs=0.005),5),by=age]
gc()
fullCJ[,q9995:=round(quantile(b,probs=0.9995),5),by=age]
fullCJ[,q0005:=round(quantile(b,probs=0.0005),5),by=age]
gc()

sub_fullCJ <- fullCJ[1:100,] #The first 100 lines contain already the information 
finalNoEffectDat <- unique(sub_fullCJ[,list(age,q25,q75,q025,q975,q995,q005,q9995,q0005)])

#This takes too much memory
#finalNoEffectDat <- unique(fullCJ[,list(age,q025,q975,q995,q005,q9995,q0005)])
fwrite(finalNoEffectDat,file="/nfs/scistore13/robingrp/sojavee/figureData/NoSignificantEffectDistribution.txt")

if(F){


p_mapping <- ggplot(data=allDatLong,aes(x=as.factor(variable),y=value,group=SNP,color=`Effect change`)) + geom_line(aes(size=`Effect change`)) +  theme_bw()  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
    theme( legend.text = element_text(size=13), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
    theme(title = element_text(color="gray20", size=12)) +
    theme(axis.title = element_text(color="gray20", size=14)) +
    theme(axis.text = element_text(color="gray45",size=13)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill="white")) +
    theme(strip.text=element_text(size = 8, color = "gray10")) +
    xlab("Age") + ylab("Effect size (log HR)")  + scale_size_manual(values = c(0.1,0.5)) + guides(size = "none")

p_mapping <- p_mapping + scale_color_nejm(name=NULL)
ggsave(plot=p_mapping,filename="/nfs/scistore13/robingrp/sojavee/figures/MenopauseSignificantEffectChange.pdf",width=7,height=5)
}