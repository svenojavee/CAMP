library(data.table)

#Read the full data with results
dat <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/Menopause_Linear_FullCov.res")
failFile <- fread("/nfs/scistore13/robingrp/human_data/pheno/timeVaryingPhenotypes/menopauseUK1.fail")
Ntot <- sum(failFile$V1)

#Calculate at which age the chi-squared statistic is the largest
dat[,maxT:=49 + (B0 * COV_B0_B1 - B1 * SE0**2)/(B1 * COV_B0_B1 - B0 * SE1**2)]
quantile(dat$maxT,probs=c(0.05,0.1,0.9,0.95))
dat[,chiMaxT:=(B0 + B1 * (maxT-49))**2 / (SE0**2 + (maxT-49)**2 * SE1**2 + 2 * (maxT-49) * COV_B0_B1)]
dat[,N:=Ntot]

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

#Read the clumps
clDat <- NULL
for(chr in 1:23){
    if(chr == 23){
        clDat <- rbind(clDat,fread(paste0("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/UKB_8.4M_",chr,"_clumping.clumped"))[,list(SNP)])
    }else{
        clDat <- rbind(clDat,fread(paste0("/nfs/scistore13/robingrp/human_data/geno/ldp08/LDblock_clumpfile_",chr,".tmp.snps")))
    }
}

dat2 <- dat[SNP %in% clDat$SNP,]
dat2[,medVal:=median(chiMaxUseT),by=CHR]
unique(dat2[,list(CHR,medVal)])
summary(dat2$chiMaxUseT)

params <- list(df=2)

p <- ggplot(data=dat2, aes(sample=chiMaxUseT)) + 
    stat_qq(distribution = qchisq, dparams = params["df"]) +
    stat_qq_line(distribution = qchisq, dparams = params["df"],col="red") + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
    theme(legend.title = element_blank(), legend.text = element_text(size=13), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
    theme(title = element_text(color="gray20", size=13)) +
    theme(axis.title = element_text(color="gray20", size=13)) +
    theme(axis.text = element_text(color="gray40",size=13)) +
    theme(legend.position = "top") +
	theme(strip.background = element_rect(fill="white")) +
	theme(strip.text=element_text(size = 13, color = "gray10"))  + xlab("Expected chi-squared statistic") + ylab("Observed chi-squared statisitc") +
    labs(title="\u03bb = 1.044" )

ggsave(p,filename="/nfs/scistore13/robingrp/sojavee/figures/QQPlot_ANM.jpeg",width=6,height=6)

