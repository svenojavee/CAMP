library(data.table)
library(ggplot2)


#Read the 244 loci having an effect on ANM

replDat <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/ReplicatedNovelAssociationsatMax_45_52.csv")
prevDat <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/PreviousAssociationsatMax_45_52.csv")
setnames(prevDat,"Variant","SNP")

dat <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/Menopause_Linear_FullCov.res")

#Let's read the slopes that were replicated in Estonia

slopeReplDat <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/ReplicatedNovelAssociationsSlope.csv")[,list(SNP)]
slopeReplDat <- merge(slopeReplDat,dat[,list(SNP,A1,A2,B0,SE0,B1,SE1,COV_B0_B1)],by="SNP",all.x=T)
slopeReplDat[,maxT:=49 + (B0 * COV_B0_B1 - B1 * SE0**2)/(B1 * COV_B0_B1 - B0 * SE1**2)]
slopeReplDat[,chiMaxT:=(B0 + B1 * (maxT-49))**2 / (SE0**2 + (maxT-49)**2 * SE1**2 + 2 * (maxT-49) * COV_B0_B1)]

slopeReplDat[,chi45:=(B0 + B1 * (45-49))**2 / (SE0**2 + (45-49)**2 * SE1**2 + 2 * (45-49) * COV_B0_B1)]
slopeReplDat[,chi52:=(B0 + B1 * (52-49))**2 / (SE0**2 + (52-49)**2 * SE1**2 + 2 * (52-49) * COV_B0_B1)]
slopeReplDat[,useMaxT:=maxT]
slopeReplDat[(maxT < 45 | maxT > 52) & chi45 > chi52,useMaxT:=45]
slopeReplDat[(maxT < 45 | maxT > 52) & chi52 > chi45,useMaxT:=52]
slopeReplDat[,chiMaxUseT:=(B0 + B1 * (useMaxT-49))**2 / (SE0**2 + (useMaxT-49)**2 * SE1**2 + 2 * (useMaxT-49) * COV_B0_B1)]
slopeReplDat[,BMaxUseT:=(B0 + B1 * (useMaxT-49)) ]
slopeReplDat[,SEMaxUseT:=sqrt(SE0**2 + (useMaxT-49)**2 * SE1**2 + 2 * (useMaxT-49) * COV_B0_B1) ]
slopeReplDat[,PMaxUseT:=pchisq(chiMaxUseT,df=2,lower.tail=F)] #Let's use df=2 instead



allDat <- rbind(replDat[,list(SNP)],prevDat[,list(SNP)])
allDat <- merge(allDat,dat[,list(SNP,B1,P1)],all.x=T,by.x="SNP",by.y="SNP")



countStrongEvidence <- sum(slopeReplDat$PMaxUseT < 5e-8)
countModerateEvidence <- sum(allDat$P1 < 5e-8) - countStrongEvidence
countWeakEvidence <- sum(allDat$P1 < 5e-2) - countStrongEvidence - countModerateEvidence
countNoEvidence <- nrow(allDat) - countWeakEvidence - countModerateEvidence - countStrongEvidence

datPlot <- data.table(count=c(countStrongEvidence,countModerateEvidence,countWeakEvidence,countNoEvidence),evidenceType=c("Strong evidence","Moderate evidence","Weak evidence","No evidence"))

datPlot[,evidenceType:=factor(evidenceType,levels=c("No evidence","Weak evidence","Moderate evidence","Strong evidence"))]


p_piePlot <- ggplot(data = datPlot,aes(x="",y=count,fill=evidenceType)) + 
    geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + 
    #geom_text(aes(y = ypos, label = count), color = "black", size=5)  + # theme_bw() +
    geom_text(aes(label = count), position = position_stack(vjust = 0.5),size=5) + theme_void() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
    theme(legend.text = element_text(size=13), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
    theme(title = element_text(color="gray20", size=12)) +
    theme(axis.title = element_text(color="gray20", size=14)) +
    theme(axis.text = element_blank()) +
    theme(legend.position = "right") +
    theme(strip.background = element_rect(fill="white")) +
    theme(strip.text=element_text(size = 12, color = "gray10"))  + 
    xlab("") + ylab("")  + scale_fill_manual(values = c("Strong evidence"="#A00000","Moderate evidence" = "#F95050", "Weak evidence" = "#FE8282","No evidence"="#DDDDDD"),name="Evidence for\nage-specific effect")  



