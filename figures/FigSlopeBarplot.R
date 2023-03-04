library(data.table)
library(ggplot2)


unreplSlope <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/UnreplicatedNovelAssociationsSlope.csv")

replSlope <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/ReplicatedNovelAssociationsSlope.csv")
datUse <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/Menopause_Linear_FullCov.res")

sum(unreplSlope$SNP %in% replSlope$SNP)

datSlope <- datUse[SNP %in% replSlope$SNP,]
tt_grid <- c(40,45,50,55)
datJoin <- NULL
for(tt in tt_grid){
    datSlope[,hazardUse:= exp(B0 + B1 * (tt - 49))]
    datSlope[,tUse:=tt]
    datJoin <- rbind(datJoin,datSlope[,list(SNP,CHR,BP,tUse,hazardUse)])
}
datJoin[,chrBP:=paste0("chr",CHR,":",format(BP, big.mark = ",", scientific = FALSE) )]
datJoin <- datJoin[order(CHR,BP),]
datJoin[,chrBP_fact:=factor(chrBP,levels=unique(datJoin$chrBP))]
datJoin[,chHaz:=(hazardUse-1)*100]
datJoin[,meanHaz:=mean(chHaz),by=SNP]

#
datJoin[,chrNum:=CHR]
datJoin[chrNum=="X",chrNum:="23"]
datJoin[,chrNum:=as.numeric(chrNum)]

datJoin <- datJoin[order(chrNum,BP)]
datJoin[,chrBP_fact:=factor(chrBP,levels=unique(datJoin$chrBP))]
textSize <- 8
p <- ggplot(data=datJoin, aes(x=tUse,y=chHaz)) + geom_bar(stat="identity",fill="coral") + theme_bw() +
    facet_wrap(.~chrBP_fact,nrow=4) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
    theme(legend.title = element_blank(), legend.text = element_text(size=textSize), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
    theme(title = element_text(color="gray20", size=textSize)) +
    theme(axis.title = element_text(color="gray20", size=12)) +
    theme(axis.text = element_text(color="gray40",size=textSize)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill="white")) +
    theme(strip.text=element_text(size = textSize, color = "gray10")) +
    xlab("Age") + ylab("Minor allele hazard difference (%)") + geom_hline(yintercept=0,lty=3) + geom_hline(aes(yintercept=meanHaz),col="red",lty=2)
   

ggsave(p,filename="/nfs/scistore13/robingrp/sojavee/figures/MP_SlopeChange.pdf",width=6,height=5)

###############
#Question: Do the novel slopes map to novel findings or previous discoveries
novelSNPs <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/ReplicatedNovelAssociationsAtMax.csv")
novelSNPs[SNP %in% replSlope$SNP]



replNovel <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/ReplicatedNovelAssociationsAtMax.csv")
replNovel2 <- replNovel[SNP %in% c("rs11638671","rs7946546"),]

datSlope <- datUse[SNP %in% replNovel2$SNP,]
tt_grid <- c(40,45,50,55)
datJoin <- NULL
for(tt in tt_grid){
    datSlope[,hazardUse:= exp(B0 + B1 * (tt - 49))]
    datSlope[,tUse:=tt]
    datJoin <- rbind(datJoin,datSlope[,list(SNP,CHR,BP,tUse,hazardUse)])
}
datJoin[,chrBP:=paste0("chr",CHR,":",format(BP, big.mark = ",", scientific = FALSE) )]
datJoin <- datJoin[order(CHR,BP),]
datJoin[,chrBP_fact:=factor(chrBP,levels=unique(datJoin$chrBP))]
datJoin[,chHaz:=(hazardUse-1)*100]
datJoin[,meanHaz:=mean(chHaz),by=SNP]

p <- ggplot(data=datJoin, aes(x=tUse,y=chHaz)) + geom_bar(stat="identity",fill="coral") + theme_bw() +
    facet_wrap(.~chrBP_fact,nrow=2,scales="free") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
    theme(legend.title = element_blank(), legend.text = element_text(size=13), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
    theme(title = element_text(color="gray20", size=13)) +
    theme(axis.title = element_text(color="gray20", size=13)) +
    theme(axis.text = element_text(color="gray40",size=13)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill="white")) +
    theme(strip.text=element_text(size = 13, color = "gray10")) +
    xlab("Age") + ylab("Minor allele hazard difference (%)") + geom_hline(yintercept=0,lty=3) + geom_hline(aes(yintercept=meanHaz),col="red",lty=2)
   
ggsave(p,filename="/nfs/scistore13/robingrp/sojavee/figures/MP_SlopeChange_example.pdf",width=4,height=4)


