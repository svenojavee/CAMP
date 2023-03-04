library(ggplot2)
library(data.table)
library(ggpubr)
library(ggsci)


dat <- fread("/nfs/scistore13/robingrp/sojavee/figData/Menopause_df2_atMax_distribution")

datSig <- dat[PMaxT < 5e-8,]

dat[,plotType:="All SNPs"]
datSig[,plotType:="Significant SNPs"]

datAll <- rbind(dat,datSig)
#Let's truncate between 30 and 70
datAll <- datAll[maxT > 30 & maxT < 70,]

p_maxDist <- ggplot(data=datAll,aes(x=maxT)) + geom_density(fill="coral") + theme_bw()  + facet_wrap(.~plotType,nrow=1,ncol=2) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
    theme( legend.text = element_text(size=13), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
    theme(title = element_text(color="gray20", size=12)) +
    theme(axis.title = element_text(color="gray20", size=14)) +
    theme(axis.text = element_text(color="gray40",size=13)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill="white")) +
    theme(strip.text=element_text(size = 12, color = "gray10")) +
    xlab("Age with strongest effect evidence") + ylab("Density") 

p_maxDist <- p_maxDist + scale_fill_nejm()
ggsave(plot=p_maxDist,filename="/nfs/scistore13/robingrp/sojavee/figures/MenopauseMaxEffectDistribution.pdf",width=5,height=2.5)

