library(data.table)
library(ggplot2)


datChange <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/colocalisation/ColocalisedGenes_ANM_eQTL_changeEvidenceSNPs.txt")
datNoChange <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/colocalisation/ColocalisedGenes_ANM_eQTL_noChangeEvidenceSNPs.txt")

datChange[,type:="Change"]
datNoChange[,type:="NoChange"]

mean(datChange$colocProb[datChange$colocProb > 0.01])
mean(datNoChange$colocProb[datChange$colocProb > 0.01])
datChange[colocProb>0.01]
datNoChange[colocProb>0.01]

datAll <- rbind(datChange[,list(colocProb,type)],datNoChange[,list(colocProb,type)])

p_ecdf <- ggplot(data=datAll[colocProb>0.01,],aes(x=colocProb,color=type)) + stat_ecdf(geom = "step") +  theme_bw()  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
    theme( legend.text = element_text(size=13), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
    theme(title = element_text(color="gray20", size=12)) +
    theme(axis.title = element_text(color="gray20", size=14)) +
    theme(axis.text = element_text(color="gray40",size=13)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill="white")) +
    theme(strip.text=element_text(size = 8, color = "gray10")) + xlim(c(0.01,1)) #+ xlab("Age") + ylab("Effect size (log HR)")

ggsave(file="/nfs/scistore13/robingrp/sojavee/figures/ECDF_colocDiff.pdf",plot=p_ecdf,width=5,height=4)
