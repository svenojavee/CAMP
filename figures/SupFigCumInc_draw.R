library(data.table)
library(survival)
library(cmprsk)
library(ggplot2)

load("/nfs/scistore13/robingrp/sojavee/figData/MP_UK_comprisk.RData")
mUK <- m
load("/nfs/scistore13/robingrp/sojavee/figData/MP_EST_comprisk.RData")
mEST <- m



UK_res <- mUK[[1]]
EST_res <- mEST[[1]]


datUK <- data.table(time=UK_res$time,est=UK_res$est,Data="UK Biobank")
datEST <- data.table(time=EST_res$time,est=EST_res$est,Data="Estonian Biobank")

datFull <- rbind(datUK,datEST)

datFull[,est:= 100 * est]

datFull[Data == "UK Biobank"]
datFull[time >40 & time<42]

p_Inc <- ggplot(data = datFull,aes(x=time,y=est,col=Data)) + 
     geom_line() + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
    theme(legend.title = element_blank(), legend.text = element_text(size=14), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
    theme(title = element_text(color="gray20", size=12)) +
    theme(axis.title = element_text(color="gray20", size=14)) +
    theme(axis.text = element_text(color="gray40",size=14)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill="white")) +
    theme(strip.text=element_text(size = 12, color = "gray10"))  + xlim(c(35,65)) + 
    xlab("Age") + ylab("Cumulative incidence (%)") 

library(ggpubr)
library(ggsci)

natColours <- pal_npg()(10)
natColours[11] <- '#778899'
natColours <- pal_npg()(10)
natColours[11] <- '#778899'
d1 <- natColours[c(1,4,7,10)]
p_Inc <- p_Inc + scale_colour_manual(values=d1) + scale_fill_manual(values=d1) 

ggsave(plot=p_Inc,filename="/nfs/scistore13/robingrp/sojavee/figures/SupFig_Menopause_Incidence.pdf",height=4,width=5.333)


