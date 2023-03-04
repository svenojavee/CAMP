library(data.table)
library(ggplot2)
library(ggpubr)
library(ggsci)

dat <- fread("/nfs/scistore13/robingrp/sojavee/figData/FigProportionAgeSpecSlope_replicated.txt")
dat[,perc:=countValSlope/sum(countValSlope),by=age]

dat[,SlopeChange:=factor(SlopeChange,levels=c("Time-varying effect","No change"))]

p1 <- ggplot(data=dat,aes(x=age,y=countValSlope,fill=SlopeChange)) + geom_bar(position="stack", stat="identity") + theme_bw()  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
    theme(legend.title = element_blank(), legend.text = element_text(size=13), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
    theme(title = element_text(color="gray20", size=12)) +
    theme(axis.title = element_text(color="gray20", size=14)) +
    theme(axis.text = element_text(color="gray40",size=13)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill="white")) +
    theme(strip.text=element_text(size = 8, color = "gray10")) +
    xlab("Age") + ylab("Number of siginificant regions")


p1 <- p1 + scale_fill_nejm()
ggsave(plot=p1,filename="/nfs/scistore13/robingrp/sojavee/figures/MenopausePropSignificant_replicated.pdf",width=7,height=5)


