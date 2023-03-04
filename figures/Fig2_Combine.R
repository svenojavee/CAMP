library(ggplot2)
library(data.table)
library(ggpubr)
library(ggsci)

#1) Curve of significance changing

source("/nfs/scistore13/robingrp/sojavee/ageSpecificMarginal/analyseLinearEffect/figures/Fig2_MapAgeToAgeReplicated.R")

#Make the max distribution as a separate plot
#source("/nfs/scistore13/robingrp/sojavee/ageSpecificMarginal/analyseLinearEffect/figures/Fig2_maxDistribution.R")


#2) Read the significant slopes

datSigSlopes <- fread("/nfs/scistore13/robingrp/sojavee/figureData/SignificantEffectPlotting.txt")

p_slopes <- ggplot(data=datSigSlopes,aes(x=as.factor(variable),y=value,group=SNP,color=(`Effect change`))) +
    #geom_line(aes(size=`Effect change`)) +
    geom_line(data=datSigSlopes[datSigSlopes$`Effect change`!="Increasing effect",],aes(as.factor(variable),y=value,color=`Effect change`),size=0.1) +  theme_bw()  +
    geom_line(data=datSigSlopes[datSigSlopes$`Effect change`=="Increasing effect",],aes(as.factor(variable),y=value,color=`Effect change`),size=0.3) +  theme_bw()  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
    theme( legend.text = element_text(size=13), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
    theme(title = element_text(color="gray20", size=12)) +
    theme(axis.title = element_text(color="gray20", size=14)) +
    theme(axis.text = element_text(color="gray40",size=13)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill="white")) +
    theme(strip.text=element_text(size = 8, color = "gray10")) +
    xlab("Age") + ylab("Effect (log HR)\n(significant variants)") + ylim(c(-0.22,0.2)) # + scale_size_manual(values = c(0.1,0.5)) + guides(size = "none") 

p_slopes <- p_slopes + scale_color_nejm(name=NULL) # scale_colour_manual(values = c("Increasing effect"="#EE4C97FF","Decreasing effect" = "#E18727FF"),name=NULL) #

#3) Pie chart to summarise significant findings
source("/nfs/scistore13/robingrp/sojavee/ageSpecificMarginal/analyseLinearEffect/figures/Fig2_ageSpecificProportion.R")

p <- ggarrange(p_mapping,ggarrange(p_slopes,p_piePlot,labels=c("b","c"),nrow=2,widths=c(1,1)) ,labels=c("a",""),nrow=1,widths=c(1,1), legend="top") 
ggsave(file="/nfs/scistore13/robingrp/sojavee/figures/FigAMNEffectDistribution.pdf",plot=p,width=12,height=7)
