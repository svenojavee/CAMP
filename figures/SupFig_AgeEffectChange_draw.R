library(data.table)
library(ggplot2)
library(ggpubr)
library(ggsci)

#1) Read the significant slopes

datSigSlopes <- fread("/nfs/scistore13/robingrp/sojavee/figureData/SignificantEffectPlotting.txt")

p_mapping <- ggplot(data=datSigSlopes,aes(x=as.factor(variable),y=value,group=SNP,color=`Effect change`)) + geom_line(aes(size=`Effect change`)) +  theme_bw()  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
    theme( legend.text = element_text(size=13), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
    theme(title = element_text(color="gray20", size=12)) +
    theme(axis.title = element_text(color="gray20", size=14)) +
    theme(axis.text = element_text(color="gray40",size=13)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill="white")) +
    theme(strip.text=element_text(size = 8, color = "gray10")) +
    xlab("Age") + ylab("Effect size (log HR)")  + scale_size_manual(values = c(0.1,0.5)) + guides(size = "none") + ylim(c(-0.22,0.2))

p_mapping <- p_mapping + scale_color_nejm(name=NULL)

#2) Read the quantiles of non-significant slopes
datNonSigSlopes <- fread("/nfs/scistore13/robingrp/sojavee/figureData/NoSignificantEffectDistribution.txt")
datNonSigSlopes_long <- melt(datNonSigSlopes,id.vars=c("age"))
datNonSigSlopes_long[variable %in% c("q25","q75"),CItype:="IQR"]
datNonSigSlopes_long[variable %in% c("q025","q975"),CItype:="95% CI"]
datNonSigSlopes_long[variable %in% c("q005","q995"),CItype:="99% CI"]
datNonSigSlopes_long[variable %in% c("q0005","q9995"),CItype:="99.9% CI"]

datNonSigSlopes_long_IQR_q25 <- datNonSigSlopes_long[CItype == "IQR" & variable == "q25",]
datNonSigSlopes_long_IQR_q75 <- datNonSigSlopes_long[CItype == "IQR" & variable == "q75",]

datNonSigSlopes_long_IQR_q75 <- datNonSigSlopes_long_IQR_q75[,list(age,value)]
datNonSigSlopes_long_IQR_q25 <- datNonSigSlopes_long_IQR_q25[,list(age,value)]
names(datNonSigSlopes_long_IQR_q75) <- c("age","hIQR")
names(datNonSigSlopes_long_IQR_q25) <- c("age","lIQR")

datNonSigSlopes_long_CI <- datNonSigSlopes_long[CItype != "IQR",]
datNonSigSlopes_long_CI <- merge(datNonSigSlopes_long_CI,datNonSigSlopes_long_IQR_q75,by="age",all.x=T)
datNonSigSlopes_long_CI <- merge(datNonSigSlopes_long_CI,datNonSigSlopes_long_IQR_q25,by="age",all.x=T)


p_nonSig <- ggplot(data=datNonSigSlopes_long,aes(x=age,y=value,group=variable,color=CItype)) + geom_line() +  theme_bw()  + #+ geom_ribbon(aes(ymin=lIQR,ymax=hIQR),colour = NA) 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
    theme( legend.text = element_text(size=13), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
    theme(title = element_text(color="gray20", size=12)) +
    theme(axis.title = element_text(color="gray20", size=14)) +
    theme(axis.text = element_text(color="gray40",size=13)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill="white")) +
    theme(strip.text=element_text(size = 8, color = "gray10")) +
    xlab("Age") + ylab("Effect size (log HR)") + scale_x_continuous(breaks=seq(41,55,2)) + ylim(c(-0.22,0.2))

p_nonSig <- p_nonSig + scale_colour_manual(values = c("IQR"="#A00000","95% CI" = "#F95050", "99% CI" = "#FE8282","99.9% CI"="#FFB0B0"),name=NULL) #scale_color_nejm(name=NULL)

p <- ggarrange(p_mapping,p_nonSig,labels=c("a","b"),nrow=1,widths=c(1,1), legend="top",align="h") 
ggsave(file="/nfs/scistore13/robingrp/sojavee/figures/FigSlope_Effect_NoEffect.pdf",plot=p,width=12,height=6)


#,fill="#CD0000"