library(data.table)
library(ggplot2)
library(survival)
library(dplyr)
library(ggpubr)
library(ggsci)

datUse <- fread("/nfs/scistore13/robingrp/sojavee/figData/MP_Illustration.txt")

#1) Illustration
nn <- 500
t0 <- 37
#fwrite(datUse,file="exampleCurvePlot.txt")
m <- coxph(Surv(y,rep(1,2*nn))~ x+ tt(x), tt=function(x, t, ...) {x * (t-t0)},data=datUse)
m
datNew <- data.table(xx=c(rep(1,1),0))
#fit <- survfit(m, newdata = datNew)  #Does not work

#Longer path
dtimes <- sort(unique(with(datUse,y[status==1])))

datExtended <- as.data.table(survSplit(Surv(y,status==1)~.,datUse,cut=dtimes))

datExtended[,ttk:= x * (y-t0)] 
#Call Cox PH

m2 <- coxph(Surv(tstart,y,event)~ttk+x ,data=datExtended)
m2

m_const <- coxph(Surv(y,rep(1,2*nn))~x,data=datUse)
m_const

datNew <- data.table(x=c(1,0),id=c(1,0))
fit <- survfit(m_const, newdata = datNew,se.fit=T)  
plot(fit,col=c(4,2),xlim=c(30,75),conf.int = T)


nd1 <- data.table(tstart=c(0,dtimes[1:(length(dtimes)-1)]),y=dtimes,trt=rep(1,length(dtimes)),prior=rep(2,length(dtimes)),x=rep(0,length(dtimes)))
nd1[,"event"] <- 0
nd1[,"ttk"] <- nd1$x * (nd1$y-t0)

nd2 <- nd1
nd2[,"x"] <- 1
nd2[,"ttk"] <- nd2$x * (nd1$y-t0)

nd1[,"id"] <- "1"
nd2[,"id"] <- "2"


#Plot the constant effect assumption
m_const_sum <- summary(m_const)$coefficients
f_const <- function(t){
  return(as.numeric(t*0+coef(m_const)))
}
f_lCI <- function(t){
  return(as.numeric(t*0 + m_const_sum[1] - 1.96 * (m_const_sum[3]+0.0085) ))
}
f_uCI <- function(t){
  return(as.numeric(t*0 + m_const_sum[1] + 1.96 * (m_const_sum[3]+0.0085)))
}



m_2_sum <- summary(m2)$coefficients
corPar <- (m2)$var[2,1]
f_2 <- function(t){
  valR <- (m_2_sum[1,1] * (t-t0) + m_2_sum[2,1])
  return(valR * 1)
}
f2_lCI <- function(t){
  return(f_2(t) - 1* 1.96 * sqrt(m_2_sum[2,3] ** 2 + 1 * (t-t0) ** 2 * m_2_sum[1,3] **2  + 1.0* 2*(t-t0)*corPar  + 0.00008*(t-53)**2) )
}
f2_uCI <- function(t){
  return(f_2(t) + 1.96 * sqrt(m_2_sum[2,3] ** 2 + 1* (t-t0) ** 2 * m_2_sum[1,3] **2  + 1.0 * 2*(t-t0)*corPar + 0.00008*(t-53)**2)  )
}


#########
#Let's do a nice ggplot
sf_change <- (survfit(m2,newdata=rbind(nd1,nd2),id=id,se.fit=T))
sf_const <- survfit(m_const, newdata = datNew,se.fit=T)  

datChange <- data.table(time=sf_change$time,strat=c(rep(1,length(sf_change$time)/2),rep(2,length(sf_change$time)/2) ),surv=sf_change$surv,lci=sf_change$lower,uci=sf_change$upper)
datConst_tmp <- data.table(time=sf_const$time,surv=sf_const$surv,lci=sf_const$lower,uci=sf_const$upper)

datConst1 <- datConst_tmp[,list(time,2,surv.1,lci.1,uci.1)]
datConst2 <- datConst_tmp[,list(time,1,surv.2,lci.2,uci.2)]

names(datConst1) <- c("time","strat","surv","lci","uci")
names(datConst2) <- c("time","strat","surv","lci","uci")
datConst <- rbind(datConst1,datConst2)

datChange[,mType:="Linear effect assumption"]
datConst[,mType:="Constant effect assumption"]
datCurve <- rbind(datChange,datConst)

#datCurve[,mType:=factor(mType,levels=c("Linear effect assumption","Constant effect assumption"))]


datCurve[strat == 1,strChar:="Risk\nvariant"]
datCurve[strat == 2,strChar:="Baseline"]
datCurve[,strChar:=factor(strChar,levels=c("Risk\nvariant","Baseline"))]

p_curves <- ggplot(data = datCurve,aes(x=time,y=surv,col=strChar,fill=strChar)) + geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, colour = NA) +
  geom_line(size = 0.8, lty=1) + facet_wrap(.~mType) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
  theme(legend.title = element_blank(), legend.text = element_text(size=12), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
  theme(title = element_text(color="gray20", size=12)) +
  theme(axis.title = element_text(color="gray20", size=12)) +
  theme(axis.text = element_text(color="gray40",size=12)) +
  theme(legend.position = "bottom") +
  theme(strip.background = element_rect(fill="white"),strip.text.x = element_blank()) + 
  theme(strip.text=element_text(size = 12, color = "gray10")) +
  xlab("Age") + ylab("Probability of no event") + xlim(c(30,70)) + scale_fill_nejm() + scale_color_nejm()


#And a nice linear change plot
linDatChange <- data.table(tt=seq(30,70,0.05))
linDatChange[,eff:=f_2(tt)]
linDatChange[,lci:=f2_lCI(tt)]
linDatChange[,uci:=f2_uCI(tt)]
linDatChange[,mType:="Linear effect assumption"]

constDat <- data.table(tt=seq(30,70,0.05))
constDat[,eff:=f_const(tt)]
constDat[,lci:=f_lCI(tt)]
constDat[,uci:=f_uCI(tt)]
constDat[,mType:="Constant effect assumption"]

dat <- rbind(linDatChange,constDat)
dat[,mType:=factor(mType,levels=c("Linear effect assumption","Constant effect assumption"))]


#Tune the plot with segments
dat[mType == "Linear effect assumption",xSeg:=52.25]
dat[mType == "Linear effect assumption",xEnd:=70]

dat[mType == "Linear effect assumption",ySeg:=0]
dat[mType == "Linear effect assumption",yEnd:=0]

#
true_fun <- function(t){
  return(-0.55+0.015*t)
}
dat[,trueEff:=true_fun(tt)]

dat1 <- copy(dat)
dat1[,eff:=trueEff]
dat[,trueEff:=NULL]
dat1[,trueEff:=NULL]
dat1[,lt:="True\neffect"]
dat[,lt:="Estimated\neffect"]

datAll <- rbind(dat,dat1)
datAll[,mType:=factor(mType,levels=c("Linear effect assumption","Constant effect assumption"))]

p_effect <- ggplot(data = datAll,aes(x=tt,y=eff)) +  geom_ribbon(fill="red",aes(ymin = lci, ymax = uci), alpha = 0.2) +
  geom_line(size = 1.1, aes(lty=lt)) + facet_wrap(.~mType) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
  theme(legend.title = element_blank(), legend.text = element_text(size=12), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
  theme(title = element_text(color="gray20", size=12),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(axis.title = element_text(color="gray20", size=12)) +
  theme(axis.text = element_text(color="gray40",size=12)) +
  theme(legend.position = "top") +
  theme(strip.background = element_rect(fill="white")) +
  theme(strip.text=element_text(size = 12, color = "gray10")) +
  xlab("Age") + ylab("Effect size") + geom_hline(yintercept=0,lty=3)  #+ geom_abline(data = data.frame(xint=-0.55,sl=0.015), aes(intercept = xint,slope=sl,linetype="True value"))




p_seg <- p_effect + geom_segment(col=2,lty=1,aes(x = xSeg, y = ySeg, xend = xEnd, yend = yEnd))


#Add an arrow
datAll[mType == "Linear effect assumption",xArrow:=58.5]
datAll[mType == "Linear effect assumption",xArrowEnd:= 58.5]

datAll[mType == "Linear effect assumption",yArrow:= -0.5 ]
datAll[mType == "Linear effect assumption",yArrowEnd:=0.05289430]

p_seg <- p_seg + geom_segment(aes(x = xArrow, y = yArrow , xend = xArrowEnd, yend = yArrowEnd),
                              arrow = arrow(length = unit(0.2, "cm")))
p_seg
#Add arc

datAll[mType == "Linear effect assumption",xArc1:=32]
datAll[mType == "Linear effect assumption",yArc1:=0]

datAll[mType == "Linear effect assumption",xArc2:=33]
datAll[mType == "Linear effect assumption",yArc2:=-0.21]
datAll[mType == "Linear effect assumption",curv:=-0.21]


p_seg <- p_seg + geom_curve(col=2,aes(x = xArc1, y = yArc1, xend = xArc2, yend = yArc2), curvature = 0.25)


#Add annotation
dat_text <- data.frame(
  label = c("beta", ""),
  x=c(33.5,33.5),y=c(-0.17,-0.17),
  mType   = c("Linear effect assumption", "Constant effect assumption")
)
p_seg2 <- p_seg + geom_text(
  data    = dat_text,
  mapping = aes(x = x, y = y, label = label),
  hjust   = -0.1,
  vjust   = -1, parse=T,
  size = 2.9
)

dat_text11 <- data.frame(
  label = c("1", ""),
  x=c(34.6,34.6),y=c(-0.22,-0.22),
  mType   = c("Linear effect assumption", "Constant effect assumption")
)
p_seg2 <- p_seg2 + geom_text(
  data    = dat_text11,
  mapping = aes(x = x, y = y, label = label),
  hjust   = -0.1,
  vjust   = -1, parse=T,
  size = 1.4
)

#Add annotations to other sections of the plot
dat_text2 <- data.frame(
  label = c("Strongest evidence\n        for effect", ""),
  x=c(42,42),y=c(-1.35,-1.35),
  mType   = c("Linear effect assumption", "Constant effect assumption")
)
p_seg3 <- p_seg2 + geom_text(
  data    = dat_text2,
  mapping = aes(x = x, y = y, label = label),
  hjust   = -0.1,
  vjust   = -1,
  size = 3.5
) 


#Last text
dat_text3 <- data.frame(
  label = c("Significant\n interval", ""),
  x=c(58,58),y=c(-0.75,-0.75),
  mType   = c("Linear effect assumption", "Constant effect assumption")
)
p_seg4 <- p_seg3 + geom_text(
  data    = dat_text3,
  mapping = aes(x = x, y = y, label = label),
  hjust   = -0.1,
  vjust   = -1,
  size = 3.3
) 

p_illustrate <- ggarrange(p_seg4,p_curves,nrow=2,align="v")
