library(data.table)
library(ggplot2)
library(survival)
library(dplyr)
library(ggpubr)
library(ggsci)
library(stringr)

source("/nfs/scistore13/robingrp/sojavee/ageSpecificMarginal/analyseLinearEffect/figures/Fig1_Illustration.R")
source("/nfs/scistore13/robingrp/sojavee/ageSpecificMarginal/analyseLinearEffect/figures/Fig1_Manhattan.R")

axis_set_mod <- as.data.table(axis_set)
axis_set_mod[,CHR_write:=as.character(CHR)]
axis_set_mod <- axis_set_mod[!(CHR %in% c(11,13,15,16,18,19,21,22)),]
axis_set_mod[CHR_write == "23",CHR_write:="X"] 

#X_all_sub <- X_all_sub[order(CHR2,decreasing=T),]
X_all_sub1 <- X_all_sub[CHR2 == "Novel discovery",]
X_all_sub2 <- X_all_sub[CHR2 != "Novel discovery",]
X_all_sub_ordered <- rbind(X_all_sub2,X_all_sub1)

my_colors <- c("#b5b4b1", "#636361","#BC3C29FF")
# add names so bar for 'c' gets fill, too
names(my_colors) <- c(" ","  ","Novel discovery")

X_all_sub_ordered[P==0,P:=10**(-249)]
X_all_sub_ordered[,log10P:=-log10(P)]

X_all_sub_ordered[log10P>249,log10P:=249]

p_manhattan <- ggplot(data=X_all_sub_ordered, aes(x=bpAdj,y=log10P,col=CHR2)) + geom_point(size=1) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
    theme(legend.title = element_blank(), legend.text = element_text(size=13), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
    theme(title = element_text(color="gray20", size=13)) +
    theme(axis.title = element_text(color="gray20", size=13)) +
    theme(axis.text = element_text(color="gray40",size=13)) +
    theme(legend.position = "top") +
	theme(strip.background = element_rect(fill="white")) +
	theme(strip.text=element_text(size = 13, color = "gray10")) +
    xlab("Chromosome") + ylab("-log10(p)") + geom_hline(col="coral",yintercept=-log10(5e-8),lty=3) +
    scale_color_manual(values = my_colors,limits=c("Novel discovery"," ") ) + #,breaks=c("Novel discovery","","") #,"col_odd","col_even"
    scale_x_continuous(label = axis_set_mod$CHR_write, breaks = axis_set_mod$center) + 
    ylim(c(0,250))#+ scale_y_continuous(trans='log2')

ggsave(p_manhattan,file="/nfs/scistore13/robingrp/sojavee/figures/Fig1_Illustrate_Manhattan.jpeg",height=7,width=6)


#3) Combine
p <- ggarrange(p_illustrate,p_manhattan,labels=c("a","b"),nrow=1,widths=c(1,1.1), legend="bottom") 

ggsave(p,file="/nfs/scistore13/robingrp/sojavee/figures/Fig1_Illustrate_Manhattan.jpeg",height=6,width=12)



#
my_colors <- c("#b5b4b1", "#636361")
# add names so bar for 'c' gets fill, too
names(my_colors) <- c(" ","  ")
p_manhattan <- ggplot(data=X_all_sub_ordered, aes(x=bpAdj,y=log10P,col=CHR2)) + geom_point(size=1) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
    theme(legend.title = element_blank(), legend.text = element_text(size=13), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
    theme(title = element_text(color="gray20", size=13)) +
    theme(axis.title = element_text(color="gray20", size=13)) +
    theme(axis.text = element_text(color="gray40",size=13)) +
    theme(legend.position = "top") +
	theme(strip.background = element_rect(fill="white")) +
	theme(strip.text=element_text(size = 13, color = "gray10")) +
    xlab("Chromosome") + ylab("-log10(p)") + #geom_hline(col="coral",yintercept=-log10(5e-8),lty=3) +
    scale_color_manual(values = my_colors) + #,breaks=c("Novel discovery","","") #,"col_odd","col_even"
    scale_x_continuous(label = axis_set_mod$CHR_write, breaks = axis_set_mod$center) + 
    ylim(c(0,250)) + guides(color = "none")
 #+ scale_y_continuous(trans='log2') 
ggsave(p_manhattan,file="/nfs/scistore13/robingrp/sojavee/figures/Fig1_Illustrate_Manhattan.jpeg",height=5,width=8)




