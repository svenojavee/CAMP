library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)


X <- NULL
clumpDat <- NULL

for(chr in 1:23){
    print(chr)
    X <- rbind(X,fread(paste0("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/Menopause_LinearEvaluated/Menopause_df2_atMax_45_52_",chr,".res")))
}


clumpDat <- NULL
for(chr in 1:23){
    dat <- fread(paste0("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/Menopause_atMax_",chr,"_clumping_p1_df2.clumped"))
    splitsNeeded <- max(dat$TOTAL)+1

    datSplit <- str_split_fixed(dat$SP2, pattern=",",n=splitsNeeded)
    datSplit[datSplit == "" | datSplit == "NONE"] <- NA
    datSplit <- as.data.table(datSplit)
    datSplit[,indexSNP:=dat$SNP]
    
    datLong <- melt(datSplit,id.vars="indexSNP")
    datLong[,value:=gsub("(1)","",value,fixed=T)]
    datLong[,variable:=NULL]

    indexSNPData <- copy(datLong)
    indexSNPData[,value:=indexSNP]
    indexSNPData <- unique(indexSNPData)

    #Add the index SNPs to the main file
    datLong <- rbind(datLong,indexSNPData)
    datLong <- datLong[!is.na(value),]

    datLong <- merge(datLong,dat[,list(SNP,CHR,BP,P)],by.x="indexSNP",by.y="SNP",all.x=T)

    clumpDat <- rbind(clumpDat,datLong)
}


#Read the novel/previously found results + clumps
replDat <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/ReplicatedNovelAssociationsatMax_45_52.csv")
prevDat <- fread("/mnt/beegfs/robingrp/sojavee/MenopauseMarginal/summaryAtAgeMenopause/outputs/PreviousAssociationsatMax_45_52.csv")

replDat <- replDat[,list(SNP)]
prevDat <- prevDat[,list(Variant)]
names(prevDat) <- "SNP"
replDat[,specType:="Novel discovery"]
prevDat[,specType:="Previous discovery"]

mergDat <- rbind(replDat,prevDat)
mergDat[specType == "Previous discovery",specType:=NA]

clumpDat2 <- merge(clumpDat,mergDat,by.x="indexSNP",by.y="SNP",all.x=T)
clumpDat3 <- clumpDat2[!is.na(specType),list(value,specType)]

X <- merge(X,clumpDat3,by.x="SNP",by.y="value",all.x=T)

X_all_sub <- X[P<0.001,]
X_all_sub <- X_all_sub[order(CHR,BP),]


X_all_sub[,max_bp:=max(BP),by="CHR"]
toAdd <- unique(X_all_sub[,list(CHR,max_bp)])
toAdd[,toAdd:=lag(max_bp)]
toAdd[is.na(toAdd),toAdd:=0]
toAdd[,cumToAdd:=cumsum(as.numeric(toAdd))]

X_all_sub <- merge(X_all_sub,toAdd[,list(CHR,cumToAdd)],by="CHR",all.x=T)
X_all_sub[,bpAdj:=BP+cumToAdd]

axis_set <- X_all_sub %>% 
  group_by(CHR) %>% 
  summarize(center = mean(bpAdj))

X_all_sub[,CHR2:=""]
X_all_sub[!is.na(specType),CHR2:=specType]


p <- ggplot(data=X_all_sub, aes(x=bpAdj,y=-log10(P),col=CHR2)) + geom_point(size=0.2) + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
    theme(legend.title = element_blank(), legend.text = element_text(size=13), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
    theme(title = element_text(color="gray20", size=13)) +
    theme(axis.title = element_text(color="gray20", size=13)) +
    theme(axis.text = element_text(color="gray40",size=13)) +
    theme(legend.position = "top") +
	theme(strip.background = element_rect(fill="white")) +
	theme(strip.text=element_text(size = 13, color = "gray10")) +
    xlab("Chromosome") + ylab("-log10(p)") + geom_hline(yintercept=-log10(5e-8),lty=3) +
    #scale_color_manual(values = c(rep(c("#276FBF", "#183059"), unique(length(axis_set$CHR))), "red","green") ) +
    scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) + scale_y_continuous(trans='log2')

ggsave(p,filename="/nfs/scistore13/robingrp/sojavee/figures/MP_atMax_8.4M.pdf",width=10,height=6)
ggsave(p,filename="/nfs/scistore13/robingrp/sojavee/figures/MP_atMax_8.4M.jpeg",width=10,height=6)
fwrite(X_all_sub,file="/nfs/scistore13/robingrp/sojavee/figData/MP_atMax_8.4M.txt")
fwrite(axis_set,file="/nfs/scistore13/robingrp/sojavee/figData/MP_atMax_8.4M_help.txt")

if(F){

########
X[is.na(specType),specType:=""]
library(lattice)
manhattan.plot<-function(chr, pos, pvalue, 
	sig.level=NA, annotate=NULL, ann.default=list(),
	should.thin=T, thin.pos.places=2, thin.logp.places=2, 
	xlab="Chromosome", ylab=expression(-log[10](p-value)),
	col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {

	if (length(chr)==0) stop("chromosome vector is empty")
	if (length(pos)==0) stop("position vector is empty")
	if (length(pvalue)==0) stop("pvalue vector is empty")

	#make sure we have an ordered factor
	if(!is.ordered(chr)) {
		chr <- ordered(chr)
	} else {
		chr <- chr[,drop=T]
	}

	#make sure positions are in kbp
	if (any(pos>1e6)) pos<-pos/1e6;

	#calculate absolute genomic position
	#from relative chromosomal positions
	posmin <- tapply(pos,chr, min);
	posmax <- tapply(pos,chr, max);
	posshift <- head(c(0,cumsum(posmax)),-1);
	names(posshift) <- levels(chr)
	genpos <- pos + posshift[chr];
	getGenPos<-function(cchr, cpos) {
		p<-posshift[as.character(cchr)]+cpos
		return(p)
	}

	#parse annotations
	grp <- NULL
	ann.settings <- list()
	label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
		col=NULL, fontface=NULL, fontsize=NULL, show=F)
	parse.label<-function(rawval, groupname) {
		r<-list(text=groupname)
		if(is.logical(rawval)) {
			if(!rawval) {r$show <- F}
		} else if (is.character(rawval) || is.expression(rawval)) {
			if(nchar(rawval)>=1) {
				r$text <- rawval
			}
		} else if (is.list(rawval)) {
			r <- modifyList(r, rawval)
		}
		return(r)
	}

	if(!is.null(annotate)) {
		if (is.list(annotate)) {
			grp <- annotate[[1]]
		} else {
			grp <- annotate
		} 
		if (!is.factor(grp)) {
			grp <- factor(grp)
		}
	} else {
		grp <- factor(rep(1, times=length(pvalue)))
	}
  
	ann.settings<-vector("list", length(levels(grp)))
	ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)

	if (length(ann.settings)>1) { 
		lcols<-trellis.par.get("superpose.symbol")$col 
		lfills<-trellis.par.get("superpose.symbol")$fill
		for(i in 2:length(levels(grp))) {
			ann.settings[[i]]<-list(pch=pch, 
				col=lcols[(i-2) %% length(lcols) +1 ], 
				fill=lfills[(i-2) %% length(lfills) +1 ], 
				cex=cex, label=label.default);
			ann.settings[[i]]$label$show <- T
		}
		names(ann.settings)<-levels(grp)
	}
	for(i in 1:length(ann.settings)) {
		if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
		ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
			parse.label(ann.settings[[i]]$label, levels(grp)[i]))
	}
	if(is.list(annotate) && length(annotate)>1) {
		user.cols <- 2:length(annotate)
		ann.cols <- c()
		if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
			ann.cols<-match(names(annotate)[-1], names(ann.settings))
		} else {
			ann.cols<-user.cols-1
		}
		for(i in seq_along(user.cols)) {
			if(!is.null(annotate[[user.cols[i]]]$label)) {
				annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
					levels(grp)[ann.cols[i]])
			}
			ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
				annotate[[user.cols[i]]])
		}
	}
 	rm(annotate)

	#reduce number of points plotted
	if(should.thin) {
		thinned <- unique(data.frame(
			logp=round(-log10(pvalue),thin.logp.places), 
			pos=round(genpos,thin.pos.places), 
			chr=chr,
			grp=grp)
		)
		logp <- thinned$logp
		genpos <- thinned$pos
		chr <- thinned$chr
		grp <- thinned$grp
		rm(thinned)
	} else {
		logp <- -log10(pvalue)
	}
	rm(pos, pvalue)
	gc()

	#custom axis to print chromosome names
	axis.chr <- function(side,...) {
		if(side=="bottom") {
			panel.axis(side=side, outside=T,
				at=((posmax+posmin)/2+posshift),
				labels=levels(chr), 
				ticks=F, rot=0,
				check.overlap=F
			)
		} else if (side=="top" || side=="right") {
			panel.axis(side=side, draw.labels=F, ticks=F);
		}
		else {
			axis.default(side=side,...);
		}
	 }

	#make sure the y-lim covers the range (plus a bit more to look nice)
	prepanel.chr<-function(x,y,...) { 
		A<-list();
		maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+.5;
		A$ylim=c(0,maxy);
		A;
	}

	xyplot(logp~genpos, chr=chr, groups=grp,
		axis=axis.chr, ann.settings=ann.settings, 
		prepanel=prepanel.chr, scales=list(axs="i"),
		panel=function(x, y, ..., getgenpos) {
			if(!is.na(sig.level)) {
				#add significance line (if requested)
				panel.abline(h=-log10(sig.level), lty=2);
			}
			panel.superpose(x, y, ..., getgenpos=getgenpos);
			if(!is.null(panel.extra)) {
				panel.extra(x,y, getgenpos, ...)
			}
		},
		panel.groups = function(x,y,..., subscripts, group.number) {
			A<-list(...)
			#allow for different annotation settings
			gs <- ann.settings[[group.number]]
			A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
			A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
			A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
			A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
			A$x <- x
			A$y <- y
			do.call("panel.xyplot", A)
			#draw labels (if requested)
			if(gs$label$show) {
				gt<-gs$label
				names(gt)[which(names(gt)=="text")]<-"labels"
				gt$show<-NULL
				if(is.character(gt$x) | is.character(gt$y)) {
					peak = which.max(y)
					center = mean(range(x))
					if (is.character(gt$x)) {
						if(gt$x=="peak") {gt$x<-x[peak]}
						if(gt$x=="center") {gt$x<-center}
					}
					if (is.character(gt$y)) {
						if(gt$y=="peak") {gt$y<-y[peak]}
					}
				}
				if(is.list(gt$x)) {
					gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
				}
				do.call("panel.text", gt)
			}
		},
		xlab=xlab, ylab=ylab, 
		panel.extra=panel.extra, getgenpos=getGenPos, ...
	);
}
pdf("/nfs/scistore13/robingrp/sojavee/figures/MP_atMax_8.4M.pdf",width=10,height=6)
manhattan.plot(X_all_sub$CHR, X_all_sub$BP, X_all_sub$P, annotate=X_all_sub$specType)
dev.off()


####
library(qqman)

manhattan2 <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray30", 
    "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05), 
    genomewideline = -log10(5e-08), highlight = NULL, highlight2 = NULL, highlight3 = NULL, logp = TRUE, 
    annotatePval = NULL, annotateTop = TRUE, ...) 
{
    CHR = BP = P = index = NULL
    if (!(chr %in% names(x))) 
        stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(x))) 
        stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(x))) 
        stop(paste("Column", p, "not found!"))
    if (!(snp %in% names(x))) 
        warning(paste("No SNP column found. OK unless you're trying to highlight."))
    if (!is.numeric(x[[chr]])) 
        stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(x[[bp]])) 
        stop(paste(bp, "column should be numeric."))
    if (!is.numeric(x[[p]])) 
        stop(paste(p, "column should be numeric."))
    if (!is.null(x[[snp]])) 
        d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
            pos = NA, index = NA, SNP = x[[snp]], stringsAsFactors = FALSE)
    else d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
        pos = NA, index = NA)
    d <- d[order(d$CHR, d$BP), ]
    if (logp) {
        d$logp <- -log10(d$P)
    }
    else {
        d$logp <- d$P
    }
    d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP, 
        d$CHR, length))
    nchr = length(unique(d$CHR))
    if (nchr == 1) {
        d$pos = d$BP
        xlabel = paste("Chromosome", unique(d$CHR), "position")
    }
    else {
        lastbase = 0
        ticks = NULL
        for (i in unique(d$index)) {
            if (i == 1) {
                d[d$index == i, ]$pos = d[d$index == i, ]$BP
            }
            else {
                lastbase = lastbase + max(d[d$index == (i - 1), 
                  "BP"])
                d[d$index == i, "BP"] = d[d$index == i, "BP"] - 
                  min(d[d$index == i, "BP"]) + 1
                d[d$index == i, "pos"] = d[d$index == i, "BP"] + 
                  lastbase
            }
        }
        ticks <- tapply(d$pos, d$index, quantile, probs = 0.5)
        xlabel = "Chromosome"
        labs <- unique(d$CHR)
    }
    xmax = ceiling(max(d$pos) * 1.03)
    xmin = floor(max(d$pos) * -0.03)
    def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
        las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
            ceiling(max(d$logp)+3)), xlab = xlabel, ylab = expression(-log[10](italic(p))))
    dotargs <- list(...)
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
        names(dotargs)]))
    if (!is.null(chrlabs)) {
        if (is.character(chrlabs)) {
            if (length(chrlabs) == length(labs)) {
                labs <- chrlabs
            }
            else {
                warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
            }
        }
        else {
            warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
        }
    }
    if (nchr == 1) {
        axis(1, ...)
    }
    else {
        axis(1, at = ticks, labels = labs, ...)
    }
    col = rep_len(col, max(d$index))
    if (nchr == 1) {
        with(d, points(pos, logp, pch = 20, col = col[1], ...))
    }
    else {
        icol = 1
        for (i in unique(d$index)) {
            points(d[d$index == i, "pos"], d[d$index == i, "logp"], 
                col = col[icol], pch = 20, ...)
            icol = icol + 1
        }
    }
    if (suggestiveline) 
        abline(h = suggestiveline, col = "blue")
    if (genomewideline) 
        abline(h = genomewideline, col = "coral",lty=2)
    if (!is.null(highlight)) {
        if (any(!(highlight %in% d$SNP))) 
            warning("You're trying to highlight SNPs that don't exist in your results.")
        d.highlight = d[which(d$SNP %in% highlight), ]
        with(d.highlight, points(pos, logp, col = "deepskyblue2", pch = 20, 
            ...))
    }
	if (!is.null(highlight2)) {
        if (any(!(highlight2 %in% d$SNP))) 
            warning("You're trying to highlight2 SNPs that don't exist in your results.")
        d.highlight2 = d[which(d$SNP %in% highlight2), ]
        with(d.highlight2, points(pos, logp, col = "firebrick1", pch = 20, 
            ...))
    }
	if (!is.null(highlight3)) {
        if (any(!(highlight3 %in% d$SNP))) 
            warning("You're trying to highlight3 SNPs that don't exist in your results.")
        d.highlight3 = d[which(d$SNP %in% highlight3), ]
        with(d.highlight3, points(pos, logp, col = "darkorchid1", pch = 20, 
            ...))
    }
    if (!is.null(annotatePval)) {
        if (logp) {
            topHits = subset(d, P <= annotatePval)
        }
        else topHits = subset(d, P >= annotatePval)
        par(xpd = TRUE)
        if (annotateTop == FALSE) {
            if (logp) {
                with(subset(d, P <= annotatePval), textxy(pos, 
                  -log10(P), offset = 0.625, labs = topHits$SNP, 
                  cex = 0.45), ...)
            }
            else with(subset(d, P >= annotatePval), textxy(pos, 
                P, offset = 0.625, labs = topHits$SNP, cex = 0.45), 
                ...)
        }
        else {
            topHits <- topHits[order(topHits$P), ]
            topSNPs <- NULL
            for (i in unique(topHits$CHR)) {
                chrSNPs <- topHits[topHits$CHR == i, ]
                topSNPs <- rbind(topSNPs, chrSNPs[1, ])
            }
            if (logp) {
                textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
                  labs = topSNPs$SNP, cex = 0.5, ...)
            }
            else textxy(topSNPs$pos, topSNPs$P, offset = 0.625, 
                labs = topSNPs$SNP, cex = 0.5, ...)
        }
    }
    par(xpd = FALSE)
}

snpsOfInterest <- X_all_sub$SNP[!is.na(X_all_sub$specType)]
X_all_sub[P==0,P:=6.645652e-270]
X_all_sub[P<1e-200,P:=1e-200]

pdf("/nfs/scistore13/robingrp/sojavee/figures/MP_atMax_8.4M_qqman.pdf",width=10,height=6)
manhattan2(X_all_sub,highlight=snpsOfInterest,suggestiveline=F)
dev.off()

jpeg("/nfs/scistore13/robingrp/sojavee/figures/MP_atMax_8.4M_qqplot.png",width=6,height=6)
qq(X$P)
dev.off()
#

X_qq <- X[,list(SNP,P)]
fwrite(X_qq,file="/nfs/scistore13/robingrp/sojavee/figData/qqPlotDat.txt")

}