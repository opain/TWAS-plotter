#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at Cardiff University under the supervision of Richard Anney.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--twas", action="store", default=NA, type='character',
		help="Path to genome-wide TWAS results [required]"),
make_option("--sig_z", action="store", default=NA, type='numeric',
		help="Absolute Z score threshold for transcriptome-wide significance [optional]"),
make_option("--sig_p", action="store", default=NA, type='numeric',
		help="p-value threshold for transcriptome-wide significance [optional]"),
make_option("--output", action="store", default=NA, type='character',
		help="Path to save output [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load libraries
suppressMessages(library(data.table))
suppressMessages(library(ggrepel))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

# Create plotting function
TWAS_manhattan = function(dataframe, title=NULL, ylimit=max(abs(dataframe$TWAS.Z),na.rm=T)+1, Sig_Z_Thresh=qnorm(1-(0.05/length(dataframe$TWAS.Z))/2)) {

	d=dataframe[order(dataframe$CHR, dataframe$P0),]
	d=d[!is.na(d$TWAS.P),]

	d$pos=NA
	ticks=NULL
	lastbase=0
	numchroms=length(unique(d$CHR))

	for (i in unique(d$CHR)) {
		if (i==1) {
			d[d$CHR==i, ]$pos=d[d$CHR==i, ]$P0
	}	else {
			lastbase=lastbase+tail(subset(d,CHR==i-1)$P0, 1)
			d[d$CHR==i, ]$pos=d[d$CHR==i, ]$P0+lastbase
			}
	ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
	}

	ticklim=c(min(d$pos),max(d$pos))

	mycols=rep(c("gray35","gray72"),60)

	d$Sig_Z_Thresh<-Sig_Z_Thresh

	d_sig<-d[which(abs(d$TWAS.Z) > d$Sig_Z_Thresh),]
	d_sig<-d_sig[rev(order(abs(d_sig$TWAS.Z))),]
	d_sig<-d_sig[!duplicated(d_sig$ID),]

	ggplot(d,aes(x=pos,y=TWAS.Z,colour=factor(CHR))) +
		geom_point(size=0.5) +
		scale_x_continuous(name="Chromosome", breaks=ticks, labels=(unique(d$CHR))) +
		scale_y_continuous(name='Z score',limits=c(-ylimit,ylimit)) +
		scale_colour_manual(values=mycols, guide=FALSE) +
		geom_hline(yintercept=0,colour="black") +
		geom_hline(yintercept=Sig_Z_Thresh,colour="blue") +
		geom_hline(yintercept=-Sig_Z_Thresh,colour="blue") +
		geom_point(data=d_sig, aes(x=pos,y=TWAS.Z), colour="red", fill='red', size=1.5) +
		geom_text_repel(data=d_sig, aes(x=pos,y=TWAS.Z, label=ID), colour='black', size=3) +
		theme(	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
				axis.line = element_line(colour = "black"),
				axis.text.x = element_text(angle=45, size=8))
}

# Read in the TWAS data
twas<-data.frame(fread(opt$twas))

# Make plot
if(!is.na(opt$sig_z)){
	png(paste0(opt$output,'.png'), unit='px', res=300, width = 2000, height = 1250)
	print(TWAS_manhattan(dataframe=twas, Sig_Z_Thresh=opt$sig_z))
	dev.off()
}
if(!is.na(opt$sig_p)){
	png(paste0(opt$output,'.png'), unit='px', res=300, width = 2000, height = 1250)
	print(TWAS_manhattan(dataframe=twas, Sig_Z_Thresh=qnorm(1 - opt$sig_p)))
	dev.off()
}
if(is.na(opt$sig_p) & is.na(opt$sig_z)){
	png(paste0(opt$output,'.png'), unit='px', res=300, width = 2000, height = 1250)
	TWAS_manhattan(dataframe=twas)
	dev.off()
}
