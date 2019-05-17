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
		help="Path to save output [required]"),
make_option("--width", action="store", default=2000, type='numeric',
		help="width of plot [optional]"),
make_option("--height", action="store", default=1250, type='numeric',
		help="height of plot [optional]"),
make_option("--res", action="store", default=300, type='numeric',
		help="height of plot [optional]")
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

	if(sum(d_sig$TWAS.Z > 0) > 0){
		d_sig_pos<-d_sig[d_sig$TWAS.Z > 0,]
	}
	
	if(sum(d_sig$TWAS.Z < 0) > 0){
		d_sig_neg<-d_sig[d_sig$TWAS.Z < 0,]
	}

	chr_labs<-as.character(unique(d$CHR))
	chr_labs[chr_labs == '19'| chr_labs == '21']<-' '

	if(dim(d_sig)[1] == 0){
		p<-ggplot(d,aes(x=pos,y=TWAS.Z,colour=factor(CHR))) +
			geom_point(size=0.5) +
			scale_x_continuous(name="Chromosome", breaks=ticks, labels=chr_labs) +
			scale_y_continuous(name='Z score',limits=c(-ylimit,ylimit)) +
			scale_colour_manual(values=mycols, guide=FALSE) +
			geom_hline(yintercept=0,colour="black") +
			geom_hline(yintercept=Sig_Z_Thresh,colour="blue") +
			geom_hline(yintercept=-Sig_Z_Thresh,colour="blue")
			
	} else {
	
	p<-ggplot(d,aes(x=pos,y=TWAS.Z,colour=factor(CHR))) +
		geom_point(size=0.5) +
		scale_x_continuous(name="Chromosome", breaks=ticks, labels=chr_labs) +
		scale_y_continuous(name='Z score',limits=c(-ylimit,ylimit)) +
		scale_colour_manual(values=mycols, guide=FALSE) +
		geom_hline(yintercept=0,colour="black") +
		geom_hline(yintercept=Sig_Z_Thresh,colour="blue") +
		geom_hline(yintercept=-Sig_Z_Thresh,colour="blue") +
		geom_point(data=d_sig, aes(x=pos,y=TWAS.Z), colour="red", fill='red', size=1.5)
		
		if(sum(d_sig$TWAS.Z > 0) > 0){
			p<-p+geom_text_repel(data=d_sig_pos, aes(x=pos,y=TWAS.Z, label=ID), colour='black', nudge_y=1, size=2.5, force=5, segment.alpha=0.25, ylim=c(Sig_Z_Thresh+0.1,NA))
		}
		
		if(sum(d_sig$TWAS.Z < 0) > 0){
			p<-p+geom_text_repel(data=d_sig_neg, aes(x=pos,y=TWAS.Z, label=ID), colour='black', nudge_y=-1, size=2.5, force=5, segment.alpha=0.25, ylim=c(NA,-Sig_Z_Thresh-0.1))
		}
	}
	
		p<-p+theme(	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
					axis.line = element_line(colour = "black"),
					axis.text.x = element_text(angle=45, size=8, hjust=1))
					
		p
}

# Read in the TWAS data
twas<-data.frame(fread(opt$twas))

# Make plot
if(!is.na(opt$sig_z)){
	png(paste0(opt$output,'.png'), unit='px', res=opt$res, width = opt$width, height = opt$height)
	print(TWAS_manhattan(dataframe=twas, Sig_Z_Thresh=opt$sig_z))
	dev.off()
}
if(!is.na(opt$sig_p)){
	png(paste0(opt$output,'.png'), unit='px', res=opt$res, width = opt$width, height = opt$height)
	print(TWAS_manhattan(dataframe=twas, Sig_Z_Thresh=qnorm(1 - opt$sig_p)))
	dev.off()
}
if(is.na(opt$sig_p) & is.na(opt$sig_z)){
	png(paste0(opt$output,'.png'), unit='px', res=opt$res, width = opt$width, height = opt$height)
	print(TWAS_manhattan(dataframe=twas))
	dev.off()
}
