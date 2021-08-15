#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at Cardiff University under the supervision of Richard Anney and Andrew Pocklington.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--twas", action="store", default=NA, type='character',
	help="File containing TWAS results [required]"),
make_option("--pos", action="store", default=NA, type='character',
	help="File containing TWAS feature locations [required]"),
make_option("--post_proc_prefix", action="store", default=NA, type='character',
	help="Prefix of files created by post process script [required]"),
make_option("--window", action="store", default=NA, type='numeric',
	help="Window for plot and defining independent associations [required]"),
make_option("--gene_loc", action="store", default=NA, type='character',
	help="File containing gene locations for genes that aren't in the TWAS [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

if(!file.exists(paste0(opt$post_proc_prefix,'.joint_included.dat'))){
  cat(paste0(opt$post_proc_prefix,'.joint_included.dat does not exist\n'))
  q()
}

suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

# Read in the data
twas<-data.frame(fread(opt$twas))

# Create WGT column to match pos file
twas$tmp<-gsub(".*/", "", twas$FILE)
for(i in 1:dim(twas)[1]){
	twas$tmp2[i]<-gsub(paste0('/',twas$tmp[i]), "", twas$FILE[i])
}
twas$tmp3<-gsub(".*/", "", twas$tmp2)
twas$WGT<-paste0(twas$tmp3,'/',twas$tmp)
twas$tmp<-NULL
twas$tmp2<-NULL
twas$tmp3<-NULL

# Update positions with pos data
twas$P0<-NULL
twas$P1<-NULL
pos<-data.frame(fread(opt$pos))
twas<-merge(twas,pos[c('WGT','P0','P1')],by='WGT')

# Retain the most significant verion of each gene ID
twas<-twas[order(twas$TWAS.P),]
twas<-twas[!duplicated(twas$ID),]

# Integrate genes from reference that are not included in the TWAS
ref<-data.frame(fread(opt$gene_loc))
ref_not_in_twas<-ref[!(ref$V4 %in% twas$ID),]
names(ref_not_in_twas)<-c('CHR_ref','P0_ref','P1_ref','ID')
twas<-merge(twas,ref_not_in_twas, by='ID', all=T) 
twas$P0[is.na(twas$P0)]<-twas$P0_ref[is.na(twas$P0)]
twas$P1[is.na(twas$P1)]<-twas$P1_ref[is.na(twas$P1)]
twas$CHR[is.na(twas$CHR)]<-twas$CHR_ref[is.na(twas$CHR)]
twas$CHR_ref<-NULL
twas$P0_ref<-NULL
twas$P1_ref<-NULL

# Read in post process files
joint<-data.frame(fread(paste0(opt$post_proc_prefix,'.joint_included.dat')),Joint=T)
n_line<-system(paste0('wc -l ',opt$post_proc_prefix,'.joint_dropped.dat'),intern=T)
n_line<-as.numeric(unlist(strsplit(n_line,' '))[1])
if(n_line > 1){
	dropped<-data.frame(fread(paste0(opt$post_proc_prefix,'.joint_dropped.dat')),Joint=F)
} else {
	dropped<-NULL
}

# Read in SNP p-values before and after conditioning on joint model
post_proc_file<-sub(".*/", '', opt$post_proc_prefix)
post_proc_dir<-sub(paste0('/',post_proc_file), '', opt$post_proc_prefix)
cond_list<-list.files(path=post_proc_dir, pattern=paste0(post_proc_file,'.loc_.*cond$'))
cond<-NULL
for(i in cond_list){
	temp<-data.frame(fread(paste0(post_proc_dir,'/',i)))
	cond<-rbind(cond,temp)
}

for(i in cond_list){
	# Read in the first chunk as determined by the post-process files
	cond_i<-data.frame(fread(paste0(post_proc_dir,'/',i)))
	
	# Skip if all NA (can happen in MHC region)
	if(sum(complete.cases(cond_i) == T) == 0){
	  next
	}
	  
	# Remove SNPs with missing values
	cond_i<-cond_i[complete.cases(cond_i),]
	
	# Identify the chromosome number
	chr_i<-twas$CHR[twas$ID == joint$ID[1]]
	
	# Extract TWAS results for features in this region
	twas_i<-twas[twas$CHR == chr_i & twas$P1 > min(cond_i$POS) & twas$P0 < max(cond_i$POS),]
	
	# Identify jointly significant genes within locus
	twas_i_joint<-twas_i[(twas_i$FILE %in% joint$FILE),]
	twas_i_joint<-twas_i_joint[order(twas_i_joint$P0),]
	
	#Expand window if there are more than 1 joint genes in locus
	if(dim(twas_i_joint)[1] > 1){
		opt$window<-(max(twas_i_joint$P1) - min(twas_i_joint$P0))+opt$window
	}
	
	# Extract TWAS features within window of jointly significant genes
	twas_i_win<-twas_i[twas_i$P0 < (max(twas_i_joint$P1) + opt$window) & twas_i$P1 > (min(twas_i_joint$P0) - opt$window),]
	twas_i_win<-twas_i_win[order(twas_i_win$P0),]
	
	# variable for right boundary to account for labels
	twas_i_win$P1_lab<-twas_i_win$P1 + (nchar(twas_i_win$ID)*1e5*((opt$window/1)/1e6))
	Right_boundary<-(max(twas_i_win$P1_lab)-opt$window/2)/1e6
	True_window<-Right_boundary*1e6-min(twas_i_joint$P0)

	# Assign genes 'lines' for gene locations to be plotted.
	for(j in 1:dim(twas_i_win)[1]){
		if(j == 1){
			twas_i_win$Line<-NA
			twas_i_win$Line[j] <- 1
		} else {
			for(k in unique(twas_i_win$Line[1:(j-1)])){
				if(twas_i_win$P0[j] < max(twas_i_win$P1_lab[which(twas_i_win$Line == k)])){
					if(k != max(unique(twas_i_win$Line[1:j-1]))){
					} else {
						twas_i_win$Line[j]<-max(unique(twas_i_win$Line[1:j-1]))+1
					}
				} else {
					twas_i_win$Line[j]<-k
					break()
				}
			}
		}
	}
	
	# Label genes as either jointly or marginally significant
	twas_i_win$Label[(twas_i_win$ID %in% twas$ID)]<-'NS'
	twas_i_win$Label[(twas_i_win$ID %in% ref_not_in_twas$ID)]<-NA
	twas_i_win$Label[(twas_i_win$ID %in% dropped$ID)]<-'Marginally'
	twas_i_win$Label[(twas_i_win$ID %in% joint$ID)]<-'Jointly'
	twas_i_win$Label<-factor(twas_i_win$Label, levels=c('NA','NS','Marginally','Jointly'))
	
	# Extract SNPs within opt$window of jointly significant genes and right boundary
	cond_i_win<-cond_i[which(cond_i$POS < (Right_boundary+((opt$window*1.1)/5e6))*1e6 & cond_i$POS > (min(twas_i_joint$P0) - (opt$window*1.3))) ,]
	cond_i_win_for_plot<-rbind(data.frame(cond_i_win[c('SNP','POS')],P=cond_i_win$GWAS.P,Type='Marginal'), data.frame(cond_i_win[c('SNP','POS')], P=cond_i_win$GWAS_cond.P,Type='Conditional'))
	
	# Italicise gene IDs
	twas_i_win$ID_ital<-paste0("italic('",twas_i_win$ID,"')")
	
	# Create plot showing SNP p-values before and after conditioning on the joint model.
	SNP_plot<-ggplot(cond_i_win_for_plot, aes(x=POS/1e6, y=-log10(P), colour=Type)) +
		geom_point() +
		xlim(c(((min(twas_i_joint$P0) - (opt$window*1.3))/1e6),Right_boundary+((opt$window*1.1)/5e6))) +
		labs(x=paste0("Position on Chromosome ",chr_i," (Mb)"), y="-log10(P-value)") +
	  theme_cowplot() +
		theme(legend.title = element_blank()) +
		scale_color_manual(values=c("#3333FF","#999999"))
	
	# Create plot showing gene locations
	Gene_plot<-	ggplot() +
			geom_rect(data=twas_i_win, mapping=aes(xmin=P0/1e6, xmax=P1/1e6, ymin=Line+0.1, ymax=Line+1-0.1, colour=Label, fill=Label)) +
			geom_text(data=twas_i_win, aes(x=P1/1e6, y=Line+0.5, label=ID_ital, colour=Label), hjust=-0.1, parse=T) +
			xlim(c(((min(twas_i_joint$P0) - (opt$window*1.3))/1e6),Right_boundary+((opt$window*1.1)/5e6))) +
			labs(x=paste0('Chromosome ',chr_i, " (Mb)")) +
	    theme_cowplot() +
			theme(axis.line=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.title = element_blank())
	
	# Combine the two plots
	out<-paste0(post_proc_dir,'/',sub('.cond','',i),'.OP.png')
	png(out, units='px', res=300, width=2500,height=1000+max(twas_i_win$Line)*50)
	print(plot_grid(Gene_plot,SNP_plot,rel_heights = c(max(twas_i_win$Line)/10,1), ncol=1, align='v', axis='l'))
	dev.off()
}

