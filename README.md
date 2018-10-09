# TWAS-plotter

This repository contains functions for plotting TWAS results:

TWAS-plotter is a tool that creates a Manhattan-style plot for TWAS results using ggplot2. It is designed to work with the output of [**FUSION.assoc_test.R** ](https://github.com/gusevlab/fusion_twas/blob/master/FUSION.assoc_test.R).

TWAS-locus-plotter is a tool that creates locus plots for chosen loci using ggplot2. The figures are very similar to those created using the post-process script but the the labels do not overlap, and genes that are available in the TWAS are highlighted. 

 Feel free to edit the plotting function to meet your needs.



## Prerequisites

* R and the required packages:

```R
install.packages(c('data.table','optparse','ggplot2','ggrepel','cowplot'))
```

* Perform TWAS using FUSION:

  - Instructions on how to perform a TWAS are available [here](http://gusevlab.org/projects/fusion/).

  - For example:

```sh
for chr in $(seq 1 22); do
    qsub -b y Rscript ~/Software/fusion_twas-master/FUSION.assoc_test.R \
        --sumstats ~/Data/GWAS_sumstats/ukbiobank-2017-1160-prePRS-fusion.tsv \
        --weights ~/Data/FUSION-models/CMC.BRAIN.RNASEQ.pos \
        --weights_dir ~/Data/FUSION-models \
        --ref_ld_chr ~/Software/fusion_twas-master/LDREF/1000G.EUR. \
        --out ~/Analyses/TWAS/ukbiobank-2017-1160-prePRS-fusion.tsv.chr${chr} \
        --chr ${chr}
done
```

  - Per chromosome results files should combined without duplicating the header.
      - For example:

```sh
# Combine per chromomsome TWAS results
head -1 ~/Analyses/TWAS/ukbiobank-2017-1160-prePRS-fusion.tsv.chr1 > ~/Analyses/TWAS/ukbiobank-2017-1160-prePRS-fusion.tsv.GW
tail -n +2 -q ~/Analyses/TWAS/ukbiobank-2017-1160-prePRS-fusion.tsv.chr* >> ~/Analyses/TWAS/ukbiobank-2017-1160-prePRS-fusion.tsv.GW
```

* Perform post-TWAS analysis using FUSION:

  * Instructions on how to perform a post-TWAS are available [here](http://gusevlab.org/projects/fusion/).

  * For example:

```sh
# Extract significant features
x=0.000001
cat ~/Analyses/TWAS/ukbiobank-2017-1160-prePRS-fusion.tsv.GW | awk -v var="${x}" 'NR == 1 || $19 < var' > ~/Analyses/TWAS/ukbiobank-2017-1160-prePRS-fusion.tsv.GW.Sig

# Change directory to where glist-hg19 is saved
cd ~/Data

# Run post-process script for chromosomes containing significant genes 
for chr in $(seq 1 22); do
    status=$(awk -v var="${chr}" '$3 == var {print "Present";exit;}' ~/Analyses/TWAS/ukbiobank-2017-1160-prePRS-fusion.tsv.GW.Sig)
    if [ "$status" == "Present" ]; then
        qsub -b y -cwd Rscript ~/Software/fusion_twas-master/FUSION.post_process.R \
            --input ~/Analyses/TWAS/ukbiobank-2017-1160-prePRS-fusion.tsv.GW.Sig \
            --sumstats ~/Data/GWAS_sumstats/ukbiobank-2017-1160-prePRS-fusion.tsv \
            --ref_ld_chr ~/Software/fusion_twas-master/LDREF/1000G.EUR. \
            --out ~/Analyses/TWAS/ukbiobank-2017-1160-prePRS-fusion.tsv.GW.Sig.PostProc.${chr} \
            --chr ${chr} \
            --save_loci \
            --plot \
            --locus_win 500000
    fi
done
```



## TWAS-plotter.V1.0.r

Creates Manhattan-style plot of TWAS results.

![](![Alt text](FUSION-locus-example.jpg?raw=true "Optional Title"))


### Parameters

| Flag     | Description                                                  | Default |
| :------- | ------------------------------------------------------------ | :-----: |
| --twas   | Path to genome-wide TWAS results                             |   NA    |
| --sig_z  | Absolute Z score threshold for transcriptome-wide significance |   NA    |
| --sig_p  | p-value threshold for transcriptome-wide significance        |   NA    |
| --output | Path to save output                                          |   NA    |

If neither sig_z or sig_p or specified, a Bonferroni correction for the number of features will be used. If there are multiple features with the same ID achieving significance, only the most significant will be labeled. 

### Output files

A .png file will be created.

### Examples

```sh
# Default Bonferroni corrected significance threshold.
Rscript TWAS-plotter.V1.0.r \
	--twas ukbiobank-2017-1160-prePRS-fusion.tsv.GW \
	--output Test
	
# Specify significance threshold as a Z-score.
Rscript TWAS-plotter.V1.0.r \
	--twas ukbiobank-2017-1160-prePRS-fusion.tsv.GW \
	--sig_z 5 \
	--output Test
	
# Specify significance threshold as a p-value.
Rscript TWAS-plotter.V1.0.r \
	--twas ukbiobank-2017-1160-prePRS-fusion.tsv.GW \
	--sig_p 1e-6 \
	--output Test
```



## TWAS-locus-plotter.V1.0.r

Creates locus plots for TWAS loci. Similar to plots created by FUSION script but it uses ggplot so is more adaptable (for me) and it highlights which genes were available in the TWAS.

###### TWAS-locus-plotter example:

![](C:\Users\mpmop\Desktop\TWAS-locus-example.png)

###### FUSION equivalent:

![](C:\Users\mpmop\Desktop\PDFtoJPG.me-1 (8).jpg)

### Parameters

| Flag               | Description                                                  | Default |
| :----------------- | ------------------------------------------------------------ | :-----: |
| --twas             | File containing TWAS results                                 |   NA    |
| --pos              | File containing TWAS feature locations                       |   NA    |
| --post_proc_prefix | Prefix of files created by post process script               |   NA    |
| --window           | Window for plot and defining independent associations        |   NA    |
| --gene_loc         | File containing gene locations for genes that aren't in the TWAS |         |

### Output files

<post_proc_prefix>.OP.png will be created.

### Examples

```R
# Run TWAS-locus-plotter for all chromomsomes with significant features
for chr in $(seq 1 22); do
    status=$(awk -v var="${chr}" '$3 == var {print "Present";exit;}' ~/Analyses/TWAS/ukbiobank-2017-1160-prePRS-fusion.tsv.GW.Sig)
    if [ "$status" == "Present" ]; then
        qsub -b y Rscript TWAS-locus-plotter.V1.0.r \
            --twas ~/Analyses/TWAS/ukbiobank-2017-1160-prePRS-fusion.tsv.GW.Sig \
            --pos ~/Data/FUSION-models/CMC.BRAIN.RNASEQ.pos \
            --post_proc_prefix ~/Analyses/TWAS/ukbiobank-2017-1160-prePRS-fusion.tsv.GW.Sig.PostProc.${chr} \
            --window 0.5e6 \
            --gene_loc ~/Data/glist-hg19
    fi
done

```



## Help

This script was written by Dr Oliver Pain under the supervision of Dr Richard Anney whilst at the MRC Centre for Neuropsychiatric Genetics and Genomics, Cardiff University.

If you have any questions or comments use the [google group](https://groups.google.com/forum/#!forum/twas-related-r-scripts).







