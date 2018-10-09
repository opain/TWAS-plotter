# TWAS-plotter

This repository contains functions for plotting TWAS results:

TWAS-plotter is a tool that creates a Manhattan-style plot for TWAS results using ggplot2. It is designed to work with the output of [**FUSION.assoc_test.R** ](https://github.com/gusevlab/fusion_twas/blob/master/FUSION.assoc_test.R).

TWAS-locus-plotter is a tool that creates locus plots for chosen loci using ggplot2. The figures are very similar to those created using the post-process script but the the labels do not overlap, and genes that are available in the TWAS are highlighted. 

 Feel free to edit the plotting function to meet your needs.



## Getting started

### Prerequisites

* R and the required packages:

```R
install.packages(c('data.table','optparse','ggplot2','ggrepel','cowplot'))
```

* Perform TWAS using FUSION:

  - Instructions on how to perform a TWAS are available [here](http://gusevlab.org/projects/fusion/).
  - Per chromosome results files should combined without duplicating the header (see below).
  - You can also combine results from multiple TWAS.

```sh
# Combine per chromomsome TWAS results
head -1 ukbiobank-2017-1160-prePRS-fusion.tsv.chr1 > ukbiobank-2017-1160-prePRS-fusion.tsv.GW
tail -n +2 -q ukbiobank-2017-1160-prePRS-fusion.tsv.chr* >> ukbiobank-2017-1160-prePRS-fusion.tsv.GW
```



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



## Examples

```shell
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



## Help

This script was written by Dr Oliver Pain under the supervision of Dr Richard Anney whilst at the MRC Centre for Neuropsychiatric Genetics and Genomics, Cardiff University.

If you have any questions or comments use the [google group](https://groups.google.com/forum/#!forum/twas-related-r-scripts).







