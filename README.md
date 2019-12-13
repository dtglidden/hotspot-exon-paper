# Scripts used in the Hotspot Exon Paper

This repository holds all R and Python scripts used to generate the figures seen in this paper. There are 6 main subdirectories:
* data: Intermediate datasets are saved/loaded here.
* example\_data: Datasets for testing and peer review purposes.
* lib: Common R functions used by several scripts.
* maxent: A mirror of the [Perl MAXENT scripts from the Burge Lab](http://genes.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html). Slight modifications to input handling were made in order to call MAXENT functions from R.
* plots: Any saved plots are stored here.
* scripts: All R/Python scripts.Run scripts with --help for description of arguments

For testing and peer review purposes the following R scripts should be invoked by running "Rscript script.R --test". Running the scripts this way utilizes test datasets to plots that resemble the figures in the manuscript.
* scripts/calc\_mapsy\_muts.R
* scripts/feature\_pie.R
* scripts/conditional\_prob.R
* scripts/compare\_ss\_usage\_and\_hgmd\_splice.R
* scripts/plotEscGroups.R
* scripts/wt\_spleff\_correlations.R

## System requirements
Code requires R >= 3.5 and Python 3. Running time for each script should be no more than a few minutes. Most will complete in a matter of seconds. Code was tested on Debian Linux 3.16.7. Non-standard hardward is not required.

## Installation
Make sure R and Python are installed with the necessary packages for each script. Install time should only take a few minutes if R and Python have not been installed.


## Python code
splice\_site\_usage.py
	inputs: GENCODE GTF annotation, splice junction read count files from STAR-mapped RNA-seq data, sample name and desired output directory
	output: two files containing the inclusion reads, exclusion reads and splice site usage of 3' and 5' splice sites within the provided splice junction data

estimate\_hotspot\_prevalence.py
	inputs: GENCODE GTF annotation file, gene haploinsufficency score file, usage data, sample name and desired output directory
	outputs: a file containing potential hotspot exons called based on whether their splice site usage is in the bottom 10% for splice sites within genes with similar number of introns and HI scores

usage\_by\_introns.py
	inputs: GENCODE GTF annotation file, intron BED annotation file, usage data, sample name and desired output directory
	outputs: a file containing the each gene's average splice site usage and number of introns

usage\_by\_HI\_score.py
	inputs: GENCODE GTF annotation file, gene haploinsufficency score file, usage data, sample name and desired output directory
	outputs: a file containing the each gene's average splice site usage and HI score
