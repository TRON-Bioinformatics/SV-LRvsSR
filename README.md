# SV-LRvsSR

A study to compare structural variation (SV) predictions from 10X Genomics linked-reads sequencing (10XWGS) and conventional Illumina short-reads sequencing (cWGS). We predict SVs from Long Ranger pipeline of 10X Genomics, NAIBR and GROC-SV for linked-reads. And we predict SV from Delly, Lumpy and SvABA and combine the predictions with linked-read tools. The MCF7 breast cancer cell line and a primary breast tumor sample was analysed for the comparison. Apart from comparing different predictions by two technologies we also propose precision enrichment logistic regression models for retaining true positive SVs amongst the high false positive calls from both the technologies.

The code for analysis was written in python (version 2.7.15) and precision enrichment models were trained and built in R (version 3.6).

# Tools and versions
Delly (v 0.7.6)
Lumpy (v 0.2.13)
SvABA (v 134)
Long Ranger (v 2.2.2)
NAIBR (v 1.0)
GROC-SV (v 0.2.5)

# Requirements
- python 2.7.15
	- pandas 0.23.4
	- numpy 1.15.4
	- pysam 0.15.4
- samtools 1.9
- bwa 0.7.17
- sambamba 0.6.8
- samblaster 0.1.24
- bcftools 1.9
- R 3.6.0
	- dplyr
	- optparse
	- tidyverse
	- caret
	- plyr

# Install python modules
```
conda install -c conda-forge pandas=0.23.4
conda install -c bioconda pysam=0.15.2 numpy=1.15.4 samtools=1.9 sambamba=0.6.8 bwa=0.7.17 bcftools=1.9 samblaster=0.1.24
```
# Install R packages
```
install.packages(c("dplyr","tidyverse","caret","plyr","optparse"))
```
# Usage

## Pipeline to predict SVs from short-reads sequencing and linked-reads sequencing tools are in Commands.org

## Integrate SV calls from short-reads sequencing tools
python Integrate.py \
-chromLengths `<length of chromosomes eg: /data/chromInfo_hg38.txt>` \
-delly `<comma separated vcf files generated from Delly in following order: Deletion.vcf,Duplication.vcf,Inversion.vcf,Translocation.vcf>` \
-lumpy `<vcf file with SV calls from Lumpy>` \
-svaba `<vcf file with SV calls from SvABA>` \
-frag `<Window size of <INT> bases for integrating overlapping SV calls from tools. Default=500 bases>` \
-outdir `<Output directory path>` \
-mode `<germline or somatic. Default=germline>` \
-finalFile `<Output file name>`\

This generates a tab delimited file with SV calls from all three tools. Example file: /example/Combined_SR.tsv

## Integrate SV calls from linked-reads sequencing tools and also combine it with short-reads calls

python Combine_SR_LR.py \
-SR `<tab delimited file containing integrated SV calls from short-reads. example: /example/Combined_SR.tsv>` \
-LR1 `<vcf file from Long Ranger containing small and mid sized deletions>` \
-LR2 `<vcf file from Long Ranger containing large sized SV>` \
-NAIBR `<bed file containing SV calls from NAIBR>` \
-GROCSV `<vcf file containing SV calls from GROC-SV>` \
-w `<Window size of <INT> bass for integrating overlapping SV calls. Default=500 bases>` \
-outdir `<Output directory path>`\

It generates a tab delimited file containing all overlapping and non overlapping SV calls from short-reads and linked-reads tools. Example file: /example/Combined_SR_LR.tsv

## Perform requantification and add features to each SV call

python Requantification.py \
-inputFile `<tab delimited file containing all SV calls. example: /example/Combined_SR_LR.tsv>` \
-out `<Ouput file name>` \
-n `<Number of cores to use. Default=1>` \
-refBit `<.2bit file of th reference genome>` \
-area `<Area in <INT> bases around the breakpoints for the rquantification genomic template. Default=500 bases>` \
-Read1_SR `<Paired-end read 1 from Illumina short-reads sequencing>` \
-Read2_SR `<Paired-end read 2 from Illumina short-reads sequencing>` \
-Read1_LR `<Paired-end read 1 from 10X Genomics linked-reads sequencing with barcodes trimmed>` \
-Read2_LR `<Paired-end read 2 from 10X Genomics linked-reads sequencing with barcodes trimmed>` \
-outdir `<Output directory path>` \
-tmpdir `<Temportary directory path>` \
-cutoff `<a read after re-alignment to genomic template has to overlap breakpoint with atleast <INT> bases. Deafult=10 bases>` \
-lengths `<length of chromosomes eg: /data/chromInfo_hg38.txt>` \

It generates a tab delimited SV file with requantification related features like junction reads and spanning pairs calculated from short-read and/or linked-reads. Example: /example/Combined_SR_LR_requant.tsv

## Include GEM counts for all SVs from Lariat linked-reads generated alignment file

python GEM_count.py \
-inputFile `<tab delimited file containing SV calls. Example: /example/Combined_SR_LR_requant.tsv>`\
-bam `<Lariat generated bam file with linked-reads>` \
-n `<Number of cores to be used. Default=1>` \
-outdir `<Output directory path>` \
-output `<Output file name>` \
-w `<Window in <INT> bases to calculate supporting GEM for SV call. Default=1000 bases>` \

This generates tab delimited file containing suppporting GEM counts for each SV. Example: /example/Combined_SR_LR_requant_GEM.tsv

## Include local coverage around breakpoints and annotate breakpoints with repetitive and poor mappability regions

python Include_annotations.py \
-Annotate `<Annotate breakpoints with repetitive and poor mappability regions. Default=Yes>` \
-repeatMasker `<Repeat masker file. example: /data/RepeatMasker_hg38.bed>` \
-mappability `<Mappability track file. example: /data/k100.umap.bed>` \
-File `<tab delimited file containing SV calls. example: /example/Combined_SR_LR_requant_GEM.tsv>` \
-BAM `<short-reads BAM alignment file>` \
-BAMLinked `<linked-reads BAM alignment file>` \
-Threads `<Number of cores. Default=1>` \
-outdir `<Output directory path>`\

This generates tab delimited file with local coverage around breakpoints calculated from short-reads and/or linked-reads alignment file. If "Annotate" argument is Yes, then each breakpoint is annotated with repetitive and porr mappability regions. Example file: /example/Combined_SR_LR_requant_GEM_annotated.tsv