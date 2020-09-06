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

1. Predict SVs using short-reads sequencing and linked-reads sequencing tools.\
Refer to Commands.org

2. Integrate SV calls from short-reads sequencing tools

python Integrate.py \
-chromLengths `<length of chromosomes eg: /data/chromInfo_hg38.txt>` \
-delly `<comma separated vcf files generated from Delly in following order: Deletion.vcf,Duplication.vcf,Inversion.vcf,Translocation.vcf>` \
-lumpy `<vcf file with SV calls from Lumpy>` \
-svaba `<vcf file with SV calls from SvABA>` \
-frag `<Window size of <INT> bases for integrating overlapping SV calls from tools. Default=500 bases>` \
-outdir `<Output directory path>` \
-mode `<germline or somatic. Default=germline>` \
-finalFile `<Output file name>`

This generates a tab delimited file with SV calls from all three tools. Example file: /example/Combined_SR.tsv

3. Combind SV calls from linked-reads and short-reads tools

python Combine_SR_LR.py \
-SR `<tab delimited file containing integrated SV calls from short-reads. example: /example/Combined_SR.tsv>` \
-LR1 `<vcf file from Long Ranger containing small and mid sized deletions>` \
-LR2 `<vcf file from Long Ranger containing large sized SV>` \
-NAIBR `<bed file containing SV calls from NAIBR>` \
-GROCSV `<vcf file containing SV calls from GROC-SV>` \
-w `<Window size of <INT> bass for integrating overlapping SV calls. Default=500 bases>` \
-outdir `<Output directory path>`

It generates a tab delimited file containing all overlapping and non overlapping SV calls from short-reads and linked-reads tools. Example file: /example/Combined_SR_LR.tsv

4. Perform requantification and add features

python Requantification.py \
-inputFile `<tab delimited file containing all SV calls. example: /example/Combined_SR_LR.tsv>` \
-out `<Ouput file name>` \
-n `<Number of cores to use. Default=1>` \
-refBit `<.2bit file of the reference genome>` \
-area `<Area in <INT> bases around the breakpoints for the rquantification genomic template. Default=500 bases>` \
-Read1_SR `<Paired-end read 1 from Illumina short-reads sequencing>` \
-Read2_SR `<Paired-end read 2 from Illumina short-reads sequencing>` \
-Read1_LR `<Paired-end read 1 from 10X Genomics linked-reads sequencing with barcodes trimmed>` \
-Read2_LR `<Paired-end read 2 from 10X Genomics linked-reads sequencing with barcodes trimmed>` \
-outdir `<Output directory path>` \
-tmpdir `<Temportary directory path>` \
-cutoff `<a read after re-alignment to genomic template has to overlap breakpoint with atleast <INT> bases. Deafult=10 bases>` \
-lengths `<length of chromosomes eg: /data/chromInfo_hg38.txt>`

It generates a tab delimited SV file with requantification related features like junction reads and spanning pairs calculated from short-read and/or linked-reads. Example: /example/Combined_SR_LR_requant.tsv. If you wish to only use short-reads or linked-reads calls then provide reads only from relevant technology.\
TwoBit file of the reference genome can be created as: [https://genome.ucsc.edu/goldenPath/help/twoBit.html]

5. Count supporting GEM from linked-reads sequencing alignment

python GEM_count.py \
-inputFile `<tab delimited file containing SV calls. Example: /example/Combined_SR_LR_requant.tsv>`\
-bam `<Lariat generated bam file with linked-reads>` \
-n `<Number of cores to be used. Default=1>` \
-outdir `<Output directory path>` \
-output `<Output file name>`
-w `<Window in <INT> bases to calculate supporting GEM for SV call. Default=1000 bases>` \

This generates tab delimited file containing suppporting GEM counts for each SV. Example: /example/Combined_SR_LR_requant_GEM.tsv

6. Include local coverage around breakpoints and annotate with repetitive and poor mappability regions

python Include_annotations.py \
-Annotate `<Annotate breakpoints with repetitive and poor mappability regions. Default=Yes>` \
-repeatMasker `<Repeat masker file. example: /data/RepeatMasker_hg38.bed>` \
-mappability `<Mappability track file. example: /data/k100.umap.bed>` \
-File `<tab delimited file containing SV calls. example: /example/Combined_SR_LR_requant_GEM.tsv>` \
-BAM `<short-reads BAM alignment file>` \
-BAMLinked `<linked-reads BAM alignment file>` \
-Threads `<Number of cores. Default=1>` \
-outdir `<Output directory path>`

This generates tab delimited file with local coverage around breakpoints calculated from short-reads and/or linked-reads alignment file. If "Annotate" argument is Yes, then each breakpoint is annotated with repetitive and poor mappability regions. Example file: /example/Combined_SR_LR_requant_GEM_annotated.tsv.\
Repeat masker and mappability files can be downloaded here: [http://genome.ucsc.edu/cgi-bin/hgTables]

7. Train and/or predict probability of true SV from short-reads technology

Rscript Model_cWGS.R\
-f `<Enter .tsv (tab delimited file) with features for each SV. example: /example/Combined_SR_LR_requant_GEM_annotated.tsv>`\
-t `<Enter whether to train a new model (yes or no). Default=no>`\
-m `<Enter the trained model for predicting true SV calls. example: /Model/SR_Model.rds>`\
-r `<Number of total paired-end reads in the sample (required for normalization of requantification features)>`\
-s `<Mark positive SV calls with probability threshold score > <FLOAT>. Default=0.6 (max=1.0, min=0.0)>`\
-o `<Output file name with 'PredictionSR' field with positive SV calls>`\
-d `<Output directory path. Default=current directory>`

This generates a tab delimitd file with a field "PredictionSR" reporting "Positive" for true SV calls. Example file: /example/cWGS_predictions.tsv \
If you wish to train a new logistic regression model then include a column "Validation" with Positive or Negative calls in input file (-f argument). And pass -t argument with "yes". This file should also contain features/columns derived from short-reads sequencing alignments. The columns are: JR_SR (junction reads), SP_SR (spanning pairs), SVType (Dels, Dups, Invs or Trans), Size (size of SV), LocalCoverage_Pos1_SR (local coverage around breakpoint 1), LocalCoverage_Pos2_SR (local coverage around breakpoint 2).\
Or else the trained model (/Model/SR_Model.rds) could be used for predicting the probability of true SV calls.

8. Train and/or predict probability of true SV from linked-reads technology

Rscript Model_10XWGS.R\
-f `<Enter .tsv (tab delimited file) with features for each SV. example: /example/Combined_SR_LR_requant_GEM_annotated.tsv>`\
-t `<Enter whether to train a new model (yes or no). Default=no>`\
-m `<Enter the trained model for predicting true SV calls. example: /Model/LR_Model.rds>`\
-r `<Number of total paired-end reads in the sample (required for normalization of requantification features)>`\
-g `<Number of GEMs detected in the linked-reads experiment (required for normalization of GEM count)>`\
-s `<Mark positive SV calls with probability threshold score > <FLOAT>. Default=0.6 (max=1.0, min=0.0)>`\
-o `<Output file name with 'PredictionLR' field with positive SV calls>`\
-d `<Output directory path. Default=current directory>`

This generates a tab delimitd file with a field "PredictionLR" reporting "Positive" for true SV calls. Example file: /example/10XWGS_predictions.tsv \
If you wish to train a new logistic regression model then include a column "Validation" with Positive or Negative calls in input file (-f argument). And pass -t argument with "yes". This file should also contain features/columns derived from linked-reads sequencing alignments. The columns are: JR_LR (junction reads), SP_LR (spanning pairs), SVType (Dels, Dups, Invs or Trans), Size (size of SV), LocalCoverage_Pos1_LR (local coverage around breakpoint 1), LocalCoverage_Pos2_LR (local coverage around breakpoint 2) and GEM (normalized GEM count supporting the SV).\
Or else the trained model (/Model/LR_Model.rds) could be used for predicting the probability of true SV calls.
