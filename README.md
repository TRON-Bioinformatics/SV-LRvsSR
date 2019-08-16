# SV-LRvsSR

A study to compare structural variation (SV) predictions from 10X Genomics linked-reads sequencing (LR) and conventional Illumina short-reads sequencing (SR). We predict SVs from Long Ranger pipeline of 10X Genomics and combine predictions from three tools (Delly, Lumpy, SvABA) for SR. The MCF7 breast cancer cell line and a primary breast tumor sample was analysed for the comparison. Apart from comparing different predictions by two technologies we also propose precision enrichment models for retaining true positive SVs amongst the high false positive calls from both the technologies.

The code for analysis was written in python (version 2.7.15) and precision enrichment models were trained and built in R (version 3.4.2). The commands for running all tools and analysis pipelines are mentioned in commands.org. The codes are present in folder Script.
