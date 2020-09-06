### training and testing of cWGS Model
require(optparse)
require(caret) 
library(dplyr)
library(plyr)
library(tidyverse)

# Parse commandline arguments --------------------------------------------------
options(stringsAsFactors = FALSE)

argument_list <- list(
  make_option(c("-f", "--SVFile"), default="", help="Enter .tsv file with list of SV"), 
  make_option(c("-t", "--training", default="no", help="Enter whether model needs to be trained or not (yes or no; default=no")),
  make_option(c("-r", "--total_reads", default="", help="Enter number of total paired-end reads in the sample")),
  make_option(c("-m", "--model_file"), default="", help="Path to .rds file containing the machine learning model (./Model/SR_Model.rds)"),
  make_option(c("-s", "--prediction_threshold"), default=0.6, 
              help="Prediction score threshold above which a SV is called positive (default=0.7, max=1, min=0)"),
  make_option(c("-o", "--output"), default="", help="Enter output file name"),
  make_option(c("-d", "--outdir", default=".", help="Enter output directory"))
)
opt <- parse_args(OptionParser(option_list=argument_list))

### check required arguments
if(is.na(opt$SVFile) | opt$SVFile == "") {
  stop("Missing input file with list of SV...")
}

if(is.na(opt$output) | opt$output == "") {
  stop("Missing output file name...")
}

if(is.na(opt$total_reads) | opt$total_reads == "") {
  stop("Missing total number of paired-end reads in the sample...")
}


### read files and arguments
Tog<- read.csv(opt$SVFile, sep = '\t', header=T, stringsAsFactors = F, na.strings = c("","NA"," ",NA,NaN))

if ("Category" %in% colnames(Tog)){
  Tog<- Tog %>% 
    dplyr::filter(Size>50|Size==0) %>% 
    dplyr::mutate(JR_SR = log2(((JR_SR/(2*opt$total_reads))*100000000) + 1),
                  SP_SR = log2(((SP_SR/(opt$total_reads))*100000000) + 1),
                  LocalCoverage_Pos1 = log2(LocalCoverage_Pos1_SR +1),
                  LocalCoverage_Pos2 = log2(LocalCoverage_Pos2_SR +1),
                  Size = log10(Size+1))
  
} else {
  Tog<- Tog %>% 
    dplyr::filter(Size>50|Size==0) %>% 
    dplyr::mutate(SVType = factor(SVType, levels=c("Dels","Dups","Invs","Trans")),
                  JR_SR = log2(((JR_SR/(2*opt$total_reads))*100000000) + 1),
                  SP_SR = log2(((SP_SR/(opt$total_reads))*100000000) + 1),
                  LocalCoverage_Pos1 = log2(LocalCoverage_Pos1_SR +1),
                  LocalCoverage_Pos2 = log2(LocalCoverage_Pos2_SR +1),
                  Size = log10(Size+1))
}

### Train a logistic regression model if training argument is yes

'%ni%' <- Negate('%in%')
featuresPresent="yes"
if (opt$training=="yes" | opt$training=="Yes"){
  ### check if all required features/predictors present in input file
  if ("JR_SR" %ni% colnames(Tog)){
    stop("Junction reads (JR_SR) missing in input file")
    featuresPresent="no"
  }
  if ("SP_SR" %ni% colnames(Tog)){
    stop("Spanning reads (SP_SR) missing in input file")
    featuresPresent="no"
  }
  if ("LocalCoverage_Pos1_SR" %ni% colnames(Tog)){
    stop("Local coverage around breakpoint 1 (LocalCoverage_Pos1_SR) missing in input file")
    featuresPresent="no"
  }
  if ("LocalCoverage_Pos2_SR" %ni% colnames(Tog)){
    stop("Local coverage around breakpoint 2 (LocalCoverage_Pos2_SR) missing in input file")
    featuresPresent="no"
  }
  if ("Size" %ni% colnames(Tog)){
    stop("Size of SV missing in input file")
    featuresPresent="no"
  }
  if ("SVType" %ni% colnames(Tog)){
    stop("Type of SV (Dels,Dups,Invs,Trans) missing in input file")
    featuresPresent="no"
  }
  if ("Validation" %ni% colnames(Tog)){
    stop("Validated data points (Positive or Negative) missing in input file")
    featuresPresent="no"
  }
  
  ### train the model
  if (featuresPresent=="yes"){
    set.seed(9867)
    train.control<- trainControl(method="repeatedcv",
                                 number=10,
                                 repeats=3,
                                 classProbs=TRUE,
                                 summaryFunction = twoClassSummary,
                                 savePredictions= TRUE)
    FinalModel<- train(Validation ~ JR_SR + SP_SR + SVType + Size + LocalCoverage_Pos1 + LocalCoverage_Pos2,
                       data = dplyr::filter(Tog, is.na(Validation)),
                       method = "glm", maximize=TRUE, metric="ROC",
                       trControl= train.control, family=binomial(link='logit'))
    write_rds(FinalModel, paste0(opt$outdir,"/SR_Model.rds")) 
    
  }
  else {
    stop("Missing features")
  }
} else {
  ### check missing features
  if ("JR_SR" %ni% colnames(Tog)){
    stop("Junction reads (JR_SR) missing in input file")
  }
  if ("SP_SR" %ni% colnames(Tog)){
    stop("Spanning reads (SP_SR) missing in input file")
  }
  if ("LocalCoverage_Pos1_SR" %ni% colnames(Tog)){
    stop("Local coverage around breakpoint 1 (LocalCoverage_Pos1_SR) missing in input file")
  }
  if ("LocalCoverage_Pos2_SR" %ni% colnames(Tog)){
    stop("Local coverage around breakpoint 2 (LocalCoverage_Pos2_SR) missing in input file")
  }
  if ("Size" %ni% colnames(Tog)){
    stop("Size of SV missing in input file")
  }
  if ("SVType" %ni% colnames(Tog)){
    stop("Type of SV (Dels,Dups,Invs,Trans) missing in input file")
  }
  FinalModel<- read_rds(opt$model_file)
}

### make predictions

ApplyModel<- function(ForPrediction, Model, threshold){
  DF_predict<- predict(Model, 
                       ForPrediction,
                       type = "prob")
  ForPrediction<- cbind(ForPrediction, pred_class = DF_predict$Positive)
  if ("Category" %in% colnames(ForPrediction) & "Common" %in% unique(ForPrediction$Category)){
    ForPrediction<- ForPrediction %>% 
      dplyr::mutate(PredictionSR = ifelse(Category=="Common", "Positive", ifelse(pred_class>threshold, "Positive", "Negative")))  
  } else {
    ForPrediction<- ForPrediction %>% 
      dplyr::mutate(PredictionSR = ifelse(pred_class>threshold, "Positive", "Negative"))  
  }
  return(ForPrediction)  
}

FinalPredictions<- ApplyModel(Tog, FinalModel, opt$prediction_threshold)

### write prediction for all SV calls
write.table(FinalPredictions, paste0(opt$outdir,"/",opt$output,".tsv"), sep='\t', quote=F, row.names = F)

