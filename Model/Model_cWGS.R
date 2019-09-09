require(caret) 
library(dplyr)
library(DMwR) 
library(purrr)
require(readr)
library(MASS)

Tog<- read.csv("MCF7_Final_call.tsv", sep='\t', header=T,
               stringsAsFactors = F, na.strings = c("","NA"," ",NA,NaN))

### prepare data for features
Tog<- Tog %>% 
  filter(!((SVType=="Invs" & Size< 10000) | (SVType=="Dups" & Size <10000))) %>% 
  filter(Size>50|Size==0) %>% 
  filter(Linked_SVType!='UNK'|is.na(Linked_SVType)) %>% 
  filter(Predicted_by %in% c("Common","Only cWGS"))
colnames(Tog)[which(colnames(Tog)=="Predicted_by")]<-"Category"
sample<-"MCF7"

prepareDF<- function(DF){
  File<- DF %>% 
    dplyr::select(Norm_JR_SR,Norm_SP_SR,PCR,Category,SVType)
  File$JRS_SR=log2(File$Norm_JR_SR + File$Norm_SP_SR +1)
  File$Norm_JR_SR<- log2(File$Norm_JR_SR+1)
  File$Norm_SP_SR<- log2(File$Norm_SP_SR+1)
  File$SVType<- as.factor(File$SVType)
  File$Category<- as.factor(File$Category)
  colnames(File)<-c("JR_SR","SP_SR","PCR","Category", "SVType","JRS")
  File$PCR[File$PCR=="FALSE"]<-"Negative"
  File$PCR[File$PCR=="TRUE"]<-"Positive"
  File$PCR<- factor(File$PCR, levels=c("Negative","Positive"))
  return(File)
}

ForTraining<- prepareDF(filter(Tog, !is.na(PCR)))
ForPrediction<- prepareDF(filter(Tog, is.na(PCR)))

### plot number of positive and negative class samples in training data

plotData<- function(DF, NameFile){
  nDF<- DF %>% 
    dplyr::count(PCR) %>% 
    mutate(percent= (n /sum(n)) *100)
  nDF$PCR<- as.character(nDF$PCR)
  nDF$PCR<- factor(nDF$PCR, levels = c("Negative","Positive"))
  col_pal<-c("#A6CEE3","#1F78B4")
  p <- ggplot(nDF, aes(
    x = PCR, 
    y = percent, 
    label = paste0(round(percent, 2), "%\nn = ", n), 
    fill = PCR)) + 
    geom_bar(stat = "identity", color = "black") + 
    geom_text(size = 5)+#, vjust = "inward") + 
    scale_fill_manual(values = col_pal)+
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          axis.text= element_text(size=10),
          axis.title = element_text(size=10))+
    ggtitle(paste("Labels for training data\n(Total=",sum(nDF$n),")", sep='')) +
    labs(x="PCR Validation Status", y="Percent")
  p + ggsave(file = NameFile, w =3, h = 6)
  
}

### funtion to train and test logistic regression model
cal_result<- function(Training_SR, Testing_SR, train.control){
  set.seed(617)
  model_fit <- train(PCR ~ .,
                     data = Training_SR,
                     method = "glm", maximize=TRUE, metric="ROC", 
                     trControl= train.control, family=binomial(link='logit'))
  print(model_fit)
  Modelpredict<-predict(model_fit, Testing_SR, type="prob")
  Modelpredict<- cbind(Modelpredict, SVType=Testing_SR$SVType)
  Model2_predict<- Modelpredict %>% 
    mutate(pred_class=ifelse(Positive>0.75,"Positive","Negative"))
  conf<-table(Predicted=Model2_predict$pred_class, Reference=Testing_SR$PCR)
  Table<- confusionMatrix(conf, positive = "Positive", mode = "everything")
  Result<- c(Table$byClass['Sensitivity'], Table$byClass['Specificity'], Table$byClass['Precision'], 
             Table$byClass['F1'], Table$byClass['Balanced Accuracy'])
  Parameter<-c("Sensitivity","Specificity","Precision","F1","Accuracy")
  DF_Result<- data.frame(Parameter=factor(Parameter),
                         Values=Result)
                         
  return(list(model_fit,DF_Result))
}

Modeling<-function(TrainingSet, samplingType){
  DF_Result<- data.frame()
  DataSR<-TrainingSet[,c("JR_SR","SP_SR","SVType","PCR")]  
  train.control<- trainControl(method="repeatedcv",
                               number=10,
                               repeats=3,
                               classProbs=TRUE,
                               summaryFunction = twoClassSummary,
                               savePredictions= TRUE)
  set.seed(156)
  for (repeats in 1:10){
    DataSR<- DataSR[sample(nrow(DataSR)),]
    intrain_SR<- createDataPartition(y=DataSR$PCR, p=0.7, list=FALSE)
    Training_SR<- DataSR[intrain_SR,]
    plotData(Training_SR, "Training.SR.plot.png")  
    Testing_SR<- DataSR[-intrain_SR,]
      if (samplingType!="none"){
        train.control$sampling<- samplingType }
    NewResult<- cal_result(Training_SR, Testing_SR, train.control)
    DF_Result<- rbind(DF_Result, NewResult[[2]])
  }
  return(list(NewResult[[1]],DF_Result))
}

### Train logistic regression model with unbalanced data set
Unbalanced<- Modeling(ForTraining, "none")

### Train logistic regression model with training data balanced with downsampling
Bal_down<- Modeling(ForTraining, "down")

### Train logistic regression model with training data balanced with upsampling
Bal_up<- Modeling(ForTraining, "up")

### Train logistic regression model with training data balanced with SMOTE
Bal_smote<- Modeling(ForTraining, "smote")

### plot performance of trained models with data set unbalanced or balanced by different sammpling type
Results<- data.frame()
Results<- rbind(Results, cbind(Unbalanced[[2]], SamplingType="Unbalanced"))
Results<- rbind(Results, cbind(Bal_down[[2]], SamplingType="downsampling"))
Results<- rbind(Results, cbind(Bal_up[[2]], SamplingType="upsampling"))
Results<- rbind(Results, cbind(Bal_smote[[2]], SamplingType="SMOTE"))

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

FinalDF<- data_summary(Results,
                       varname = "Values",
                       groupnames = c("Parameter","SamplingType"))
write.csv(FinalDF, "Performance_SR.csv", quote=F, row.names = F)

p1<- ggplot(FinalDF, aes(x=Parameter, y=Values, fill=SamplingType))+
  geom_bar(stat = "identity", color="black", position = position_dodge())+
  geom_errorbar(aes(ymin=Values-sd, ymax=Values+sd), width=0.2, position = position_dodge(0.9))+
  geom_text(aes(label=round(Values,2)), position = position_dodge(width=0.9), vjust=-0.4,hjust=1.1, size=2, fontface= "bold")+
  theme_bw()+
  scale_fill_brewer(palette = "Set3")+
  theme(legend.text = element_text(size=8),
        # legend.position = "bottom",
        axis.text =  element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.title.x=element_blank()) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1))
p1 + ggsave("MCF7_performance.png", width=10, height = 4)

### final model trained on complete data set (unbalanced data)
set.seed(9867)
train.control<- trainControl(method="repeatedcv",
                             number=10,
                             repeats=3,
                             classProbs=TRUE,
                             summaryFunction = twoClassSummary,
                             savePredictions= TRUE)
FinalModel<- train(PCR ~ JR_SR+SP_SR+SVType,
                   data = ForTraining,
                   method = "glm", maximize=TRUE, metric="ROC",
                   trControl= train.control, family=binomial(link='logit'))
write_rds(FinalModel,"SR_Model.rds") # the trained model for cWGS technology

### predict all labels for remaining SVs
ApplyModel<- function(ForPrediction, Model, sample, tech, samplingType){
  DF_predict<- predict(Model, 
                       ForPrediction,
                       type = "prob")
  ForPrediction<- cbind(ForPrediction, pred_class=DF_predict$Positive)
  ForPrediction<- mutate(ForPrediction, PredictionSR=ifelse(pred_class>0.75,"Positive","Negative"))
  ForPrediction$SVType<- as.character(ForPrediction$SVType)
  ForPrediction$SVType[ForPrediction$SVType=="Dels"]<-"Deletion"
  ForPrediction$SVType[ForPrediction$SVType=="Dups"]<-"Duplication"
  ForPrediction$SVType[ForPrediction$SVType=="Invs"]<-"Inversion"
  ForPrediction$SVType[ForPrediction$SVType=="Trans"]<-"Translocation"
  ForPrediction$SVType<- factor(ForPrediction$SVType, levels = c("Deletion","Duplication","Inversion", "Translocation"))
  Count_Table<- ForPrediction %>%
    dplyr::group_by(SVType, PredictionSR) %>%
    dplyr::summarise(Count=n()) %>%
    dplyr::mutate(Percent=Count/sum(Count) * 100)
  Table2<- ForPrediction %>%
    dplyr::group_by(SVType) %>%
    dplyr::summarise(Count=n())
  p1<- ggplot(Count_Table, aes(x= SVType, y=Count, fill=PredictionSR))+
    geom_bar(stat="identity")+
    scale_fill_brewer(palette = "Paired")+
    theme_bw()+
    theme(legend.position = "bottom",
          legend.text = element_text(size=10),
          plot.title = element_text(size=10,hjust = 0.5, face = "bold"),
          axis.text = element_text(size=10),
          axis.title = element_text(size=10,face="bold")) +
    labs(x="Predicted SV", y="Counts")+ 
    geom_text(data=filter(Count_Table, PredictionSR=="Positive"),
              aes(x= SVType, y= Count, label=paste(round(Percent,2),"%",sep='')),
              size=3, hjust=0.5, vjust=-0.2 ) +
    annotate("text", x=c(1,2,3,4), y=max(Table2$Count)+500, label=paste("N=",Table2$Count,sep=''))+
    ggtitle(paste(sample,": ", tech," decision tree predictions", sep=''))
  p1<-p1 + ggsave(paste(sample,"_",tech,"_",samplingType,"_predictions.png",sep=''), width=4.5, height = 5)
  return(ForPrediction$PredictionSR)  
}

UnbalancedSR<-ApplyModel(prepareDF(Tog), FinalModel, sample = "MCF7", tech = "cWGS", samplingType = "Unbalanced")
DownSR<-ApplyModel(prepareDF(Tog), Bal_down[[1]], sample = "MCF7", tech = "cWGS", samplingType = "DownSampling")
UpSR<-ApplyModel(prepareDF(Tog), Bal_up[[1]], sample = "MCF7", tech = "cWGS", samplingType = "UpSampling")
SMOTESR<-ApplyModel(prepareDF(Tog), Bal_smote[[1]], sample = "MCF7", tech = "cWGS", samplingType = "SMOTE")

write.table(cbind(Tog, UnbalancedPrediction=UnbalancedSR, DownSamplingPrediction=DownSR,
                UpSamplingPrediction=UpSR, SMOTESamplingPrediction=SMOTESR), "cWGS_Predictions.tsv",sep='\t',
          quote=F, row.names = F)

### Prediction of labels for primary tumor's structural variations
TogTumor<- read.csv("Primary_tumor_Final_calls.tsv", sep='\t', header=T,
                    stringsAsFactors = F, na.strings = c("","NA"," ",NA,NaN))
FinalModel<- read_rds("SR_Model.rds")
TogTumor<- TogTumor %>% 
  filter(!((SVType=="Invs" & Size< 10000) | (SVType=="Dups" & Size <10000))) %>% 
  filter(Size>50|Size==0) %>% 
  filter(Linked_SVType!='UNK'|is.na(Linked_SVType)) %>% 
  filter(Predicted_by %in% c("Common","Only cWGS"))
TogTumor$JRS_SR=log2(TogTumor$Norm_JR_SR + TogTumor$Norm_SP_SR +1)
TogTumor$Norm_JR_SR<- log2(TogTumor$Norm_JR_SR+1)
TogTumor$Norm_SP_SR<- log2(TogTumor$Norm_SP_SR+1)
TogTumor$SVType<- as.factor(TogTumor$SVType)
colnames(TogTumor)[which(colnames(TogTumor)=="Predicted_by")]<-"Category"
TogTumor$Category<- as.factor(TogTumor$Category)

UnbalancedSR_Tumor<-ApplyModel(TogTumor, FinalModel, sample = "Primary tumor", 
                               tech = "cWGS", samplingType = "Unbalanced")
write.table(cbind(TogTumor, UnbalancedPrediction=UnbalancedSR_Tumor), "Tumor_cWGS_Predictions.tsv",sep='\t',
            quote=F, row.names = F)
