require(caret)     
library(dplyr)
library(DMwR) 
library(purrr)
require(readr)

Tog<- read.csv("MCF7_Final_call.tsv", sep='\t', header=T,
               stringsAsFactors = F, na.strings = c("","NA"," ",NA,NaN))

### prepare data for features
Tog<- Tog %>% 
  filter(!((SVType=="Invs" & Size< 10000) | (SVType=="Dups" & Size <10000))) %>% 
  filter(Size>50|Size==0) %>% 
  filter(Linked_SVType!='UNK'|is.na(Linked_SVType)) %>% 
  filter(Predicted_by %in% c("Common","Only 10XWGS"))
colnames(Tog)[which(colnames(Tog)=="Predicted_by")]<-"Category"
sample<-"MCF7"

prepareDF<- function(Tog){
File<- Tog %>% 
  dplyr::select("GEM","Norm_JR_LR","Norm_SP_LR","Size","PCR",
         "Category","SVType")
File$JRS_LR=log2(File$Norm_JR_LR + File$Norm_SP_LR +1)
File$Norm_JR_LR<- log2(File$Norm_JR_LR+1)
File$Norm_SP_LR<- log2(File$Norm_SP_LR+1)
File$SVType<- as.factor(File$SVType)
File$Category<- as.factor(File$Category)
File$GEM<- log2(File$GEM+1)
colnames(File)<-c("GEM","JR_LR","SP_LR","Size","PCR","Category", "SVType", 
                  "JRS_LR")
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
cal_result<- function(Training_LR, Testing_LR, train.control){
  set.seed(5627)
  model_fit <- train(PCR ~ .,
                     data = Training_LR,
                     method = "glm", maximize=TRUE,metric="ROC",
                     trControl= train.control)
  summary(model_fit)
  Modelpredict<-predict(model_fit, Testing_LR, type="prob")
  Modelpredict<- cbind(Modelpredict, SVType=Testing_LR$SVType)
  Model2_predict<- Modelpredict %>% 
    mutate(pred_class=ifelse(Positive>0.75,"Positive","Negative"))
  conf<-table(Predicted=Model2_predict$pred_class, Reference=Testing_LR$PCR)
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
  DataLR<-TrainingSet[,c("GEM","JR_LR","SP_LR","SVType", "PCR")]  
  train.control<- trainControl(method="repeatedcv",
                               number=10,
                               repeats=3,
                               classProbs=TRUE,
                               summaryFunction = twoClassSummary,
                               savePredictions= TRUE)
  set.seed(1256)
  for (repeats in 1:10){
    DataLR<- DataLR[sample(nrow(DataLR)),]
    intrain_LR<- createDataPartition(y=DataLR$PCR, p=0.7, list=FALSE)
    Training_LR<- DataLR[intrain_LR,]
    plotData(Training_LR, "Training.LR.plot.png")  
    Testing_LR<- DataLR[-intrain_LR,]
    if (samplingType!="none"){
        train.control$sampling<- samplingType  
      }
    NewResult<- cal_result(Training_LR, Testing_LR, train.control)
    DF_Result<- rbind(DF_Result, NewResult[[2]])
  }
  return(list(NewResult[[1]], DF_Result))
  }

### Train logistic regression model
Normal<- Modeling(ForTraining, "none")

### plot performance of trained model
Results<- as.data.frame(Normal[[2]])

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
                       groupnames = c("Parameter"))

write.csv(FinalDF, "Performance_LR.csv", quote=F, row.names = F)

p1<- ggplot(FinalDF, aes(x=Parameter, y=Values, fill=Parameter))+
  geom_bar(stat = "identity", color="black", position = position_dodge())+
  geom_text(aes(label=round(Values,2)), position = position_dodge(width=0.9), vjust=-0.2,hjust=1.1, size=2.8)+
  geom_errorbar(aes(ymin=Values-sd, ymax=Values+sd), width=0.2, position = position_dodge(0.9))+
  theme_bw()+
  scale_fill_brewer(palette = "Pastel1")+
  theme(axis.text =  element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.title.x=element_blank(),
        legend.text = element_blank(),
        legend.position = "none") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1))
p1 + ggsave("MCF7_performance_LR.png", width=5, height = 4)

### Final model trained on complete data set
train.control<- trainControl(method="repeatedcv",
                             number=10,
                             repeats=3,
                             classProbs=TRUE,
                             summaryFunction = twoClassSummary,
                             savePredictions= TRUE)
FinalModel <- train(PCR ~ GEM+JR_LR+SP_LR+SVType,
                   data = ForTraining,
                   method = "glm", maximize=TRUE,metric="ROC",
                   trControl= train.control)
write_rds(FinalModel,"LR_Model.rds") # the traned model for 10XWGS technology
FinalModel<- read_rds("LR_Model.rds")

### predict all labels for remaining SVs
ApplyModel<- function(ForPrediction, Model, sample, tech, samplingType){
  DF_predict<- predict(Model, 
                       ForPrediction,
                       type = "prob")
  ForPrediction<- cbind(ForPrediction, pred_class=DF_predict$Positive)
  ForPrediction<- mutate(ForPrediction, PredictionLR=ifelse(pred_class>0.75,"Positive","Negative"))
  ForPrediction$SVType<- as.character(ForPrediction$SVType)
  ForPrediction$SVType[ForPrediction$SVType=="Dels"]<-"Deletion"
  ForPrediction$SVType[ForPrediction$SVType=="Dups"]<-"Duplication"
  ForPrediction$SVType[ForPrediction$SVType=="Invs"]<-"Inversion"
  ForPrediction$SVType[ForPrediction$SVType=="Trans"]<-"Translocation"
  ForPrediction$SVType<- factor(ForPrediction$SVType, levels = c("Deletion","Duplication","Inversion", "Translocation"))
  Count_Table<- ForPrediction %>%
    dplyr::group_by(SVType, PredictionLR) %>%
    dplyr::summarise(Count=n()) %>%
    dplyr::mutate(Percent=Count/sum(Count) * 100)
  Table2<- ForPrediction %>%
    dplyr::group_by(SVType) %>%
    dplyr::summarise(Count=n())
  
  p1<- ggplot(Count_Table, aes(x= SVType, y=Count, fill=PredictionLR))+
    geom_bar(stat="identity")+
    scale_fill_brewer(palette = "Paired")+
    theme_bw()+
    theme(legend.position = "bottom",
          legend.text = element_text(size=10),
          plot.title = element_text(size=10,hjust = 0.5, face = "bold"),
          axis.text = element_text(size=10),
          axis.title = element_text(size=10,face="bold")) +
    labs(x="Predicted SV", y="Counts")+  
    geom_text(data=filter(Count_Table, PredictionLR=="Positive"),
              aes(x= SVType, y= Count, label=paste(round(Percent,2),"%",sep='')),
              size=3, hjust=0.5, vjust=-0.2 ) +
    annotate("text", x=c(1,2,3,4), y=max(Table2$Count)+500, label=paste("N=",Table2$Count,sep=''))+
    ggtitle(paste(sample,": ", tech," decision tree predictions", sep=''))
  p1<-p1 + ggsave(paste(sample,"_",tech,"_",samplingType,"_predictions.png",sep=''), width=4.5, height = 5)
  return(ForPrediction$PredictionLR)  
}

NormalLR<-ApplyModel(prepareDF(Tog), FinalModel, sample = "MCF7", tech = "10XWGS", samplingType = "Normal")

### Prediction of labels for primary tumor's structural variations
TogTumor<- read.csv("Primary_tumor_Final_calls.tsv", sep='\t', header=T,
                    stringsAsFactors = F, na.strings = c("","NA"," ",NA,NaN))
TogTumor<- TogTumor %>% 
  filter(!((SVType=="Invs" & Size< 10000) | (SVType=="Dups" & Size <10000))) %>% 
  filter(Size>50|Size==0) %>% 
  filter(Linked_SVType!='UNK'|is.na(Linked_SVType)) %>% 
  filter(Predicted_by %in% c("Common","Only 10XWGS"))
TogTumor$JRS_LR=log2(TogTumor$Norm_JR_LR + TogTumor$Norm_SP_LR +1)
TogTumor$Norm_JR_LR<- log2(TogTumor$Norm_JR_LR+1)
TogTumor$Norm_SP_LR<- log2(TogTumor$Norm_SP_LR+1)
TogTumor$SVType<- as.factor(TogTumor$SVType)
colnames(TogTumor)[which(colnames(TogTumor)=="Predicted_by")]<-"Category"
TogTumor$Category<- as.factor(TogTumor$Category)
TogTumor$GEM<- log2(TogTumor$GEM+1)

NormalLR_Tumor<-ApplyModel(TogTumor, FinalModel, sample = "Primary tumor", tech = "10XWGS", samplingType = "Normal")

write.table(cbind(Tog, Prediction=NormalLR), "10XWGS_Predictions.tsv",sep='\t',
          quote=F, row.names = F)
write.table(cbind(TogTumor, Prediction=NormalLR_Tumor), "Tumor_10XWGS_Predictions.tsv",sep='\t',
            quote=F, row.names = F)
