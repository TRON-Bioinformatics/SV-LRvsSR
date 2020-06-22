### Training and testing of 10XWGS model
require(caret) 
library(dplyr)
library(plyr)
library(tidyverse)

### read file containing all the features (SupplementaryTable3 )
Tog<- read.csv("SupplementaryTable3.tsv", sep='\t', header=T,
               stringsAsFactors = F, na.strings = c("","NA"," ",NA,NaN))

### prepare data for features
Tog<- Tog %>% 
  filter(Size>50|Size==0) %>% 
  filter(Predicted_by %in% c("Common","Only 10XWGS"))
colnames(Tog)[which(colnames(Tog)=="Predicted_by")]<-"Category"
Tog$Internal<-NA
Tog$Internal[grep("Internal", Tog$Comment)]<- "Internal"
Tog$Study1<-NA
Tog$Study1[grep("Weaver", Tog$Comment)]<- "Study1"
Tog$Study2<-NA
Tog$Study2[grep("Hillmer", Tog$Comment)]<- "Study2"
sample<-"MCF7"

### prepare data for training
ForTraining<- Tog %>% 
  dplyr::filter(!is.na(Validation) & !is.na(Internal)) %>% 
  dplyr::mutate(SVType = factor(SVType, levels=c("Dels","Dups","Invs","Trans")),
                Category = factor(Category),
                Norm_JR_LR = log2(Norm_JR_LR +1),
                Norm_SP_LR = log2(Norm_SP_LR +1),
                JRS= log2(Norm_JR_LR + Norm_SP_LR +1),
                LocalCoverage_Pos1 = log2(LocalCoverage_Pos1_LR +1),
                LocalCoverage_Pos2 = log2(LocalCoverage_Pos2_LR +1),
                PCR = ifelse(Validation=="TRUE", "Positive", "Negative"),
                GEM = log2(GEM +1),
                Size = log10(Size+1)) %>% 
  dplyr::mutate(PCR = factor(PCR, levels = c("Negative", "Positive"))) %>% 
  dplyr::select(Norm_JR_LR, Norm_SP_LR, JRS, PCR, SVType, Category, LocalCoverage_Pos1, LocalCoverage_Pos2, GEM, Size) %>% 
  dplyr::rename(JR_LR = Norm_JR_LR, SP_LR = Norm_SP_LR)

### prepare data for prediction
ForPrediction<- Tog %>% 
  dplyr::filter(is.na(Internal)) %>% 
  dplyr::mutate(SVType = factor(SVType, levels=c("Dels","Dups","Invs","Trans")),
                Category = factor(Category),
                Norm_JR_LR = log2(Norm_JR_LR +1),
                Norm_SP_LR = log2(Norm_SP_LR +1),
                JRS= log2(Norm_JR_LR + Norm_SP_LR +1),
                LocalCoverage_Pos1 = log2(LocalCoverage_Pos1_LR +1),
                LocalCoverage_Pos2 = log2(LocalCoverage_Pos2_LR +1),
                PCR = ifelse(Validation=="TRUE", "Positive", "Negative"),
                GEM = log2(GEM +1),
                Size = log10(Size+1)) %>% 
  dplyr::mutate(PCR = factor(PCR, levels = c("Negative", "Positive"))) %>% 
  dplyr::select(Norm_JR_LR, Norm_SP_LR, JRS, PCR, SVType, Category, LocalCoverage_Pos1, LocalCoverage_Pos2, GEM, Size) %>% 
  dplyr::rename(JR_LR = Norm_JR_LR, SP_LR = Norm_SP_LR)

prepareDF<- function(DF){
  DF<- DF %>% 
    dplyr::mutate(SVType = factor(SVType, levels=c("Dels","Dups","Invs","Trans")),
                  Category = factor(Category),
                  Norm_JR_LR = log2(Norm_JR_LR +1),
                  Norm_SP_LR = log2(Norm_SP_LR +1),
                  JRS= log2(Norm_JR_LR + Norm_SP_LR +1),
                  LocalCoverage_Pos1 = log2(LocalCoverage_Pos1_LR +1),
                  LocalCoverage_Pos2 = log2(LocalCoverage_Pos2_LR +1),
                  PCR = ifelse(Validation=="TRUE", "Positive", "Negative"),
                  GEM = log2(GEM +1),
                  Size = log10(Size +1)) %>% 
    dplyr::mutate(PCR = factor(PCR, levels = c("Negative", "Positive"))) %>% 
    dplyr::select(Norm_JR_LR, Norm_SP_LR, JRS, PCR, SVType, Category, LocalCoverage_Pos1, LocalCoverage_Pos2, GEM, Size) %>% 
    dplyr::rename(JR_LR = Norm_JR_LR, SP_LR = Norm_SP_LR)
  return(DF)
}

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
  Imp <- varImp(model_fit)$importance %>% as.data.frame() %>% rownames_to_column(var="Features")
  
  Modelpredict<-predict(model_fit, Testing_LR, type="prob")
  Modelpredict<- cbind(Modelpredict, SVType=Testing_LR$SVType)
  Model2_predict<- Modelpredict %>% 
    mutate(pred_class=ifelse(Positive>0.6,"Positive","Negative"))
  conf<-table(Predicted=Model2_predict$pred_class, Reference=Testing_LR$PCR)
  Table<- confusionMatrix(conf, positive = "Positive", mode = "everything")
  Result<- c(Table$byClass['Sensitivity'], Table$byClass['Specificity'], Table$byClass['Precision'], 
             Table$byClass['F1'], Table$byClass['Balanced Accuracy'])
  Parameter<-c("Sensitivity","Specificity","Precision","F1","Accuracy")
  DF_Result<- data.frame(Parameter=factor(Parameter),
                         Values=Result)
  return(list(model_fit, DF_Result, Imp))
}

Modeling<-function(TrainingSet, samplingType){
  DF_Result<- data.frame()
  DataLR<-TrainingSet[,c("GEM","JR_LR","SP_LR","SVType", "PCR", "Size", "LocalCoverage_Pos1", "LocalCoverage_Pos2")]  
  train.control<- trainControl(method="repeatedcv",
                               number=10,
                               repeats=3,
                               classProbs=TRUE,
                               summaryFunction = twoClassSummary,
                               savePredictions= TRUE)
  set.seed(1256)
  Importance<-data.frame()
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
    Importance<- rbind.fill(Importance, NewResult[[3]])
  }
  return(list(NewResult[[1]], DF_Result, Importance))
  }

### Train logistic regression model
Normal<- Modeling(ForTraining, "none")

### plot feature importance plot
df2<- Normal[[3]] %>% 
  dplyr::group_by(Features) %>% 
  dplyr::summarise(Mean = mean(Overall),
                   SD = sd(Overall),
                   Count = n(),
                   SE = sd(Overall)/sqrt(n()))

df2$Features<- gsub("JR_LR", "Junction\nreads", df2$Features)
df2$Features<- gsub("LocalCoverage_Pos1", "Local\ncoverage\nPos1", df2$Features)
df2$Features<- gsub("LocalCoverage_Pos2", "Local\ncoverage\nPos2", df2$Features)
df2$Features<- gsub("SP_LR", "Spanning\npairs", df2$Features)
df2$Features<- gsub("SVTypeDups","Dups", df2$Features)
df2$Features<- gsub("SVTypeInvs","Invs", df2$Features)
df2$Features<- gsub("SVTypeTrans","Trans", df2$Features)
p<- ggplot(df2, aes(x = Features, y = Mean))+
  geom_col()+
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width=0.2) +
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) +
  labs(y="% Feature importance")
p+ ggsave("LR_Feature_Imp.png", dpi=600, width = 6, height = 4)

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
p1 + ggsave("10XWGSModel_performance.png", width=5, height = 4)

### Final model trained on complete data set
train.control<- trainControl(method="repeatedcv",
                             number=10,
                             repeats=3,
                             classProbs=TRUE,
                             summaryFunction = twoClassSummary,
                             savePredictions= TRUE)
FinalModel <- train(PCR ~ GEM+JR_LR+SP_LR+SVType+Size+LocalCoverage_Pos1+LocalCoverage_Pos2,
                   data = ForTraining,
                   method = "glm", maximize=TRUE,metric="ROC",
                   trControl= train.control)
write_rds(FinalModel,"LR_Model.rds") # the traned model for 10XWGS technology

### predict all labels for remaining SVs
ApplyModel<- function(ForPrediction, Model, sample, tech, samplingType){
  DF_predict<- predict(Model, 
                       ForPrediction,
                       type = "prob")
  ForPrediction<- cbind(ForPrediction, pred_class=DF_predict$Positive)
  ForPrediction<- mutate(ForPrediction, PredictionLR=ifelse(Category=="Common", "Positive",
                                                                   ifelse(pred_class>0.6,"Positive","Negative")))
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
    theme(legend.position = "none",
          plot.title = element_text(size=10,hjust = 0.5, face = "bold"),
          axis.text = element_text(size=10),
          axis.title = element_text(size=12,face="bold")) +
    labs(x="Predicted SV", y="Counts")+  
    geom_text(data=filter(Count_Table, PredictionLR=="Positive"),
              aes(x= SVType, y= Count, label=paste(round(Percent,2),"%",sep='')),
              size=4, hjust=0.5, vjust=-0.2 ) +
    annotate("text", x=c(1,2,3,4), y=max(Table2$Count)+500, label=paste("N=",Table2$Count,sep=''))+
    ggtitle(paste(sample,":Combined model prediction for 10XWGS SVs", sep=''))
  p1<-p1 + ggsave(paste(sample,"_",tech,"_",samplingType,"_predictions_combined.png",sep=''), width=4.5, height = 5)
  return(ForPrediction$PredictionLR)  
}

NormalLR<-ApplyModel(prepareDF(Tog), FinalModel, sample = "MCF7", tech = "10XWGS", 
                     samplingType = "Normal")

### write prediction for all SV calls
write.table(cbind(Tog, `10XWGS_Prediction`=NormalLR), "10XWGS_Predictions.tsv",sep='\t',
                      quote=F, row.names = F)

