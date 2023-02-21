### Linear regression
#!/usr/bin Rscript
library(openxlsx)
library(mltools)
library(data.table)
load("./example/SBS_signature.RData")
load("./Kassandra_TME.RData")
head(SBS_signature)
head(infi)

SBS_signature$rowsum<- rowSums(SBS_signature[,c(3:ncol(SBS_signature))])
SBS_signature$rowsum[SBS_signature$rowsum==0]<-1
SBS_signature[,c(3:ncol(SBS_signature))] <- SBS_signature[,c(3:ncol(SBS_signature))]/SBS_signature$rowsum
SBS_signature$rowsum<-NULL
SBS_signature<-na.omit(SBS_signature)
mus<- SBS_signature
mus$cancer<-NULL
infi_mus<- left_join( infi,mus,by = "sample")
total <- na.omit(infi_mus)
set.seed(555)
total$cancer <- as.factor(total$cancer)
newdata <- one_hot(as.data.table(total))
newdata <-as.data.frame(newdata)
infi_id<- colnames(infi)[3:ncol(infi)]
mus_id<- colnames(mus)[2:ncol(mus)]
cancer_id <-colnames(newdata)[grep("cancer",colnames(newdata))]
lm_result<-c()
i=1
for (i in 1:length(infi_id)) {
  j=1
  for (j in 1:length(mus_id)) {
    k=1
    for (k in 1:length(cancer_id)) {
      data <- newdata[,c(infi_id[i],mus_id[j],cancer_id[k])]
      colnames(data)<- c("infi","mus","cancer_id")
      head(data)
      cancer_sample<- nrow(data[data$cancer_id==1,])
      data$cancer_id<-factor(data$cancer_id)
      data$cancer_id <- relevel(data$cancer_id, ref = "0")
      E<- data$mus[data$cancer_id==1]
      non_zero_sample =length(E[E!=0])
      fit <- lm(infi ~ mus + cancer_id + mus:cancer_id, data=data )
      h1<-summary(fit)
      h2<-nrow(h1$coefficients)
      estimate_mus <-    h1$coefficients[2,1]
      mus_Std.Error <- h1$coefficients[2,2]
      mus_pvalue <- h1$coefficients[2,4]
      sign <- paste(mus_id[j],infi_id[i],sep = ":")
      oo<-  gsub("cancer_","",cancer_id[k])
      label<- paste(infi_id[i],oo,sep = ":")
      final<- data.frame(cancer_id=oo,mus_id=mus_id[j],infi_id=infi_id[i],sign,label, cancer_sample,non_zero_sample,estimate_mus,mus_pvalue)
      lm_result<- rbind(lm_result,final)
    }
    j=j+1
  }
  i=i+1
}
lm_result <-as.data.frame(lm_result)
lm_result<- lm_result[lm_result$non_zero_sample>0.15*lm_result$cancer_sample,]
write.xlsx(lm_result,file="./lm_result.xlsx", colNames = TRUE, borders = "columns")





