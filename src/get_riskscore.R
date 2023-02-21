#!/usr/bin Rscript
library(survival)
library(survminer)
library(glmnet)
library(dplyr)
load("./example/SBS_signature.RData")
load("./Kassandra_TME.RData")
load("./example/OS.RData")
load("./example/actual_lm_regression_result.RData")
load("./example/TME_significant_survival.RData")

head(SBS_signature)
head(infi)
head(survdata)
head(final_result)
head(TME_significant_survival)

## Riskscore
#step1
cancer_id <- unique(final_result$cancer_id)
q=1
c2<-c()
for (q in 1:length(cancer_id)) {
  c1<- final_result[final_result$cancer_id%in%cancer_id[q],]
  c1$mus_pvalue[is.na(c1$mus_pvalue)]<- 1
  c1$mus_pvalue<-p.adjust(c1$mus_pvalue,"BH")
  c2<-rbind(c2,c1)
}
final_result<-c2
final_result<-final_result[final_result$mus_pvalue<0.1,]
target_infi<- unique(TME_significant_survival$infi_id)
target_mus<- final_result$mus_id[final_result$infi_id%in%target_infi]
target_mus<- unique(target_mus)
mus<-SBS_signature[,c("sample",target_mus)]
surv <- left_join(survdata,mus, by = "sample")
rownames(surv)<- surv$sample
names(surv)[2:3]<- c("status","time")
surv$time<-surv$time/365
surv$time<-surv$time*12
surv<-na.omit(surv)

#step2
risk_score_data <-list()
cancer_id <- unique(surv$cancer)
i=1
for (i in 1:length(cancer_id)) {
  surv_1 <- surv[surv$cancer==cancer_id[i],]
  final_result_1 <- final_result[final_result$cancer_id==cancer_id[i],]
  mus_target<- unique(final_result_1$mus_id)
  y <- Surv(surv_1$time,surv_1$status)
  dat<-surv_1[, colnames(surv_1)%in%mus_target]
  x<- as.matrix(dat)
  library(glmnet)
  if(length(mus_target)>1){
    cvfit<- cv.glmnet(x, y, family = "cox", type.measure = "deviance", nfolds = 10)
    coef.lambda.min <- coef(cvfit, s =cvfit$lambda.min)
    cvfit.min.out <- names(coef.lambda.min[which(coef.lambda.min !=0),])}else{
      cvfit.min.out <- mus_target
    }
  if(length(cvfit.min.out)!=0){
    target <- cvfit.min.out
  }else{
    target <-mus_target
  }
  a<- surv_1[,c("time","status",target)]
  fit <- coxph(Surv(time,status)~.,a)
  ho <- fit$linear.predictors
  cancer=cancer_id[i]
  N=nrow(a)
  h1<- data.frame(sample=surv_1$sample,cancer=cancer,N,time=surv_1$time,status=surv_1$status,risk_score= ho)
  rownames(h1)<-h1$sample
  res.cut <- surv_cutpoint(h1,
                           time = "time",
                           event = "status",
                           variables = c("risk_score"))
  h2<-summary(res.cut)
  id<- h2[1,1]
  h1$group[h1$risk_score>id]<- "high"
  h1$group[h1$risk_score<=id]<- "low"
  h1$surv_cutpoint<-id
  h1$mus<- paste(target,collapse = ",")
  risk_score_data[[cancer]]<-h1
}
# The risk score of each cancer is stored in the list "risk_score_data"
save(risk_score_data, file = "./risk_score_data.RData")
