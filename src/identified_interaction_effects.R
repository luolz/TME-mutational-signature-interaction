#!/usr/bin Rscript
library(survival)
library(survminer)
library(dplyr)
load("./example/SBS_signature.RData")
load("./Kassandra_TME.RData")
load("./example/OS.RData")
head(SBS_signature)
head(infi)
head(survdata)

# step1 TME-mutational signature interaction model
mus<-SBS_signature
mus$cancer<-NULL
data <- left_join(survdata,mus, by = "sample")
TME<-infi
TME$cancer<-NULL
newdata <- left_join(data,TME, by = "sample")
newdata<-na.omit(newdata)
rownames(newdata)<- newdata$sample
newdata$time<-newdata$time/365
newdata$time<-newdata$time*12
infi_id <- colnames(infi)[3:ncol(infi)]
mus_id <- colnames(mus)[2:ncol(mus)]
cancer_id <- unique(newdata$cancer)
colnames(newdata)
result_interaction<-c()
j=1
for (j in 1:length(infi_id)) {
  k=1
  for (k in 1:length(mus_id)) {
    surv_0 <-newdata[,c("sample", "OS", "OS.time", "cancer", "Age", "Gender", "Stage",infi_id[j],mus_id[k])]
    names(surv_0)[2:3]<- c("status","time")
    surv_0$sample<-NULL
    surv_0$Gender<-as.numeric(surv_0$Gender)
    surv_0$Stage <-as.numeric(surv_0$Stage)
    i=1
    for (i in 1:length(cancer_id)) {
      surv <- surv_0[surv_0$cancer==cancer_id[i],]
      surv$cancer<-NULL
      surv<-na.omit(surv)
      head(surv)
      colnames(surv)[c(ncol(surv)-1,ncol(surv))]<- c("pivot","partner")
      surv$Interaction<- surv$partner*surv$pivot
      if(nrow(surv) >= 20) {
        death_rate = sum(surv[,1])/dim(surv)[1]
        if(death_rate >= 0.1){
          E<- surv$partner
          non_zero_sample =length(E[E!=0])
          if(non_zero_sample > 0.15*nrow(surv)){
            surv_1 = Surv(surv[,2], surv[,1])
            fit <- coxph(surv_1~.,data = surv[,-c(1,2)])
            summary(fit)$coef
            cancer=cancer_id[i]
            n=nrow(surv)
            reg.summary = summary(fit)$coef
            z = reg.summary["Interaction", c("z")]
            p = reg.summary["Interaction", c("Pr(>|z|)")]
            result<- data.frame(cancer_id=cancer_id[i],n,non_zero_sample,death_rate,infi_id=infi_id[j],mus_id=mus_id[k],z,p)
            result_interaction =  rbind(result_interaction,result)
          }
        }
      }
    }
  }
}

head(result_interaction)
# step2 log-rank test for significant interaction effects
sig<- result_interaction[result_interaction$p<0.05,]
sig$label<- paste0(sig$cancer_id,":",sig$infi_id,":",sig$mus_id)
infi_id <- unique(sig$infi_id )
mus_id <-  unique(sig$mus_id)
cancer_id <- unique(sig$cancer_id)
labels<- unique(sig$label)
std_error <- function(x) sd(x) / sqrt(length(x))
result_interaction_logrank <-c()
j=1
for (j in 1:length(infi_id)) {
  k=1
  for (k in 1:length(mus_id)) {
    i=1
    for (i in 1:length(cancer_id)) {
      surv_0 <-newdata[,c("sample", "OS", "OS.time", "cancer",infi_id[j],mus_id[k])]
      names(surv_0)[2:3]<- c("status","time")
      surv_0<- surv_0[surv_0$time>0,]
      surv_0$sample<-NULL
      surv <- surv_0[surv_0$cancer==cancer_id[i],]
      surv$cancer<-NULL
      surv<-na.omit(surv)
      head(surv)
      names(surv)[3:4]<- c("infi","mus")
      label <- paste0(cancer_id[i],":",infi_id[j],":",mus_id[k])
      num<- intersect(label,labels)
      if(length(num)>0){
        id<- mean(surv$infi)
        surv$infi_group[surv$infi>id]<- "High"
        surv$infi_group[surv$infi<=id]<- "Low"
        data <- surv$mus
        n1<-std_error(data)
        n2<-mean(surv$mus)
        id1<- n2+n1
        surv$mus_group[surv$mus>id1]<- "High"
        surv$mus_group[surv$mus<=id1]<- "Low"
        n<- nrow(surv)
        h1<-surv[surv$mus_group=="High",]
        head(h1)
        high_signature_n<- nrow(h1)
        high_signature_high_TME_n<- nrow(h1[h1$infi_group=="High",])
        high_signature_low_TME_n<- nrow(h1[h1$infi_group=="Low",])
        surv_diff  <- survdiff(Surv(time, status) ~ infi_group, data =  h1)
        high_signature_pvalue <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)

        L1<-surv[surv$mus_group=="Low",]
        low_signature_n<- nrow(L1)
        low_signature_high_TME_n<- nrow(L1[L1$infi_group=="High",])
        low_signature_low_TME_n<- nrow(L1[L1$infi_group=="Low",])
        surv_diff  <- survdiff(Surv(time, status) ~ infi_group, data =  L1)
        low_signature_pvalue <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
        result<- data.frame(i,j,k,label ,cancer=cancer_id[i],infi_id=infi_id[j],mus_id=mus_id[k],n,high_signature_n,high_signature_high_TME_n,high_signature_low_TME_n,low_signature_n,low_signature_high_TME_n,low_signature_low_TME_n,  high_signature_pvalue,low_signature_pvalue)
        result_interaction_logrank<-rbind(result_interaction_logrank,result)
      }
    }
  }
}
sig_high <- result_interaction_logrank[result_interaction_logrank$high_signature_pvalue<0.05 & result_interaction_logrank$low_signature_pvalue>0.05,]
sig_low<- result_interaction_logrank[result_interaction_logrank$high_signature_pvalue>0.05 & result_interaction_logrank$low_signature_pvalue<0.05,]
sig_final <- rbind(sig_high,sig_low)%>%as.data.frame()

write.table(result_interaction,file = "result_interaction.txt",row.names = F,col.names = T,quote = F,sep = "\t")
write.table(sig_final,file = "significant_result_interaction_logrank.txt",row.names = F,col.names = T,quote = F,sep = "\t")








