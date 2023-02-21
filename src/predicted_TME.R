#!/usr/bin Rscript
library(glmnet)
library(dplyr)
load("./example/SBS_signature.RData")
load("./Kassandra_TME.RData")
load("./example/actual_lm_regression_result.RData")
head(SBS_signature)
head(infi)
head(final_result)
## Predicted TME
#step1 adjusted pvalue: fdr
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

#step2 lasso linear regression
lasso_infi<-infi
lasso_data<-list()
lasso_result<-c()
c2<-na.omit(final_result)
names(c2)[8:9]<-c("coef", "pvalue")
c2$pvalue[is.infinite(c2$pvalue)] <-1
head(c2)
u=1
for (u in 1:length(cancer_id)) {
  q1<- c2[c2$cancer_id%in%cancer_id[u],]
  q2<- q1[q1$pvalue<0.1,]
  infi_id_vo <- unique(q2$infi_id)
  w=1
  for (w in 1:length(infi_id_vo)) {
    qq<- q2[q2$infi_id%in%infi_id_vo[w],]
    mus_id <- unique(qq$mus_id)
    mus_num <- length(mus_id)
    if(mus_num>0){
      pre<- SBS_signature
      pre<-pre[pre$cancer%in%cancer_id[u],c("sample",mus_id)]
      pre_infi <- lasso_infi
      pre_infi <- pre_infi[pre$sample,c("sample",infi_id_vo[w])]
      infi_mus <- left_join(pre_infi,pre, by = "sample")
      infi_mus <- na.omit(infi_mus)
      rownames(infi_mus)<- infi_mus$sample
      U<- infi_mus[,3:ncol(infi_mus)]
      U<- as.matrix(U)
      if(length(U[U!=0])>10){
        if(mus_num==1){
          lm_fit <- lm(infi_mus[,2] ~ infi_mus[,3], data=infi_mus )
          h1<-summary(lm_fit)
          estimate_mus <-    h1$coefficients[2,1]
          Intercept <-    h1$coefficients[1,1]
          y_predicted <- estimate_mus*infi_mus[,3]+Intercept
          y<- infi_mus[,2]
          h<- cor.test(y,y_predicted)
          cor_p<- h$p.value
          cor<- h$estimate
          cvfit.min.out3<-mus_id
          lambda= ""
          model_type<-"lm"
        }

        if(mus_num>1){
          set.seed(1995)
          model_type<-"glmnet"
          cvfit <- cv.glmnet( as.matrix(infi_mus[,3:ncol(infi_mus)]),infi_mus[,2],family = "gaussian",type.measure="mse", alpha=1,nfolds = 10)
          fit<- glmnet(as.matrix(infi_mus[,3:ncol(infi_mus)]),infi_mus[,2], family = "gaussian",lambda = cvfit$lambda.min,alpha=1)
          y_predicted <- predict(fit, newx=as.matrix(infi_mus[,3:ncol(infi_mus)]),type="link", s=cvfit$lambda.min)
          y<- infi_mus[,2]
          h<- cor.test(y,y_predicted)
          cor_p<- h$p.value
          cor<- h$estimate
          coef.lambda.min <- coef(fit, s =cvfit$lambda.min)
          cvfit.min.out<- coef.lambda.min[which(coef.lambda.min !=0),]
          if(length(cvfit.min.out)>1){
            cvfit.min.out2 <- matrix(cvfit.min.out,length(cvfit.min.out),1)
            rownames(cvfit.min.out2)<- names(cvfit.min.out)
            cvfit.min.out3<- cvfit.min.out2[-1,]
            cvfit.min.out3<-cvfit.min.out3[order(abs(cvfit.min.out3),decreasing = T)]
            cvfit.min.out3<- names(cvfit.min.out3)
            cvfit.min.out3 <- paste0(cvfit.min.out3,collapse = ":")
          }else{
            cvfit.min.out3<-""
          }
          lambda=cvfit$lambda.min
        }
        cancer_sample<- nrow(infi_mus)
        name<- paste(infi_id_vo[w],cancer_id[u],sep=":")
        name1 <-paste(file_id[j],name,sep=":")
        result<- data.frame(signature_type="SBS",cancer_id=cancer_id[u],cancer_sample,infi_id=infi_id_vo[w],name1,name,mus_num,initial_features=paste0(mus_id,collapse = ":"),preserved_features=cvfit.min.out3,model_type, lambda=lambda,  cor, cor_p)
        lasso_result <- rbind(lasso_result,result)
        lasso_result <-as.data.frame(lasso_result)
        y_predicted<-as.numeric(y_predicted)
        W<- data.frame(sample=infi_mus$sample,  predict=y_predicted)
        head(W)
        Q<- left_join(infi_mus,W,by = "sample")
        lasso_data[[name1]] <- Q
      }
    }
  }
}
head(lasso_result)
# The prediction result of each cancer is stored in the list "lasso_data"
save(lasso_result, file = "./lasso_result.RData")

