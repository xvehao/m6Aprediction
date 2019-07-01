rm(list=ls())
set.seed(10086)
num = 1000
train_size = 25
nflds = 5
library(readr)
library(glmnet)
library(e1071)
library(caret)
FPKM_methylation_input <- read.csv("FPKM_methylation_input.csv")
FPKM_methylation_IP <- read.csv("FPKM_methylation_IP.csv")
FPKM_methylation_IP = FPKM_methylation_IP[,-1]
FPKM_genes_input <- read.csv("FPKM_genes_input.csv")
FPKM_genes_input = FPKM_genes_input[,-1]
X = t(FPKM_methylation_input[, -1])
Y = t(FPKM_methylation_IP[, -1])
idx = colSums(X) != 0
X = log(X[1:30, idx]+.1)
Y = log(Y[1:30, idx]+.1)
M = Y-X
Y_scale = scale(Y)
Yvar = apply(Y_scale, 2, var)
selected_site = which(order(Yvar) > (ncol(Y) - num))
features = matrix(0,nrow=num, ncol=num)
for (i in 1:length(selected_site)){
  site = selected_site[i]
  if(site<num/2){
    features[i,]=1:num
  }
  else if(site>(ncol(X)-num/2)){
    features[i,]=(ncol(X)-num+1):ncol(X)
  }
  else{
    features[i,]=(site-num/2+1):(site+num/2) 
  }
}
Y_select = Y[,selected_site]
M_select = M[,selected_site]
#X = X[,1:num]
colnames(X) = paste0("site", c(seq(ncol(X))))
X_train = X[1:train_size,]
X_test = X[(train_size+1):30,]
lab = ifelse(M_select > 0, yes="yes", no="no")
#lab = ifelse(M_select > 0, yes=1, no=0)
lab_train = lab[1:train_size,]
lab_test = lab[(train_size+1):30,]
logit = list()
logit.auc = list()
SVM = list()
SVM.auc = list()
rf = list()
rf.auc = list()
ind=list()
i=1
k=0
tc <- trainControl(method="cv", number =5, summaryFunction=twoClassSummary, classProb=T)
while (i<=100) {
  print(paste0("now priting result of: ",i))
  k = k+1
  try({
    logit[[i]] = train(X_train[,features[k,]], factor(lab_train[,k]), method="glmnet", trControl=tc, metric =  "ROC")
    logit.auc[[i]] = max(logit[[i]]$results$ROC)
    print(paste0("auc of logit is: ", logit.auc[[i]]))
    SVM[[i]] = train(X_train[,features[k,]], factor(lab_train[,k]), method="svmRadial", trControl=tc, metric =  "ROC")
    SVM.auc[[i]] = max(SVM[[i]]$results$ROC)
    print(paste0("auc of svm is: ", SVM.auc[[i]]))
    rf[[i]] = train(X_train[,features[k,]], factor(lab_train[,k]), method="rf", trControl=tc, metric =  "ROC")
    rf.auc[[i]] = max(rf[[i]]$results$ROC)
    print(paste0("auc of rf is: ", rf.auc[[i]]))
    ind[[i]] = k
    i=i+1
  }, silent=TRUE)
}
AUC = data.frame("logistic"=unlist(logit.auc),
                 "SVM"=unlist(SVM.auc),
                 "RF"=unlist(rf.auc))
ggplot(melt(AUC), aes(x=variable, y =value, color=variable, fill=variable)) + 
  theme(legend.position="none") +
  geom_boxplot(alpha=0.2) +
  labs(y = "AUC of each classifier", x=NULL)
save.image(file = "0701.RData")

logit.pre.auc=list()
logit.pre=matrix(0,nrow=(30-train_size),ncol=100)
for (i in 1:length(ind)) {
  logit.pre[,i] = predict(logit[[i]], newdata = X_test[,features[ind[[i]],]], type = "prob")$yes
}
for (i in 1:length(ind)) {
  logit.pre.auc[[i]]=roc(lab_test[i,unlist(ind)], logit.pre[i,])$auc
}
#roc(lab_test[i,unlist(ind)], logit.pre[i,], plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE)

SVM.pre.auc=list()
SVM.pre=matrix(0,nrow=(30-train_size),ncol=100)
for (i in 1:length(ind)) {
  SVM.pre[,i] = predict(SVM[[i]], newdata = X_test[,features[ind[[i]],]], type = "prob")$yes
}
for (i in 1:length(ind)) {
  SVM.pre.auc[[i]]=roc(lab_test[i,unlist(ind)], SVM.pre[i,])$auc
}
rf.pre.auc=list()
rf.pre=matrix(0,nrow=(30-train_size),ncol=100)
for (i in 1:length(ind)) {
  rf.pre[,i] = predict(rf[[i]], newdata = X_test[,features[ind[[i]],]], type = "prob")$yes
}
for (i in 1:(30-train_size)) {
  rf.pre.auc[[i]]=roc(lab_test[i,unlist(ind)], rf.pre[i,])$auc
}
results = data.frame("logistic test auc"=unlist(logit.pre.auc),
                    "SVM test auc"=unlist(SVM.pre.auc),
                    "RF test auc"=unlist(rf.pre.auc))
