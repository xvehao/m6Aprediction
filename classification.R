rm(list=ls())
library(readr)
library(glmnet)
library(e1071)
library(caret)
library(kernlab)
library(randomForest)
library(pROC)
library(reshape2)

###initialize parameters and load data
set.seed(10086)
num = 1000 #specify window size
num_y = 100 #specify number of sites of interet
nflds = 5
ptm <- proc.time()
mouse_RPKM_input <- read.csv("mouse_RPKM_input2.csv")
mouse_RPKM_IP <- read.csv("mouse_RPKM_IP2.csv")
X = t(mouse_RPKM_input[, -1])
Y = t(mouse_RPKM_IP[, -1])
###
#use the code below to read human data instead
#input <- read.csv("FPKM_methylation_input.csv")
#IP <- read.csv("FPKM_methylation_IP.csv")
#IP = IP[,-1]
#X = t(input[, 2:31])
#Y = t(IP[, 2:31])
###
test_size = 15
train_size = nrow(X)-test_size 
idx = colSums(X) != 0 #exclude sites with identical zero vaues
X = log2(X[, idx]+.1) #log2 transformation
Y = log2(Y[, idx]+.1) 
M = Y-X
Y_scale = scale(Y)
Yvar = apply(Y_scale, 2, var) #select sites with highest variance in m6a level
selected_site = which(order(Yvar) > (ncol(Y) - 2*num_y))
features = matrix(0,nrow=num, ncol=num) #construct the corresponding gene expression level profile at the selected sites
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
colnames(X) = paste0("site", c(seq(ncol(X))))
X_train = X[1:train_size,]
X_test = X[(train_size+1):(train_size+test_size),]
lab = ifelse(M_select > 0, yes="yes", no="no") #assign label
lab_train = lab[1:train_size,]
lab_test = lab[(train_size+1):(train_size+test_size),]
logit = logit.auc = list()
SVM = SVM.auc = list()
rf = rf.auc = list()
ind=list()
i=1
k=0
tc <- trainControl(method="cv", number =5, summaryFunction=twoClassSummary, classProb=T) #5-fold cross-validation
while (i<=num_y) {
  print(paste0("now priting result of: ",i))
  k = k+1 #skip the sites with unbalanced labels
  ### training
  try({
    logit[[i]] = train(X_train[,features[k,]], factor(lab_train[,k]), 
                       method="glmnet", trControl=tc, metric =  "ROC")
    logit.auc[[i]] = max(logit[[i]]$results$ROC)
    print(paste0("auc of logit is: ", logit.auc[[i]]))
    SVM[[i]] = train(X_train[,features[k,]], factor(lab_train[,k]),
                     method="svmRadial", trControl=tc, metric =  "ROC")
    SVM.auc[[i]] = max(SVM[[i]]$results$ROC)
    print(paste0("auc of svm is: ", SVM.auc[[i]]))
    rf[[i]] = train(X_train[,features[k,]], factor(lab_train[,k]), 
                    method="rf", trControl=tc, metric =  "ROC")
    rf.auc[[i]] = max(rf[[i]]$results$ROC)
    print(paste0("auc of rf is: ", rf.auc[[i]]))
    ind[[i]] = k
    i=i+1
  }, silent=TRUE)
}
### compute in-sampel AUC
AUC = data.frame("logistic"=unlist(logit.auc),
                 "SVM"=unlist(SVM.auc),
                 "RF"=unlist(rf.auc))
ggplot(melt(AUC), aes(x=variable, y =value, color=variable, fill=variable)) + 
  theme(legend.position="none", text = element_text(size=25)) +
  geom_boxplot(alpha=0.2) +
  labs(y = "AUC of each classifier", x=NULL)
#ggsave("classificationboxplot.pdf", width = 5, height = 4, units = "in")

### compuete out-sampel AUC
logit.pre.auc=list()
logit.pre=matrix(0,nrow=test_size,ncol=num_y)
for (i in 1:num_y) {
  logit.pre[,i] = predict(logit[[i]], newdata = X_test[,features[ind[[i]],]], type = "prob")$yes
}
for (i in 1:(test_size)) {
  logit.pre.auc[[i]]=roc(lab_test[i,unlist(ind)], logit.pre[i,])$auc
}
SVM.pre.auc=list()
SVM.pre=matrix(0,nrow=(test_size),ncol=num_y)
for (i in 1:num_y) {
  SVM.pre[,i] = predict(SVM[[i]], newdata = X_test[,features[ind[[i]],]], type = "prob")$yes
}
for (i in 1:(test_size)) {
  SVM.pre.auc[[i]]=roc(lab_test[i,unlist(ind)], SVM.pre[i,])$auc
}
rf.pre.auc=list()
rf.pre=matrix(0,nrow=(test_size),ncol=num_y)
for (i in 1:num_y) {
  rf.pre[,i] = predict(rf[[i]], newdata = X_test[,features[ind[[i]],]], type = "prob")$yes
}
for (i in 1:(test_size)) {
  rf.pre.auc[[i]]=roc(lab_test[i,unlist(ind)], rf.pre[i,])$auc
}
AUC.pre = data.frame("ENLR"=unlist(logit.pre.auc),
                     "SVM"=unlist(SVM.pre.auc),
                     "RF"=unlist(rf.pre.auc))

### plot the result
ggplot(melt(AUC.pre), aes(x=variable, y =value, color=variable, fill=variable)) + 
  theme(legend.position="none", text = element_text(size=21.5)) +
  geom_boxplot(alpha=0.2) +
  labs(y = "AUC", x=NULL)+
  ggtitle("a")+
  theme(plot.title = element_text(hjust = 0.5))


# save.image(file = "0101classification.RData")
