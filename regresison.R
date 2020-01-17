rm(list=ls())
library(readr)
library(glmnet)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(caret)

###initialize parameters and load data
set.seed(10086)
num = 1000 #specify window size
num_y = 100 #specify number of sites of interet
nflds = 5
input <- read.csv("mouse_RPKM_input2.csv") 
IP <- read.csv("mouse_RPKM_IP2.csv")
X = t(input[, -1])
Y = t(IP[, -1])
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
idx2 = colSums(X) != 0 #exclude sites with identical zero vaues
idx1 = sample(1:(nrow(X))) #randomize samples
X = X[idx1, idx2]
colnames(X) = paste0("site", c(seq(ncol(X))))
Y = Y[idx1, idx2]
Yscaled = scale(log2(Y+.1)) #scale Y after log2 transformation
Yvar = apply(Y, 2, var)
selected_site = which(order(Yvar) > (ncol(Y) - num_y)) #select sites with highest variance in m6a level
Yselected = Yscaled[,selected_site]
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
Xscaled = scale(log2(X+.1)) #scale X after log2 transformation
Xtrain = Xscaled[1:train_size,]
Xtest = Xscaled[(train_size+1):(train_size+test_size),]
Ytrain = Yselected[1:train_size,]
Ytest = Yselected[(train_size+1):(train_size+test_size),]
flds <- createFolds(Ytrain[,1], k = nflds, list = TRUE, returnTrain = FALSE)
Yvar = apply(Y, 2, var)
lasso = lasso_pre_out = lasso_pre_in = list()
ridge = ridge_pre_in = ridge_pre_out = list()
en = en_pre_in = en_pre_out = list()
enTune = enTune_pre_in = enTune_pre_out = list()
error2 = list()

### train model
tc <- trainControl(method="cv", number =5) #5-fold cross validation
for (i in 1:num_y){
  penaltyfactor=rep(1,length(features[i,])) # do not penalize the coefficient of target site
  penaltyfactor[i]=0
  lasso[[i]] = cv.glmnet(Xtrain[,features[i,]], Ytrain[,i], 
                         alpha = 1, nfolds = nflds, grouped = FALSE, 
                         type.measure = "mse", penalty.factor=penaltyfactor)
  ridge[[i]] = cv.glmnet(Xtrain[,features[i,]], Ytrain[,i], 
                         alpha = 0, nfolds = nflds, grouped = FALSE, 
                         type.measure = "mse", penalty.factor=penaltyfactor)
  en[[i]] = train(Xtrain[,features[i,]], Ytrain[,i], #use default parameter searching grid
                  method = "glmnet",
                  trControl = tc,
                  metric = "RMSE")
  enetGrid <- expand.grid(lambda = ridge[[i]]$lambda[seq(1, length(ridge[[i]]$lambda), by = 5)], #use customized parameter searching grid
                          alpha = seq(0, 1, length = 10))
  enTune[[i]] <- train(Xtrain[,features[i,]], Ytrain[,i],
                       method = "glmnet",
                       tuneGrid = enetGrid,
                       trControl = tc,
                       metric = "RMSE")
### test model
  lasso_pre_in[[i]] = predict(lasso[[i]], newx = Xtrain[,features[i,]]) #predict with traning data
  lasso_pre_out[[i]] = predict(lasso[[i]], newx = Xtest[,features[i,]]) # predict with test data
  ridge_pre_in[[i]] = predict(ridge[[i]], newx = Xtrain[,features[i,]])
  ridge_pre_out[[i]] = predict(ridge[[i]], newx = Xtest[,features[i,]])
  en_pre_out[[i]] = predict(en[[i]], newdata = Xtest[,features[i,]])
  en_pre_in[[i]] = predict(en[[i]], newdata = Xtrain[,features[i,]])
  enTune_pre_out[[i]] = predict(enTune[[i]], newdata = Xtest[,features[i,]])
  enTune_pre_in[[i]] = predict(enTune[[i]], newdata = Xtrain[,features[i,]])
}
lasso_error_in= (data.frame(lasso_pre_in) - Ytrain)^2 #compute in-sample error
lasso_error_out= (data.frame(lasso_pre_out) - Ytest)^2 #compute out-sample error
ridge_error_in= (data.frame(ridge_pre_in) - Ytrain)^2
ridge_error_out= (data.frame(ridge_pre_out) - Ytest)^2
en_error_in= (data.frame(en_pre_in) - Ytest)^2
en_error_out= (data.frame(en_pre_out) - Ytest)^2
enTune_error_in= (data.frame(enTune_pre_in) - Ytest)^2
enTune_error_out= (data.frame(enTune_pre_out) - Ytest)^2

lasso_mse_in = rowMeans(lasso_error_in) #compute in-sample mse
lasso_mse_out = rowMeans(lasso_error_out) #compute out-sample mse
ridge_mse_in = rowMeans(ridge_error_in)
ridge_mse_out = rowMeans(ridge_error_out)
en_mse_in = rowMeans(en_error_in)
en_mse_out = rowMeans(en_error_out)
enTune_mse_in = rowMeans(enTune_error_in)
enTune_mse_out = rowMeans(enTune_error_out)

mse_out = data.frame("lasso_outsample_mse" = lasso_mse_out,
                     "ridge_outsample_mse" = ridge_mse_out,
                     "en_outsample_mse" = en_mse_out,
                     "enTune_outsample_mse" = enTune_mse_out)

print(mse_out)

### compute PCC and SCC between predicted value and atucal value
pearson_en=spearman_en=array()
for (i in 1:test_size){
  pearson_en[i] = cor(t(data.frame(en_pre_out)[i,]),Ytest[i,])
  spearman_en[i] = cor(t(data.frame(en_pre_out)[i,]),Ytest[i,],method = "spearman")
}
correlation = data.frame("PCC"=pearson_en,
                         "SCC"=spearman_en)

### plot the result
rownames(correlation)=rownames(mse_out)
ggplot(melt(correlation), aes(x=variable, y=value, color=variable, fill=variable)) + 
  geom_boxplot(alpha=0.2) +
  labs(y = NULL, x=NULL) +
  theme(legend.position="none", text = element_text(size=25)) +
  ylim(0.3,0.85)+
  ggtitle("b")+
  theme(plot.title = element_text(hjust = 0.5))

#save.image(file = "mouseregression.RData")
