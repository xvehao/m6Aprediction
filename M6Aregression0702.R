rm(list=ls())
set.seed(10086)
num = 1000
train_size = 25
nflds = 5
library(readr)
library(glmnet)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(caret)
source(file = "mylasso.R")
ptm <- proc.time()

FPKM_methylation_input <- read.csv("FPKM_methylation_input.csv")
FPKM_methylation_IP <- read.csv("FPKM_methylation_IP.csv")
FPKM_methylation_IP = FPKM_methylation_IP[,-1]
FPKM_genes_input <- read.csv("FPKM_genes_input.csv")
FPKM_genes_input = FPKM_genes_input[,-1]
A = readRDS("AjacencyMatrix.rds")
A = forceSymmetric(A)
X = t(FPKM_methylation_input[, -1])
Y = t(FPKM_methylation_IP[, -1])
idx = colSums(X) != 0
X = X[, idx]
Y = scale(log2(Y[, idx]+.1))
A = A[idx,idx]
Yvar = apply(Y, 2, var)
selected_site = which(order(Yvar) > (ncol(Y) - num))
Y_s = Y[,selected_site]
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
X = scale(log2(X+.1))
X_train = X[1:train_size,]
X_test = X[(train_size+1):30,]
#A = A[1:num, 1:num]
#A[which(abs(X_cor)>0.8)] = 1
Y_train = Y_s[1:train_size,]
Y_test = Y_s[(train_size+1):30,]
flds <- createFolds(Y_train[,1], k = nflds, list = TRUE, returnTrain = FALSE)
Yvar = apply(Y, 2, var)
lasso = list()
lasso_pre_out = list()
lasso_pre_in = list()
ridge = list()
ridge_pre_in = list()
ridge_pre_out = list()
enet = list()
enet_pre_in = list()
enet_pre_out = list()
g = list()
g_cv = list()
g_cv_pre = list()
g_cv_error2 = list()
g_pre_in = list()
g_pre_out = list()
error2 = list()

for (i in 1:100){#ncol(Y_s)
  X_cor = cor(X_train[,features[i,]])
  ppnet= A[features[i,],features[i,]]
  ppnet[which(abs(X_cor)>0.8)] = 1
  diag(ppnet)=0
  ppnet = ppnet
  print(i)
  g_cv_min_mse = Inf
  myg_cv_min_mse = Inf
  lasso[[i]] = cv.glmnet(X_train[,features[i,]], Y_train[,i], alpha = 1, nfolds = nflds, grouped = FALSE, type.measure = "mse")
  ridge[[i]] = cv.glmnet(X_train[,features[i,]], Y_train[,i], alpha = 0, nfolds = nflds, grouped = FALSE, type.measure = "mse")
  enet[[i]] = cv.glmnet(X_train[,features[i,]], Y_train[,i], alpha = 0.5, nfolds = nflds, grouped = FALSE, type.measure = "mse")
  
  lasso_pre_in[[i]] = predict(lasso[[i]], newx = X_train[,features[i,]])
  lasso_pre_out[[i]] = predict(lasso[[i]], newx = X_test[,features[i,]])
  ridge_pre_in[[i]] = predict(ridge[[i]], newx = X_train[,features[i,]])
  ridge_pre_out[[i]] = predict(ridge[[i]], newx = X_test[,features[i,]])
  enet_pre_in[[i]] = predict(enet[[i]], newx = X_train[,features[i,]])
  enet_pre_out[[i]] = predict(enet[[i]], newx = X_test[,features[i,]])
  for (lambda in c(.01, .1, .5)){#c(.1, .2, .3, .4, .5, .6, .7, .8, .9)
    for (rho in c(.01, .1, .3)){#c(.01, .02, .03, .05, .1, .2, .3, .4, .5)
      print(paste0("rho:",rho))
      print(paste0("lambda:",lambda))
      
      for (alpha in seq(1,.1,-.3)){
        print(paste0("alpha:",alpha))
        carry_on = TRUE
        for (j in  1:nflds){
          if (carry_on == TRUE){
            cv_test_ind = flds[[j]]
            X_cv_train = X_train[-cv_test_ind, features[i,]]
            X_cv_test = X_train[cv_test_ind, features[i,]]
            Y_cv_train = Y_train[-cv_test_ind,i]
            Y_cv_test = Y_train[cv_test_ind,i]
            g_cv = myGIREN(X_cv_train, Y_cv_train, ppnet, lambda, rho*lambda, alpha, iter=500)
            if (sum(abs(g_cv)) >= 1e+6){
              carry_on = FALSE
              next
            }
            g_cv_pre[[j]] = X_cv_test%*%g_cv[-1] + g_cv[1]
            g_cv_error2[[j]] = (g_cv_pre[[j]] - Y_cv_test)^2
          }
        }
        if(carry_on==TRUE){
          g_cv_mse = mean(unlist(g_cv_error2))
          if (g_cv_mse < g_cv_min_mse){
            best_g_lambda = lambda
            best_g_rho = rho
            best_g_alpha = alpha
            g_cv_min_mse = g_cv_mse 
          }
        }
      }
    }
  }
  g[[i]] = myGIREN(X_train[,features[i,]], Y_train[,i], ppnet, best_g_lambda, best_g_rho*best_g_lambda, best_g_alpha)
  g_pre_in[[i]] = X_train[,features[i,]]%*%g[[i]][-1] + g[[i]][1]
  g_pre_out[[i]] = X_test[,features[i,]]%*%g[[i]][-1] + g[[i]][1]
}
lasso_error_in= (data.frame(lasso_pre_in) - Y_train[,1:length(lasso_pre_in)])^2
lasso_error_out= (data.frame(lasso_pre_out) - Y_test[,1:length(lasso_pre_out)])^2
ridge_error_in= (data.frame(ridge_pre_in) - Y_train[,1:length(ridge_pre_in)])^2
ridge_error_out= (data.frame(ridge_pre_out) - Y_test[,1:length(ridge_pre_out)])^2
enet_error_in= (data.frame(enet_pre_in) - Y_train[,1:length(enet_pre_in)])^2
enet_error_out= (data.frame(enet_pre_out) - Y_test[,1:length(enet_pre_out)])^2
g_error_in= (data.frame(g_pre_in) - Y_train[,1:length(g_pre_in)])^2
g_error_out= (data.frame(g_pre_out) - Y_test[,1:length(g_pre_out)])^2

lasso_mse_in = rowMeans(lasso_error_in)
lasso_mse_out = rowMeans(lasso_error_out)
ridge_mse_in = rowMeans(ridge_error_in)
ridge_mse_out = rowMeans(ridge_error_out)
enet_mse_in = rowMeans(enet_error_in)
enet_mse_out = rowMeans(enet_error_out)
g_mse_in = rowMeans(g_error_in)
g_mse_out = rowMeans(g_error_out)
mse_out = data.frame("lasso_outsample_mse" = lasso_mse_out,
                     "ridge_outsample_mse" = ridge_mse_out,
                     "enet_outsample_mse" = enet_mse_out, 
                     "g_outsample_mse" = g_mse_out)
medianse_out = data.frame("lasso_outsample_medianse" = apply(lasso_error_out,1, median),
                          "ridge_outsample_medianse" = apply(ridge_error_out,1, median),
                          "enet_outsample_medianse" = apply(enet_error_out,1, median), 
                          "g_outsample_medianse" = apply(g_error_out,1, median))
mse_in = data.frame("lasso_insample_mse" = lasso_mse_in,
                    "ridge_insample_mse" = ridge_mse_in,
                    "enet_insample_mse" = enet_mse_in,
                    "g_insample_mse" = g_mse_in)
medianse_in = data.frame("lasso_outsample_medianse" = apply(lasso_error_in,1, median),
                          "ridge_outsample_medianse" = apply(ridge_error_in,1, median),
                          "enet_outsample_medianse" = apply(enet_error_in,1, median), 
                          "g_outsample_medianse" = apply(g_error_in,1, median))
print(mse_out)
print(mse_in)
proc.time() - ptm

hist(as.numeric(lasso_error_out[1,]),breaks=25)
hist(as.numeric(g_error_out[1,]), breaks = 25, add=T,col=rgb(1,0,0,1/4))
lasso_n_coef=list()
g_n_coef=list()
myg_n_coef=list()
enet_n_coef=list()
for (i in 1:length(g)){
  lasso_n_coef[[i]] = sum(coef(lasso[[i]])!=0)
  g_n_coef[[i]] = sum(g[[i]]!=0)
  enet_n_coef[[i]] = sum(coef(enet[[i]])!=0)
}
myData = list()
myData$mse = mse
myData$lasso_n_coef = lasso_n_coef
myData$g_n_coef = g_n_coef
for(i in 1:nrow(X_test)){
  error2[[i]] = data.frame("lasso" = as.numeric(lasso_error_out[i,]),
                           "ridge" = as.numeric(ridge_error_out[i,]),
                           "EN" = as.numeric(enet_error_out[i,]),
                           "GIREN" = as.numeric(g_error_out[i,]))#"myg_error" = as.numeric(myg_error[1,])
  ggplot(melt(error2[[i]]), aes(x=variable, y=value, color=variable, fill=variable)) + 
    geom_boxplot(alpha=0.2) +
    labs(y = "MSE of each model", x=NULL) +
    theme(legend.position="none", text = element_text(size=25)) +
    ylim(0,1.5)
  #ggsave(paste0("msebox", i, ".pdf"), width = 5, height = 4, units = "in")
}
error2[[1]] = data.frame("lasso_error" = as.numeric(lasso_error_out[1,]),
                         "g_error" = as.numeric(g_error_out[1,]))#"myg_error" = as.numeric(myg_error[1,])
ggplot(melt(error2[[1]]), aes(x=variable, y=value)) + 
  geom_boxplot() +
  ylim(0,1)
ggplot(melt(error2[[2]]), aes(x=variable, y=value)) + 
  geom_boxplot()

#save.image(file = "0702.RData")

#write.csv(myData, file = "MyData0501.csv")

#Z = colMeans(X[1:25,])
#W = colMeans(Y[1:25,])
#linearMod = lm(W~Z)
#beta0 = coef(linearMod)[1]
#beta1 = coef(linearMod)[2]
#linearMod.pre = X[26:30,]*beta1 + beta0
#plot(X[26,], Y[26,])
#abline(linearMod)