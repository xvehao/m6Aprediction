# m6Aprediction
This directory provides the R code and data for predicting RNA methylation status from gene expression data using classification and regression methods. 

For more technical details about the scripts and data, please refer to:

Predicting RNA methylation status from gene expression data using classification and regression methods by Hao Xue, Zhen Wei, Kunqi Chen, Yujiao Tang, Xiangyu Wu, Jionglong Su and Jia Meng

Abstract 
Background: RNA N6-methyladenosine methylation (m6A) has emerged as an important epigenetic modification for its role in regulating the stability, structure, processing and translation of RNA. Instability of m6A homeostasis would result in flaws in stem cell regulation, decrease in fertility and risk of cancer. Till this day, experimental detection and quantification of RNA m6A modification are still time-consuming and labor-intensive. There are only limited number of epitranscriptome samples accumulated in existing databases, and a matched RNA methylation profile is often not available for a biological problem of interests. Because gene expression data is usually readily available for most biology problems, it could be appealing if it is possible to predict RNA methylation status from gene expression data using in silico methods. 

Results: In this study, we explored the possibility of computational prediction of RNA methylation status from gene expression data using classification and regression methods. For classification methods, we considered Elastic Net-regularized Logistic Regression (EN-LR), Support Vector Machine (SVM) and Random Forests (RF). For regression analysis, we considered Lasso, Ridge, Elastic Net (EN) and Gene Interaction-Regularized Elastic Net (GIREN). The RNA methylation data from 30 experimental conditions were collected and used for training (25 samples) and testing (5 samples), respectively. For methylation site prediction (classification analysis), the best performance was achieved using Random Forests (AUC = 0.82) While for methylation level prediction (regression analysis), the best performance was achieved using GIREN (median squared error = 0.28). Our exploratory study suggested that gene expression data could be used as to construct predictor for m6A methylation status. The scripts and the data used in this project are publically available from GitHub at: https://github.com/xvehao/m6Aprediction . 

The reads of input control and immunoprecipitated samples were stored, respectivly, in FPKM_methylation_input.csv and FPRKM_methylation_IP.csv. Scripts for classification and regression are contained, respectively, in M6AClassfication.R and M6ARegression.R. As for other files, one will find protein-protein interaction network in AjacencyMatrix.rds and functions implementing Gene Interaction-Regularized Elastic Net in mylasso.R  

Please contact Hao Xue haoxue9705@163.com if you have any questions regarding this package.
