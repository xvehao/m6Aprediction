# m6Aprediction
This directory provides the R code and data for predicting RNA methylation status from gene expression data using classification and regression methods. 

For more technical details about the scripts and data, please refer to:

Predicting RNA methylation status from gene expression data using classification and regression methods by Hao Xue, Zhen Wei, Kunqi Chen, Yujiao Tang, Xiangyu Wu, Jionglong Su and Jia Meng

ABSTRACT
Background: RNA N6-methyladenosine (m6A) has emerged as an important epigenetic modification for its role in regulating the stability, structure, processing and translation of RNA. Instability of m6A homeostasis may result in flaws in stem cell regulation, decrease in fertility and risk of cancer. To this day, experimental detection and quantification of RNA m6A modification are still time-consuming and labor-intensive. There is only a limited number of epitranscriptome samples in existing databases, and a matched RNA methylation profile is often not available for a biological problem of interests. Since gene expression data is usually readily available for most biological problems; it could be appealing if we can estimate the RNA methylation status from gene expression data using in silico methods. 

Results: In this study, we explored the possibility of computational prediction of RNA methylation status from the gene expression data using classification and regression methods based on mouse RNA methylation data collected from 73 experimental conditions. Elastic Net-regularized Logistic Regression (ENLR), Support Vector Machine (SVM) and Random Forests (RF) were constructed for classification. Both SVM and RF achieved the best performance with the mean area under the curve (AUC) 0.84 across samples, but SVM had a narrower AUC spread. Gene Site Enrichment Analysis was conducted on those sites selected by ENLR as predictors to access the biological significance of the model. Three functional annotation terms were found statistically significant: phosphoprotein, SRC Homology 3 (SH3) domain and endoplasmic reticulum (ER). All three terms were found to be closely related to m6A pathway. For regression analysis, Elastic Net (EN) was implemented, which yielded a mean Pearson Correlation Coefficient (PCC) 0.68 and a mean Spearman Correlation Coefficient (SCC) 0.64. Our exploratory study suggested that gene expression data could be used to construct predictors for m6A methylation status with adequate accuracy. 

Conclusion: Our work showed for the first time that RNA methylation status may be predicted from the matched gene expression data. This finding may facilitate RNA modification research in various biological contexts when a matched RNA methylation profile is not available, especially in the very early stage of the study.

Please contact Hao Xue haoxue9705@163.com if you have any questions regarding this package.
