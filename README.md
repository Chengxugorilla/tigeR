## 1.Built-in Model

```
library(tigeR)
library(e1071)
library(pROC)
data(bayesmodel, package = "tigeR")
data(Stem.Sig, package = "tigeR")

test_Expr <- extract_mtr('MEL_GSE78220_exp')
test_Expr <- dataPreprocess(test_Expr, Stem.Sig)
test_response <- extract_label('MEL_GSE78220_meta')
predict_response_R <- predict(bayesmodel, t(test_Expr), type = 'class') == 'R'

#Obtaining the meta informations of patients whose prediction results are 'Response'.
data("MEL_GSE78220_meta")

rc <- MEL_GSE78220_meta[predict_response_R,]
rc$response <- sub('CR|MR|PR|SD|CRPR', 'R', rc$response)
rc$response <- sub('PD', 'NR', rc$response)
rc$response <- as.factor(rc$response)

#Drawing roc curve and calculating the AUC of roc curver.
roc1 <- roc(rc$overall.survival..days., response = rc$vital.status)
plot(roc1)
auc(roc1)

```

## 2.Naive Bayes Model

```
#Please load data set in Baidu cloud before running the code！
data(Stem.Sig,package = 'tigeR')

library(SummarizedExperiment)
SE1 <- MEL_GSE91061
SE2 <- MEL_phs000452
SE3 <- RCC_Braun_2020

SElist <- list(SE1, SE2, SE3)

#building model
library(tigeR)
mymodel <- build_model('NB', SElist, Stem.Sig, response_NR = TRUE)

##testing model
library(pROC)

#read tigeR Built-in datasets
library(magrittr)
extract_mtr('MEL_GSE78220_exp') %>% dataPreprocess(Stem.Sig, turn2HL = TRUE) -> test_Expr1
extract_mtr('MEL_PRJEB23709_exp') %>% dataPreprocess(Stem.Sig, turn2HL = TRUE) -> test_Expr2
feature <- intersect(rownames(test_Expr1),rownames(test_Expr2))
test_Expr <- cbind(test_Expr1[rownames(test_Expr1) %in% feature,], test_Expr2[rownames(test_Expr2) %in% feature,])

#Obtaining the meta informations
rc <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)

value <- as.numeric(predict(mymodel, t(test_Expr), type = 'raw')[,1])
#Drawing roc curve and calculating the AUC of roc curver.
roc1 <- roc(rc$response_NR, value)
plot(roc1)
auc(roc1)

```
## 3.Random Forest Model

```
#Please load data set in Baidu cloud before running the code！
data(Stem.Sig,package = 'tigeR')

library(SummarizedExperiment)
SE1 <- MEL_GSE91061
SE2 <- MEL_phs000452
SE3 <- RCC_Braun_2020

SElist <- list(SE1, SE2, SE3)

#building model
library(tigeR)
mymodel <- build_model('RF', SElist, Stem.Sig, rmBE = TRUE,response_NR = TRUE)

##testing model
library(pROC)

#read tigeR Built-in datasets
library(magrittr)
selected_gene <- rownames(mymodel$importance)
test_Expr1 <- extract_mtr('MEL_GSE78220_exp')
test_Expr1 <- test_Expr1[rownames(test_Expr1) %in% selected_gene,]
test_Expr2 <- extract_mtr('MEL_PRJEB23709_exp')
test_Expr2 <- test_Expr2[rownames(test_Expr2) %in% selected_gene,]
test_Expr <- cbind(test_Expr1, test_Expr2[rownames(test_Expr2) %in% rownames(test_Expr1),])

#Obtaining the meta informations
rc <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)

value <- as.numeric(predict(mymodel, t(test_Expr), type = 'vote')[,1])
#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
roc1 <- roc(rc$response_NR, value)
plot(roc1)
auc(roc1)

```
## 4.SVM Model

```
#Please load data set in Baidu cloud before running the code！
data(Stem.Sig, package = "tigeR")

library(SummarizedExperiment)
SE1 <- MEL_GSE91061
SE2 <- MEL_phs000452
SE3 <- RCC_Braun_2020

SElist <- list(SE1, SE2)

#building model
library(tigeR)
mymodel <- build_model('SVM', SElist, Stem.Sig, rmBE = FALSE,response_NR = TRUE)

##testing model
library(pROC)

#read tigeR Built-in datasets
library(magrittr)
selected_gene <- colnames(mymodel$SV)
test_Expr1 <- extract_mtr('MEL_GSE78220_exp')
test_Expr1 <- test_Expr1[rownames(test_Expr1) %in% selected_gene,]
test_Expr2 <- extract_mtr('MEL_PRJEB23709_exp')
test_Expr2 <- test_Expr2[rownames(test_Expr2) %in% selected_gene,]
test_Expr <- cbind(test_Expr1, test_Expr2[rownames(test_Expr2) %in% rownames(test_Expr1),])

#Obtaining the meta informations
rc <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)

value <- predict(mymodel, t(test_Expr),type='eps-regression')
#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
roc1 <- roc(rc$response, value)
plot(roc1)
auc(roc1)

```
## 5.Cancerclass Model

```
#Please load data set in Baidu cloud before running the code！
data(Stem.Sig, package = "tigeR")

library(SummarizedExperiment)
SE1 <- MEL_GSE91061
SE2 <- MEL_phs000452
SE3 <- RCC_Braun_2020

SElist <- list(SE1, SE2)

#building model
library(tigeR)
mymodel <- build_model('CC', SElist, Stem.Sig, rmBE = TRUE)

exprs <- cbind(extract_mtr('MEL_GSE78220_exp'), extract_mtr('MEL_PRJEB23709_exp'))
exprs <- dataPreprocess(exprs, Stem.Sig, turn2HL = FALSE)
pData <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)
rownames(pData) <- pData$sample_id
pData$class <- c(extract_label('MEL_GSE78220_meta'), extract_label('MEL_PRJEB23709_meta'))
identical(rownames(pData),colnames(exprs))
metadata <- data.frame(labelDescription = colnames(pData), row.names = colnames(pData))
adf <- new("AnnotatedDataFrame", data = as.data.frame(pData), varMetadata = metadata)

exampleSet2 <- new("ExpressionSet", exprs = exprs, phenoData = adf)
prediction <- cancerclass::predict(mymodel, exampleSet2, positive = "R", dist = "cor")
value <- as.numeric(prediction@prediction[,3])

library(pROC)
#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
roc <- roc(rc$response_NR, value)
plot(roc)
auc(roc)

```

## 6. Adaboost Model

```
#Please load data set in Baidu cloud before running the code！
data(Stem.Sig, package = "tigeR")

library(SummarizedExperiment)
SE1 <- MEL_GSE91061
SE2 <- MEL_phs000452
SE3 <- RCC_Braun_2020

SElist <- list(SE1, SE2)

#building model
library(tigeR)
mymodel <- build_model('ADB', SElist, Stem.Sig, rmBE = FALSE)

##testing model
library(pROC)

#read tigeR Built-in datasets
library(magrittr)
selected_gene <- names(mymodel$importance)
test_Expr1 <- extract_mtr('MEL_GSE78220_exp')
test_Expr1 <- test_Expr1[rownames(test_Expr1) %in% selected_gene,]
test_Expr2 <- extract_mtr('MEL_PRJEB23709_exp')
test_Expr2 <- test_Expr2[rownames(test_Expr2) %in% selected_gene,]
test_Expr <- cbind(test_Expr1, test_Expr2[rownames(test_Expr2) %in% rownames(test_Expr1),])

#Obtaining the meta informations
rc <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)

pred <- predict.boosting(mymodel, as.data.frame(t(test_Expr)))
value <- pred$votes[,1]

#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
roc1 <- roc(rc$response_NR, value)
plot(roc1)
auc(roc1)

```

## 7. Logitboost Model

```
#Please load data set in Baidu cloud before running the code！
data(Stem.Sig, package = "tigeR")

library(SummarizedExperiment)
SE1 <- MEL_GSE91061
SE2 <- MEL_phs000452
SE3 <- RCC_Braun_2020

SElist <- list(SE1, SE2)

#building model
library(tigeR)
mymodel <- build_model('LGB', SElist, Stem.Sig, rmBE = TRUE)

##testing model
library(pROC)

#read tigeR Built-in datasets
library(magrittr)
extract_mtr('MEL_GSE78220_exp') %>% dataPreprocess(Stem.Sig, turn2HL = TRUE) -> test_Expr1
extract_mtr('MEL_PRJEB23709_exp') %>% dataPreprocess(Stem.Sig, turn2HL = TRUE) -> test_Expr2
feature <- intersect(rownames(test_Expr1),rownames(test_Expr2))
test_Expr <- cbind(test_Expr1[rownames(test_Expr1) %in% feature,], test_Expr2[rownames(test_Expr2) %in% feature,])

#Obtaining the meta informations of patients whose prediction results are 'Response'.
rc <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)

#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
roc1 <- roc(rc$response, value)
plot(roc1)
auc(roc1)

```
## 8. Logistics Model

```
#Please load data set in Baidu cloud before running the code！
data(Stem.Sig, package = "tigeR")

library(SummarizedExperiment)
SE1 <- MEL_GSE91061
SE2 <- MEL_phs000452
SE3 <- RCC_Braun_2020

SElist <- list(SE1, SE2)

#building model
library(tigeR)
mymodel <- build_model('LGT', SElist, Stem.Sig[1:10], rmBE = FALSE, response_NR = TRUE)

##testing model
library(pROC)

#read tigeR Built-in datasets
library(magrittr)
extract_mtr('MEL_GSE78220_exp') %>% dataPreprocess(Stem.Sig[1:10], turn2HL = FALSE) -> test_Expr1
extract_mtr('MEL_PRJEB23709_exp') %>% dataPreprocess(Stem.Sig[1:10], turn2HL = FALSE) -> test_Expr2
test_Expr <- cbind(test_Expr1, test_Expr2[rownames(test_Expr2) %in% rownames(test_Expr1),])

#Obtaining the meta informations
rc <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)

df <- data.frame(t(max_min_normalization(test_Expr)))
value <- as.numeric(predict(mymodel, df, type = 'response'))
#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
roc1 <- roc(rc$response_NR, value)
plot(roc1)
auc(roc1)

```
## 9. Immunotherapy Response
```
library(tigeR)
SE <- MEL_GSE91061
Immunotherapy_Response(gene='CD274', SE)

```
## 10. Visualization
```
library(tigeR)

plt_diff('CD274',MEL_GSE91061,'Treatment') # Treatment vs UnTreatment
plt_diff('CD274',MEL_GSE91061,'Response') # Responder vs Non-Responder
plt_surv('CD274',MEL_GSE91061) # Survival analysis

```
## 11. Cibersort
```
library(tigeR)

data("MEL_GSE78220_exp")
mixture <- as.matrix(MEL_GSE78220_exp[,-1])
rownames(mixture) <- unlist(MEL_GSE78220_exp[,1])
data("LM22",package = "tigeR")

result <- CIBERSORT(LM22,mixture,perm=10, QN=T)

```