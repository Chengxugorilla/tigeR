# tigeR
`tigeR` is a R package for Tumor Immunotherapy Gene Expression Analysis
![Screenshot](https://github.com/Chengxugorilla/tigeR/raw/main/man/figures/logo.png)

TIGER web server(http://tiger.canceromics.org/#/)
## Install
devtools::install_github("Chengxugorilla/tigeR")

## Quick Start
`tigeR` generally supports the quantification and visualization of tumor immunotherapy data. 

### 1. Load packages and demo data
The demo data include 4 Melanoma RNA-seq datasets and 1 Renal Cell Carcinoma dataset. Use `Dataloader()` to get an overview of all datasets.
```
library(tigeR)
Dataloader(c(4,5,13,14,18))
```
### 2. Calculate signature scores of existing immunotherapy prognosis signatures
```{r}
Sig_scores <- calculate_Signature_Score(assay(MEL_GSE78220))
View(Sig_scores)
```
### 3. Build machine learning model for immunotherapy prognosis prediction
    tigeR allows users to build machine learning prediction models for immunotherapy prognosis. There 7 model including Naive-Bayes Model, Random Forest Model, SVM Model, Cancerclass Model, Adaboost Model, Logitboost Model, Logistics Model.
    We use Naive-Bayes Model to make an example:
```
SElist <- list(MEL_GSE91061, MEL_phs000452, RCC_Braun_2020)

#building model
mymodel <- build_Model('NB', SElist, Stem.Sig, response_NR = TRUE)

#read tigeR Built-in datasets
library(magrittr)
extract_mtr('MEL_GSE78220') %>% dataPreprocess(Stem.Sig, turn2HL = TRUE) -> test_Expr1
extract_mtr('MEL_PRJEB23709') %>% dataPreprocess(Stem.Sig, turn2HL = TRUE) -> test_Expr2
feature <- intersect(rownames(test_Expr1),rownames(test_Expr2))

test_Expr <- cbind(test_Expr1[rownames(test_Expr1) %in% feature,], test_Expr2[rownames(test_Expr2) %in% feature,])
response <- extract_label(c('MEL_GSE78220','MEL_PRJEB23709'))
value <- as.numeric(predict(mymodel, t(test_Expr), type = 'raw')[,1])

#Drawing roc curve and calculating the AUC of roc curver.
library(pROC)
ROC <- roc(response, value)
plot(ROC)
auc(ROC)
```
### 4. Immunotherapy Response
```
library(tigeR)

DataLoader(5)
Immunotherapy_Response(gene='CD274', MEL_GSE91061)

```
### 5. Visualization
```
library(tigeR)

DataLoader(5)
plt_diff('CD274',MEL_GSE91061,'Treatment') # Treatment vs UnTreatment
plt_diff('CD274',MEL_GSE91061,'Response') # Responder vs Non-Responder
plt_surv('CD274',MEL_GSE91061) # Survival analysis

```
### 6. Cibersort
```
library(tigeR)

Dataloader(4)
mixture <- as.matrix(MEL_GSE78220_exp[,-1])
rownames(mixture) <- unlist(MEL_GSE78220_exp[,1])
data("LM22",package = "tigeR")

result <- CIBERSORT(LM22,mixture,perm=10, QN=T)

```
