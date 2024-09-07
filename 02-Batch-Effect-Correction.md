# 🖎🏻 Batch Effect Correction

## Batch Correction Method
|Name|Algorithm|Data type|
|:--:|:-------:|:-------:|
|Combat|empirical Bayes frameworks|cleaned and normalized data|
|Combat_seq|negative binomial regression|count data|
|limma|linear model|log transfered data|


```r
library(tigeR)
library(tigeR.data)
library(SummarizedExperiment)
SE <- cbind(MEL_GSE145996,MEL_GSE78220) 
SE_Combat <- remove_BE(SE,method = "Combat")
SE_Combat_seq <- remove_BE(SE,method = "Combat_seq")
SE_limma <- remove_BE(SE,method = "limma")
SE_DWD <- remove_BE(SE,method = "DWD")

SE_batch <- browse_BE(SE) + ggtitle("Origin")
Combat_batch <- browse_BE(SE_Combat) + ggtitle("Combat")
CombatSQ_batch <- browse_BE(SE_Combat_seq) + ggtitle("Combat_seq")
limma_batch <- browse_BE(SE_limma) + ggtitle("limma")
gridExtra::grid.arrange(grobs=list(SE_batch,Combat_batch,
                                   CombatSQ_batch,limma_batch))
```
<p align="center">
<img src="./figs/BatchEffect.svg" alt="SVG Image">
</p>