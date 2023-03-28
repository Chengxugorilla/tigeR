# tigeR
## 1.Built-in Model

```
library(tigeR)
library(e1071)
library(pROC)
data(bayesmodel, package = "tigeR")
Stem.Sig <- c("ENO1","RBM17","RHCG","H1FX","UPF2","C19orf12","TXNDC5","SERPINA1","HNRNPH2","PHB2","LHX2","NNAT","CHCHD5","FAM204A","IFT20","COPA","UQCRH","KRT1","S100A7","RYK","RPP38","NUDT19","RNF11","VPS35","CEACAM6","PTMA","FOXG1","AC018890.6","CISD2","PCBP2","SMIM26","ACKR1","SERBP1","MYL6","GSDMC","HIST1H2AE","MSRB2","LRP3","EEF2","MTRNR2L1","C18ORF32","SELENOW","CTDNEP1","KCTD8","SPINK1","LECT1","MMP24OS","MADCAM1","BCL10","SLC25A3","LY6D","TPD52L1","EIF5AL1","LSM14A","SUMO3","SPHK1","TMED5","SLC25A6","ALKBH5","OOEP","RENBP","VRK1","ABCA7","FABP3","DENND2C","CNIH1","CALML5","PDCD2","MARVELD1","GPI","SHC4","HM13","DERL2","UBB","PIPOX","STEAP1B","MRPL23","CIB2","PVR","GPR137B","S100A6","SERF2","DSG1","NDUFA4","LDHB","GRAMD1A","RRS1","LBP","VPS29","YBX1","TRIM47","SLC26A3","NTS","PKMYT1","SELENOM","HTR2B","JTB","FEM1B","KRTDAP","NUPR2","PFDN5","SYNE4","SND1","MYDGF","ATP5G2","C11orf31","SMOX","MAMDC2","CLN3","RPA1","DESI1","ZNF330","LRRFIP2","ALDOA","CITED4","SHFM1","GTSF1","FAM98C","PRDX1","C19orf70","C9ORF78","ATP5B","C20orf24","PLEKHA7","PPP4C","ALYREF","PIGP","EIF4EBP3","TMA7","PFN1","ANKRD35","ATP5J2","PRSS27","DYRK1B","CFH","UBA52","POLR2E","MID1IP1","RNF126","SIK3","VCX2","SF3A2","UFD1L","SNX8","TMEM41A","PSMB3","NBPF14","C7orf73","TCEB2","EGLN2","C4BPA","TMPRSS3","ZFAND6","KNG1","SH3GL1","GGACT","ADAMTSL1","ASF1B","UBE2L3","RGS20","ATP5I","KRT10","MRPS21","ARHGEF35","BCL7C","AP2S1","AGT","RACK1","COMMD8","OMG","RAD23A","C14orf166","DKK2","DDA1","EIF5A2","MLANA","H2AFZ","H3F3B","C1orf43","FAM104B","IRX3","PIH1D1","KYNU","EIF3M","EMC3","SLC25A29","SSBP4","ALKBH7","EIF2S3","FUZ","PNCK","LRRN4CL","UBE2D2","ACTG1","APOA1BP","NGFRAP1","IRX5","UBE2S","ATP5G3","AGR2","C5ORF15","SNF8","KXD1","ZNF90","EMX2","PDXP","PTHLH","TM7SF3","NPM1","RIOK3","TAGLN2","VGLL1","CALB2","RBX1","UGT2B15","MTHFD2","RBM3","INAFM1","CCDC106","CDHR1","EPHA3","LMF2","AC007325.4","GNPTAB","TUBB","RALY","PIGC","SFRP1","KIF1C","DSCR8","GC","SUB1","SET","LINC00844","IGSF21","LMF1","EPHA7","NCAPH2","AL353997.6","ZNF106","PERP","TICAM1","SNRPE","ENY2","C1QBP","HNRNPA1","PPIA","C11ORF58","SUMO2","TXNL4A","PLPPR5","SHC3","NTN5","AC011043.1","MIR205HG","TSPAN10","ACTB","C19orf43","PPP2R5A","AARD","TMEM256","CHCHD2","IGFBP1","ERP29","DERL1","TMPRSS5","FAM89A","EHD2","WIF1","IGKC","NECTIN1","BIRC7","AC090498.1","TIAM1","H3F3A","SLURP1","EIF4A1","TSPYL5","PON3","TMED10","GRTP1","GLRX5","GRID2","POLE4","CENPU","SNORD100","ECM1","CLPTM1","PLP2","GRHL3","C1orf131","C8orf82","KRT19","AL365205.1","MTRNR2L10","RER1","C6ORF62","ATP5J","GNB2L1","POM121","SMC2","HOMER1","C1orf56","BAIAP2L2","EFNB1","PKP1","CEBPZOS","TPM2","ATP5G1","WDR77","DIAPH2","AP3S1","CASP6","ANKRD6","KLRC2","SERPINA3","GTSE1","CPEB4","SPRR1B","SLC25A5","GPR87","COX7A2L","PRKCDBP","NME1","NDC1","FGL1","COPB1","BECN1","COX4I1","NKD1","STEAP3","PCP4","CDK13","SPRR2D","PDP1","ADGRF4","PEX13","FTH1","HN1","FZD8","TXN","EEF1A1","PSMF1","PFKFB2","CYP27B1","GNB1","HPCA","TAC1","SFTPB","CDKN2A","FAM83B","TCF7L1","EEF1G","FAM110A","TP53BP2","ORM1","FAU","ANP32B","CRYGS","C1ORF61","PARL","TMEM18","BEX5","FABP6","CDKN2B","CALML3","RPIA","FAM89B","SNRPB","RSL1D1","ORM2","CCNDBP1","ATP5F1B","NUDT4","CMTM5","G6PD","POU4F3","TPH1","MUC21","UBAP1","HCAR2","FAHD2A","CCS","COMMD7","EIF2S2","C5","AATF","HINT1","SPOCD1","ATP23","FAAP20","NRM","AGRN","SLC22A31","SEC61B","ZNF750","DBI","DHCR7","BPI","DDOST","GAPDH","NOX1","HNRNPD","KCNF1","TTC9B","NR0B2","UBE2A","CGN","MRPL54","ZBTB43","SERPINB2","NAA50","DGAT2","PSMA7","MKI67","TPI1","CHMP4A","LYPLA1","OSBPL6","C1ORF122","APH1A","MCTS1","SFTA2","NAPSA","PPP1CA","TINCR","ROPN1","BARX2","MRGBP","AC012640.4","COX6A1","GSTA1","MTRNR2L2","MTHFD2L","C17ORF89","CRP","IDH3G","GCLC","AL031727.1","GSTP1","IL36G","ROPN1B","ZMYND11","POP4","PRDX3","RAN","GTF3A","NACA","ARSJ","LINC01268","APOA2","TIMM23","FAM83A","QARS")

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
Stem.Sig <- c("ENO1","RBM17","RHCG","H1FX","UPF2","C19orf12","TXNDC5","SERPINA1","HNRNPH2","PHB2","LHX2","NNAT","CHCHD5","FAM204A","IFT20","COPA","UQCRH","KRT1","S100A7","RYK","RPP38","NUDT19","RNF11","VPS35","CEACAM6","PTMA","FOXG1","AC018890.6","CISD2","PCBP2","SMIM26","ACKR1","SERBP1","MYL6","GSDMC","HIST1H2AE","MSRB2","LRP3","EEF2","MTRNR2L1","C18ORF32","SELENOW","CTDNEP1","KCTD8","SPINK1","LECT1","MMP24OS","MADCAM1","BCL10","SLC25A3","LY6D","TPD52L1","EIF5AL1","LSM14A","SUMO3","SPHK1","TMED5","SLC25A6","ALKBH5","OOEP","RENBP","VRK1","ABCA7","FABP3","DENND2C","CNIH1","CALML5","PDCD2","MARVELD1","GPI","SHC4","HM13","DERL2","UBB","PIPOX","STEAP1B","MRPL23","CIB2","PVR","GPR137B","S100A6","SERF2","DSG1","NDUFA4","LDHB","GRAMD1A","RRS1","LBP","VPS29","YBX1","TRIM47","SLC26A3","NTS","PKMYT1","SELENOM","HTR2B","JTB","FEM1B","KRTDAP","NUPR2","PFDN5","SYNE4","SND1","MYDGF","ATP5G2","C11orf31","SMOX","MAMDC2","CLN3","RPA1","DESI1","ZNF330","LRRFIP2","ALDOA","CITED4","SHFM1","GTSF1","FAM98C","PRDX1","C19orf70","C9ORF78","ATP5B","C20orf24","PLEKHA7","PPP4C","ALYREF","PIGP","EIF4EBP3","TMA7","PFN1","ANKRD35","ATP5J2","PRSS27","DYRK1B","CFH","UBA52","POLR2E","MID1IP1","RNF126","SIK3","VCX2","SF3A2","UFD1L","SNX8","TMEM41A","PSMB3","NBPF14","C7orf73","TCEB2","EGLN2","C4BPA","TMPRSS3","ZFAND6","KNG1","SH3GL1","GGACT","ADAMTSL1","ASF1B","UBE2L3","RGS20","ATP5I","KRT10","MRPS21","ARHGEF35","BCL7C","AP2S1","AGT","RACK1","COMMD8","OMG","RAD23A","C14orf166","DKK2","DDA1","EIF5A2","MLANA","H2AFZ","H3F3B","C1orf43","FAM104B","IRX3","PIH1D1","KYNU","EIF3M","EMC3","SLC25A29","SSBP4","ALKBH7","EIF2S3","FUZ","PNCK","LRRN4CL","UBE2D2","ACTG1","APOA1BP","NGFRAP1","IRX5","UBE2S","ATP5G3","AGR2","C5ORF15","SNF8","KXD1","ZNF90","EMX2","PDXP","PTHLH","TM7SF3","NPM1","RIOK3","TAGLN2","VGLL1","CALB2","RBX1","UGT2B15","MTHFD2","RBM3","INAFM1","CCDC106","CDHR1","EPHA3","LMF2","AC007325.4","GNPTAB","TUBB","RALY","PIGC","SFRP1","KIF1C","DSCR8","GC","SUB1","SET","LINC00844","IGSF21","LMF1","EPHA7","NCAPH2","AL353997.6","ZNF106","PERP","TICAM1","SNRPE","ENY2","C1QBP","HNRNPA1","PPIA","C11ORF58","SUMO2","TXNL4A","PLPPR5","SHC3","NTN5","AC011043.1","MIR205HG","TSPAN10","ACTB","C19orf43","PPP2R5A","AARD","TMEM256","CHCHD2","IGFBP1","ERP29","DERL1","TMPRSS5","FAM89A","EHD2","WIF1","IGKC","NECTIN1","BIRC7","AC090498.1","TIAM1","H3F3A","SLURP1","EIF4A1","TSPYL5","PON3","TMED10","GRTP1","GLRX5","GRID2","POLE4","CENPU","SNORD100","ECM1","CLPTM1","PLP2","GRHL3","C1orf131","C8orf82","KRT19","AL365205.1","MTRNR2L10","RER1","C6ORF62","ATP5J","GNB2L1","POM121","SMC2","HOMER1","C1orf56","BAIAP2L2","EFNB1","PKP1","CEBPZOS","TPM2","ATP5G1","WDR77","DIAPH2","AP3S1","CASP6","ANKRD6","KLRC2","SERPINA3","GTSE1","CPEB4","SPRR1B","SLC25A5","GPR87","COX7A2L","PRKCDBP","NME1","NDC1","FGL1","COPB1","BECN1","COX4I1","NKD1","STEAP3","PCP4","CDK13","SPRR2D","PDP1","ADGRF4","PEX13","FTH1","HN1","FZD8","TXN","EEF1A1","PSMF1","PFKFB2","CYP27B1","GNB1","HPCA","TAC1","SFTPB","CDKN2A","FAM83B","TCF7L1","EEF1G","FAM110A","TP53BP2","ORM1","FAU","ANP32B","CRYGS","C1ORF61","PARL","TMEM18","BEX5","FABP6","CDKN2B","CALML3","RPIA","FAM89B","SNRPB","RSL1D1","ORM2","CCNDBP1","ATP5F1B","NUDT4","CMTM5","G6PD","POU4F3","TPH1","MUC21","UBAP1","HCAR2","FAHD2A","CCS","COMMD7","EIF2S2","C5","AATF","HINT1","SPOCD1","ATP23","FAAP20","NRM","AGRN","SLC22A31","SEC61B","ZNF750","DBI","DHCR7","BPI","DDOST","GAPDH","NOX1","HNRNPD","KCNF1","TTC9B","NR0B2","UBE2A","CGN","MRPL54","ZBTB43","SERPINB2","NAA50","DGAT2","PSMA7","MKI67","TPI1","CHMP4A","LYPLA1","OSBPL6","C1ORF122","APH1A","MCTS1","SFTA2","NAPSA","PPP1CA","TINCR","ROPN1","BARX2","MRGBP","AC012640.4","COX6A1","GSTA1","MTRNR2L2","MTHFD2L","C17ORF89","CRP","IDH3G","GCLC","AL031727.1","GSTP1","IL36G","ROPN1B","ZMYND11","POP4","PRDX3","RAN","GTF3A","NACA","ARSJ","LINC01268","APOA2","TIMM23","FAM83A","QARS")

#standardization of expression matrix
exp1 <- as.matrix(MEL_GSE91061_exp[,-1])
exp1 <- apply(exp1, 2, as.numeric)
rownames(exp1) <- MEL_GSE91061_exp[,1]

exp2 <- as.matrix(MEL_phs000452_exp[,-1])
exp2 <- apply(exp2, 2, as.numeric)
rownames(exp2) <- MEL_phs000452_exp[,1]

meta1 <- MEL_GSE91061_meta
meta2 <- MEL_phs000452_meta

library(tigeR)
response1 <- meta1$response
response1 <- response_standardize(response1)
response2 <- meta2$response
response2 <- response_standardize(response2)

#obtain index of NA samples
filt1 <- grep('UNK',response1)
filt2 <- grep('UNK',response2)

#remove UNK. exp2 has no UNK
exp1 <- exp1[,-filt1]
response1 <- response1[-filt1]
#exp2 <- exp2[,-filt1]
#response2 <- response2[-filt2]

library(SummarizedExperiment)
colData1 <- DataFrame(response = response1)
SE1 <- SummarizedExperiment(assays = exp1,
                            colData = colData1)

colData2 <- DataFrame(response = response2)
SE2 <- SummarizedExperiment(assays = exp2,
                            colData = colData2)

SElist <- list(SE1, SE2)

#building model
mymodel <- build_NB_model(SElist, Stem.Sig, rmBE = TRUE)

##testing model
library(pROC)

#read tigeR Built-in datasets
test_Expr <- extract_mtr('MEL_GSE78220_exp')
test_Expr <- dataPreprocess(test_Expr, Stem.Sig)
test_response <- extract_label('MEL_GSE78220_meta')

#the index of sample which prediction result is Responder
predict_response_R <- predict(mymodel, t(test_Expr), type = 'class') == 'R'

#Obtaining the meta informations of patients whose prediction results are 'Response'.
data("MEL_GSE78220_meta")

rc <- MEL_GSE78220_meta[predict_response_R,]
rc$response <- sub('CR|MR|PR|CRPR', 'R', rc$response)
rc$response <- sub('PD|SD', 'NR', rc$response)
rc$response <- as.factor(rc$response)

#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
roc1 <- roc(rc$overall.survival..days., response = rc$vital.status)
plot(roc1)
auc(roc1)

```
## 3.Random Forest Model

```
#Please load data set in Baidu cloud before running the code！
Stem.Sig <- c("ENO1","RBM17","RHCG","H1FX","UPF2","C19orf12","TXNDC5","SERPINA1","HNRNPH2","PHB2","LHX2","NNAT","CHCHD5","FAM204A","IFT20","COPA","UQCRH","KRT1","S100A7","RYK","RPP38","NUDT19","RNF11","VPS35","CEACAM6","PTMA","FOXG1","AC018890.6","CISD2","PCBP2","SMIM26","ACKR1","SERBP1","MYL6","GSDMC","HIST1H2AE","MSRB2","LRP3","EEF2","MTRNR2L1","C18ORF32","SELENOW","CTDNEP1","KCTD8","SPINK1","LECT1","MMP24OS","MADCAM1","BCL10","SLC25A3","LY6D","TPD52L1","EIF5AL1","LSM14A","SUMO3","SPHK1","TMED5","SLC25A6","ALKBH5","OOEP","RENBP","VRK1","ABCA7","FABP3","DENND2C","CNIH1","CALML5","PDCD2","MARVELD1","GPI","SHC4","HM13","DERL2","UBB","PIPOX","STEAP1B","MRPL23","CIB2","PVR","GPR137B","S100A6","SERF2","DSG1","NDUFA4","LDHB","GRAMD1A","RRS1","LBP","VPS29","YBX1","TRIM47","SLC26A3","NTS","PKMYT1","SELENOM","HTR2B","JTB","FEM1B","KRTDAP","NUPR2","PFDN5","SYNE4","SND1","MYDGF","ATP5G2","C11orf31","SMOX","MAMDC2","CLN3","RPA1","DESI1","ZNF330","LRRFIP2","ALDOA","CITED4","SHFM1","GTSF1","FAM98C","PRDX1","C19orf70","C9ORF78","ATP5B","C20orf24","PLEKHA7","PPP4C","ALYREF","PIGP","EIF4EBP3","TMA7","PFN1","ANKRD35","ATP5J2","PRSS27","DYRK1B","CFH","UBA52","POLR2E","MID1IP1","RNF126","SIK3","VCX2","SF3A2","UFD1L","SNX8","TMEM41A","PSMB3","NBPF14","C7orf73","TCEB2","EGLN2","C4BPA","TMPRSS3","ZFAND6","KNG1","SH3GL1","GGACT","ADAMTSL1","ASF1B","UBE2L3","RGS20","ATP5I","KRT10","MRPS21","ARHGEF35","BCL7C","AP2S1","AGT","RACK1","COMMD8","OMG","RAD23A","C14orf166","DKK2","DDA1","EIF5A2","MLANA","H2AFZ","H3F3B","C1orf43","FAM104B","IRX3","PIH1D1","KYNU","EIF3M","EMC3","SLC25A29","SSBP4","ALKBH7","EIF2S3","FUZ","PNCK","LRRN4CL","UBE2D2","ACTG1","APOA1BP","NGFRAP1","IRX5","UBE2S","ATP5G3","AGR2","C5ORF15","SNF8","KXD1","ZNF90","EMX2","PDXP","PTHLH","TM7SF3","NPM1","RIOK3","TAGLN2","VGLL1","CALB2","RBX1","UGT2B15","MTHFD2","RBM3","INAFM1","CCDC106","CDHR1","EPHA3","LMF2","AC007325.4","GNPTAB","TUBB","RALY","PIGC","SFRP1","KIF1C","DSCR8","GC","SUB1","SET","LINC00844","IGSF21","LMF1","EPHA7","NCAPH2","AL353997.6","ZNF106","PERP","TICAM1","SNRPE","ENY2","C1QBP","HNRNPA1","PPIA","C11ORF58","SUMO2","TXNL4A","PLPPR5","SHC3","NTN5","AC011043.1","MIR205HG","TSPAN10","ACTB","C19orf43","PPP2R5A","AARD","TMEM256","CHCHD2","IGFBP1","ERP29","DERL1","TMPRSS5","FAM89A","EHD2","WIF1","IGKC","NECTIN1","BIRC7","AC090498.1","TIAM1","H3F3A","SLURP1","EIF4A1","TSPYL5","PON3","TMED10","GRTP1","GLRX5","GRID2","POLE4","CENPU","SNORD100","ECM1","CLPTM1","PLP2","GRHL3","C1orf131","C8orf82","KRT19","AL365205.1","MTRNR2L10","RER1","C6ORF62","ATP5J","GNB2L1","POM121","SMC2","HOMER1","C1orf56","BAIAP2L2","EFNB1","PKP1","CEBPZOS","TPM2","ATP5G1","WDR77","DIAPH2","AP3S1","CASP6","ANKRD6","KLRC2","SERPINA3","GTSE1","CPEB4","SPRR1B","SLC25A5","GPR87","COX7A2L","PRKCDBP","NME1","NDC1","FGL1","COPB1","BECN1","COX4I1","NKD1","STEAP3","PCP4","CDK13","SPRR2D","PDP1","ADGRF4","PEX13","FTH1","HN1","FZD8","TXN","EEF1A1","PSMF1","PFKFB2","CYP27B1","GNB1","HPCA","TAC1","SFTPB","CDKN2A","FAM83B","TCF7L1","EEF1G","FAM110A","TP53BP2","ORM1","FAU","ANP32B","CRYGS","C1ORF61","PARL","TMEM18","BEX5","FABP6","CDKN2B","CALML3","RPIA","FAM89B","SNRPB","RSL1D1","ORM2","CCNDBP1","ATP5F1B","NUDT4","CMTM5","G6PD","POU4F3","TPH1","MUC21","UBAP1","HCAR2","FAHD2A","CCS","COMMD7","EIF2S2","C5","AATF","HINT1","SPOCD1","ATP23","FAAP20","NRM","AGRN","SLC22A31","SEC61B","ZNF750","DBI","DHCR7","BPI","DDOST","GAPDH","NOX1","HNRNPD","KCNF1","TTC9B","NR0B2","UBE2A","CGN","MRPL54","ZBTB43","SERPINB2","NAA50","DGAT2","PSMA7","MKI67","TPI1","CHMP4A","LYPLA1","OSBPL6","C1ORF122","APH1A","MCTS1","SFTA2","NAPSA","PPP1CA","TINCR","ROPN1","BARX2","MRGBP","AC012640.4","COX6A1","GSTA1","MTRNR2L2","MTHFD2L","C17ORF89","CRP","IDH3G","GCLC","AL031727.1","GSTP1","IL36G","ROPN1B","ZMYND11","POP4","PRDX3","RAN","GTF3A","NACA","ARSJ","LINC01268","APOA2","TIMM23","FAM83A","QARS")

#standardization of expression matrix
exp1 <- as.matrix(MEL_GSE91061_exp[,-1])
exp1 <- apply(exp1, 2, as.numeric)
rownames(exp1) <- MEL_GSE91061_exp[,1]

exp2 <- as.matrix(MEL_phs000452_exp[,-1])
exp2 <- apply(exp2, 2, as.numeric)
rownames(exp2) <- MEL_phs000452_exp[,1]

meta1 <- MEL_GSE91061_meta
meta2 <- MEL_phs000452_meta

library(tigeR)
response1 <- meta1$response
response1 <- response_standardize(response1)
response2 <- meta2$response
response2 <- response_standardize(response2)

#obtain index of NA samples
filt1 <- grep('UNK',response1)
filt2 <- grep('UNK',response2)

#remove UNK. exp2 has no UNK
exp1 <- exp1[,-filt1]
response1 <- response1[-filt1]
#exp2 <- exp2[,-filt1]
#response2 <- response2[-filt2]

library(SummarizedExperiment)
colData1 <- DataFrame(response = response1)
SE1 <- SummarizedExperiment(assays = exp1,
                            colData = colData1)

colData2 <- DataFrame(response = response2)
SE2 <- SummarizedExperiment(assays = exp2,
                            colData = colData2)

SElist <- list(SE1, SE2)

#building model
mymodel <- build_RF_model(SElist, Stem.Sig, rmBE = TRUE)

##testing model
library(pROC)

#read tigeR Built-in datasets
test_Expr <- cbind(extract_mtr('MEL_GSE78220_exp'), extract_mtr('MEL_PRJEB23709_exp'))
test_Expr <- dataPreprocess(test_Expr, Stem.Sig, turn2HL = FALSE)
test_response <- rbind(extract_label('MEL_GSE78220_meta'), extract_label('MEL_PRJEB23709_meta'))

#the index of sample which prediction result is Responder
predict_response_R <- predict(mymodel, t(test_Expr), type = 'class') == 'R'

#Obtaining the meta informations of patients whose prediction results are 'Response'.
data("MEL_GSE78220_meta")
data("MEL_PRJEB23709_meta")

rc <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)[predict_response_R,]
rc$response <- sub('CR|MR|PR|CRPR', 'R', rc$response)
rc$response <- sub('PD|SD', 'NR', rc$response)
rc$response <- as.factor(rc$response)

#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
roc1 <- roc(rc$overall.survival..days., response = rc$vital.status)

rc2 <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)
rc2$response <- sub('CR|MR|PR|CRPR', 'R', rc2$response)
rc2$response <- sub('PD|SD', 'NR', rc2$response)
rc2$response <- as.factor(rc2$response)

roc2 <- roc(rc2$overall.survival..days., response = rc2$vital.status)

par(mfrow = c(1,2))
plot(roc1)
plot(roc2)

print("The AUC of selected patient is:")
auc(roc1)
print("The AUC of all patient is:")
auc(roc2)

```
## 4.SVM Model

```
#Please load data set in Baidu cloud before running the code！
Stem.Sig <- c("ENO1","RBM17","RHCG","H1FX","UPF2","C19orf12","TXNDC5","SERPINA1","HNRNPH2","PHB2","LHX2","NNAT","CHCHD5","FAM204A","IFT20","COPA","UQCRH","KRT1","S100A7","RYK","RPP38","NUDT19","RNF11","VPS35","CEACAM6","PTMA","FOXG1","AC018890.6","CISD2","PCBP2","SMIM26","ACKR1","SERBP1","MYL6","GSDMC","HIST1H2AE","MSRB2","LRP3","EEF2","MTRNR2L1","C18ORF32","SELENOW","CTDNEP1","KCTD8","SPINK1","LECT1","MMP24OS","MADCAM1","BCL10","SLC25A3","LY6D","TPD52L1","EIF5AL1","LSM14A","SUMO3","SPHK1","TMED5","SLC25A6","ALKBH5","OOEP","RENBP","VRK1","ABCA7","FABP3","DENND2C","CNIH1","CALML5","PDCD2","MARVELD1","GPI","SHC4","HM13","DERL2","UBB","PIPOX","STEAP1B","MRPL23","CIB2","PVR","GPR137B","S100A6","SERF2","DSG1","NDUFA4","LDHB","GRAMD1A","RRS1","LBP","VPS29","YBX1","TRIM47","SLC26A3","NTS","PKMYT1","SELENOM","HTR2B","JTB","FEM1B","KRTDAP","NUPR2","PFDN5","SYNE4","SND1","MYDGF","ATP5G2","C11orf31","SMOX","MAMDC2","CLN3","RPA1","DESI1","ZNF330","LRRFIP2","ALDOA","CITED4","SHFM1","GTSF1","FAM98C","PRDX1","C19orf70","C9ORF78","ATP5B","C20orf24","PLEKHA7","PPP4C","ALYREF","PIGP","EIF4EBP3","TMA7","PFN1","ANKRD35","ATP5J2","PRSS27","DYRK1B","CFH","UBA52","POLR2E","MID1IP1","RNF126","SIK3","VCX2","SF3A2","UFD1L","SNX8","TMEM41A","PSMB3","NBPF14","C7orf73","TCEB2","EGLN2","C4BPA","TMPRSS3","ZFAND6","KNG1","SH3GL1","GGACT","ADAMTSL1","ASF1B","UBE2L3","RGS20","ATP5I","KRT10","MRPS21","ARHGEF35","BCL7C","AP2S1","AGT","RACK1","COMMD8","OMG","RAD23A","C14orf166","DKK2","DDA1","EIF5A2","MLANA","H2AFZ","H3F3B","C1orf43","FAM104B","IRX3","PIH1D1","KYNU","EIF3M","EMC3","SLC25A29","SSBP4","ALKBH7","EIF2S3","FUZ","PNCK","LRRN4CL","UBE2D2","ACTG1","APOA1BP","NGFRAP1","IRX5","UBE2S","ATP5G3","AGR2","C5ORF15","SNF8","KXD1","ZNF90","EMX2","PDXP","PTHLH","TM7SF3","NPM1","RIOK3","TAGLN2","VGLL1","CALB2","RBX1","UGT2B15","MTHFD2","RBM3","INAFM1","CCDC106","CDHR1","EPHA3","LMF2","AC007325.4","GNPTAB","TUBB","RALY","PIGC","SFRP1","KIF1C","DSCR8","GC","SUB1","SET","LINC00844","IGSF21","LMF1","EPHA7","NCAPH2","AL353997.6","ZNF106","PERP","TICAM1","SNRPE","ENY2","C1QBP","HNRNPA1","PPIA","C11ORF58","SUMO2","TXNL4A","PLPPR5","SHC3","NTN5","AC011043.1","MIR205HG","TSPAN10","ACTB","C19orf43","PPP2R5A","AARD","TMEM256","CHCHD2","IGFBP1","ERP29","DERL1","TMPRSS5","FAM89A","EHD2","WIF1","IGKC","NECTIN1","BIRC7","AC090498.1","TIAM1","H3F3A","SLURP1","EIF4A1","TSPYL5","PON3","TMED10","GRTP1","GLRX5","GRID2","POLE4","CENPU","SNORD100","ECM1","CLPTM1","PLP2","GRHL3","C1orf131","C8orf82","KRT19","AL365205.1","MTRNR2L10","RER1","C6ORF62","ATP5J","GNB2L1","POM121","SMC2","HOMER1","C1orf56","BAIAP2L2","EFNB1","PKP1","CEBPZOS","TPM2","ATP5G1","WDR77","DIAPH2","AP3S1","CASP6","ANKRD6","KLRC2","SERPINA3","GTSE1","CPEB4","SPRR1B","SLC25A5","GPR87","COX7A2L","PRKCDBP","NME1","NDC1","FGL1","COPB1","BECN1","COX4I1","NKD1","STEAP3","PCP4","CDK13","SPRR2D","PDP1","ADGRF4","PEX13","FTH1","HN1","FZD8","TXN","EEF1A1","PSMF1","PFKFB2","CYP27B1","GNB1","HPCA","TAC1","SFTPB","CDKN2A","FAM83B","TCF7L1","EEF1G","FAM110A","TP53BP2","ORM1","FAU","ANP32B","CRYGS","C1ORF61","PARL","TMEM18","BEX5","FABP6","CDKN2B","CALML3","RPIA","FAM89B","SNRPB","RSL1D1","ORM2","CCNDBP1","ATP5F1B","NUDT4","CMTM5","G6PD","POU4F3","TPH1","MUC21","UBAP1","HCAR2","FAHD2A","CCS","COMMD7","EIF2S2","C5","AATF","HINT1","SPOCD1","ATP23","FAAP20","NRM","AGRN","SLC22A31","SEC61B","ZNF750","DBI","DHCR7","BPI","DDOST","GAPDH","NOX1","HNRNPD","KCNF1","TTC9B","NR0B2","UBE2A","CGN","MRPL54","ZBTB43","SERPINB2","NAA50","DGAT2","PSMA7","MKI67","TPI1","CHMP4A","LYPLA1","OSBPL6","C1ORF122","APH1A","MCTS1","SFTA2","NAPSA","PPP1CA","TINCR","ROPN1","BARX2","MRGBP","AC012640.4","COX6A1","GSTA1","MTRNR2L2","MTHFD2L","C17ORF89","CRP","IDH3G","GCLC","AL031727.1","GSTP1","IL36G","ROPN1B","ZMYND11","POP4","PRDX3","RAN","GTF3A","NACA","ARSJ","LINC01268","APOA2","TIMM23","FAM83A","QARS")

#standardization of expression matrix
exp1 <- as.matrix(MEL_GSE91061_exp[,-1])
exp1 <- apply(exp1, 2, as.numeric)
rownames(exp1) <- MEL_GSE91061_exp[,1]

exp2 <- as.matrix(MEL_phs000452_exp[,-1])
exp2 <- apply(exp2, 2, as.numeric)
rownames(exp2) <- MEL_phs000452_exp[,1]

meta1 <- MEL_GSE91061_meta
meta2 <- MEL_phs000452_meta

library(tigeR)
response1 <- meta1$response
response1 <- response_standardize(response1)
response2 <- meta2$response
response2 <- response_standardize(response2)

#obtain index of NA samples
filt1 <- grep('UNK',response1)
filt2 <- grep('UNK',response2)

#remove UNK. exp2 has no UNK
exp1 <- exp1[,-filt1]
response1 <- response1[-filt1]
#exp2 <- exp2[,-filt1]
#response2 <- response2[-filt2]

library(SummarizedExperiment)
colData1 <- DataFrame(response = response1)
SE1 <- SummarizedExperiment(assays = exp1,
                            colData = colData1)

colData2 <- DataFrame(response = response2)
SE2 <- SummarizedExperiment(assays = exp2,
                            colData = colData2)

SElist <- list(SE1, SE2)

#building model
mymodel <- build_SVM_model(SElist, Stem.Sig, rmBE = TRUE)

##testing model
library(pROC)

#read tigeR Built-in datasets
test_Expr <- extract_mtr('MEL_GSE78220_exp')
test_Expr <- dataPreprocess(test_Expr, Stem.Sig, turn2HL = FALSE)
test_response <- extract_label('MEL_GSE78220_meta')

#the index of sample which prediction result is Responder
predict_response_R <- predict(mymodel, t(test_Expr), type = 'class') == 'R'

#Obtaining the meta informations of patients whose prediction results are 'Response'.
data("MEL_GSE78220_meta")

rc <- MEL_GSE78220_meta[predict_response_R,]
rc$response <- sub('CR|MR|PR|CRPR', 'R', rc$response)
rc$response <- sub('PD|SD', 'NR', rc$response)
rc$response <- as.factor(rc$response)

#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
roc1 <- roc(rc$overall.survival..days., response = rc$vital.status)
plot(roc1)
auc(roc1)

```

## 5.Cancerclass Model

```
#Please load data set in Baidu cloud before running the code！
Stem.Sig <- c("ENO1","RBM17","RHCG","H1FX","UPF2","C19orf12","TXNDC5","SERPINA1","HNRNPH2","PHB2","LHX2","NNAT","CHCHD5","FAM204A","IFT20","COPA","UQCRH","KRT1","S100A7","RYK","RPP38","NUDT19","RNF11","VPS35","CEACAM6","PTMA","FOXG1","AC018890.6","CISD2","PCBP2","SMIM26","ACKR1","SERBP1","MYL6","GSDMC","HIST1H2AE","MSRB2","LRP3","EEF2","MTRNR2L1","C18ORF32","SELENOW","CTDNEP1","KCTD8","SPINK1","LECT1","MMP24OS","MADCAM1","BCL10","SLC25A3","LY6D","TPD52L1","EIF5AL1","LSM14A","SUMO3","SPHK1","TMED5","SLC25A6","ALKBH5","OOEP","RENBP","VRK1","ABCA7","FABP3","DENND2C","CNIH1","CALML5","PDCD2","MARVELD1","GPI","SHC4","HM13","DERL2","UBB","PIPOX","STEAP1B","MRPL23","CIB2","PVR","GPR137B","S100A6","SERF2","DSG1","NDUFA4","LDHB","GRAMD1A","RRS1","LBP","VPS29","YBX1","TRIM47","SLC26A3","NTS","PKMYT1","SELENOM","HTR2B","JTB","FEM1B","KRTDAP","NUPR2","PFDN5","SYNE4","SND1","MYDGF","ATP5G2","C11orf31","SMOX","MAMDC2","CLN3","RPA1","DESI1","ZNF330","LRRFIP2","ALDOA","CITED4","SHFM1","GTSF1","FAM98C","PRDX1","C19orf70","C9ORF78","ATP5B","C20orf24","PLEKHA7","PPP4C","ALYREF","PIGP","EIF4EBP3","TMA7","PFN1","ANKRD35","ATP5J2","PRSS27","DYRK1B","CFH","UBA52","POLR2E","MID1IP1","RNF126","SIK3","VCX2","SF3A2","UFD1L","SNX8","TMEM41A","PSMB3","NBPF14","C7orf73","TCEB2","EGLN2","C4BPA","TMPRSS3","ZFAND6","KNG1","SH3GL1","GGACT","ADAMTSL1","ASF1B","UBE2L3","RGS20","ATP5I","KRT10","MRPS21","ARHGEF35","BCL7C","AP2S1","AGT","RACK1","COMMD8","OMG","RAD23A","C14orf166","DKK2","DDA1","EIF5A2","MLANA","H2AFZ","H3F3B","C1orf43","FAM104B","IRX3","PIH1D1","KYNU","EIF3M","EMC3","SLC25A29","SSBP4","ALKBH7","EIF2S3","FUZ","PNCK","LRRN4CL","UBE2D2","ACTG1","APOA1BP","NGFRAP1","IRX5","UBE2S","ATP5G3","AGR2","C5ORF15","SNF8","KXD1","ZNF90","EMX2","PDXP","PTHLH","TM7SF3","NPM1","RIOK3","TAGLN2","VGLL1","CALB2","RBX1","UGT2B15","MTHFD2","RBM3","INAFM1","CCDC106","CDHR1","EPHA3","LMF2","AC007325.4","GNPTAB","TUBB","RALY","PIGC","SFRP1","KIF1C","DSCR8","GC","SUB1","SET","LINC00844","IGSF21","LMF1","EPHA7","NCAPH2","AL353997.6","ZNF106","PERP","TICAM1","SNRPE","ENY2","C1QBP","HNRNPA1","PPIA","C11ORF58","SUMO2","TXNL4A","PLPPR5","SHC3","NTN5","AC011043.1","MIR205HG","TSPAN10","ACTB","C19orf43","PPP2R5A","AARD","TMEM256","CHCHD2","IGFBP1","ERP29","DERL1","TMPRSS5","FAM89A","EHD2","WIF1","IGKC","NECTIN1","BIRC7","AC090498.1","TIAM1","H3F3A","SLURP1","EIF4A1","TSPYL5","PON3","TMED10","GRTP1","GLRX5","GRID2","POLE4","CENPU","SNORD100","ECM1","CLPTM1","PLP2","GRHL3","C1orf131","C8orf82","KRT19","AL365205.1","MTRNR2L10","RER1","C6ORF62","ATP5J","GNB2L1","POM121","SMC2","HOMER1","C1orf56","BAIAP2L2","EFNB1","PKP1","CEBPZOS","TPM2","ATP5G1","WDR77","DIAPH2","AP3S1","CASP6","ANKRD6","KLRC2","SERPINA3","GTSE1","CPEB4","SPRR1B","SLC25A5","GPR87","COX7A2L","PRKCDBP","NME1","NDC1","FGL1","COPB1","BECN1","COX4I1","NKD1","STEAP3","PCP4","CDK13","SPRR2D","PDP1","ADGRF4","PEX13","FTH1","HN1","FZD8","TXN","EEF1A1","PSMF1","PFKFB2","CYP27B1","GNB1","HPCA","TAC1","SFTPB","CDKN2A","FAM83B","TCF7L1","EEF1G","FAM110A","TP53BP2","ORM1","FAU","ANP32B","CRYGS","C1ORF61","PARL","TMEM18","BEX5","FABP6","CDKN2B","CALML3","RPIA","FAM89B","SNRPB","RSL1D1","ORM2","CCNDBP1","ATP5F1B","NUDT4","CMTM5","G6PD","POU4F3","TPH1","MUC21","UBAP1","HCAR2","FAHD2A","CCS","COMMD7","EIF2S2","C5","AATF","HINT1","SPOCD1","ATP23","FAAP20","NRM","AGRN","SLC22A31","SEC61B","ZNF750","DBI","DHCR7","BPI","DDOST","GAPDH","NOX1","HNRNPD","KCNF1","TTC9B","NR0B2","UBE2A","CGN","MRPL54","ZBTB43","SERPINB2","NAA50","DGAT2","PSMA7","MKI67","TPI1","CHMP4A","LYPLA1","OSBPL6","C1ORF122","APH1A","MCTS1","SFTA2","NAPSA","PPP1CA","TINCR","ROPN1","BARX2","MRGBP","AC012640.4","COX6A1","GSTA1","MTRNR2L2","MTHFD2L","C17ORF89","CRP","IDH3G","GCLC","AL031727.1","GSTP1","IL36G","ROPN1B","ZMYND11","POP4","PRDX3","RAN","GTF3A","NACA","ARSJ","LINC01268","APOA2","TIMM23","FAM83A","QARS")

#standardization of expression matrix
exp1 <- as.matrix(MEL_GSE91061_exp[,-1])
exp1 <- apply(exp1, 2, as.numeric)
rownames(exp1) <- MEL_GSE91061_exp[,1]

exp2 <- as.matrix(MEL_phs000452_exp[,-1])
exp2 <- apply(exp2, 2, as.numeric)
rownames(exp2) <- MEL_phs000452_exp[,1]

meta1 <- MEL_GSE91061_meta
meta2 <- MEL_phs000452_meta

library(tigeR)
response1 <- meta1$response
response1 <- response_standardize(response1)
response2 <- meta2$response
response2 <- response_standardize(response2)

#obtain index of NA samples
filt1 <- grep('UNK',response1)
filt2 <- grep('UNK',response2)

#remove UNK. exp2 has no UNK
exp1 <- exp1[,-filt1]
response1 <- response1[-filt1]
#exp2 <- exp2[,-filt1]
#response2 <- response2[-filt2]

library(SummarizedExperiment)
colData1 <- DataFrame(response = response1)
SE1 <- SummarizedExperiment(assays = exp1,
                            colData = colData1)

colData2 <- DataFrame(response = response2)
SE2 <- SummarizedExperiment(assays = exp2,
                            colData = colData2)

SElist <- list(SE1, SE2)


data("MEL_GSE78220_meta")
data("MEL_PRJEB23709_meta")
#building model
mymodel <- build_CC_model(SElist, Stem.Sig, rmBE = TRUE)

exprs <- cbind(extract_mtr('MEL_GSE78220_exp'), extract_mtr('MEL_PRJEB23709_exp'))
exprs <- dataPreprocess(exprs, Stem.Sig, turn2HL = FALSE)
pData <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)
rownames(pData) <- pData$sample_id
pData$class <- c(extract_label('MEL_GSE78220_meta'), extract_label('MEL_PRJEB23709_meta'))
identical(rownames(pData),colnames(exprs))
metadata <- data.frame(labelDescription = colnames(pData), row.names = colnames(pData))
adf <- new("AnnotatedDataFrame", data = as.data.frame(pData), varMetadata = metadata)

exampleSet2 <- new("ExpressionSet", exprs = exprs, phenoData = adf)
prediction <- predict(mymodel, exampleSet2, positive = "R", dist = "cor")
result <- prediction@prediction

##testing model
library(pROC)

#read tigeR Built-in datasets
test_Expr <- cbind(extract_mtr('MEL_GSE78220_exp'), extract_mtr('MEL_PRJEB23709_exp'))
test_Expr <- dataPreprocess(test_Expr, Stem.Sig, turn2HL = FALSE)
test_response <- c(extract_label('MEL_GSE78220_meta'), extract_label('MEL_PRJEB23709_meta'))

#the index of sample which prediction result is Responder
predict_response_R <- result[,2] == 'R'

#Obtaining the meta informations of patients whose prediction results are 'Response'.
data("MEL_GSE78220_meta")
data("MEL_PRJEB23709_meta")

rc <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)[predict_response_R,]
rc$response <- sub('CR|MR|PR|CRPR', 'R', rc$response)
rc$response <- sub('PD|SD', 'NR', rc$response)
rc$response <- as.factor(rc$response)

#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
roc1 <- roc(rc$overall.survival..days., response = rc$vital.status)

rc2 <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)
rc2$response <- sub('CR|MR|PR|CRPR', 'R', rc2$response)
rc2$response <- sub('PD|SD', 'NR', rc2$response)
rc2$response <- as.factor(rc2$response)

roc2 <- roc(rc2$overall.survival..days., response = rc2$vital.status)

par(mfrow = c(1,2))
plot(roc1)
plot(roc2)

print("The AUC of selected patient is:")
auc(roc1)
print("The AUC of all patient is:")
auc(roc2)

```

## 6. Adaboost Model

```
#Please load data set in Baidu cloud before running the code！
Stem.Sig <- c("ENO1","RBM17","RHCG","H1FX","UPF2","C19orf12","TXNDC5","SERPINA1","HNRNPH2","PHB2","LHX2","NNAT","CHCHD5","FAM204A","IFT20","COPA","UQCRH","KRT1","S100A7","RYK","RPP38","NUDT19","RNF11","VPS35","CEACAM6","PTMA","FOXG1","AC018890.6","CISD2","PCBP2","SMIM26","ACKR1","SERBP1","MYL6","GSDMC","HIST1H2AE","MSRB2","LRP3","EEF2","MTRNR2L1","C18ORF32","SELENOW","CTDNEP1","KCTD8","SPINK1","LECT1","MMP24OS","MADCAM1","BCL10","SLC25A3","LY6D","TPD52L1","EIF5AL1","LSM14A","SUMO3","SPHK1","TMED5","SLC25A6","ALKBH5","OOEP","RENBP","VRK1","ABCA7","FABP3","DENND2C","CNIH1","CALML5","PDCD2","MARVELD1","GPI","SHC4","HM13","DERL2","UBB","PIPOX","STEAP1B","MRPL23","CIB2","PVR","GPR137B","S100A6","SERF2","DSG1","NDUFA4","LDHB","GRAMD1A","RRS1","LBP","VPS29","YBX1","TRIM47","SLC26A3","NTS","PKMYT1","SELENOM","HTR2B","JTB","FEM1B","KRTDAP","NUPR2","PFDN5","SYNE4","SND1","MYDGF","ATP5G2","C11orf31","SMOX","MAMDC2","CLN3","RPA1","DESI1","ZNF330","LRRFIP2","ALDOA","CITED4","SHFM1","GTSF1","FAM98C","PRDX1","C19orf70","C9ORF78","ATP5B","C20orf24","PLEKHA7","PPP4C","ALYREF","PIGP","EIF4EBP3","TMA7","PFN1","ANKRD35","ATP5J2","PRSS27","DYRK1B","CFH","UBA52","POLR2E","MID1IP1","RNF126","SIK3","VCX2","SF3A2","UFD1L","SNX8","TMEM41A","PSMB3","NBPF14","C7orf73","TCEB2","EGLN2","C4BPA","TMPRSS3","ZFAND6","KNG1","SH3GL1","GGACT","ADAMTSL1","ASF1B","UBE2L3","RGS20","ATP5I","KRT10","MRPS21","ARHGEF35","BCL7C","AP2S1","AGT","RACK1","COMMD8","OMG","RAD23A","C14orf166","DKK2","DDA1","EIF5A2","MLANA","H2AFZ","H3F3B","C1orf43","FAM104B","IRX3","PIH1D1","KYNU","EIF3M","EMC3","SLC25A29","SSBP4","ALKBH7","EIF2S3","FUZ","PNCK","LRRN4CL","UBE2D2","ACTG1","APOA1BP","NGFRAP1","IRX5","UBE2S","ATP5G3","AGR2","C5ORF15","SNF8","KXD1","ZNF90","EMX2","PDXP","PTHLH","TM7SF3","NPM1","RIOK3","TAGLN2","VGLL1","CALB2","RBX1","UGT2B15","MTHFD2","RBM3","INAFM1","CCDC106","CDHR1","EPHA3","LMF2","AC007325.4","GNPTAB","TUBB","RALY","PIGC","SFRP1","KIF1C","DSCR8","GC","SUB1","SET","LINC00844","IGSF21","LMF1","EPHA7","NCAPH2","AL353997.6","ZNF106","PERP","TICAM1","SNRPE","ENY2","C1QBP","HNRNPA1","PPIA","C11ORF58","SUMO2","TXNL4A","PLPPR5","SHC3","NTN5","AC011043.1","MIR205HG","TSPAN10","ACTB","C19orf43","PPP2R5A","AARD","TMEM256","CHCHD2","IGFBP1","ERP29","DERL1","TMPRSS5","FAM89A","EHD2","WIF1","IGKC","NECTIN1","BIRC7","AC090498.1","TIAM1","H3F3A","SLURP1","EIF4A1","TSPYL5","PON3","TMED10","GRTP1","GLRX5","GRID2","POLE4","CENPU","SNORD100","ECM1","CLPTM1","PLP2","GRHL3","C1orf131","C8orf82","KRT19","AL365205.1","MTRNR2L10","RER1","C6ORF62","ATP5J","GNB2L1","POM121","SMC2","HOMER1","C1orf56","BAIAP2L2","EFNB1","PKP1","CEBPZOS","TPM2","ATP5G1","WDR77","DIAPH2","AP3S1","CASP6","ANKRD6","KLRC2","SERPINA3","GTSE1","CPEB4","SPRR1B","SLC25A5","GPR87","COX7A2L","PRKCDBP","NME1","NDC1","FGL1","COPB1","BECN1","COX4I1","NKD1","STEAP3","PCP4","CDK13","SPRR2D","PDP1","ADGRF4","PEX13","FTH1","HN1","FZD8","TXN","EEF1A1","PSMF1","PFKFB2","CYP27B1","GNB1","HPCA","TAC1","SFTPB","CDKN2A","FAM83B","TCF7L1","EEF1G","FAM110A","TP53BP2","ORM1","FAU","ANP32B","CRYGS","C1ORF61","PARL","TMEM18","BEX5","FABP6","CDKN2B","CALML3","RPIA","FAM89B","SNRPB","RSL1D1","ORM2","CCNDBP1","ATP5F1B","NUDT4","CMTM5","G6PD","POU4F3","TPH1","MUC21","UBAP1","HCAR2","FAHD2A","CCS","COMMD7","EIF2S2","C5","AATF","HINT1","SPOCD1","ATP23","FAAP20","NRM","AGRN","SLC22A31","SEC61B","ZNF750","DBI","DHCR7","BPI","DDOST","GAPDH","NOX1","HNRNPD","KCNF1","TTC9B","NR0B2","UBE2A","CGN","MRPL54","ZBTB43","SERPINB2","NAA50","DGAT2","PSMA7","MKI67","TPI1","CHMP4A","LYPLA1","OSBPL6","C1ORF122","APH1A","MCTS1","SFTA2","NAPSA","PPP1CA","TINCR","ROPN1","BARX2","MRGBP","AC012640.4","COX6A1","GSTA1","MTRNR2L2","MTHFD2L","C17ORF89","CRP","IDH3G","GCLC","AL031727.1","GSTP1","IL36G","ROPN1B","ZMYND11","POP4","PRDX3","RAN","GTF3A","NACA","ARSJ","LINC01268","APOA2","TIMM23","FAM83A","QARS")

#standardization of expression matrix
exp1 <- as.matrix(MEL_GSE91061_exp[,-1])
exp1 <- apply(exp1, 2, as.numeric)
rownames(exp1) <- MEL_GSE91061_exp[,1]

exp2 <- as.matrix(MEL_phs000452_exp[,-1])
exp2 <- apply(exp2, 2, as.numeric)
rownames(exp2) <- MEL_phs000452_exp[,1]

meta1 <- MEL_GSE91061_meta
meta2 <- MEL_phs000452_meta

library(tigeR)
response1 <- meta1$response
response1 <- response_standardize(response1)
response2 <- meta2$response
response2 <- response_standardize(response2)

#obtain index of NA samples
filt1 <- grep('UNK',response1)
filt2 <- grep('UNK',response2)

#remove UNK. exp2 has no UNK
exp1 <- exp1[,-filt1]
response1 <- response1[-filt1]
#exp2 <- exp2[,-filt1]
#response2 <- response2[-filt2]

library(SummarizedExperiment)
colData1 <- DataFrame(response = response1)
SE1 <- SummarizedExperiment(assays = exp1,
                            colData = colData1)

colData2 <- DataFrame(response = response2)
SE2 <- SummarizedExperiment(assays = exp2,
                            colData = colData2)

SElist <- list(SE1, SE2)


data("MEL_GSE78220_meta")
data("MEL_PRJEB23709_meta")
#building model
library(adabag)
mymodel <- build_Adaboost_model(SElist, Stem.Sig, rmBE = FALSE)

exp_test <- cbind(extract_mtr('MEL_GSE78220_exp'), extract_mtr('MEL_PRJEB23709_exp'))
exp_test <- dataPreprocess(exp_test, Stem.Sig, turn2HL = FALSE)

pred <- predict.boosting(mymodel, as.data.frame(t(exp_test)))

result <- pred$class

##testing model
library(pROC)

#read tigeR Built-in datasets
test_Expr <- cbind(extract_mtr('MEL_GSE78220_exp'), extract_mtr('MEL_PRJEB23709_exp'))
test_Expr <- dataPreprocess(test_Expr, Stem.Sig, turn2HL = FALSE)
test_response <- c(extract_label('MEL_GSE78220_meta'), extract_label('MEL_PRJEB23709_meta'))

#the index of sample which prediction result is Responder
predict_response_R <- result == 'R'

#Obtaining the meta informations of patients whose prediction results are 'Response'.
data("MEL_GSE78220_meta")
data("MEL_PRJEB23709_meta")

rc <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)[predict_response_R,]
rc$response <- sub('CR|MR|PR|CRPR', 'R', rc$response)
rc$response <- sub('PD|SD', 'NR', rc$response)
rc$response <- as.factor(rc$response)

#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
roc1 <- roc(rc$overall.survival..days., response = rc$vital.status)

rc2 <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)
rc2$response <- sub('CR|MR|PR|CRPR', 'R', rc2$response)
rc2$response <- sub('PD|SD', 'NR', rc2$response)
rc2$response <- as.factor(rc2$response)

roc2 <- roc(rc2$overall.survival..days., response = rc2$vital.status)

par(mfrow = c(1,2))
plot(roc1)
plot(roc2)

print("The AUC of selected patient is:")
auc(roc1)
print("The AUC of all patient is:")
auc(roc2)

```

## 7. Logitboost Model

```
#Please load data set in Baidu cloud before running the code！
Stem.Sig <- c("ENO1","RBM17","RHCG","H1FX","UPF2","C19orf12","TXNDC5","SERPINA1","HNRNPH2","PHB2","LHX2","NNAT","CHCHD5","FAM204A","IFT20","COPA","UQCRH","KRT1","S100A7","RYK","RPP38","NUDT19","RNF11","VPS35","CEACAM6","PTMA","FOXG1","AC018890.6","CISD2","PCBP2","SMIM26","ACKR1","SERBP1","MYL6","GSDMC","HIST1H2AE","MSRB2","LRP3","EEF2","MTRNR2L1","C18ORF32","SELENOW","CTDNEP1","KCTD8","SPINK1","LECT1","MMP24OS","MADCAM1","BCL10","SLC25A3","LY6D","TPD52L1","EIF5AL1","LSM14A","SUMO3","SPHK1","TMED5","SLC25A6","ALKBH5","OOEP","RENBP","VRK1","ABCA7","FABP3","DENND2C","CNIH1","CALML5","PDCD2","MARVELD1","GPI","SHC4","HM13","DERL2","UBB","PIPOX","STEAP1B","MRPL23","CIB2","PVR","GPR137B","S100A6","SERF2","DSG1","NDUFA4","LDHB","GRAMD1A","RRS1","LBP","VPS29","YBX1","TRIM47","SLC26A3","NTS","PKMYT1","SELENOM","HTR2B","JTB","FEM1B","KRTDAP","NUPR2","PFDN5","SYNE4","SND1","MYDGF","ATP5G2","C11orf31","SMOX","MAMDC2","CLN3","RPA1","DESI1","ZNF330","LRRFIP2","ALDOA","CITED4","SHFM1","GTSF1","FAM98C","PRDX1","C19orf70","C9ORF78","ATP5B","C20orf24","PLEKHA7","PPP4C","ALYREF","PIGP","EIF4EBP3","TMA7","PFN1","ANKRD35","ATP5J2","PRSS27","DYRK1B","CFH","UBA52","POLR2E","MID1IP1","RNF126","SIK3","VCX2","SF3A2","UFD1L","SNX8","TMEM41A","PSMB3","NBPF14","C7orf73","TCEB2","EGLN2","C4BPA","TMPRSS3","ZFAND6","KNG1","SH3GL1","GGACT","ADAMTSL1","ASF1B","UBE2L3","RGS20","ATP5I","KRT10","MRPS21","ARHGEF35","BCL7C","AP2S1","AGT","RACK1","COMMD8","OMG","RAD23A","C14orf166","DKK2","DDA1","EIF5A2","MLANA","H2AFZ","H3F3B","C1orf43","FAM104B","IRX3","PIH1D1","KYNU","EIF3M","EMC3","SLC25A29","SSBP4","ALKBH7","EIF2S3","FUZ","PNCK","LRRN4CL","UBE2D2","ACTG1","APOA1BP","NGFRAP1","IRX5","UBE2S","ATP5G3","AGR2","C5ORF15","SNF8","KXD1","ZNF90","EMX2","PDXP","PTHLH","TM7SF3","NPM1","RIOK3","TAGLN2","VGLL1","CALB2","RBX1","UGT2B15","MTHFD2","RBM3","INAFM1","CCDC106","CDHR1","EPHA3","LMF2","AC007325.4","GNPTAB","TUBB","RALY","PIGC","SFRP1","KIF1C","DSCR8","GC","SUB1","SET","LINC00844","IGSF21","LMF1","EPHA7","NCAPH2","AL353997.6","ZNF106","PERP","TICAM1","SNRPE","ENY2","C1QBP","HNRNPA1","PPIA","C11ORF58","SUMO2","TXNL4A","PLPPR5","SHC3","NTN5","AC011043.1","MIR205HG","TSPAN10","ACTB","C19orf43","PPP2R5A","AARD","TMEM256","CHCHD2","IGFBP1","ERP29","DERL1","TMPRSS5","FAM89A","EHD2","WIF1","IGKC","NECTIN1","BIRC7","AC090498.1","TIAM1","H3F3A","SLURP1","EIF4A1","TSPYL5","PON3","TMED10","GRTP1","GLRX5","GRID2","POLE4","CENPU","SNORD100","ECM1","CLPTM1","PLP2","GRHL3","C1orf131","C8orf82","KRT19","AL365205.1","MTRNR2L10","RER1","C6ORF62","ATP5J","GNB2L1","POM121","SMC2","HOMER1","C1orf56","BAIAP2L2","EFNB1","PKP1","CEBPZOS","TPM2","ATP5G1","WDR77","DIAPH2","AP3S1","CASP6","ANKRD6","KLRC2","SERPINA3","GTSE1","CPEB4","SPRR1B","SLC25A5","GPR87","COX7A2L","PRKCDBP","NME1","NDC1","FGL1","COPB1","BECN1","COX4I1","NKD1","STEAP3","PCP4","CDK13","SPRR2D","PDP1","ADGRF4","PEX13","FTH1","HN1","FZD8","TXN","EEF1A1","PSMF1","PFKFB2","CYP27B1","GNB1","HPCA","TAC1","SFTPB","CDKN2A","FAM83B","TCF7L1","EEF1G","FAM110A","TP53BP2","ORM1","FAU","ANP32B","CRYGS","C1ORF61","PARL","TMEM18","BEX5","FABP6","CDKN2B","CALML3","RPIA","FAM89B","SNRPB","RSL1D1","ORM2","CCNDBP1","ATP5F1B","NUDT4","CMTM5","G6PD","POU4F3","TPH1","MUC21","UBAP1","HCAR2","FAHD2A","CCS","COMMD7","EIF2S2","C5","AATF","HINT1","SPOCD1","ATP23","FAAP20","NRM","AGRN","SLC22A31","SEC61B","ZNF750","DBI","DHCR7","BPI","DDOST","GAPDH","NOX1","HNRNPD","KCNF1","TTC9B","NR0B2","UBE2A","CGN","MRPL54","ZBTB43","SERPINB2","NAA50","DGAT2","PSMA7","MKI67","TPI1","CHMP4A","LYPLA1","OSBPL6","C1ORF122","APH1A","MCTS1","SFTA2","NAPSA","PPP1CA","TINCR","ROPN1","BARX2","MRGBP","AC012640.4","COX6A1","GSTA1","MTRNR2L2","MTHFD2L","C17ORF89","CRP","IDH3G","GCLC","AL031727.1","GSTP1","IL36G","ROPN1B","ZMYND11","POP4","PRDX3","RAN","GTF3A","NACA","ARSJ","LINC01268","APOA2","TIMM23","FAM83A","QARS")

#standardization of expression matrix
exp1 <- as.matrix(MEL_GSE91061_exp[,-1])
exp1 <- apply(exp1, 2, as.numeric)
rownames(exp1) <- MEL_GSE91061_exp[,1]

exp2 <- as.matrix(MEL_phs000452_exp[,-1])
exp2 <- apply(exp2, 2, as.numeric)
rownames(exp2) <- MEL_phs000452_exp[,1]

meta1 <- MEL_GSE91061_meta
meta2 <- MEL_phs000452_meta

library(tigeR)
response1 <- meta1$response
response1 <- response_standardize(response1)
response2 <- meta2$response
response2 <- response_standardize(response2)

#obtain index of NA samples
filt1 <- grep('UNK',response1)
filt2 <- grep('UNK',response2)

#remove UNK. exp2 has no UNK
exp1 <- exp1[,-filt1]
response1 <- response1[-filt1]
#exp2 <- exp2[,-filt1]
#response2 <- response2[-filt2]

library(SummarizedExperiment)
colData1 <- DataFrame(response = response1)
SE1 <- SummarizedExperiment(assays = exp1,
                            colData = colData1)

colData2 <- DataFrame(response = response2)
SE2 <- SummarizedExperiment(assays = exp2,
                            colData = colData2)

SElist <- list(SE1, SE2)

#building model
library(tigeR)
mymodel <- build_Logitboost_model(SElist, Stem.Sig, rmBE = FALSE)

##testing model
library(pROC)

#read tigeR Built-in datasets
test_Expr <- cbind(extract_mtr('MEL_GSE78220_exp'), extract_mtr('MEL_PRJEB23709_exp'))
test_Expr <- dataPreprocess(test_Expr, Stem.Sig, turn2HL = FALSE)
test_response <- rbind(extract_label('MEL_GSE78220_meta'), extract_label('MEL_PRJEB23709_meta'))

#the index of sample which prediction result is Responder
predict_response_R <- caTools::predict.LogitBoost(mymodel, t(test_Expr), type = 'class') == 'R'

#Obtaining the meta informations of patients whose prediction results are 'Response'.
data("MEL_GSE78220_meta")
data("MEL_PRJEB23709_meta")

rc <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)[predict_response_R,]
rc$response <- sub('CR|MR|PR|CRPR', 'R', rc$response)
rc$response <- sub('PD|SD', 'NR', rc$response)
rc$response <- as.factor(rc$response)

#Drawing roc curve of patients whose prediction results are 'Responder' and calculating the AUC of roc curver.
roc1 <- roc(rc$overall.survival..days., response = rc$vital.status)

rc2 <- rbind(MEL_GSE78220_meta, MEL_PRJEB23709_meta)
rc2$response <- sub('CR|MR|PR|CRPR', 'R', rc2$response)
rc2$response <- sub('PD|SD', 'NR', rc2$response)
rc2$response <- as.factor(rc2$response)

roc2 <- roc(rc2$overall.survival..days., response = rc2$vital.status)

par(mfrow = c(1,2))
plot(roc1)
plot(roc2)

print("The AUC of selected patient is:")
auc(roc1)
print("The AUC of all patient is:")
auc(roc2)

```
