#' ---	
#' title: "PAMR main processing"	
#' output:	
#'   html_document:	
#'     self_contained: false	
#' knit: (function(inputFile, encoding) {	
#'   rmarkdown::render(inputFile, encoding = encoding, output_dir = here::here("output/"))})	
#' ---	
#' 	
#' 	
knitr::opts_chunk$set(echo = TRUE)	
#' 	
#' 	
#' ## Load data and Preprocessing	
#' 	
#' Can be examined on data_load_and_preprocess.Rmd/html auxillary notebook, here it is only loaded:	
#' 	
#' 	
source(here::here("scripts","Rmds","data_load_and_preprocess.R"))	
#' 	
#' 	
#' # Fitting	
#' 	
#' 	
library(pamr)	
	
pamr=pamr.train(traindata)	
adaptive_thresholding=TRUE	
	
  if (adaptive_thresholding) {	
    thres=pamr.adaptthresh(pamr)	
    pamr=pamr.train(traindata, threshold.scale = thres)	
  }	
	
pamrcv=pamr.cv(pamr, traindata)	
#' 	
#' 	
#' The $\theta$ scaling factors:	
#' 	
#' 	
thres	
#' 	
#' 	
#' 	
pamrcv
dev.off()
pamr.plotcv(pamrcv)	
#' 	
#' 	
#' 	
#' ## Threshold λ selection	
#' 	
#' Selection of two threshold value to evaluate	
#' 	
#' $To$ - threshold "optimal", lowest cross-validated error rate <br>	
#' $Th$ - threshold "higher" , suboptimal (but still errors do not explode), with fewer genes selected	
#' 	
#' 	
To=17.38  #  Optimal threshold	
Th=24.33  #   Higher threshold	
#' 	
#' 	
#' Visualizing shrunken deltas (differences from shrunken centroid to overall centroids) and saving also the index of the retained genes for classification at the two thresholds (that can be used to run another classification algorithm that does not perform also feature selection):	
#' 	
#' <details>	
#' <summary>Click to expand</summary>	
#' 	
#' 	
selectedgenes=pamr.listgenes(pamr, traindata, threshold = To)	
selectedgenes=as.vector(selectedgenes[,1])	
print(paste("Selected genes at optimal threshold to", To, ":", length(selectedgenes)))	
sgenes_indexTo <- rownames(traindata$x) %in% selectedgenes	
selectedgenes=pamr.listgenes(pamr, traindata, threshold = Th)	
selectedgenes=as.vector(selectedgenes[,1])	
print(paste("Selected genes at optimal threshold to", To, ":", length(selectedgenes)))	
sgenes_indexTh<- rownames(traindata$x) %in% selectedgenes	
sgenes_index=data.frame(sgenes_indexTo,sgenes_indexTh)	
write.csv(sgenes_index, here::here("output/selected_genes.csv"), row.names = FALSE)	
#' 	
#' 	
#' </details>	
#' 	
#' 	
#' # Visualizing classifications	
#' 	
#' Loading required libraries:	
#' 	
#' 	
library(classmap) 	
library(classmapExt) #The novel auxillary package to classmap	
#' 	
#' 	
#' ### Training Set	
#' 	
#' Using vcr.pamr.train function to produce all the quantities needed for the visualization of the two models with the two thresholds:	
#' 	
#' 	
vcrtrainTo=vcr.pamr.train(data=traindata, pamrfit = pamr, threshold=To)	
vcrtrainTh=vcr.pamr.train(data=traindata, pamrfit = pamr, threshold=Th)	
#' 	
#' 	
#' #### Silhouttes	
#' 	
#' 	
cfm=caret::confusionMatrix(factor(vcrtrainTo$pred),factor(vcrtrainTo$y))	
cfm$table	
silplot(vcrtrainTo, main="")	
#' 	
#' 	
#' 	
cfm=caret::confusionMatrix(factor(vcrtrainTh$pred),factor(vcrtrainTh$y))	
cfm$table	
silplot(vcrtrainTh, main="")	
#' 	
#' 	
#' #### Classmaps	
#' 	
#' 	
par(mfrow = c(3, 1))
classmap(vcrtrainTo, whichclass = 1, main = "Pred of class LUAD" )	
classmap(vcrtrainTo, whichclass = 2, main = "Pred of class LUSC")	
classmap(vcrtrainTo, whichclass = 3, main = "Pred of class NORM")	
mtext(paste0("Classmaps for training at threshold: ",To), line=-5, side=3, outer=TRUE, cex=1)	
#' 	
#' 	
#' 	
par(mfrow=c(3, 1))	
classmap(vcrtrainTh, whichclass = 1, main = "Pred of class LUAD" )	
classmap(vcrtrainTh, whichclass = 2, main = "Pred of class LUSC")	
classmap(vcrtrainTh, whichclass = 3, main = "Pred of class NORM")	
mtext(paste0("Classmaps for training at threshold: ",Th), line=-5, side=3, outer=TRUE, cex=1)	
#' 	
#' 	
#' #### MDS color-scaled plots	
#' 	
#' 	
mdscolorscale(vcrtrainTo, diss=vcrtrainTo$pwd, main=	
                paste0("(Train) MDScolorscale at thresh: ",To))	
#' 	
#' 	
#' 	
mdscolorscale(vcrtrainTh, diss=vcrtrainTh$pwd, main=	
                paste0("(Train) MDScolorscale at thresh: ",Th))	
#' 	
#' 	
#' #### Quasi Residual Plots	
#' 	
#' Drawn only for higher threshold $Th$	
#' 	
#' For some continuos covariates:	
#' 	
#' 	
	
par(mfrow=c(3,3))	
	
qresplot(vcrtrainTh$PAC,traindata$pancovariates$MSI.MANTIS.Score, plotErrorBars = TRUE, 	
         main="Mantis Score", xlim=c(0.2,0.4))	
qresplot(vcrtrainTh$PAC,traindata$covariates$years_smoked, plotErrorBars = TRUE, 	
         main="Years smoked")	
qresplot(vcrtrainTh$PAC,traindata$covariates$longest_dimension, plotErrorBars = TRUE, 	
         main="Tumor dimension (longest)", xlim=c(0,3))	
qresplot(vcrtrainTh$PAC,traindata$pancovariates$Fraction.Genome.Altered, plotErrorBars = TRUE, 	
         main="Fraction of Genome Altered")	
qresplot(vcrtrainTh$PAC,traindata$pancovariates$Diagnosis.Age, plotErrorBars = TRUE, 	
         main="Age")	
qresplot(vcrtrainTh$PAC,traindata$covariates$`Alignment Quality Index`, plotErrorBars = TRUE,	
         main="Alignment Quality Index", xlim = c(0.045,0.11), grid=seq(0.04, 0.11, length.out=10))	
#' 	
#' 	
#' 	
	
par(mfrow=c(1,3))	
	
qresplot(vcrtrainTh$PAC[traindata$y=="LUAD"],	
         traindata$covariates$`Alignment Quality Index`[traindata$y=="LUAD"],plotErrorBars = TRUE,	
         main="Alignment Quality Index in LUAD", xlim = c(0.045,0.11), grid=seq(0.04, 0.11, length.out=10))	
qresplot(vcrtrainTh$PAC[traindata$y=="LUSC"],	
         traindata$covariates$`Alignment Quality Index`[traindata$y=="LUSC"],plotErrorBars = TRUE,	
         main="Alignment Quality Index in LUSC", xlim = c(0.045,0.11), grid=seq(0.04, 0.11, length.out=10))	
qresplot(vcrtrainTh$PAC[traindata$y=="NORM"],	
         traindata$covariates$`Alignment Quality Index`[traindata$y=="NORM"],plotErrorBars = TRUE,	
         main="Alignment Quality Index in NORM", xlim = c(0.045,0.11), grid=seq(0.04, 0.11, length.out=10))	
#' 	
#' 	
#' For some discrete covariates oppurtunely integer-coded ( a table with corresponding values is provided)	
#' 	
#' 	
	
legend=data.frame(matrix(nrow = 20, ncol = 0))	
legend$Code=1:20	
	
par(mfrow=c(2,2))	
	
#ajcc stage	
discrete_covariate=factor(unlist(traindata$covariates$ajcc_pathologic_stage, use.names = FALSE)) 	
levels=levels(discrete_covariate)	
discrete_cov_length=length(levels)	
if (length(levels) < 20) {	
    levels <- c(levels, rep(NA, 20 - length(levels)))	
}	
legend$`Ajcc Stage`=levels	
qresplot(vcrtrainTh$PAC,as.numeric(discrete_covariate), plotErrorBars = TRUE,	
         grid=seq(1:(1+discrete_cov_length))-0.5 , main = "Ajcc Stage")	
#histologic type	
discrete_covariate=factor(unlist(traindata$covariates$morphology, use.names = FALSE)) 	
levels=levels(discrete_covariate)	
discrete_cov_length=length(levels)	
if (length(levels) < 20) {	
    levels <- c(levels, rep(NA, 20 - length(levels)))	
}	
legend$`Histologic Type`=levels	
qresplot(vcrtrainTh$PAC,as.numeric(discrete_covariate), plotErrorBars = TRUE,	
         grid=seq(1:(1+discrete_cov_length))-0.5 , main = "Histologic Type")	
#site	
discrete_covariate=factor(unlist(traindata$covariates$site_of_resection_or_biopsy, use.names = FALSE)) 	
levels=levels(discrete_covariate)	
discrete_cov_length=length(levels)	
if (length(levels) < 20) {	
    levels <- c(levels, rep(NA, 20 - length(levels)))	
}	
legend$`Site`=levels	
qresplot(vcrtrainTh$PAC,as.numeric(discrete_covariate), plotErrorBars = TRUE,	
         grid=seq(1:(1+discrete_cov_length))-0.5 , main = "Site")	
#gender	
discrete_covariate=factor(unlist(traindata$covariates$gender, use.names = FALSE)) 	
levels=levels(discrete_covariate)	
discrete_cov_length=length(levels)	
if (length(levels) < 20) {	
    levels <- c(levels, rep(NA, 20 - length(levels)))	
}	
legend$`Gender`=levels	
qresplot(vcrtrainTh$PAC,as.numeric(discrete_covariate), plotErrorBars = TRUE,	
         grid=seq(1:(1+discrete_cov_length))-0.5 , main = "Gender", xlim=c(0,3))	
	
	
library(kableExtra)	
knitr::kable(legend, "html", align = "c") %>% 	
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = F)	
#' 	
#' 	
#' On selected genes for the classification in dataset:	
#' 	
#' 	
	
selectedgenes=which(sgenes_index$sgenes_indexTh)	
	
#in all the genes	
coeff=rep(0,nrow(traindata$x[selectedgenes,]))	
	
#getting indexes only of genes that take on more that 7 unique values (<7 values it is rather discrete)	
#also lefting out genes that are more than 35% observation that take 0 values	
index_nonsparse <- which(apply(traindata$x[selectedgenes,], 1, function(x) {	
    unique_count <- length(unique(x))	
    zero_count <- sum(x == 0)	
    return(unique_count >= 7 & zero_count / length(x) < 0.35)}))	
	
for (i in index_nonsparse) {	
  lm=lm(vcrtrainTh$PAC~traindata$x[i,])	
  coeff[i]=lm$coefficients[2]	
}	
	
sorted_indices <- order(coeff, decreasing = TRUE)	
sorted_coeff <- sort(coeff, decreasing = TRUE)	
sorted_indices[1:9]	
sorted_coeff[1:9]	
	
par(mfrow=c(3,3))	
	
for (i in sorted_indices[1:9]) {	
  gene=traindata$x[i,]	
  qresplot(vcrtrainTh$PAC,gene,plotErrorBars = TRUE, main = rownames(traindata$x)[i])	
}	
#' 	
#' 	
#' On all genes in dataset:	
#' 	
#' 	
#in all the genes	
coeff=rep(0,nrow(traindata$x))	
	
#getting indexes only of genes that take on more that 7 unique values (<7 values it is rather discrete)	
index_nonsparse <- which(apply(traindata$x, 1, function(x) length(unique(x)) >= 8)) 	
	
#getting indexes only of genes that take on more that 7 unique values (<7 values it is rather discrete)	
#also lefting out genes that are more than 35% observation that take 0 values	
index_nonsparse <- which(apply(traindata$x, 1, function(x) {	
    unique_count <- length(unique(x))	
    zero_count <- sum(x == 0)	
    return(unique_count >= 7 & zero_count / length(x) < 0.35)}))	
	
for (i in index_nonsparse) {	
  lm=lm(vcrtrainTh$PAC~traindata$x[i,])	
  coeff[i]=lm$coefficients[2]	
}	
	
sorted_indices <- order(coeff, decreasing = TRUE)	
sorted_coeff <- sort(coeff, decreasing = TRUE)	
sorted_indices[1:9]	
sorted_coeff[1:9]	
	
	
par(mfrow=c(3,3))	
	
for (i in sorted_indices[1:9]) {	
  gene=traindata$x[i,]	
  qresplot(vcrtrainTh$PAC,gene,plotErrorBars = TRUE, main = rownames(traindata$x)[i], 	
           grid=seq(0.0, 0.3, length.out=10), xlim = c(0,0.25))	
}	
#' 	
#' 	
#' ### Testing Set	
#' 	
#' Using vcr.pamr.newdata function to produce all the quantities needed for the visualization of the classification performed on test set by the two models with the two thresholds:	
#' 	
#' 	
vcrtestTo=vcr.pamr.newdata(newdata=testdata, vcr.pamr.train.out = vcrtrainTo)	
vcrtestTh=vcr.pamr.newdata(newdata=testdata, vcr.pamr.train.out = vcrtrainTh)	
#' 	
#' 	
#' #### Silhouttes	
#' 	
#' 	
cfm=caret::confusionMatrix(factor(vcrtestTo$pred),factor(vcrtestTo$ynew))	
cfm$table	
silplot(vcrtestTo, main="")	
#' 	
#' 	
#' 	
cfm=caret::confusionMatrix(factor(vcrtestTh$pred),factor(vcrtestTh$ynew))	
cfm$table	
silplot(vcrtestTh, main="")	
#' 	
#' 	
#' #### Classmaps	
#' 	
#' 	
par(mfrow=c(3,1))	
classmap(vcrtestTo, whichclass = 1, main = "Pred of class LUAD" )	
classmap(vcrtestTo, whichclass = 2, main = "Pred of class LUSC")	
classmap(vcrtestTo, whichclass = 3, main = "Pred of class NORM")	
mtext(paste0("Classmaps for testing at threshold: ",To), line=-5, side=3, outer=TRUE, cex=1)	
#' 	
#' 	
#' 	
par(mfrow=c(3,1))	
classmap(vcrtestTh, whichclass = 1, main = "Pred of class LUAD" )	
classmap(vcrtestTh, whichclass = 2, main = "Pred of class LUSC")	
classmap(vcrtestTh, whichclass = 3, main = "Pred of class NORM")	
mtext(paste0("Classmaps for testing at threshold: ",Th), line=-5, side=3, outer=TRUE, cex=1)	
#' 	
#' 	
#' #### MDS color-scaled plots	
#' 	
#' 	
mdscolorscale(vcrtestTo, diss=vcrtestTo$pwd, main=	
                paste0("(test) MDScolorscale at thresh: ",To))	
#' 	
#' 	
#' 	
mdscolorscale(vcrtestTh, diss=vcrtestTh$pwd, main=	
                paste0("(test) MDScolorscale at thresh: ",Th))	
#' 	
#' 	
#' #### Quasi Residual Plots	
#' 	
#' Drawn only for higher threshold $Th$	
#' 	
#' For some continuos covariates:	
#' 	
#' 	
	
par(mfrow=c(2,3))	
	
qresplot(vcrtestTh$PAC,testdata$pancovariates$MSI.MANTIS.Score, plotErrorBars = TRUE, 	
         main="Mantis Score", xlim=c(0.2,0.4))	
qresplot(vcrtestTh$PAC,testdata$covariates$years_smoked, plotErrorBars = TRUE, 	
         main="Years smoked")	
qresplot(vcrtestTh$PAC,testdata$covariates$longest_dimension, plotErrorBars = TRUE, 	
         main="Tumor dimension (longest)", xlim=c(0,3))	
qresplot(vcrtestTh$PAC,testdata$pancovariates$Fraction.Genome.Altered, plotErrorBars = TRUE, 	
         main="Fraction of Genome Altered")	
qresplot(vcrtestTh$PAC,testdata$pancovariates$Diagnosis.Age, plotErrorBars = TRUE, 	
         main="Age")	
qresplot(vcrtestTh$PAC,testdata$covariates$`Alignment Quality Index`, plotErrorBars = TRUE,	
         main="Alignment Quality Index", xlim = c(0.045,0.11), grid=seq(0.04, 0.11, length.out=10))	
#' 	
#' 	
#' 	
	
par(mfrow=c(1,3))	
	
qresplot(vcrtestTh$PAC[testdata$y=="LUAD"],	
         testdata$covariates$`Alignment Quality Index`[testdata$y=="LUAD"],plotErrorBars = TRUE,	
         main="Alignment Quality Index in LUAD", xlim = c(0.045,0.11), grid=seq(0.04, 0.11, length.out=10))	
qresplot(vcrtestTh$PAC[testdata$y=="LUSC"],	
         testdata$covariates$`Alignment Quality Index`[testdata$y=="LUSC"],plotErrorBars = TRUE,	
         main="Alignment Quality Index in LUSC", xlim = c(0.045,0.11), grid=seq(0.04, 0.11, length.out=10))	
qresplot(vcrtestTh$PAC[testdata$y=="NORM"],	
         testdata$covariates$`Alignment Quality Index`[testdata$y=="NORM"],plotErrorBars = TRUE,	
         main="Alignment Quality Index in NORM", xlim = c(0.045,0.11), grid=seq(0.04, 0.11, length.out=10))	
#' 	
#' 	
#' For some discrete covariates oppurtunely integer-coded ( a table with corresponding values is provided)	
#' 	
#' 	
	
legend=data.frame(matrix(nrow = 20, ncol = 0))	
legend$Code=1:20	
	
par(mfrow=c(2,2))	
	
#ajcc stage	
discrete_covariate=factor(unlist(testdata$covariates$ajcc_pathologic_stage, use.names = FALSE)) 	
levels=levels(discrete_covariate)	
discrete_cov_length=length(levels)	
if (length(levels) < 20) {	
    levels <- c(levels, rep(NA, 20 - length(levels)))	
}	
legend$`Ajcc Stage`=levels	
qresplot(vcrtestTh$PAC,as.numeric(discrete_covariate), plotErrorBars = TRUE,	
         grid=seq(1:(1+discrete_cov_length))-0.5 , main = "Ajcc Stage")	
#histologic type	
discrete_covariate=factor(unlist(testdata$covariates$morphology, use.names = FALSE)) 	
levels=levels(discrete_covariate)	
discrete_cov_length=length(levels)	
if (length(levels) < 20) {	
    levels <- c(levels, rep(NA, 20 - length(levels)))	
}	
legend$`Histologic Type`=levels	
qresplot(vcrtestTh$PAC,as.numeric(discrete_covariate), plotErrorBars = TRUE,	
         grid=seq(1:(1+discrete_cov_length))-0.5 , main = "Histologic Type")	
#site	
discrete_covariate=factor(unlist(testdata$covariates$site_of_resection_or_biopsy, use.names = FALSE)) 	
levels=levels(discrete_covariate)	
discrete_cov_length=length(levels)	
if (length(levels) < 20) {	
    levels <- c(levels, rep(NA, 20 - length(levels)))	
}	
legend$`Site`=levels	
qresplot(vcrtestTh$PAC,as.numeric(discrete_covariate), plotErrorBars = TRUE,	
         grid=seq(1:(1+discrete_cov_length))-0.5 , main = "Site")	
#gender	
discrete_covariate=factor(unlist(testdata$covariates$gender, use.names = FALSE)) 	
levels=levels(discrete_covariate)	
discrete_cov_length=length(levels)	
if (length(levels) < 20) {	
    levels <- c(levels, rep(NA, 20 - length(levels)))	
}	
legend$`Gender`=levels	
qresplot(vcrtestTh$PAC,as.numeric(discrete_covariate), plotErrorBars = TRUE,	
         grid=seq(1:(1+discrete_cov_length))-0.5 , main = "Gender", xlim=c(0,3))	
	
	
library(kableExtra)	
knitr::kable(legend, "html", align = "c") %>% 	
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = F)	
#' 	
#' 	
#' On selected genes for the classification in dataset:	
#' 	
#' 	
	
selectedgenes=which(sgenes_index$sgenes_indexTh)	
	
#in all the genes	
coeff=rep(0,nrow(testdata$x[selectedgenes,]))	
	
#getting indexes only of genes that take on more that 7 unique values (<7 values it is rather discrete)	
#also lefting out genes that are more than 35% observation that take 0 values	
index_nonsparse <- which(apply(testdata$x[selectedgenes,], 1, function(x) {	
    unique_count <- length(unique(x))	
    zero_count <- sum(x == 0)	
    return(unique_count >= 7 & zero_count / length(x) < 0.35)}))	
	
for (i in index_nonsparse) {	
  lm=lm(vcrtestTh$PAC~testdata$x[i,])	
  coeff[i]=lm$coefficients[2]	
}	
	
sorted_indices <- order(coeff, decreasing = TRUE)	
sorted_coeff <- sort(coeff, decreasing = TRUE)	
sorted_indices[1:9]	
sorted_coeff[1:9]	
	
par(mfrow=c(3,3))	
	
for (i in sorted_indices[1:9]) {	
  gene=testdata$x[i,]	
  qresplot(vcrtestTh$PAC,gene,plotErrorBars = TRUE, main = rownames(testdata$x)[i])	
}	
#' 	
#' 	
#' On all genes in dataset:	
#' 	
#' 	
#in all the genes	
coeff=rep(0,nrow(testdata$x))	
	
#getting indexes only of genes that take on more that 7 unique values (<7 values it is rather discrete)	
index_nonsparse <- which(apply(testdata$x, 1, function(x) length(unique(x)) >= 8)) 	
	
#getting indexes only of genes that take on more that 7 unique values (<7 values it is rather discrete)	
#also lefting out genes that are more than 35% observation that take 0 values	
index_nonsparse <- which(apply(testdata$x, 1, function(x) {	
    unique_count <- length(unique(x))	
    zero_count <- sum(x == 0)	
    return(unique_count >= 7 & zero_count / length(x) < 0.35)}))	
	
for (i in index_nonsparse) {	
  lm=lm(vcrtestTh$PAC~testdata$x[i,])	
  coeff[i]=lm$coefficients[2]	
}	
	
sorted_indices <- order(coeff, decreasing = TRUE)	
sorted_coeff <- sort(coeff, decreasing = TRUE)	
sorted_indices[1:9]	
sorted_coeff[1:9]	
	
	
par(mfrow=c(3,3))	
	
for (i in sorted_indices[1:9]) {	
  gene=testdata$x[i,]	
  qresplot(vcrtestTh$PAC,gene,plotErrorBars = TRUE, main = rownames(testdata$x)[i], 	
           grid=seq(0.0, 0.3, length.out=10), xlim = c(0,0.25))	
}	
#' 	
