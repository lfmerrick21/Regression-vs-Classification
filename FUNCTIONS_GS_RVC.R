calc_missrate <- function(gt_mat)
{
  col_func <- function(gt_col)
  {
    missrate <- sum(is.na(gt_col)) / length(gt_col)
    return(missrate)
  }

  missrate_vect <- apply(gt_mat, 2, col_func)

  return(missrate_vect)
}

# Calculates the minor allele frequency for every marker in a genotype matrix (coded as c(-1,0,1))
calc_maf_apply <- function(gt_mat, encoding = c(-1, 0, 1))
{
  col_func1 <- function(gt_col)
  {
    allele1_ct <- (sum(gt_col == -1, na.rm = T) * 2) + sum(gt_col == 0, na.rm = T)
    allele2_ct <- (sum(gt_col == 1, na.rm = T) * 2) + sum(gt_col == 0, na.rm = T)

    maf <- min(c(allele1_ct, allele2_ct)) / (sum(!is.na(gt_col))*2)
  }

  col_func2 <- function(gt_col)
  {
    allele1_ct <- (sum(gt_col == 0, na.rm = T) * 2) + sum(gt_col == 1, na.rm = T)
    allele2_ct <- (sum(gt_col == 2, na.rm = T) * 2) + sum(gt_col == 1, na.rm = T)

    maf <- min(c(allele1_ct, allele2_ct)) / (sum(!is.na(gt_col))*2)
  }

  if (all(encoding == c(-1, 0, 1)))
  {
    maf_vect <- apply(gt_mat, 2, col_func1)
  } else if (all(encoding == c(0, 1, 2)))
  {
    maf_vect <- apply(gt_mat, 2, col_func2)
  } else{
    print('Encoding not recognized, returning NULL')
    maf_vect <- NULL
  }

  return(maf_vect)
}

# This is a function that will split data into a list of k-folds
make_CV_sets <- function(list_length, k = 5)
{
  rand_values <- rnorm(list_length)
  k_quantiles <- quantile(rand_values, 0:k/k)
  k_assign <- cut(rand_values, k_quantiles, include.lowest = T, labels = F)

  cv_list <- list()
  for (i in 1:k)
  {
    fold_assignment <- k_assign != i
    cv_list[[i]] <- fold_assignment
  }
  return(cv_list)
}


test_all_models_BLUP_pc_mean_recode <- function(genotypes, phenotype,PCA=NULL,CV=NULL,nIter = 5000, burnIn = 2000, folds = 5)
{

  library(BGLR)
  library(rrBLUP)
  library(caret)
  library(tidyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)

  BGLR_acc_results <- list()
  BGLR_acc_predictions<- list()
  Predictions_ALL<-c()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype[fold_indices,]
    pheno_test=phenotype[-fold_indices,]

    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_train=apply(myGD_train,2,as.numeric)
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_test=apply(myGD_test,2,as.numeric)
    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]
    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train)
    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=predictions,obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(pheno_test,GEBV=predictions,RE=pred_effects,FE=fix_effects)

    rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                   Z = myGD_train,
                                   X = myPCA_train)
    pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
    fix_effects_PC <- myPCA_test  %*% rrBLUP_model_PC$beta
    predictions_PC <- c(pred_effects_PC) + c(fix_effects_PC)
    acc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete")
    sacc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics_PC=postResample(pred=predictions_PC,obs=myY_test)
    results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
    prediction_PC=data.frame(pheno_test,GEBV=predictions_PC,RE=pred_effects_PC,FE=fix_effects_PC)
    Predictions<-cbind(prediction,prediction_PC[,3:5])
    Predictions_ALL=rbind(Predictions_ALL,Predictions)
    BGLR_acc_results[[i]] <- list(results,results_PC)
    #BGLR_acc_predictions[[i]] <- list(results,results_PC,prediction,prediction_PC)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold, na.rm = TRUE)[2:11]
  names(Predictions)[6:8]<-c("GEBV_PC","RE_PC","FE_PC")
  results_ALL=list(results,Predictions_ALL)
  return(results_ALL)
}

test_all_models_BLUP_vs_pc_mean_recode <- function(train_genotypes, train_phenotype,train_PCA=NULL,test_genotypes, test_phenotype,test_PCA=NULL)
{
  library(BGLR)
  library(tidyr)
  library(rrBLUP)
  library(caret)
  library(dplyr)

  # Calculate the GS model using rrBLUP
  myGD_train <- as.matrix(train_genotypes)
  myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_train=apply(myGD_train,2,as.numeric)

  myGD_test <- as.matrix(test_genotypes)
  myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_test=apply(myGD_test,2,as.numeric)

  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]
  myPCA_train <- train_PCA
  myPCA_test <- test_PCA

  gc()
  rrBLUP_model <- mixed.solve(y = myY_train,
                              Z = myGD_train)
  pred_effects <- myGD_test %*% rrBLUP_model$u
  fix_effects <- rrBLUP_model$beta
  predictions <- c(pred_effects) + c(fix_effects)
  acc <- cor(predictions, myY_test, use = "pairwise.complete")
  sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
  metrics=postResample(pred=predictions,obs=myY_test)
  results=c(ACC=acc,SACC=sacc,metrics)
  prediction=data.frame(test_phenotype,GEBV=predictions,RE=pred_effects,FE=fix_effects)

  rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                 Z = myGD_train,
                                 X = myPCA_train)
  pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
  fix_effects_PC <- myPCA_test  %*% rrBLUP_model_PC$beta
  predictions_PC <- c(pred_effects_PC) + c(fix_effects_PC)
  acc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete")
  sacc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete", method = c("spearman"))
  metrics_PC=postResample(pred=predictions_PC,obs=myY_test)
  results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
  prediction_PC=data.frame(test_phenotype,GEBV=predictions_PC,RE=pred_effects_PC,FE=fix_effects_PC)
  Accuracy=c(results,results_PC)
  names(Accuracy)<-c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  Predictions<-cbind(prediction,prediction_PC[,3:5])
  names(Predictions)[6:8]<-c("GEBV_PC","RE_PC","FE_PC")
  results_ALL=list(Accuracy,Predictions)
  return(results_ALL)
}

BGLR_cv_Ordinal_pc_Mean <- function(genotypes, phenotype,PCA=NULL, nIter = 80000, burnIn = 10000, folds = 5)
{
  library(BGLR)
  library(tidyr)
  library(caret)
  # Make the CV list


  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)

  BGLR_acc_results <- list()
  BGLR_acc_metrics<- list()
  Predictions_ALL<-list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype[,2]
    pheno_train[pheno_train=="NaN"]<-NA
    pheno_train=droplevels(pheno_train)
    pheno_train[-fold_indices] <- NA
    # Calculate the GS model using BGLR
    ##Ordinal
    #Without PCs
    BO_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
    BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
    BO_predictions <- predict(BO_results)
    BO_acc <- cor(as.numeric(phenotype[-fold_indices,2]), BO_predictions[-fold_indices],use = "complete.obs")
    BO_sacc <- cor(as.numeric(phenotype[-fold_indices,2]), BO_predictions[-fold_indices],use = "complete.obs", method = c("spearman"))
    DF=BO_results$probs[-fold_indices,]
    probs=colnames(DF)[max.col(replace(DF, cbind(seq_len(nrow(DF)), max.col(DF,ties.method="first")), -Inf), "first")]
    #probs=colnames(DF)[max.col(DF,ties.method="first")]
    BO_acc_cat=cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(probs),use = "complete.obs")
    metrics=postResample(as.factor(probs),phenotype[-fold_indices,2])
    tests=data.frame(obs=phenotype[-fold_indices,2],BO_results$probs[-fold_indices,],pred=factor(probs))
    mets=confusionMatrix(data = tests$pred, reference = tests$obs)
    results=c(ACC=BO_acc,SACC=BO_sacc,C_ACC=BO_acc_cat,metrics)
    #With PCs
    gc()
    BO_ETA_PC<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BL"))
    BO_results_PC <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA_PC, nIter=nIter, burnIn=burnIn)
    BO_predictions_PC <- predict(BO_results_PC)

    BO_acc_PC <- cor(as.numeric(phenotype[-fold_indices,2]), BO_predictions_PC[-fold_indices],use = "complete.obs")
    BO_sacc_PC <- cor(as.numeric(phenotype[-fold_indices,2]), BO_predictions_PC[-fold_indices],use = "complete.obs", method = c("spearman"))
    DF_PC=BO_results_PC$probs[-fold_indices,]
    probs_PC=colnames(DF)[max.col(replace(DF, cbind(seq_len(nrow(DF)), max.col(DF,ties.method="first")), -Inf), "first")]
    #probs_PC=colnames(DF_PC)[max.col(DF_PC,ties.method="first")]
    BO_acc_cat_PC=cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(probs_PC),use = "complete.obs")
    metrics_PC=postResample(as.factor(probs_PC),phenotype[-fold_indices,2])
    tests_PC=data.frame(obs=phenotype[-fold_indices,2],BO_results_PC$probs[-fold_indices,],pred=factor(probs_PC))
    mets_PC=confusionMatrix(data = tests_PC$pred, reference = tests_PC$obs)
    results_PC=c(ACC_PC=BO_acc_PC,SACC_PC=BO_sacc_PC,C_ACC_PC= BO_acc_cat_PC,  metrics_PC)
    #ALL
    Predictions<-data.frame(phenotype[-fold_indices,1],tests,tests_PC[,-1])
    Predictions_ALL[[i]]=list(Predictions)
    BGLR_acc_results[[i]] <- list(results,results_PC)
    BGLR_acc_metrics[[i]] <- list(mets,mets_PC)

  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("R_acc","S_acc","Cat_acc","R2","Kappa","R_acc_PC","S_acc_PC","Cat_acc_PC","R2_PC","Kappa_PC")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[,2:11], na.rm = TRUE)
  results_ALL=list(results,BGLR_acc_metrics,Predictions_ALL)
  return(results_ALL)
}


BGLR_vs_Ordinal_pc_Mean <- function(train_genotypes, train_phenotype,train_PCA=NULL,test_genotypes, test_phenotype,test_PCA=NULL,nIter = 80000, burnIn = 10000)
{
  library(caret)
  library(BGLR)
  library(tidyr)
  phenotype <- c(train_phenotype[,2],test_phenotype[,2])
  pheno_train <- c(train=train_phenotype[,2],test=test_phenotype[,2])
  pheno_train[grep("test",names(pheno_train),value=FALSE)] <- NA
  genotypes<-rbind(train_genotypes,test_genotypes)
  PCA<-rbind(train_PCA,test_PCA)
  # Split into training and testing data
  # Calculate the GS model using BGLR
  ##Ordinal
  #Without PCs
  BO_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
  BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
  BO_predictions <- predict(BO_results)

  BO_acc <- cor(phenotype[-c(1:length(train_phenotype[,2]))], BO_predictions[-c(1:length(train_phenotype[,2]))],use = "complete.obs")
  BO_sacc <- cor(phenotype[-c(1:length(train_phenotype[,2]))], BO_predictions[-c(1:length(train_phenotype[,2]))],use = "complete.obs", method = c("spearman"))

  DF=BO_results$probs[-c(1:length(train_phenotype[,2])),]
  #probs=colnames(DF)[max.col(replace(DF, cbind(seq_len(nrow(DF)), max.col(DF,ties.method="first")), -Inf), "first")]
  probs=colnames(DF)[max.col(DF,ties.method="first")]
  BO_acc_cat=cor(as.numeric(phenotype[-c(1:length(train_phenotype[,2]))]), as.numeric(probs),use = "complete.obs")

  tests=data.frame(obs=as.factor(phenotype[-c(1:length(train_phenotype[,2]))]),BO_results$probs[-c(1:length(train_phenotype[,2])),],pred=factor(probs))



  mynames <- unique(c(levels(tests$pred), levels(tests$obs)))
  tests$pred <- factor(tests$pred, levels = mynames)
  tests$obs <- factor(tests$obs, levels = mynames)

  metrics=postResample(tests$pred, tests$obs)
  mets=confusionMatrix(data = tests$pred, reference = tests$obs)
  results=c(ACC=BO_acc,SACC=BO_sacc,C_ACC=BO_acc_cat,metrics)
  #With PCs
  #gc()
  #BO_ETA_PC<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BL"))
  #BO_results_PC <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA_PC, nIter=nIter, burnIn=burnIn)
  #BO_predictions_PC <- predict(BO_results_PC)

  #BO_acc_PC <- cor(phenotype[-c(1:length(train_phenotype[,2]))], BO_predictions_PC[-c(1:length(train_phenotype[,2]))],use = "complete.obs")
  #BO_sacc_PC <- cor(phenotype[-c(1:length(train_phenotype[,2]))], BO_predictions_PC[-c(1:length(train_phenotype[,2]))],use = "complete.obs", method = c("spearman"))

  #DF=BO_results_PC$probs[-c(1:length(train_phenotype[,2])),]
  #probs=colnames(DF)[max.col(replace(DF, cbind(seq_len(nrow(DF)), max.col(DF,ties.method="first")), -Inf), "first")]
  #probs=colnames(DF)[max.col(DF,ties.method="first")]
  #BO_acc_cat_PC=cor(as.numeric(phenotype[-c(1:length(train_phenotype[,2]))]), as.numeric(probs),use = "complete.obs")

  #tests_PC=data.frame(obs=as.factor(phenotype[-c(1:length(train_phenotype[,2]))]),BO_results_PC$probs[-c(1:length(train_phenotype[,2])),],pred=factor(probs))

  #mynames_PC <- unique(c(levels(tests_PC$pred), levels(tests_PC$obs)))
  #tests_PC$pred <- factor(tests_PC$pred, levels = mynames_PC)
  #tests_PC$obs <- factor(tests_PC$obs, levels = mynames_PC)

  #metrics_PC=postResample(tests_PC$pred,tests_PC$obs)
  #mets_PC=confusionMatrix(data = tests_PC$pred, reference = tests_PC$obs)
  #results_PC=c(ACC_PC=BO_acc_PC,SACC_PC=BO_sacc_PC,C_ACC_PC=BO_acc_cat_PC,metrics_PC)

  #ALL
  #Predictions<-data.frame(test_phenotype[,2],tests,tests_PC[,-1])
  #Accuracy=c(results,results_PC)
  Predictions<-data.frame(test_phenotype[,1],tests)
  Accuracy=c(results)
  names(Predictions)[1]<-c("Genotype")
  #results_ALL=list(Accuracy,Predictions,mets,mets_PC)
  results_ALL=list(Accuracy,Predictions,mets)
  return(results_ALL)
}

Caret_Models_Mean_M <- function(genotypes, phenotype,CV=NULL,model="rf", folds = 5,markers=5000){
  library(caret)
  library(tidyr)
  library(caret)    # for classification and regression training
  library(kernlab)
  library(doParallel)# for fitting SVMs
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)
  svm_results <- list()
  svm_results_metrics<- list()
  Predictions_ALL<-list()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    phenotype[phenotype=="NaN"]<-NA
    phenotype=droplevels(phenotype)

    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]
    myY_train=droplevels(myY_train)
    #samp=sample(1:length(genotypes), markers)
    #m_samp=genotypes[,samp]
    #myGD_train <- m_samp[fold_indices,]
    #myGD_test <- m_samp[-fold_indices,]

    myGD_train <- genotypes[fold_indices,]
    myGD_test <- genotypes[-fold_indices,]

    maf <- calc_maf_apply(myGD_train, encoding = c(0, 1, 2))
    mono_indices <- which(maf == 0)
    #taxa <- row.names(myGD)
    if(length(mono_indices)==0){
      myGD_train = myGD_train
      myGD_test = myGD_test

    }else{

      myGD_train = myGD_train[,-mono_indices]
      myGD_test = myGD_test[,-mono_indices]
    }
    mydata <- data.frame(myY_train, myGD_train)
    mydata$myY_train=as.factor(mydata$myY_train)
    myY_test=as.factor(myY_test)
    colnames(mydata)[1]<-c("Y")
    mydata=mydata
    gc()

    svmFit1 <- train(Y ~ ., data = mydata,
                     method = model
    ,
                     preProcess = c("center", "scale"),
                     trControl = trainControl(method = "cv", number = 10),
                    tuneLength = 10)
    svm.linear_pred1 <- predict(svmFit1, myGD_test)
    svm.linear_pred1_probs <- predict(svmFit1, myGD_test, type = "prob")

    mynames <- unique(c(levels(svm.linear_pred1), levels(myY_test)))
    svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
    myY_test <- factor(myY_test, levels = mynames)

    acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
    sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
    mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
    mets1=confusionMatrix(data = svm.linear_pred1, reference = myY_test, mode = "prec_recall")

    results=c(ACC=acc,SACC=sacc,metrics)
    Predictions_ALL[[i]]=list(Genotypes=phenotype[-fold_indices,1],obs=myY_test,pred=svm.linear_pred1,svm.linear_pred1_probs)
    svm_results[[i]] <- list(results)
    svm_results_metrics[[i]] <- list(mets,mets1)

    #cm$overall

  }
  model_vect <- c("R_acc","S_acc","R2","Kappa")
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(svm_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(svm_results[[i]]))
    SVM_acc_table <- rbind(SVM_acc_table, results_long)
  }
  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold, na.rm = TRUE)[2:5]
  results_ALL=list(results,svm_results_metrics,Predictions_ALL)
  return(results_ALL)
}

Caret_Models_Mean_Matrix <- function(genotypes, phenotype,CV=NULL,model="svmRadial", folds = 5,markers=5000,Matrix="VanRaden",sampling="up"){
  library(caret)
  library(tidyr)   # for classification and regression training
  library(kernlab)
  library(doParallel)
  library(GAPIT3)
  # for fitting SVMs
  #library(DMwR)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)
  svm_results <- list()
  svm_results_metrics<- list()
  Predictions_ALL<-list()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    phenotype[phenotype=="NaN"]<-NA
    phenotype=droplevels(phenotype)

    if(Matrix=="Markers"){
    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]
    myY_train=droplevels(myY_train)
    #samp=sample(1:length(genotypes), markers)
    #m_samp=genotypes[,samp]
    #myGD_train <- m_samp[fold_indices,]
    #myGD_test <- m_samp[-fold_indices,]

    myGD_train <- genotypes[fold_indices,]
    myGD_test <- genotypes[-fold_indices,]

    #maf <- calc_maf_apply(myGD_train, encoding = c(0, 1, 2))
    #mono_indices <- which(maf == 0)
    #taxa <- row.names(myGD)
    #if(length(mono_indices)==0){
      #myGD_train2 = myGD_train
      #myGD_test1 = myGD_test

    #}else{

      #myGD_train2 = myGD_train[,-mono_indices]
      #myGD_test1 = myGD_test[,-mono_indices]
    #}

    myGD_train2=as.matrix(sapply(myGD_train, as.numeric))
    myGD_test1=as.matrix(sapply(myGD_test, as.numeric))

    #set.seed(112233)
    #library(parallel)
    # Calculate the number of cores
    #no_cores <- detectCores()

    #library(doParallel)
    # create the cluster for caret to use
    #cl <- makePSOCKcluster(no_cores)
    #registerDoParallel(cl)

    #gc()
    #,sampling = "down"
    ctrl=trainControl(method = "repeatedcv", repeats = 5, savePredictions="none",sampling = sampling)
    #ctrl=trainControl(method = "cv", number = 10, savePredictions="none",sampling = sampling)
    svmFit1 <- train(myGD_train2,myY_train,
                     method = model,
                     preProcess = c("center", "scale","nzv"),
                     trControl = ctrl,
                     metric="Kappa",
                     tuneLength = 10,
                   allowParallel=TRUE)
    svm.linear_pred1 <- predict(svmFit1, myGD_test1)
    #svm.linear_pred1_probs <- predict(svmFit1, myGD_test, type = "prob")

    mynames <- unique(c(levels(svm.linear_pred1), levels(myY_test)))
    svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
    myY_test <- factor(myY_test, levels = mynames)

    acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
    sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))

    metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
    mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
    #mets1=confusionMatrix(data = svm.linear_pred1, reference = myY_test, mode = "prec_recall")
    results=c(ACC=acc,SACC=sacc,metrics)
    #Predictions_ALL[[i]]=list(Genotypes=phenotype[-fold_indices,1],obs=myY_test,pred=svm.linear_pred1,svm.linear_pred1_probs)
    Predictions_ALL[[i]]=list(Genotypes=phenotype[-fold_indices,1],obs=myY_test,pred=svm.linear_pred1)
    svm_results[[i]] <- list(results)
    #svm_results_metrics[[i]] <- list(mets,mets1)
    svm_results_metrics[[i]] <- list(mets)
}

if(Matrix=="VanRaden"){
  #pheno_train <- phenotype
  #pheno_train[-fold_indices,2] <- NA
  #myY_train <- pheno_train[,2]
  #myY_train=droplevels(myY_train)

  VR=GAPIT.kinship.VanRaden(genotypes)

  myY_train <- phenotype[fold_indices,2]
  myY_test <- phenotype[-fold_indices,2]
  myY_train=droplevels(myY_train)

  myGD_train2 <- VR[fold_indices,]
  myGD_test1 <- VR[-fold_indices,]

  ctrl=trainControl(method = "repeatedcv", repeats = 5, savePredictions="none",sampling = sampling)
  #ctrl=trainControl(method = "cv", number = 10, savePredictions="none",sampling = sampling)
  svmFit1 <- train(myGD_train2,myY_train,
                     method = model,
                     preProcess = c("center", "scale","nzv"),
                     trControl = ctrl,
                     metric="Kappa",
                     tuneLength = 10,
                     allowParallel=TRUE)
svm.linear_pred1 <- predict(svmFit1,myGD_test1)

mynames <- unique(c(levels(svm.linear_pred1), levels(phenotype[,2])))
      svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
      phenotype[,2] <- factor(phenotype[,2], levels = mynames)

acc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete")
sacc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete", method = c("spearman"))
metrics=postResample(pred=svm.linear_pred1,obs=phenotype[-fold_indices,2])
mets=confusionMatrix(data = svm.linear_pred1, reference = phenotype[-fold_indices,2])
results=c(ACC=acc,SACC=sacc,metrics)
#Predictions_ALL[[i]]=list(Genotypes=phenotype[-fold_indices,1],obs=myY_test,pred=svm.linear_pred1,svm.linear_pred1_probs)
Predictions_ALL[[i]]=list(Genotypes=phenotype[-fold_indices,1],obs=phenotype[-fold_indices,2],pred=svm.linear_pred1)
svm_results[[i]] <- list(results)
#svm_results_metrics[[i]] <- list(mets,mets1)
svm_results_metrics[[i]] <- list(mets)
}

if(Matrix=="Zhang"){
  #pheno_train <- phenotype
  #pheno_train[-fold_indices,2] <- NA
  #myY_train <- pheno_train[,2]
  #myY_train=droplevels(myY_train)
  ZZ=GAPIT.kinship.Zhang(genotypes)

  myY_train <- phenotype[fold_indices,2]
  myY_test <- phenotype[-fold_indices,2]
  myY_train=droplevels(myY_train)

  myGD_train2 <- ZZ[fold_indices,]
  myGD_test1 <- ZZ[-fold_indices,]

  ctrl=trainControl(method = "repeatedcv", repeats = 5, savePredictions="none",sampling = sampling)
  svmFit1 <- train(myGD_train2,myY_train,
                     method = model,
                     preProcess = c("center", "scale","nzv"),
                     trControl = ctrl,
                     metric="Kappa",
                     tuneLength = 10,
                     allowParallel=TRUE)
  svm.linear_pred1 <- predict(svmFit1,myGD_test1)

  mynames <- unique(c(levels(svm.linear_pred1), levels(phenotype[,2])))
        svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
        phenotype[,2] <- factor(phenotype[,2], levels = mynames)

  acc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete")
  sacc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete", method = c("spearman"))
  metrics=postResample(pred=svm.linear_pred1,obs=phenotype[-fold_indices,2])
  mets=confusionMatrix(data = svm.linear_pred1, reference = phenotype[-fold_indices,2])
  results=c(ACC=acc,SACC=sacc,metrics)
  #Predictions_ALL[[i]]=list(Genotypes=phenotype[-fold_indices,1],obs=myY_test,pred=svm.linear_pred1,svm.linear_pred1_probs)
  Predictions_ALL[[i]]=list(Genotypes=phenotype[-fold_indices,1],obs=phenotype[-fold_indices,2],pred=svm.linear_pred1)
  svm_results[[i]] <- list(results)
  #svm_results_metrics[[i]] <- list(mets,mets1)
  svm_results_metrics[[i]] <- list(mets)
}

if(Matrix=="Endelman"){
  #pheno_train <- phenotype
  #pheno_train[-fold_indices,2] <- NA
  #myY_train <- pheno_train[,2]
  #myY_train=droplevels(myY_train)
  G=A.mat(genotypes)

  myY_train <- phenotype[fold_indices,2]
  myY_test <- phenotype[-fold_indices,2]
  myY_train=droplevels(myY_train)

  myGD_train2 <- G[fold_indices,]
  myGD_test1 <- G[-fold_indices,]

  ctrl=trainControl(method = "repeatedcv", repeats = 5, savePredictions="none",sampling = sampling)
  svmFit1 <- train(myGD_train2,myY_train,
                     method = model,
                     preProcess = c("center", "scale","nzv"),
                     trControl = ctrl,
                     metric="Kappa",
                     tuneLength = 10,
                     allowParallel=TRUE)
  svm.linear_pred1 <- predict(svmFit1,myGD_test1)

  mynames <- unique(c(levels(svm.linear_pred1), levels(phenotype[,2])))
        svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
        phenotype[,2] <- factor(phenotype[,2], levels = mynames)

  acc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete")
  sacc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete", method = c("spearman"))
  metrics=postResample(pred=svm.linear_pred1,obs=phenotype[-fold_indices,2])
  mets=confusionMatrix(data = svm.linear_pred1, reference = phenotype[-fold_indices,2])
  results=c(ACC=acc,SACC=sacc,metrics)
  #Predictions_ALL[[i]]=list(Genotypes=phenotype[-fold_indices,1],obs=myY_test,pred=svm.linear_pred1,svm.linear_pred1_probs)
  Predictions_ALL[[i]]=list(Genotypes=phenotype[-fold_indices,1],obs=phenotype[-fold_indices,2],pred=svm.linear_pred1)
  svm_results[[i]] <- list(results)
  #svm_results_metrics[[i]] <- list(mets,mets1)
  svm_results_metrics[[i]] <- list(mets)
}



    #stopCluster(cl)
    #registerDoSEQ()
    #cm$overall
    #rm(myGD_train,myGD_test,myY_train,myY_train,svmFit1,svm.linear_pred1,metrics,mets,results)
  }
  model_vect <- c("R_acc","S_acc","R2","Kappa")
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(svm_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(svm_results[[i]]))
    SVM_acc_table <- rbind(SVM_acc_table, results_long)
  }
  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold, na.rm = TRUE)[2:5]
  results_ALL=list(Accuracy=results,CM=svm_results_metrics,Predictions=Predictions_ALL)
  return(results_ALL)
}

Caret_Models_Mean_VS_M <- function(train_genotypes, train_phenotype,test_genotypes, test_phenotype,CV=NULL,model="svmRadial", folds = 5,markers=3000)
{
  gc()
  library(caret)
  library(tidyr)
  library(caret)    # for classification and regression training
  library(kernlab)  # for fitting SVMs
  # Make the CV list
  # Split into training and testing data
  # Calculate the GS model using rrBLUP
  train_phenotype[train_phenotype=="NaN"]<-NA
  train_phenotype=droplevels(train_phenotype)

  test_phenotype[test_phenotype=="NaN"]<-NA
  test_phenotype=droplevels(test_phenotype)

  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]
  samp=sample(1:length(train_genotypes), markers)

  myGD_train <- train_genotypes[,samp]
  myGD_test <- test_genotypes[,samp]

  maf <- calc_maf_apply(myGD_train, encoding = c(0, 1, 2))
  mono_indices <- which(maf == 0)

  if(length(mono_indices)==0){
    myGD_train = myGD_train
    myGD_test = myGD_test
  }else{
    myGD_train = myGD_train[,-mono_indices]
    myGD_test = myGD_test[,-mono_indices]}

  mydata <- data.frame(myY_train, myGD_train)
  colnames(mydata)[1]<-c("Y")
  gc()
  svmFit1 <- train(Y ~ ., data = mydata,
                   method = model,
                   preProcess = c("center", "scale"),
                   trControl = trainControl(method = "cv", number = 10),
                   tuneLength = 10)
  svm.linear_pred1 <- predict(svmFit1, myGD_test)
  svm.linear_pred1_probs <- predict(svmFit1, myGD_test, type = "prob")

  mynames <- unique(c(levels(svm.linear_pred1), levels(myY_test)))
  svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
  myY_test <- factor(myY_test, levels = mynames)

  acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
  sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))
  metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
  mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
  mets1=confusionMatrix(data = svm.linear_pred1, reference = myY_test, mode = "prec_recall")

  results=c(ACC=acc,SACC=sacc,metrics)
  Predictions<-cbind(Genotypes=test_phenotype[,1],obs=myY_test,pred=svm.linear_pred1,svm.linear_pred1_probs)

  results_ALL=list(results,Predictions,mets,mets1)
  results_ALL=list(results,Predictions)
  return(results_ALL)
}

Caret_Models_Mean_VS_Matrix <- function(train_genotypes, train_phenotype,test_genotypes, test_phenotype,CV=NULL,model="svmRadial", folds = 5,markers=5000,Matrix="VanRaden",sampling="up")
{
  library(tidyr)
  library(caret)    # for classification and regression training
  library(kernlab)  # for fitting SVMs
  library(GAPIT3)
  #library(DMwR)
  # Make the CV list
  # Split into training and testing data
  # Calculate the GS model using rrBLUP
  train_phenotype[train_phenotype=="NaN"]<-NA
  train_phenotype=droplevels(train_phenotype)

  test_phenotype[test_phenotype=="NaN"]<-NA
  test_phenotype=droplevels(test_phenotype)

  if(Matrix=="Markers"){
  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]
  #samp=sample(1:length(train_genotypes), markers)

  #myGD_train <- train_genotypes[,samp]
  #myGD_test <- test_genotypes[,samp]

  #maf <- calc_maf_apply(myGD_train, encoding = c(0, 1, 2))
  #mono_indices <- which(maf == 0)

  #if(length(mono_indices)==0){
    #myGD_train = myGD_train
    #myGD_test = myGD_test
  #}else{
    #myGD_train = myGD_train[,-mono_indices]
    #myGD_test = myGD_test[,-mono_indices]}

  myGD_train2=as.matrix(sapply(train_genotypes, as.numeric))
  myGD_test1=as.matrix(sapply(test_genotypes, as.numeric))

  #ctrl=trainControl(method = "repeatedcv", repeats = 5, savePredictions="none",sampling = sampling)
  ctrl=trainControl(method = "repeatedcv", repeats = 5, savePredictions="none")
  #ctrl=trainControl(method = "cv", number = 10,sampling = "up")
  svmFit1 <- train(myGD_train2,myY_train,
                   method = model,
                   preProcess = c("center", "scale","nzv"),
                   trControl = ctrl,
                   metric="Kappa",
                   tuneLength = 10,
                   allowParallel=TRUE)
  svm.linear_pred1 <- predict(svmFit1, myGD_test1)
  #svm.linear_pred1_probs <- predict(svmFit1, myGD_test1, type = "prob")

  mynames <- unique(c(levels(svm.linear_pred1), levels(myY_test)))
  svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
  myY_test <- factor(myY_test, levels = mynames)

  acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
  sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))

  metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
  mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
  #mets1=confusionMatrix(data = svm.linear_pred1, reference = myY_test, mode = "prec_recall")

  results=c(ACC=acc,SACC=sacc,metrics)
  #Predictions<-cbind(Genotypes=test_phenotype[,1],obs=myY_test,pred=svm.linear_pred1,svm.linear_pred1_probs)
  Predictions<-cbind(Genotypes=test_phenotype[,1],obs=myY_test,pred=svm.linear_pred1)
  #results_ALL=list(Accuracy=results,Predictions=Predictions,CM=mets,mets1)
  results_ALL=list(Accuracy=results,CM=mets,Predictions=Predictions)
  }

  if(Matrix=="VanRaden"){

  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]
  myY_train=droplevels(myY_train)

  genotypes=rbind(train_genotypes,test_genotypes)
  VR=GAPIT.kinship.VanRaden(genotypes)
  myGD_train2=VR[c(1:length(train_phenotype[,2])),]
  myGD_test1=VR[-c(1:length(train_phenotype[,2])),]

  #ctrl=trainControl(method = "repeatedcv", repeats = 5, savePredictions="none",sampling = sampling)
  #ctrl=trainControl(method = "cv", number = 10,sampling = "up")
  svmFit1 <- train(myGD_train2,myY_train,
                   method = model,
                   preProcess = c("center", "scale","nzv"),
                   #trControl = ctrl,
                   metric="Kappa",
                   tuneLength = 10,
                   allowParallel=TRUE)
  svm.linear_pred1 <- predict(svmFit1, myGD_test1)
  #svm.linear_pred1_probs <- predict(svmFit1, myGD_test1, type = "prob")

  mynames <- unique(c(levels(svm.linear_pred1), levels(myY_test)))
  svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
  myY_test <- factor(myY_test, levels = mynames)

  acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
  sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))

  metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
  mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
  #mets1=confusionMatrix(data = svm.linear_pred1[-c(1:length(train_phenotype[,2]))], reference = myY_test, mode = "prec_recall")

  results=c(ACC=acc,SACC=sacc,metrics)
  #Predictions<-cbind(Genotypes=test_phenotype[,1],obs=myY_test,pred=svm.linear_pred1[-c(1:length(train_phenotype[,2]))],svm.linear_pred1_probs[-c(1:length(train_phenotype[,2]))])
  Predictions<-cbind(Genotypes=test_phenotype[,1],obs=myY_test,pred=svm.linear_pred1)
  #results_ALL=list(Accuracy=results,Predictions=Predictions,CM=mets,mets1)
  results_ALL=list(Accuracy=results,CM=mets,Predictions=Predictions)
  }

  if(Matrix=="Zhang"){

  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]
  myY_train=droplevels(myY_train)

  genotypes=rbind(train_genotypes,test_genotypes)
  ZZ=GAPIT.kinship.Zhang(genotypes)
  myGD_train2=ZZ[c(1:length(train_phenotype[,2])),]
  myGD_test1=ZZ[-c(1:length(train_phenotype[,2])),]

  ctrl=trainControl(method = "repeatedcv", repeats = 5, savePredictions="none",sampling = sampling)
  #ctrl=trainControl(method = "cv", number = 10,sampling = "up")
  svmFit1 <- train(myGD_train2,myY_train,
                   method = model,
                   preProcess = c("center", "scale","nzv"),
                   trControl = ctrl,
                   metric="Kappa",
                   tuneLength = 10,
                   allowParallel=TRUE)
  svm.linear_pred1 <- predict(svmFit1, myGD_test1)
  #svm.linear_pred1_probs <- predict(svmFit1, myGD_test1, type = "prob")

  mynames <- unique(c(levels(svm.linear_pred1), levels(myY_test)))
  svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
  myY_test <- factor(myY_test, levels = mynames)

  acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
  sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))

  metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
  mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
  #mets1=confusionMatrix(data = svm.linear_pred1[-c(1:length(train_phenotype[,2]))], reference = myY_test, mode = "prec_recall")

  results=c(ACC=acc,SACC=sacc,metrics)
  #Predictions<-cbind(Genotypes=test_phenotype[,1],obs=myY_test,pred=svm.linear_pred1[-c(1:length(train_phenotype[,2]))],svm.linear_pred1_probs[-c(1:length(train_phenotype[,2]))])
  Predictions<-cbind(Genotypes=test_phenotype[,1],obs=myY_test,pred=svm.linear_pred1)
  #results_ALL=list(Accuracy=results,Predictions=Predictions,CM=mets,mets1)
  results_ALL=list(Accuracy=results,CM=mets,Predictions=Predictions)
  }

  if(Matrix=="Endelman"){

  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]
  myY_train=droplevels(myY_train)

  genotypes=rbind(train_genotypes,test_genotypes)
  G=A.mat(genotypes)
  myGD_train2=G[c(1:length(train_phenotype[,2])),]
  myGD_test1=G[-c(1:length(train_phenotype[,2])),]

  ctrl=trainControl(method = "repeatedcv", repeats = 5, savePredictions="none",sampling = sampling)
  #ctrl=trainControl(method = "cv", number = 10,sampling = "up")
  svmFit1 <- train(myGD_train2,myY_train,
                   method = model,
                   preProcess = c("center", "scale","nzv"),
                   trControl = ctrl,
                   metric="Kappa",
                   tuneLength = 10,
                   allowParallel=TRUE)
  svm.linear_pred1 <- predict(svmFit1, myGD_test1)
  #svm.linear_pred1_probs <- predict(svmFit1, myGD_test1, type = "prob")

  mynames <- unique(c(levels(svm.linear_pred1), levels(myY_test)))
  svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
  myY_test <- factor(myY_test, levels = mynames)

  acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
  sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))

  metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
  mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
  #mets1=confusionMatrix(data = svm.linear_pred1[-c(1:length(train_phenotype[,2]))], reference = myY_test, mode = "prec_recall")

  results=c(ACC=acc,SACC=sacc,metrics)
  #Predictions<-cbind(Genotypes=test_phenotype[,1],obs=myY_test,pred=svm.linear_pred1[-c(1:length(train_phenotype[,2]))],svm.linear_pred1_probs[-c(1:length(train_phenotype[,2]))])
  Predictions<-cbind(Genotypes=test_phenotype[,1],obs=myY_test,pred=svm.linear_pred1)
  #results_ALL=list(Accuracy=results,Predictions=Predictions,CM=mets,mets1)
  results_ALL=list(Accuracy=results,CM=mets,Predictions=Predictions)
  }

  return(results_ALL)
}

Caret_Models_Reg_Mean_M <- function(genotypes, phenotype,CV=NULL,model="rf", folds = 5,markers=5000){
  library(caret)
  library(tidyr)
  library(caret)    # for classification and regression training
  library(kernlab)
  library(doParallel)# for fitting SVMs
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)
  svm_results <- list()
  Predictions_ALL<-list()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    phenotype[phenotype=="NaN"]<-NA

    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]

    samp=sample(1:length(genotypes), markers)
    m_samp=genotypes[,samp]
    myGD_train <- m_samp[fold_indices,]
    myGD_test <- m_samp[-fold_indices,]

    maf <- calc_maf_apply(myGD_train, encoding = c(0, 1, 2))
    mono_indices <- which(maf == 0)
    #taxa <- row.names(myGD)
    if(length(mono_indices)==0){
      myGD_train = myGD_train
      myGD_test = myGD_test

    }else{

      myGD_train = myGD_train[,-mono_indices]
      myGD_test = myGD_test[,-mono_indices]
    }
    mydata <- data.frame(myY_train, myGD_train)
    colnames(mydata)[1]<-c("Y")
    gc()
    svmFit1 <- train(Y ~ ., data = mydata,
                     method = model
                     ,
                     preProcess = c("center", "scale"),
                     trControl = trainControl(method = "cv", number = 10),
                     tuneLength = 10)
    svm.linear_pred1 <- predict(svmFit1, myGD_test)
    acc <- cor(svm.linear_pred1, myY_test, use = "pairwise.complete")
    sacc <- cor(svm.linear_pred1, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=svm.linear_pred1,obs=myY_test)

    results=c(ACC=acc,SACC=sacc,metrics)
    Predictions_ALL[[i]]=list(Genotypes=phenotype[-fold_indices,1],obs=myY_test,pred=svm.linear_pred1)
    svm_results[[i]] <- list(results)
  }
  model_vect <- c("R_acc","S_acc","RMSE","R2","MAE")
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(svm_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(svm_results[[i]]))
    SVM_acc_table <- rbind(SVM_acc_table, results_long)
  }
  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold, na.rm = TRUE)[2:6]
  results_ALL=list(results,Predictions_ALL)
  return(results_ALL)
}

Caret_Models_Mean_Reg_VS_M <- function(train_genotypes, train_phenotype,test_genotypes, test_phenotype,CV=NULL,model="rf", folds = 5,markers=5000)
{
  gc()
  library(caret)
  library(tidyr)
  library(caret)    # for classification and regression training
  library(kernlab)  # for fitting SVMs
  # Make the CV list
  # Split into training and testing data
  # Calculate the GS model using rrBLUP
  train_phenotype[train_phenotype=="NaN"]<-NA
  train_phenotype=droplevels(train_phenotype)

  test_phenotype[test_phenotype=="NaN"]<-NA
  test_phenotype=droplevels(test_phenotype)

  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]
  samp=sample(1:length(train_genotypes), markers)

  myGD_train <- train_genotypes[,samp]
  myGD_test <- test_genotypes[,samp]

  maf <- calc_maf_apply(myGD_train, encoding = c(0, 1, 2))
  mono_indices <- which(maf == 0)

  if(length(mono_indices)==0){
    myGD_train = myGD_train
    myGD_test = myGD_test
  }else{
    myGD_train = myGD_train[,-mono_indices]
    myGD_test = myGD_test[,-mono_indices]}

  mydata <- data.frame(myY_train, myGD_train)
  colnames(mydata)[1]<-c("Y")
  gc()
  svmFit1 <- train(Y ~ ., data = mydata,
                   method = model,
                   preProcess = c("center", "scale"),
                   trControl = trainControl(method = "cv", number = 10),
                   tuneLength = 10)
  svm.linear_pred1 <- predict(svmFit1, myGD_test)
  svm.linear_pred1_probs <- predict(svmFit1, myGD_test, type = "prob")
  acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
  sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))
  metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
  mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
  mets1=confusionMatrix(data = svm.linear_pred1, reference = myY_test, mode = "prec_recall")

  results=c(ACC=acc,SACC=sacc,metrics)
  Predictions<-cbind(Genotypes=test_phenotype[,1],obs=myY_test,pred=svm.linear_pred1,svm.linear_pred1_probs)

  results_ALL=list(results,Predictions,mets,mets1)
  return(results_ALL)
}



GLM_cv_Mean <- function(genotypes, phenotype,PCA=NULL,model="GLM",fam="poisson", folds = 5){
  library(mpath)
  library(glmnetUtils)
  library(glmnet)
  library(MASS)
  library(tidyr)
  library(caret)
  library(Metrics)
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)

  svm_results <- list()
  Predictions_ALL<-list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data

    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]

    myGD_train <- genotypes[fold_indices,]
    myGD_test <- genotypes[-fold_indices,]

    myGD_train2=as.matrix(sapply(myGD_train, as.numeric))
    myGD_test1=as.matrix(sapply(myGD_test, as.numeric))

    #glmnet.control(mxitnr = 50)
    if(fam=="poisson"){
      A1_RR=glmnetUtils::cv.glmnet(myGD_train2, myY_train, family = poisson(),
                                   alpha=1,type.measure="mse", standardize = FALSE,
                                   intercept = FALSE)
    }
    if(fam=="quasipoisson"){
      A1_RR=glmnetUtils::cv.glmnet(myGD_train2, myY_train, family = quasipoisson(),
                                   alpha=1,type.measure="mse", standardize = FALSE,
                                   intercept = FALSE)
    }
    if(fam=="negative.binomial"){
      A1_RR=glmnetUtils::cv.glmnet(myGD_train2, myY_train, family = negative.binomial(theta = 5),
                                   alpha=1,type.measure="mse", standardize = FALSE,
                                   intercept = FALSE)
    }
    pred= as.numeric(predict(A1_RR,newx=myGD_test1,s='lambda.min',type='response'))

    acc <- cor(pred, myY_test, use = "pairwise.complete")
    sacc <- cor(pred, myY_test, use = "pairwise.complete", method = c("spearman"))
    rmse=rmse(myY_test,pred)
    r2<- cor(pred, myY_test, use = "pairwise.complete")^2
    #acc_max=max(acc,na.rm = TRUE)
    #sacc_max=max(sacc,na.rm = TRUE)
    #r2_max=max(r2,na.rm = TRUE)
    results=c(ACC=acc,SACC=sacc,RMSE=rmse,rsq=r2)
    Predictions_ALL[[i]]=list(Genotypes=phenotype[-fold_indices,1],obs=myY_test,pred=pred)
    svm_results[[i]] <- list(results)
  }
  model_vect <- c("R_acc","S_acc","RMSE","R2")
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(svm_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(svm_results[[i]]))
    SVM_acc_table <- rbind(SVM_acc_table, results_long)
  }
  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold, na.rm = TRUE)[2:5]
  results_ALL=list(Accurracy=results,Predictions=Predictions_ALL)
  return(results_ALL)
}

GLMP_cv_Mean <- function(genotypes, phenotype,PCA=NULL,model="GLM",fam="poisson", folds = 5){
  library(glmnetUtils)
  library(glmnet)
  library(MASS)
  library(tidyr)
  library(caret)
  library(Metrics)
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)

  svm_results <- list()
  Predictions_ALL<-list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data

    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]

    myGD_train <- genotypes[fold_indices,]
    myGD_test <- genotypes[-fold_indices,]

    myGD_train2=as.matrix(sapply(myGD_train, as.numeric))
    myGD_test1=as.matrix(sapply(myGD_test, as.numeric))


      A1_RR=glmnetUtils::cv.glmnet(myGD_train2, myY_train, family = "poisson",
                                   alpha=1,type.measure="mse")

    pred= as.numeric(predict(A1_RR,newx=myGD_test1,s='lambda.min',type='response'))

    acc <- cor(pred, myY_test, use = "pairwise.complete")
    sacc <- cor(pred, myY_test, use = "pairwise.complete", method = c("spearman"))
    rmse=rmse(myY_test,pred)
    r2<- cor(pred, myY_test, use = "pairwise.complete")^2
    #acc_max=max(acc,na.rm = TRUE)
    #sacc_max=max(sacc,na.rm = TRUE)
    #r2_max=max(r2,na.rm = TRUE)
    results=c(ACC=acc,SACC=sacc,RMSE=rmse,rsq=r2)
    Predictions_ALL[[i]]=list(Genotypes=phenotype[-fold_indices,1],obs=myY_test,pred=pred)
    svm_results[[i]] <- list(results)
  }
  model_vect <- c("R_acc","S_acc","RMSE","R2")
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(svm_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(svm_results[[i]]))
    SVM_acc_table <- rbind(SVM_acc_table, results_long)
  }
  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold, na.rm = TRUE)[2:5]
  results_ALL=list(Accurracy=results,Predictions=Predictions_ALL)
  return(results_ALL)
}

GLM_vs_Mean <- function(train_genotypes, train_phenotype,test_genotypes, test_phenotype,PCA=NULL,model="GLM",fam="poisson", folds = 5){
  library(glmnetUtils)
  library(glmnet)
  library(MASS)
  library(tidyr)
  library(caret)
  library(Metrics)
  # Split into training and testing data

  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]

  myGD_train <- train_genotypes
  myGD_test <- test_genotypes

  myGD_train2=as.matrix(sapply(myGD_train, as.numeric))
  myGD_test1=as.matrix(sapply(myGD_test, as.numeric))

  #glmnet.control(mxitnr = 50)
  if(fam=="poisson"){
    A1_RR=glmnetUtils::cv.glmnet(myGD_train2, myY_train, family = poisson(),
                                 alpha=1,type.measure="mse", standardize = FALSE,
                                 intercept = FALSE)
  }
  if(fam=="quasipoisson"){
    A1_RR=glmnetUtils::cv.glmnet(myGD_train2, myY_train, family = quasipoisson(),
                                 alpha=1,type.measure="mse", standardize = FALSE,
                                 intercept = FALSE)
  }
  if(fam=="negative.binomial"){
    A1_RR=glmnetUtils::cv.glmnet(myGD_train2, myY_train, family = negative.binomial(theta = 5),
                                 alpha=1,type.measure="mse", standardize = FALSE,
                                 intercept = FALSE)
  }
  pred= as.numeric(predict(A1_RR,newx=myGD_test1,s='lambda.min',type='response'))

  acc <- cor(pred, myY_test, use = "pairwise.complete")
  sacc <- cor(pred, myY_test, use = "pairwise.complete", method = c("spearman"))
  rmse=rmse(myY_test,pred)
  r2<- cor(pred, myY_test, use = "pairwise.complete")^2
  #acc_max=max(acc,na.rm = TRUE)
  #sacc_max=max(sacc,na.rm = TRUE)
  #r2_max=max(r2,na.rm = TRUE)
  results=c(ACC=acc,SACC=sacc,RMSE=rmse,rsq=r2)
  Predictions_ALL=cbind(Genotypes=test_phenotype[,1],obs=myY_test,pred=pred)
  results_ALL=list(Accurracy=results,Predictions=Predictions_ALL)
  return(results_ALL)
}

GLMP_vs_Mean <- function(train_genotypes, train_phenotype,test_genotypes, test_phenotype,PCA=NULL,model="GLM",fam="poisson", folds = 5){
  library(glmnetUtils)
  library(glmnet)
  library(MASS)
  library(tidyr)
  library(caret)
  library(Metrics)
  # Split into training and testing data

  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]

  myGD_train <- train_genotypes
  myGD_test <- test_genotypes

  myGD_train2=as.matrix(sapply(myGD_train, as.numeric))
  myGD_test1=as.matrix(sapply(myGD_test, as.numeric))

  A1_RR=glmnetUtils::cv.glmnet(myGD_train2, myY_train, family = "poisson",
                               alpha=1,type.measure="mse")

  pred= as.numeric(predict(A1_RR,newx=myGD_test1,s='lambda.min',type='response'))

  acc <- cor(pred, myY_test, use = "pairwise.complete")
  sacc <- cor(pred, myY_test, use = "pairwise.complete", method = c("spearman"))
  rmse=rmse(myY_test,pred)
  r2<- cor(pred, myY_test, use = "pairwise.complete")^2
  #acc_max=max(acc,na.rm = TRUE)
  #sacc_max=max(sacc,na.rm = TRUE)
  #r2_max=max(r2,na.rm = TRUE)
  results=c(ACC=acc,SACC=sacc,RMSE=rmse,rsq=r2)
  Predictions_ALL=cbind(Genotypes=test_phenotype[,1],obs=myY_test,pred=pred)
  results_ALL=list(Accurracy=results,Predictions=Predictions_ALL)
  return(results_ALL)
}

test_all_models_BLUP_pc_mean_recode_sqrt <- function(genotypes, phenotype,PCA=NULL,CV=NULL,nIter = 5000, burnIn = 2000, folds = 5)
{

  library(BGLR)
  library(rrBLUP)
  library(caret)
  library(tidyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)

  BGLR_acc_results <- list()
  BGLR_acc_predictions<- list()
  Predictions_ALL<-c()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype[fold_indices,]
    pheno_test=phenotype[-fold_indices,]

    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_train=apply(myGD_train,2,as.numeric)
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_test=apply(myGD_test,2,as.numeric)
    myY_train <- sqrt(phenotype[fold_indices,2])
    myY_test <- sqrt(phenotype[-fold_indices,2])

    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train)
    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=predictions,obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(Genotypes=phenotype[-fold_indices,1],Y=myY_test,GEBV=predictions,RE=pred_effects,FE=fix_effects)

    rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                   Z = myGD_train,
                                   X = myPCA_train)
    pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
    fix_effects_PC <- myPCA_test  %*% rrBLUP_model_PC$beta
    predictions_PC <- c(pred_effects_PC) + c(fix_effects_PC)
    acc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete")
    sacc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics_PC=postResample(pred=predictions_PC,obs=myY_test)
    results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
    prediction_PC=data.frame(Genotypes=phenotype[-fold_indices,1],Y=myY_test,GEBV=predictions_PC,RE=pred_effects_PC,FE=fix_effects_PC)
    Predictions<-cbind(prediction,prediction_PC[,3:5])
    Predictions_ALL=rbind(Predictions_ALL,Predictions)
    BGLR_acc_results[[i]] <- list(results,results_PC)
    #BGLR_acc_predictions[[i]] <- list(results,results_PC,prediction,prediction_PC)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold, na.rm = TRUE)[2:11]
  names(Predictions)[6:8]<-c("GEBV_PC","RE_PC","FE_PC")
  results_ALL=list(results,Predictions_ALL)
  return(results_ALL)
}

test_all_models_BLUP_pc_mean_recode_log <- function(genotypes, phenotype,PCA=NULL,CV=NULL,nIter = 5000, burnIn = 2000, folds = 5)
{

  library(BGLR)
  library(rrBLUP)
  library(caret)
  library(tidyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)

  BGLR_acc_results <- list()
  BGLR_acc_predictions<- list()
  Predictions_ALL<-c()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype[fold_indices,]
    pheno_test=phenotype[-fold_indices,]

    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_train=apply(myGD_train,2,as.numeric)
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_test=apply(myGD_test,2,as.numeric)
    myY_train <- log(phenotype[fold_indices,2]+0.000001)
    myY_test <- log(phenotype[-fold_indices,2]+0.000001)
    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train)
    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=predictions,obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(Genotypes=phenotype[-fold_indices,1],Y=myY_test,GEBV=predictions,RE=pred_effects,FE=fix_effects)

    rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                   Z = myGD_train,
                                   X = myPCA_train)
    pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
    fix_effects_PC <- myPCA_test  %*% rrBLUP_model_PC$beta
    predictions_PC <- c(pred_effects_PC) + c(fix_effects_PC)
    acc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete")
    sacc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics_PC=postResample(pred=predictions_PC,obs=myY_test)
    results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
    prediction_PC=data.frame(Genotypes=phenotype[-fold_indices,1],Y=myY_test,GEBV=predictions_PC,RE=pred_effects_PC,FE=fix_effects_PC)
    Predictions<-cbind(prediction,prediction_PC[,3:5])
    Predictions_ALL=rbind(Predictions_ALL,Predictions)
    BGLR_acc_results[[i]] <- list(results,results_PC)
    #BGLR_acc_predictions[[i]] <- list(results,results_PC,prediction,prediction_PC)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold, na.rm = TRUE)[2:11]
  names(Predictions)[6:8]<-c("GEBV_PC","RE_PC","FE_PC")
  results_ALL=list(results,Predictions_ALL)
  return(results_ALL)
}

test_all_models_BLUP_pc_mean_recode_boxcox <- function(genotypes, phenotype,PCA=NULL,CV=NULL,nIter = 5000, burnIn = 2000, folds = 5)
{

  library(BGLR)
  library(rrBLUP)
  library(caret)
  library(tidyr)
  library(forecast)
  boxcox_t=function(vector){
    # to find optimal lambda
    lambda = BoxCox.lambda(vector )
    # now to transform vector
    T_box = BoxCox(vector, lambda)
    return(T_box)
  }
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)

  BGLR_acc_results <- list()
  BGLR_acc_predictions<- list()
  Predictions_ALL<-c()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype[fold_indices,]
    pheno_test=phenotype[-fold_indices,]

    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_train=apply(myGD_train,2,as.numeric)
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_test=apply(myGD_test,2,as.numeric)
    myY_train <- boxcox_t(phenotype[fold_indices,2])
    myY_test <- boxcox_t(phenotype[-fold_indices,2])

    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train)
    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=predictions,obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(Genotypes=phenotype[-fold_indices,1],Y=myY_test,GEBV=predictions,RE=pred_effects,FE=fix_effects)

    rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                   Z = myGD_train,
                                   X = myPCA_train)
    pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
    fix_effects_PC <- myPCA_test  %*% rrBLUP_model_PC$beta
    predictions_PC <- c(pred_effects_PC) + c(fix_effects_PC)
    acc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete")
    sacc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics_PC=postResample(pred=predictions_PC,obs=myY_test)
    results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
    prediction_PC=data.frame(Genotypes=phenotype[-fold_indices,1],Y=myY_test,GEBV=predictions_PC,RE=pred_effects_PC,FE=fix_effects_PC)
    Predictions<-cbind(prediction,prediction_PC[,3:5])
    Predictions_ALL=rbind(Predictions_ALL,Predictions)
    BGLR_acc_results[[i]] <- list(results,results_PC)
    #BGLR_acc_predictions[[i]] <- list(results,results_PC,prediction,prediction_PC)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold, na.rm = TRUE)[2:11]
  names(Predictions)[6:8]<-c("GEBV_PC","RE_PC","FE_PC")
  results_ALL=list(results,Predictions_ALL)
  return(results_ALL)
}

Caret_Models_Regression_Mean_Matrix <- function(genotypes, phenotype,CV=NULL,model="svmRadial", folds = 5,markers=5000,Matrix="Markers",sampling="up"){
  library(caret)
  library(tidyr)   # for classification and regression training
  library(kernlab)
  library(doParallel)
  library(GAPIT3)
  # for fitting SVMs
  #library(DMwR)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)
  svm_results <- list()
  Predictions_ALL<-list()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    if(Matrix=="Markers"){
    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]
    myY_train=droplevels(myY_train)
    #samp=sample(1:length(genotypes), markers)
    #m_samp=genotypes[,samp]
    #myGD_train <- m_samp[fold_indices,]
    #myGD_test <- m_samp[-fold_indices,]

    myGD_train <- genotypes[fold_indices,]
    myGD_test <- genotypes[-fold_indices,]

    #maf <- calc_maf_apply(myGD_train, encoding = c(0, 1, 2))
    #mono_indices <- which(maf == 0)
    #taxa <- row.names(myGD)
    #if(length(mono_indices)==0){
      #myGD_train2 = myGD_train
      #myGD_test1 = myGD_test

    #}else{

      #myGD_train2 = myGD_train[,-mono_indices]
      #myGD_test1 = myGD_test[,-mono_indices]
    #}

    myGD_train2=as.matrix(sapply(myGD_train, as.numeric))
    myGD_test1=as.matrix(sapply(myGD_test, as.numeric))

    #set.seed(112233)
    #library(parallel)
    # Calculate the number of cores
    #no_cores <- detectCores()

    #library(doParallel)
    # create the cluster for caret to use
    #cl <- makePSOCKcluster(no_cores)
    #registerDoParallel(cl)

    #gc()
    #,sampling = "down"
    ctrl=trainControl(method = "repeatedcv", repeats = 5, savePredictions="none")
    svmFit1 <- train(myGD_train2,myY_train,
                     method = model,
                     preProcess = c("center", "scale","nzv"),
                     trControl = ctrl,
                     tuneLength = 10,
                   allowParallel=TRUE)
    svm.linear_pred1 <- predict(svmFit1, myGD_test1)
    #svm.linear_pred1_probs <- predict(svmFit1, myGD_test, type = "prob")
    acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
    sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))

    metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    #Predictions_ALL[[i]]=list(Genotypes=phenotype[-fold_indices,1],obs=myY_test,pred=svm.linear_pred1,svm.linear_pred1_probs)
    Predictions_ALL[[i]]=list(Genotypes=phenotype[-fold_indices,1],obs=myY_test,pred=svm.linear_pred1)
    svm_results[[i]] <- list(results)
}

if(Matrix=="VanRaden"){
  #pheno_train <- phenotype
  #pheno_train[-fold_indices,2] <- NA
  #myY_train <- pheno_train[,2]
  #myY_train=droplevels(myY_train)

  VR=GAPIT.kinship.VanRaden(genotypes)

  myY_train <- phenotype[fold_indices,2]
  myY_test <- phenotype[-fold_indices,2]

  myGD_train2 <- VR[fold_indices,]
  myGD_test1 <- VR[-fold_indices,]

  ctrl=trainControl(method = "repeatedcv", repeats = 5, savePredictions="none")
  svmFit1 <- train(myGD_train2,myY_train,
                     method = model,
                     preProcess = c("center", "scale","nzv"),
                     trControl = ctrl,
                     tuneLength = 10,
                     allowParallel=TRUE)
svm.linear_pred1 <- predict(svmFit1,myGD_test1)

acc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete")
sacc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete", method = c("spearman"))
metrics=postResample(pred=svm.linear_pred1,obs=phenotype[-fold_indices,2])
results=c(ACC=acc,SACC=sacc,metrics)

Predictions_ALL[[i]]=list(Genotypes=phenotype[-fold_indices,1],obs=phenotype[-fold_indices,2],pred=svm.linear_pred1)
svm_results[[i]] <- list(results)

}

if(Matrix=="Zhang"){
  #pheno_train <- phenotype
  #pheno_train[-fold_indices,2] <- NA
  #myY_train <- pheno_train[,2]
  #myY_train=droplevels(myY_train)
  ZZ=GAPIT.kinship.Zhang(genotypes)

  myY_train <- phenotype[fold_indices,2]
  myY_test <- phenotype[-fold_indices,2]

  myGD_train2 <- ZZ[fold_indices,]
  myGD_test1 <- ZZ[-fold_indices,]

  ctrl=trainControl(method = "repeatedcv", repeats = 5, savePredictions="none")
  svmFit1 <- train(myGD_train2,myY_train,
                     method = model,
                     preProcess = c("center", "scale","nzv"),
                     trControl = ctrl,
                     tuneLength = 10,
                     allowParallel=TRUE)
  svm.linear_pred1 <- predict(svmFit1,myGD_test1)


  acc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete")
  sacc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete", method = c("spearman"))
  metrics=postResample(pred=svm.linear_pred1,obs=phenotype[-fold_indices,2])
  results=c(ACC=acc,SACC=sacc,metrics)

  Predictions_ALL[[i]]=list(Genotypes=phenotype[-fold_indices,1],obs=phenotype[-fold_indices,2],pred=svm.linear_pred1)
  svm_results[[i]] <- list(results)
}

if(Matrix=="Endelman"){
  #pheno_train <- phenotype
  #pheno_train[-fold_indices,2] <- NA
  #myY_train <- pheno_train[,2]
  #myY_train=droplevels(myY_train)
  G=A.mat(genotypes)

  myY_train <- phenotype[fold_indices,2]
  myY_test <- phenotype[-fold_indices,2]


  myGD_train2 <- G[fold_indices,]
  myGD_test1 <- G[-fold_indices,]

  ctrl=trainControl(method = "repeatedcv", repeats = 5, savePredictions="none")
  svmFit1 <- train(myGD_train2,myY_train,
                     method = model,
                     preProcess = c("center", "scale","nzv"),
                     trControl = ctrl,
                     tuneLength = 10,
                     allowParallel=TRUE)
  svm.linear_pred1 <- predict(svmFit1,myGD_test1)


  acc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete")
  sacc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete", method = c("spearman"))
  metrics=postResample(pred=svm.linear_pred1,obs=phenotype[-fold_indices,2])
  results=c(ACC=acc,SACC=sacc,metrics)
  #Predictions_ALL[[i]]=list(Genotypes=phenotype[-fold_indices,1],obs=myY_test,pred=svm.linear_pred1,svm.linear_pred1_probs)
  Predictions_ALL[[i]]=list(Genotypes=phenotype[-fold_indices,1],obs=phenotype[-fold_indices,2],pred=svm.linear_pred1)
  svm_results[[i]] <- list(results)

}



    #stopCluster(cl)
    #registerDoSEQ()
    #cm$overall
    #rm(myGD_train,myGD_test,myY_train,myY_train,svmFit1,svm.linear_pred1,metrics,mets,results)
  }
  model_vect <- c("R_acc","S_acc","RMSE","R2","MAE")
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(svm_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(svm_results[[i]]))
    SVM_acc_table <- rbind(SVM_acc_table, results_long)
  }
  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold, na.rm = TRUE)[2:6]
  results_ALL=list(Accuracy=results,Predictions=Predictions_ALL)
  return(results_ALL)
}

Caret_Models_Regression_Mean_VS_Matrix <- function(train_genotypes, train_phenotype,test_genotypes, test_phenotype,CV=NULL,model="svmRadial", folds = 5,markers=5000,Matrix="Markers",sampling="up")
{
  library(tidyr)
  library(caret)    # for classification and regression training
  library(kernlab)  # for fitting SVMs
  library(GAPIT3)
  #library(DMwR)
  # Make the CV list
  # Split into training and testing data
  # Calculate the GS model using rrBLUP

  if(Matrix=="Markers"){
  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]
  #samp=sample(1:length(train_genotypes), markers)

  #myGD_train <- train_genotypes[,samp]
  #myGD_test <- test_genotypes[,samp]

  #maf <- calc_maf_apply(myGD_train, encoding = c(0, 1, 2))
  #mono_indices <- which(maf == 0)

  #if(length(mono_indices)==0){
    #myGD_train = myGD_train
    #myGD_test = myGD_test
  #}else{
    #myGD_train = myGD_train[,-mono_indices]
    #myGD_test = myGD_test[,-mono_indices]}

  myGD_train2=as.matrix(sapply(train_genotypes, as.numeric))
  myGD_test1=as.matrix(sapply(test_genotypes, as.numeric))

  ctrl=trainControl(method = "repeatedcv", repeats = 5, savePredictions="none")
  #ctrl=trainControl(method = "cv", number = 10,sampling = "up")
  svmFit1 <- train(myGD_train2,myY_train,
                   method = model,
                   preProcess = c("center", "scale","nzv"),
                   trControl = ctrl,
                   tuneLength = 10,
                   allowParallel=TRUE)
  svm.linear_pred1 <- predict(svmFit1, myGD_test1)
  #svm.linear_pred1_probs <- predict(svmFit1, myGD_test1, type = "prob")

  acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
  sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))

  metrics=postResample(pred=svm.linear_pred1,obs=myY_test)

  results=c(ACC=acc,SACC=sacc,metrics)
  Predictions<-cbind(Genotypes=test_phenotype[,1],obs=myY_test,pred=svm.linear_pred1)
  #results_ALL=list(Accuracy=results,Predictions=Predictions,CM=mets,mets1)
  results_ALL=list(Accuracy=results,Predictions=Predictions)
  }

  if(Matrix=="VanRaden"){

  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]


  genotypes=rbind(train_genotypes,test_genotypes)
  VR=GAPIT.kinship.VanRaden(genotypes)
  myGD_train2=VR[c(1:length(train_phenotype[,2])),]
  myGD_test1=VR[-c(1:length(train_phenotype[,2])),]

  ctrl=trainControl(method = "repeatedcv", repeats = 5, savePredictions="none")
  #ctrl=trainControl(method = "cv", number = 10,sampling = "up")
  svmFit1 <- train(myGD_train2,myY_train,
                   method = model,
                   preProcess = c("center", "scale","nzv"),
                   trControl = ctrl,
                   tuneLength = 10,
                   allowParallel=TRUE)
  svm.linear_pred1 <- predict(svmFit1, myGD_test1)

  acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
  sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))

  metrics=postResample(pred=svm.linear_pred1,obs=myY_test)

  results=c(ACC=acc,SACC=sacc,metrics)
  Predictions<-cbind(Genotypes=test_phenotype[,1],obs=myY_test,pred=svm.linear_pred1)
  results_ALL=list(Accuracy=results,Predictions=Predictions)
  }

  if(Matrix=="Zhang"){

  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]

  genotypes=rbind(train_genotypes,test_genotypes)
  ZZ=GAPIT.kinship.Zhang(genotypes)
  myGD_train2=ZZ[c(1:length(train_phenotype[,2])),]
  myGD_test1=ZZ[-c(1:length(train_phenotype[,2])),]

  ctrl=trainControl(method = "repeatedcv", repeats = 5, savePredictions="none")
  #ctrl=trainControl(method = "cv", number = 10,sampling = "up")
  svmFit1 <- train(myGD_train2,myY_train,
                   method = model,
                   preProcess = c("center", "scale","nzv"),
                   trControl = ctrl,
                   tuneLength = 10,
                   allowParallel=TRUE)
  svm.linear_pred1 <- predict(svmFit1, myGD_test1)

  acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
  sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))

  metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
  results=c(ACC=acc,SACC=sacc,metrics)
  #Predictions<-cbind(Genotypes=test_phenotype[,1],obs=myY_test,pred=svm.linear_pred1[-c(1:length(train_phenotype[,2]))],svm.linear_pred1_probs[-c(1:length(train_phenotype[,2]))])
  Predictions<-cbind(Genotypes=test_phenotype[,1],obs=myY_test,pred=svm.linear_pred1)
  #results_ALL=list(Accuracy=results,Predictions=Predictions,CM=mets,mets1)
  results_ALL=list(Accuracy=results,Predictions=Predictions)
  }

  if(Matrix=="Endelman"){

  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]

  genotypes=rbind(train_genotypes,test_genotypes)
  G=A.mat(genotypes)
  myGD_train2=G[c(1:length(train_phenotype[,2])),]
  myGD_test1=G[-c(1:length(train_phenotype[,2])),]

  ctrl=trainControl(method = "repeatedcv", repeats = 5, savePredictions="none",sampling = sampling)
  #ctrl=trainControl(method = "cv", number = 10,sampling = "up")
  svmFit1 <- train(myGD_train2,myY_train,
                   method = model,
                   preProcess = c("center", "scale","nzv"),
                   trControl = ctrl,
                   tuneLength = 10,
                   allowParallel=TRUE)
  svm.linear_pred1 <- predict(svmFit1, myGD_test1)

  acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
  sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))

  metrics=postResample(pred=svm.linear_pred1,obs=myY_test)

  results=c(ACC=acc,SACC=sacc,metrics)
  Predictions<-cbind(Genotypes=test_phenotype[,1],obs=myY_test,pred=svm.linear_pred1)

  results_ALL=list(Accuracy=results,Predictions=Predictions)
  }

  return(results_ALL)
}

RE=function(results,results1,results2,si=0.15,nrep,pop,year,loc,model,trait,pheno){
  Pred=list()
  RE_Pred=list()
  RE_Pred1=list()
  RE_Pred2=list()
  RE_PS_Pred=list()
  RE_GS_Pred=list()
  
  for (i in 1:length(results[2,])){
    res=results[2,i][[1]]
    library(janitor)
    res=clean_names(res)
    colnames(res)[1:3]=c("Genotype","Y","GEBV")
    res$GEBV=as.numeric(res$GEBV)
    res$Y=as.numeric(res$Y)
    #means1=results %>%group_by(Genotype) %>% summarise_each(funs(mean))
    means=res %>%group_by(Genotype) %>% summarise_each(funs(mean))
    #Info=results[!duplicated(results[,c(1:7)]),c(1:7)]
    #Predicted_means=left_join(Info,means,by="Genotype")
    
    #means[order(abs(means$l17_IT_R)),]
    ps=means[order(means$Y),]
    #means[order(abs(means$GEBV)),]
    gs=means[order(means$GEBV),]
    differential=round(nrow(ps)*si)
    ps_select=ps[1:differential,]
    gs_select=gs[1:differential,]
    
    REs=(mean(gs_select$Y)-mean(ps$Y))/(mean(ps_select$Y)-mean(ps$Y))
    
    res1=results1[2,i][[1]]
    library(janitor)
    res1=clean_names(res1)
    colnames(res1)[1:3]=c("Genotype","Y1","GEBV1")
    res1$GEBV1=as.numeric(res1$GEBV1)
    res1$Y1=as.numeric(res1$Y1)
    
    res2=results2[2,i][[1]]
    library(janitor)
    res2=clean_names(res2)
    colnames(res2)[1:3]=c("Genotype","Y2","GEBV2")
    res2$GEBV2=as.numeric(res2$GEBV2)
    res2$Y2=as.numeric(res2$Y2)
    
    
    res3=left_join(res,res1,by="Genotype")
    res4=left_join(res3,res2,by="Genotype")
    
    #means1=results %>%group_by(Genotype) %>% summarise_each(funs(mean))
    means1=res4 %>%group_by(Genotype) %>% summarise_each(funs(mean))
    #Info=results[!duplicated(results[,c(1:7)]),c(1:7)]
    #Predicted_means=left_join(Info,means,by="Genotype")
    
    #means[order(abs(means$l17_IT_R)),]
    ps1=means1[order(means1$Y),]
    #means[order(abs(means$GEBV)),]
    gs1=means1[order(means1$GEBV),]
    differential1=round(nrow(ps1)*si)
    ps_select1=ps1[1:differential1,]
    gs_select1=gs1[1:differential1,]
    
    
    
    
    REs1=(mean(gs_select1$Y1)-mean(ps1$Y1))/(mean(ps_select1$Y1)-mean(ps1$Y1))
    
    #means[order(abs(means$l17_IT_R)),]
    ps2=means1[order(means1$Y2),]
    #means[order(abs(means$GEBV)),]
    gs2=means1[order(means1$GEBV),]
    differential2=round(nrow(ps2)*si)
    ps_select2=ps2[1:differential2,]
    gs_select2=gs2[1:differential2,]
    
    REs2=(mean(gs_select2$Y2)-mean(ps2$Y2))/(mean(ps_select2$Y2)-mean(ps2$Y2))
    
    Pred[[i]]=res4
    
    RE_Pred[[i]]=REs
    RE_Pred1[[i]]=REs1
    RE_Pred2[[i]]=REs2
    
    RE_PS_Pred[[i]]=list(raw=ps_select1,adj=ps_select2)
    RE_GS_Pred[[i]]=list(raw=gs_select1,adj=gs_select2)
  }
  
  
  #BGLR_acc_table=c()
  #for (i in 1:length(RE_Pred)){
  #BGLR_acc_table <- rbind(BGLR_acc_table, RE_Pred[[i]])
  #}
  #colnames(BGLR_acc_table)="RE"
  model_vect <- c("RE")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(RE_Pred)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(RE_Pred[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  
  BGLR_acc_table1 <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(RE_Pred1)){
    results_long1 <- data.frame(rep(i, length(model_vect)), model_vect, unlist(RE_Pred1[[i]]))
    BGLR_acc_table1 <- rbind(BGLR_acc_table1, results_long1)
  }
  
  BGLR_acc_table2 <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(RE_Pred2)){
    results_long2 <- data.frame(rep(i, length(model_vect)), model_vect, unlist(RE_Pred2[[i]]))
    BGLR_acc_table2 <- rbind(BGLR_acc_table2, results_long2)
  }
  
  colnames(BGLR_acc_table)<-c("Rep","Metric","Accuracy")
  colnames(BGLR_acc_table1)<-c("Rep","Metric","Accuracy")
  colnames(BGLR_acc_table2)<-c("Rep","Metric","Accuracy")
  po=rep(pop,nrep)
  yr=rep(year,nrep)
  loc=rep(loc,nrep)
  mo=rep(model,nrep)
  tr=rep(trait,nrep)
  ph=rep(pheno,nrep)
  ac=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table)
  ac1=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table1)
  ac2=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table2)
  RE_Pred_ALL=list(REs=ac,Prediction=Pred,PS=RE_PS_Pred,GS=RE_GS_Pred,REs_Raw=ac1,REs_Adj=ac2)
  return(RE_Pred_ALL)
}


RE_Factor=function(results,results1,results2,si=0.15,nrep,pop,year,loc,model,trait,pheno){
  
  BGLR_acc_table_all=list()
  for (i in 1:length(results[3,])){
    BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
    for(j in 1:length(results[3,][[i]])){
      try=data.frame(results[3,][[i]][[j]]$Genotypes,results[3,][[i]][[j]]$obs,results[3,][[i]][[j]]$pred)
      colnames(try)=c("Genotype","Y","GEBV")
      BGLR_acc_table<-rbind(BGLR_acc_table,try)
    }
    BGLR_acc_table_all[[i]] <- BGLR_acc_table
  }
  
  Pred=list()
  RE_Pred=list()
  RE_Pred1=list()
  RE_Pred2=list()
  RE_PS_Pred=list()
  RE_GS_Pred=list()
  
  for (i in 1:length(BGLR_acc_table_all)){
    res=BGLR_acc_table_all[[i]]
    res$GEBV=as.numeric(res$GEBV)
    res$Y=as.numeric(res$Y)
    #means1=results %>%group_by(Genotype) %>% summarise_each(funs(mean))
    means=res %>%group_by(Genotype) %>% summarise_each(funs(mean))
    #Info=results[!duplicated(results[,c(1:7)]),c(1:7)]
    #Predicted_means=left_join(Info,means,by="Genotype")
    
    #means[order(abs(means$l17_IT_R)),]
    ps=means[order(means$Y),]
    #means[order(abs(means$GEBV)),]
    gs=means[order(means$GEBV),]
    differential=round(nrow(ps)*si)
    ps_select=ps[1:differential,]
    gs_select=gs[1:differential,]
    REs=(mean(gs_select$Y)-mean(ps$Y))/(mean(ps_select$Y)-mean(ps$Y))
    
    
    
    res1=results1[2,1][[1]]
    library(janitor)
    res1=clean_names(res1)
    colnames(res1)[1:3]=c("Genotype","Y1","GEBV1")
    res1$GEBV1=as.numeric(res1$GEBV1)
    res1$Y1=as.numeric(res1$Y1)
    
    res2=results2[2,1][[1]]
    library(janitor)
    res2=clean_names(res2)
    colnames(res2)[1:3]=c("Genotype","Y2","GEBV2")
    res2$GEBV2=as.numeric(res2$GEBV2)
    res2$Y2=as.numeric(res2$Y2)
    
    
    res3=left_join(res,res1,by="Genotype")
    res4=left_join(res3,res2,by="Genotype")
    
    #means1=results %>%group_by(Genotype) %>% summarise_each(funs(mean))
    means1=res4 %>%group_by(Genotype) %>% summarise_each(funs(mean))
    #Info=results[!duplicated(results[,c(1:7)]),c(1:7)]
    #Predicted_means=left_join(Info,means,by="Genotype")
    
    #means[order(abs(means$l17_IT_R)),]
    ps1=means1[order(means1$Y),]
    #means[order(abs(means$GEBV)),]
    gs1=means1[order(means1$GEBV),]
    differential1=round(nrow(ps1)*si)
    ps_select1=ps1[1:differential1,]
    gs_select1=gs1[1:differential1,]
    
    
    
    
    REs1=(mean(gs_select1$Y1)-mean(ps1$Y1))/(mean(ps_select1$Y1)-mean(ps1$Y1))
    
    #means[order(abs(means$l17_IT_R)),]
    ps2=means1[order(means1$Y2),]
    #means[order(abs(means$GEBV)),]
    gs2=means1[order(means1$GEBV),]
    differential2=round(nrow(ps2)*si)
    ps_select2=ps2[1:differential2,]
    gs_select2=gs2[1:differential2,]
    
    REs2=(mean(gs_select2$Y2)-mean(ps2$Y2))/(mean(ps_select2$Y2)-mean(ps2$Y2))
    
    Pred[[i]]=res4
    
    RE_Pred[[i]]=REs
    RE_Pred1[[i]]=REs1
    RE_Pred2[[i]]=REs2
    
    RE_PS_Pred[[i]]=list(raw=ps_select1,adj=ps_select2)
    RE_GS_Pred[[i]]=list(raw=gs_select1,adj=gs_select2)
  }
  
  #BGLR_acc_table=c()
  #for (i in 1:length(RE_Pred)){
  #BGLR_acc_table <- rbind(BGLR_acc_table, RE_Pred[[i]])
  #}
  #colnames(BGLR_acc_table)="RE"
  model_vect <- c("RE")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(RE_Pred)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(RE_Pred[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  BGLR_acc_table1 <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(RE_Pred1)){
    results_long1 <- data.frame(rep(i, length(model_vect)), model_vect, unlist(RE_Pred1[[i]]))
    BGLR_acc_table1 <- rbind(BGLR_acc_table1, results_long1)
  }
  
  BGLR_acc_table2 <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(RE_Pred2)){
    results_long2 <- data.frame(rep(i, length(model_vect)), model_vect, unlist(RE_Pred2[[i]]))
    BGLR_acc_table2 <- rbind(BGLR_acc_table2, results_long2)
  }
  
  colnames(BGLR_acc_table)<-c("Rep","Metric","Accuracy")
  colnames(BGLR_acc_table1)<-c("Rep","Metric","Accuracy")
  colnames(BGLR_acc_table2)<-c("Rep","Metric","Accuracy")
  po=rep(pop,nrep)
  yr=rep(year,nrep)
  loc=rep(loc,nrep)
  mo=rep(model,nrep)
  tr=rep(trait,nrep)
  ph=rep(pheno,nrep)
  ac=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table)
  ac1=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table1)
  ac2=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table2)
  RE_Pred_ALL=list(REs=ac,Prediction=Pred,PS=RE_PS_Pred,GS=RE_GS_Pred,REs_Raw=ac1,REs_Adj=ac2)
  return(RE_Pred_ALL)
}

RE_GLM=function(results,results1,results2,si=0.15,nrep,pop,year,loc,model,trait,pheno){
  
  BGLR_acc_table_all=list()
  for (i in 1:length(results[2,])){
    BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
    for(j in 1:length(results[2,][[i]])){
      try=data.frame(results[2,][[i]][[j]]$Genotypes,results[2,][[i]][[j]]$obs,results[2,][[i]][[j]]$pred)
      colnames(try)=c("Genotype","Y","GEBV")
      BGLR_acc_table<-rbind(BGLR_acc_table,try)
    }
    BGLR_acc_table_all[[i]] <- BGLR_acc_table
  }
  
  Pred=list()
  RE_Pred=list()
  RE_Pred1=list()
  RE_Pred2=list()
  RE_PS_Pred=list()
  RE_GS_Pred=list()
  
  for (i in 1:length(BGLR_acc_table_all)){
    res=BGLR_acc_table_all[[i]]
    res$GEBV=as.numeric(res$GEBV)
    res$Y=as.numeric(res$Y)
    #means1=results %>%group_by(Genotype) %>% summarise_each(funs(mean))
    means=res %>%group_by(Genotype) %>% summarise_each(funs(mean))
    #Info=results[!duplicated(results[,c(1:7)]),c(1:7)]
    #Predicted_means=left_join(Info,means,by="Genotype")
    
    #means[order(abs(means$l17_IT_R)),]
    ps=means[order(means$Y),]
    #means[order(abs(means$GEBV)),]
    gs=means[order(means$GEBV),]
    differential=round(nrow(ps)*si)
    ps_select=ps[1:differential,]
    gs_select=gs[1:differential,]
    REs=(mean(gs_select$Y)-mean(ps$Y))/(mean(ps_select$Y)-mean(ps$Y))
    res1= results1[2,i][[1]]
    library(janitor)
    res1=clean_names(res1)
    colnames(res1)[1:3]=c("Genotype","Y1","GEBV1")
    res1$GEBV1=as.numeric(res1$GEBV1)
    res1$Y1=as.numeric(res1$Y1)
    
    res2= results2[2,i][[1]]
    library(janitor)
    res2=clean_names(res2)
    colnames(res2)[1:3]=c("Genotype","Y2","GEBV2")
    res2$GEBV2=as.numeric(res2$GEBV2)
    res2$Y2=as.numeric(res2$Y2)
    
    
    res3=left_join(res,res1,by="Genotype")
    res4=left_join(res3,res2,by="Genotype")
    
    #means1=results %>%group_by(Genotype) %>% summarise_each(funs(mean))
    means1=res4 %>%group_by(Genotype) %>% summarise_each(funs(mean))
    #Info=results[!duplicated(results[,c(1:7)]),c(1:7)]
    #Predicted_means=left_join(Info,means,by="Genotype")
    
    #means[order(abs(means$l17_IT_R)),]
    ps1=means1[order(means1$Y),]
    #means[order(abs(means$GEBV)),]
    gs1=means1[order(means1$GEBV),]
    differential1=round(nrow(ps1)*si)
    ps_select1=ps1[1:differential1,]
    gs_select1=gs1[1:differential1,]
    
    
    
    
    REs1=(mean(gs_select1$Y1)-mean(ps1$Y1))/(mean(ps_select1$Y1)-mean(ps1$Y1))
    
    #means[order(abs(means$l17_IT_R)),]
    ps2=means1[order(means1$Y2),]
    #means[order(abs(means$GEBV)),]
    gs2=means1[order(means1$GEBV),]
    differential2=round(nrow(ps2)*si)
    ps_select2=ps2[1:differential2,]
    gs_select2=gs2[1:differential2,]
    
    REs2=(mean(gs_select2$Y2)-mean(ps2$Y2))/(mean(ps_select2$Y2)-mean(ps2$Y2))
    
    Pred[[i]]=res4
    
    RE_Pred[[i]]=REs
    RE_Pred1[[i]]=REs1
    RE_Pred2[[i]]=REs2
    
    RE_PS_Pred[[i]]=list(raw=ps_select1,adj=ps_select2)
    RE_GS_Pred[[i]]=list(raw=gs_select1,adj=gs_select2)
  }
  
  #BGLR_acc_table=c()
  #for (i in 1:length(RE_Pred)){
  #BGLR_acc_table <- rbind(BGLR_acc_table, RE_Pred[[i]])
  #}
  #colnames(BGLR_acc_table)="RE"
  model_vect <- c("RE")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(RE_Pred)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(RE_Pred[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  BGLR_acc_table1 <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(RE_Pred1)){
    results_long1 <- data.frame(rep(i, length(model_vect)), model_vect, unlist(RE_Pred1[[i]]))
    BGLR_acc_table1 <- rbind(BGLR_acc_table1, results_long1)
  }
  
  BGLR_acc_table2 <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(RE_Pred2)){
    results_long2 <- data.frame(rep(i, length(model_vect)), model_vect, unlist(RE_Pred2[[i]]))
    BGLR_acc_table2 <- rbind(BGLR_acc_table2, results_long2)
  }
  
  colnames(BGLR_acc_table)<-c("Rep","Metric","Accuracy")
  colnames(BGLR_acc_table1)<-c("Rep","Metric","Accuracy")
  colnames(BGLR_acc_table2)<-c("Rep","Metric","Accuracy")
  po=rep(pop,nrep)
  yr=rep(year,nrep)
  loc=rep(loc,nrep)
  mo=rep(model,nrep)
  tr=rep(trait,nrep)
  ph=rep(pheno,nrep)
  ac=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table)
  ac1=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table1)
  ac2=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table2)
  RE_Pred_ALL=list(REs=ac,Prediction=Pred,PS=RE_PS_Pred,GS=RE_GS_Pred,REs_Raw=ac1,REs_Adj=ac2)
  return(RE_Pred_ALL)
}


RE_Factor_BGLR=function(results,results1,results2,si=0.15,nrep,pop,year,loc,model,trait,pheno){
  BGLR_acc_table_all=list()
  for (i in 1:length(results[3,])){
    BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
    for(j in 1:length(results[3,][[i]])){
      try=data.frame(results[3,][[i]][[j]][[1]][,1],results[3,][[i]][[j]][[1]]$obs,results[3,][[i]][[j]][[1]]$pred)
      colnames(try)=c("Genotype","Y","GEBV")
      BGLR_acc_table<-rbind(BGLR_acc_table,try)
    }
    BGLR_acc_table_all[[i]] <- BGLR_acc_table
  }
  
  Pred=list()
  RE_Pred=list()
  RE_Pred1=list()
  RE_Pred2=list()
  RE_PS_Pred=list()
  RE_GS_Pred=list()
  
  for (i in 1:length(BGLR_acc_table_all)){
    res=BGLR_acc_table_all[[i]]
    res$GEBV=as.numeric(as.character(res$GEBV))
    res$Y=as.numeric(as.character(res$Y))
    #means1=results %>%group_by(Genotype) %>% summarise_each(funs(mean))
    means=res %>%group_by(Genotype) %>% summarise_each(funs(mean))
    #Info=results[!duplicated(results[,c(1:7)]),c(1:7)]
    #Predicted_means=left_join(Info,means,by="Genotype")
    
    #means[order(abs(means$l17_IT_R)),]
    ps=means[order(means$Y),]
    #means[order(abs(means$GEBV)),]
    gs=means[order(means$GEBV),]
    differential=round(nrow(ps)*si)
    ps_select=ps[1:differential,]
    gs_select=gs[1:differential,]
    REs=(mean(gs_select$Y)-mean(ps$Y))/(mean(ps_select$Y)-mean(ps$Y))
    res1=results1[2,i][[1]]
    library(janitor)
    res1=clean_names(res1)
    colnames(res1)[1:3]=c("Genotype","Y1","GEBV1")
    res1$GEBV1=as.numeric(res1$GEBV1)
    res1$Y1=as.numeric(res1$Y1)
    
    res2=results2[2,i][[1]]
    library(janitor)
    res2=clean_names(res2)
    colnames(res2)[1:3]=c("Genotype","Y2","GEBV2")
    res2$GEBV2=as.numeric(res2$GEBV2)
    res2$Y2=as.numeric(res2$Y2)
    
    
    res3=left_join(res,res1,by="Genotype")
    res4=left_join(res3,res2,by="Genotype")
    
    #means1=results %>%group_by(Genotype) %>% summarise_each(funs(mean))
    means1=res4 %>%group_by(Genotype) %>% summarise_each(funs(mean))
    #Info=results[!duplicated(results[,c(1:7)]),c(1:7)]
    #Predicted_means=left_join(Info,means,by="Genotype")
    
    #means[order(abs(means$l17_IT_R)),]
    ps1=means1[order(means1$Y),]
    #means[order(abs(means$GEBV)),]
    gs1=means1[order(means1$GEBV),]
    differential1=round(nrow(ps1)*si)
    ps_select1=ps1[1:differential1,]
    gs_select1=gs1[1:differential1,]
    
    
    
    
    REs1=(mean(gs_select1$Y1)-mean(ps1$Y1))/(mean(ps_select1$Y1)-mean(ps1$Y1))
    
    #means[order(abs(means$l17_IT_R)),]
    ps2=means1[order(means1$Y2),]
    #means[order(abs(means$GEBV)),]
    gs2=means1[order(means1$GEBV),]
    differential2=round(nrow(ps2)*si)
    ps_select2=ps2[1:differential2,]
    gs_select2=gs2[1:differential2,]
    
    REs2=(mean(gs_select2$Y2)-mean(ps2$Y2))/(mean(ps_select2$Y2)-mean(ps2$Y2))
    
    Pred[[i]]=res4
    
    RE_Pred[[i]]=REs
    RE_Pred1[[i]]=REs1
    RE_Pred2[[i]]=REs2
    
    RE_PS_Pred[[i]]=list(raw=ps_select1,adj=ps_select2)
    RE_GS_Pred[[i]]=list(raw=gs_select1,adj=gs_select2)
  }
  
  #BGLR_acc_table=c()
  #for (i in 1:length(RE_Pred)){
  #BGLR_acc_table <- rbind(BGLR_acc_table, RE_Pred[[i]])
  #}
  #colnames(BGLR_acc_table)="RE"
  model_vect <- c("RE")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(RE_Pred)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(RE_Pred[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  BGLR_acc_table1 <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(RE_Pred1)){
    results_long1 <- data.frame(rep(i, length(model_vect)), model_vect, unlist(RE_Pred1[[i]]))
    BGLR_acc_table1 <- rbind(BGLR_acc_table1, results_long1)
  }
  
  BGLR_acc_table2 <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(RE_Pred2)){
    results_long2 <- data.frame(rep(i, length(model_vect)), model_vect, unlist(RE_Pred2[[i]]))
    BGLR_acc_table2 <- rbind(BGLR_acc_table2, results_long2)
  }
  
  colnames(BGLR_acc_table)<-c("Rep","Metric","Accuracy")
  colnames(BGLR_acc_table1)<-c("Rep","Metric","Accuracy")
  colnames(BGLR_acc_table2)<-c("Rep","Metric","Accuracy")
  po=rep(pop,nrep)
  yr=rep(year,nrep)
  loc=rep(loc,nrep)
  mo=rep(model,nrep)
  tr=rep(trait,nrep)
  ph=rep(pheno,nrep)
  ac=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table)
  ac1=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table1)
  ac2=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table2)
  RE_Pred_ALL=list(REs=ac,Prediction=Pred,PS=RE_PS_Pred,GS=RE_GS_Pred,REs_Raw=ac1,REs_Adj=ac2)
  return(RE_Pred_ALL)
}


RE_GLM_VS=function(results,results1,results2,si=0.15,nrep,pop,year,loc,model,trait,pheno,test_genotype){
  BGLR_acc_table_all=list()
  for (i in 1:length(results[2,])){
    try=data.frame(results[2,][[i]])
    colnames(try)=c("Genotype","Y","GEBV")
    BGLR_acc_table_all[[i]] <- try
  }
  
  Pred=list()
  RE_Pred=list()
  RE_Pred1=list()
  RE_Pred2=list()
  RE_PS_Pred=list()
  RE_GS_Pred=list()
  
  for (i in 1:length(BGLR_acc_table_all)){
    res=BGLR_acc_table_all[[i]]
    res$GEBV=as.numeric(res$GEBV)
    res$Y=as.numeric(res$Y)
    res$Genotype=as.factor(res$Genotype)
    #means1=results %>%group_by(Genotype) %>% summarise_each(funs(mean))
    means=res %>%group_by(Genotype) %>% summarise_each(funs(mean))
    #Info=results[!duplicated(results[,c(1:7)]),c(1:7)]
    #Predicted_means=left_join(Info,means,by="Genotype")
    
    #means[order(abs(means$l17_IT_R)),]
    ps=means[order(means$Y),]
    #means[order(abs(means$GEBV)),]
    gs=means[order(means$GEBV),]
    differential=round(nrow(ps)*si)
    ps_select=ps[1:differential,]
    gs_select=gs[1:differential,]
    REs=(mean(gs_select$Y)-mean(ps$Y))/(mean(ps_select$Y)-mean(ps$Y))
    
    res1=results1[2,i][[1]]
    library(janitor)
    res1=clean_names(res1)
    colnames(res1)[1:3]=c("Genotype","Y1","GEBV1")
    res1$GEBV1=as.numeric(res1$GEBV1)
    res1$Y1=as.numeric(res1$Y1)
    res1$Genotype=as.factor(res1$Genotype)
    res1=res1[,c(1:3)]
    res2=results2[2,i][[1]]
    library(janitor)
    res2=clean_names(res2)
    colnames(res2)[1:3]=c("Genotype","Y2","GEBV2")
    res2$GEBV2=as.numeric(res2$GEBV2)
    res2$Y2=as.numeric(res2$Y2)
    res2$Genotype=as.factor(res2$Genotype)
    res2=res2[,c(1:3)]
    res$Genotype=test_genotype$Genotype
    res3=left_join(res,res1,by="Genotype")
    res4=left_join(res3,res2,by="Genotype")
    
    #means1=results %>%group_by(Genotype) %>% summarise_each(funs(mean))
    means1=res4 %>%group_by(Genotype) %>% summarise_each(funs(mean))
    #Info=results[!duplicated(results[,c(1:7)]),c(1:7)]
    #Predicted_means=left_join(Info,means,by="Genotype")
    
    #means[order(abs(means$l17_IT_R)),]
    ps1=means1[order(means1$Y),]
    #means[order(abs(means$GEBV)),]
    gs1=means1[order(means1$GEBV),]
    differential1=round(nrow(ps1)*si)
    ps_select1=ps1[1:differential1,]
    gs_select1=gs1[1:differential1,]
    
    
    
    
    REs1=(mean(gs_select1$Y1)-mean(ps1$Y1))/(mean(ps_select1$Y1)-mean(ps1$Y1))
    
    #means[order(abs(means$l17_IT_R)),]
    ps2=means1[order(means1$Y2),]
    #means[order(abs(means$GEBV)),]
    gs2=means1[order(means1$GEBV),]
    differential2=round(nrow(ps2)*si)
    ps_select2=ps2[1:differential2,]
    gs_select2=gs2[1:differential2,]
    
    REs2=(mean(gs_select2$Y2)-mean(ps2$Y2,na.rm=TRUE))/(mean(ps_select2$Y2)-mean(ps2$Y2,na.rm=TRUE))
    
    Pred[[i]]=res4
    
    RE_Pred[[i]]=REs
    RE_Pred1[[i]]=REs1
    RE_Pred2[[i]]=REs2
    
    RE_PS_Pred[[i]]=list(raw=ps_select1,adj=ps_select2)
    RE_GS_Pred[[i]]=list(raw=gs_select1,adj=gs_select2)
  }
  
  #BGLR_acc_table=c()
  #for (i in 1:length(RE_Pred)){
  #BGLR_acc_table <- rbind(BGLR_acc_table, RE_Pred[[i]])
  #}
  #colnames(BGLR_acc_table)="RE"
  model_vect <- c("RE")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(RE_Pred)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(RE_Pred[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  BGLR_acc_table1 <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(RE_Pred1)){
    results_long1 <- data.frame(rep(i, length(model_vect)), model_vect, unlist(RE_Pred1[[i]]))
    BGLR_acc_table1 <- rbind(BGLR_acc_table1, results_long1)
  }
  
  BGLR_acc_table2 <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(RE_Pred2)){
    results_long2 <- data.frame(rep(i, length(model_vect)), model_vect, unlist(RE_Pred2[[i]]))
    BGLR_acc_table2 <- rbind(BGLR_acc_table2, results_long2)
  }
  
  colnames(BGLR_acc_table)<-c("Rep","Metric","Accuracy")
  colnames(BGLR_acc_table1)<-c("Rep","Metric","Accuracy")
  colnames(BGLR_acc_table2)<-c("Rep","Metric","Accuracy")
  po=rep(pop,nrep)
  yr=rep(year,nrep)
  loc=rep(loc,nrep)
  mo=rep(model,nrep)
  tr=rep(trait,nrep)
  ph=rep(pheno,nrep)
  ac=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table)
  ac1=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table1)
  ac2=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table2)
  RE_Pred_ALL=list(REs=ac,Prediction=Pred,PS=RE_PS_Pred,GS=RE_GS_Pred,REs_Raw=ac1,REs_Adj=ac2)
  return(RE_Pred_ALL)
}


RE_Factor_BGLR_VS=function(results,results1,results2,si=0.15,nrep,pop,year,loc,model,trait,pheno){
  BGLR_acc_table_all=list()
  for (i in 1:length(results[2,])){
    try=data.frame(results[2,][[i]]$Genotype,results[2,][[i]]$obs,results[2,][[i]]$pred)
    colnames(try)=c("Genotype","Y","GEBV")
    BGLR_acc_table_all[[i]] <- try
  }
  
  Pred=list()
  RE_Pred=list()
  RE_Pred1=list()
  RE_Pred2=list()
  RE_PS_Pred=list()
  RE_GS_Pred=list()
  
  for (i in 1:length(BGLR_acc_table_all)){
    res=BGLR_acc_table_all[[i]]
    res$GEBV=as.numeric(res$GEBV)
    res$Y=as.numeric(res$Y)
    res$Genotype=as.factor(res$Genotype)
    #means1=results %>%group_by(Genotype) %>% summarise_each(funs(mean))
    means=res %>%group_by(Genotype) %>% summarise_each(funs(mean))
    #Info=results[!duplicated(results[,c(1:7)]),c(1:7)]
    #Predicted_means=left_join(Info,means,by="Genotype")
    
    #means[order(abs(means$l17_IT_R)),]
    ps=means[order(means$Y),]
    #means[order(abs(means$GEBV)),]
    gs=means[order(means$GEBV),]
    differential=round(nrow(ps)*si)
    ps_select=ps[1:differential,]
    gs_select=gs[1:differential,]
    REs=(mean(gs_select$Y)-mean(ps$Y))/(mean(ps_select$Y)-mean(ps$Y))
    
    res1=results1[2,i][[1]]
    library(janitor)
    res1=clean_names(res1)
    colnames(res1)[1:3]=c("Genotype","Y1","GEBV1")
    res1$GEBV1=as.numeric(res1$GEBV1)
    res1$Y1=as.numeric(res1$Y1)
    res1$Genotype=as.factor(res1$Genotype)
    
    res2=results2[2,i][[1]]
    library(janitor)
    res2=clean_names(res2)
    colnames(res2)[1:3]=c("Genotype","Y2","GEBV2")
    res2$GEBV2=as.numeric(res2$GEBV2)
    res2$Y2=as.numeric(res2$Y2)
    res2$Genotype=as.factor(res2$Genotype)
    res$Genotype=res1$Genotype
    res3=left_join(res,res1,by="Genotype")
    res4=left_join(res3,res2,by="Genotype")
    
    #means1=results %>%group_by(Genotype) %>% summarise_each(funs(mean))
    means1=res4 %>%group_by(Genotype) %>% summarise_each(funs(mean))
    #Info=results[!duplicated(results[,c(1:7)]),c(1:7)]
    #Predicted_means=left_join(Info,means,by="Genotype")
    
    #means[order(abs(means$l17_IT_R)),]
    ps1=means1[order(means1$Y),]
    #means[order(abs(means$GEBV)),]
    gs1=means1[order(means1$GEBV),]
    differential1=round(nrow(ps1)*si)
    ps_select1=ps1[1:differential1,]
    gs_select1=gs1[1:differential1,]
    
    
    
    
    REs1=(mean(gs_select1$Y1)-mean(ps1$Y1))/(mean(ps_select1$Y1)-mean(ps1$Y1))
    
    #means[order(abs(means$l17_IT_R)),]
    ps2=means1[order(means1$Y2),]
    #means[order(abs(means$GEBV)),]
    gs2=means1[order(means1$GEBV),]
    differential2=round(nrow(ps2)*si)
    ps_select2=ps2[1:differential2,]
    gs_select2=gs2[1:differential2,]
    
    REs2=(mean(gs_select2$Y2)-mean(ps2$Y2))/(mean(ps_select2$Y2)-mean(ps2$Y2))
    
    Pred[[i]]=res4
    
    RE_Pred[[i]]=REs
    RE_Pred1[[i]]=REs1
    RE_Pred2[[i]]=REs2
    
    RE_PS_Pred[[i]]=list(raw=ps_select1,adj=ps_select2)
    RE_GS_Pred[[i]]=list(raw=gs_select1,adj=gs_select2)
  }
  
  #BGLR_acc_table=c()
  #for (i in 1:length(RE_Pred)){
  #BGLR_acc_table <- rbind(BGLR_acc_table, RE_Pred[[i]])
  #}
  #colnames(BGLR_acc_table)="RE"
  model_vect <- c("RE")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(RE_Pred)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(RE_Pred[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  BGLR_acc_table1 <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(RE_Pred1)){
    results_long1 <- data.frame(rep(i, length(model_vect)), model_vect, unlist(RE_Pred1[[i]]))
    BGLR_acc_table1 <- rbind(BGLR_acc_table1, results_long1)
  }
  
  BGLR_acc_table2 <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(RE_Pred2)){
    results_long2 <- data.frame(rep(i, length(model_vect)), model_vect, unlist(RE_Pred2[[i]]))
    BGLR_acc_table2 <- rbind(BGLR_acc_table2, results_long2)
  }
  
  colnames(BGLR_acc_table)<-c("Rep","Metric","Accuracy")
  colnames(BGLR_acc_table1)<-c("Rep","Metric","Accuracy")
  colnames(BGLR_acc_table2)<-c("Rep","Metric","Accuracy")
  po=rep(pop,nrep)
  yr=rep(year,nrep)
  loc=rep(loc,nrep)
  mo=rep(model,nrep)
  tr=rep(trait,nrep)
  ph=rep(pheno,nrep)
  ac=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table)
  ac1=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table1)
  ac2=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table2)
  RE_Pred_ALL=list(REs=ac,Prediction=Pred,PS=RE_PS_Pred,GS=RE_GS_Pred,REs_Raw=ac1,REs_Adj=ac2)
  return(RE_Pred_ALL)
}
