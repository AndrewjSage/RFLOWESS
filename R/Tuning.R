#'Partition training data
#'
#' Function to partition training cases into folds for a given number of repetitions of cross validation
#' Will typically be called within PartitionData, which is used in parameter tuning step.

#' @param ntrain number of training cases
#' @param cvfolds number of folds for cross validation (must evenly divide number of training cases)
#' @return returns an array containing
#' 1- indices of training cases
#' 2- fold numbers
#' @export

AssignFolds <- function(ntrain, cvfolds){
  if (ntrain%%cvfolds!=0) {stop('Number of folds must evenly divide number of training cases.')}
  foldsize <- ntrain/cvfolds
  #randomly permute training cases
  samp <- sample(1:ntrain, ntrain, replace=FALSE)
  #assign first foldsize permuted cases to test data for cv 1
  TestCases <- sapply(X=1:cvfolds, FUN=function(foldsize, fold){return(samp[(foldsize*(fold-1)+1):(fold*foldsize)])}, foldsize=foldsize)
  return(TestCases)
}

#' Parameter Tuning Weighting Function
#'
#' function to perform cross validation for tuning parameter alpha. Is setup to return OOB errors using all training cases,
#' as well as only non-outliers, for weighted and unweighted cross validation
#' The purpose of this function is to determine how much weight to give to training cases in the validation (or test) set
#' used in parameter tuning for alpha.
#' For a single rep and fold take training cases from training set (TRAIN1) and test (or validation) set(TRAIN 2).
#' Grow forests on TRAIN1 and TRAIN2. Get OOB weights for cases in TRAIN1 using RF1, OOB weights for TRAIN2 using RF2
#' And prediction weights of cases in TRAIN1 on predictions for cases in TRAIN2 using RF1.
#' Get OOB residuals for cases in TRAIN2 using RF2, using both RF and LOWESSRF
#' Compute weights for cases in TRAIN2 using RF and RFL OOB residuals. These are the weights they'll be given in the cross validation used to assess alpha.

#' @param TRAIN set of training data
#' @param TestInd Matrix listing test cases in each fold of cross validation
#' @param fold the fold being considered
#' @param ndsize nodesize parameter value for random forest used in cross validation
#' @param ntreestune number of trees to use for forests involved in parameter tuning
#' @return a list of 6 elements
#'  samp - list of indexes of test (or validation) cases (i.e. TRAIN2)
#'  TRAINY - responses for all training cases
#'  OOBWeights1 - OOB prediction weights for training cases in TRAIN1
#'  PredWeights1 - Prediction weights for predicting test (validation) cases using training cases
#'  BisqwtRF - Weight to be applied to each test or validation case when using cross validation to set tuning parameter, using RF outliers for downweighting
#'  BisqwtRFL- Weight to be applied to each test or validation case when using cross validation to set tuning parameter, using RFLOWESS outliers for downweighting
#' @export

      Get_CVWeights <- function(TRAIN, TestInd, fold, ndsize, ntreestune){
      samp <- TestInd[,fold] #samp is a vector with indices of test cases
      TRAINY <- TRAIN[,ncol(TRAIN)]
      TRAIN1 <- TRAIN[-samp,]   #Subset to grow forest on in CV
      TRAIN2 <- TRAIN[samp,]  #Subset to test on for CV
      #Grow RF on just TRAIN1 subset
      RF1 <- randomForestSRC::rfsrc(data=TRAIN1, Y~., ntree=ntreestune, nodesize=ndsize, forest.wt="oob", membership = T)
      RF1pred <- predict(RF1, newdata=TRAIN2, forest.wt=TRUE, membership = T)
      OOBWeights1 <- RF1$forest.wt
      PredWeights1 <- RF1pred$forest.wt
      #For TRAIN2, determine how to downweight outliers effect on error for purpose of CV
      RF2 <- randomForestSRC::rfsrc(data=TRAIN2, Y~., ntree=ntreestune, nodesize=ndsize, forest.wt="oob", membership = T)
      OOBWeights2 <- RF2$forest.wt
      #RF and RF-LOWESS OOB predictions for cases in TRAIN2 based on cases in TRAIN2
      TRAIN2RFPred <- predict(RF2)$predicted.oob
      TRAIN2RFLPred <- c(LOWESSPred(OOBWeights2, OOBWeights2, TRAINY[samp], alpha=6, tol=10^-4, method="Tukey"))[[1]] #Only want first item in list, which is predicitons
      #Determine how to downweight for purpose of CV
      RFTRAIN2Residsc <- (TRAINY[samp]-TRAIN2RFPred)/(6*median(abs((TRAINY[samp]-TRAIN2RFPred))))  #This is the only place where we use the Y's for TRAIN2. To downweight likely outliers
      BisqwtRF <- c(apply(data.frame(RFTRAIN2Residsc), 1, Tukey))
      #RFLTRAIN2Residsc <- (TRAINY[samp]-TRAIN2RFLPred)/(6*median(abs((TRAINY[samp]-TRAIN2RFLPred))))  #This is the only place where we use the Y's for TRAIN2. To downweight likely outliers
      #BisqwtRFL <- apply(data.frame(RFLTRAIN2Residsc), 1, Tukey)
      #return(list(samp, TRAINY, OOBWeights1, PredWeights1, BisqwtRF, BisqwtRFL))
      return(list(samp, TRAINY, OOBWeights1, PredWeights1, BisqwtRF))
            }


#' Parameter tuning function
#' #'
#' For each candidate parameter value, calculate the CV error resulting from growing a forest on TRAIN1,
#' predicting TRAIN2, and downweighting contributions in TRAIN2 based on either RF residuals or RFL residuals
#' using OOB predictions when TRAIN2 was predicted using RF grown on TRAIN2. In cross validation, this is
#' done for each candidate value, within each fold, within each rep.

#' @param  OOBWeights1 OOB prediction weights for training cases in TRAIN1
#' @param  PredWeights1 Prediction weights for predicting test (validation) cases using training cases
#' @param  TRAINY responses for all training cases
#' @param  samp list of indexes of test (or validation) cases (i.e. TRAIN2)
#' @param  parvec vector of candidate values for tuning parameter alpha
#' @param  ind index of parameter tuning vector to work on
#' @param  tol maximal change in interation for LOWESSRF weights in cross validation
#' @param  BisqwtRF Weight to be applied to each test or validation case when using cross validation to set tuning parameter, using RF outliers for downweighting
#' @param  OutlierInd Vector of zeros and ones indicating whether training cases came from contaminating distribution
#' @return Returns a 6 by 3 matrix with errors. Rows indicate type of weighting applied to errors in TRAIN2.
#'          1-MSE without downweighting outliers in CV error
#'          2-MAPE without downweighting outliers in CV error
#'          3-MSE downweighting outliers according to BisqwtRF
#'          4-MAPE downweighting outliers according to BisqwtRF
#'          5-MSE downweighting outliers according to BisqwtRFL
#'          6-MAPE downweighting outliers according to BisqwtRFL
#'         Columns represent error for all cases (1), only outliers (2), only nonoutliers(3)
#' @export

 #  Calculate_CV_Error <- function(OOBWeights1, PredWeights1, TRAINY, samp, parvec, ind, tol=10^-4, BisqwtRF, BisqwtRFL, OutlierInd){
    Calculate_CV_Error <- function(OOBWeights1, PredWeights1, TRAINY, samp, parvec, ind, tol=10^-4, BisqwtRF, BisqwtRFL, OutlierInd){
      CVERR <- array(NA, dim=c(4,3))
      #RF-LOWESS predictions for cases in TRAIN2 based on cases in TRAIN1 for different alpha
      #PredWeights1 are weights from TRAIN1 for predictions on TRAIN2
      LCVPred <- c(LOWESSPred(OOBWeights1, PredWeights1, TRAINY[-samp], alpha=parvec[ind], method="Tukey", tol=tol)[[1]])
      #Evaluate on only TRAIN2  i.e. cases in samp
      CVERR[1, 1] <- mean((LCVPred-TRAINY[samp])^2)   #MSE without downweighting outliers in CV error
      CVERR[2, 1] <- mean(abs(LCVPred-TRAINY[samp]))  #MAPE without downweighting outliers in CV error
      CVERR[3, 1] <- mean((LCVPred-TRAINY[samp])^2*BisqwtRF)  #MSE downweighting outliers according to BisqwtRF
      CVERR[4, 1] <- mean(abs(LCVPred-TRAINY[samp])*BisqwtRF) #MAPE downweighting outliers according to BisqwtRF
#      CVERR[5, 1] <- mean((LCVPred-TRAINY[samp])^2*BisqwtRFL) #MSE downweighting outliers according to BisqwtRFL
#      CVERR[6, 1] <- mean(abs(LCVPred-TRAINY[samp])*BisqwtRFL) #MAPE downweighting outliers according to BisqwtRFL
      #error on only outlier cases
      if(sum(OutlierInd[samp])>0){
        CVERR[1, 2] <- mean(((LCVPred-TRAINY[samp])^2)[OutlierInd[samp]==1])
        CVERR[2, 2] <- mean(abs(LCVPred-TRAINY[samp])[OutlierInd[samp]==1])
        CVERR[3, 2] <- mean(((LCVPred-TRAINY[samp])^2*(BisqwtRF))[OutlierInd[samp]==1])
        CVERR[4, 2] <- mean((abs(LCVPred-TRAINY[samp])*(BisqwtRF))[OutlierInd[samp]==1])
#        CVERR[5, 2] <- mean(((LCVPred-TRAINY[samp])^2*(BisqwtRFL))[OutlierInd[samp]==1])
#        CVERR[6, 2] <- mean((abs(LCVPred-TRAINY[samp])*(BisqwtRFL))[OutlierInd[samp]==1])
        }
      #error on only non-outliers
      CVERR[1, 3] <- mean(((LCVPred-TRAINY[samp])^2)[OutlierInd[samp]==0])
      CVERR[2, 3] <- mean(abs(LCVPred-TRAINY[samp])[OutlierInd[samp]==0])
      CVERR[3, 3] <- mean(((LCVPred-TRAINY[samp])^2*(BisqwtRF))[OutlierInd[samp]==0])
      CVERR[4, 3] <- mean((abs(LCVPred-TRAINY[samp])*(BisqwtRF))[OutlierInd[samp]==0])
#      CVERR[5, 3] <- mean(((LCVPred-TRAINY[samp])^2*(BisqwtRFL))[OutlierInd[samp]==0])
#      CVERR[6, 3] <- mean(abs((LCVPred-TRAINY[samp])*(BisqwtRFL))[OutlierInd[samp]==0])
  return(CVERR)
}


    #' Apply Parameter Tuning
    #'
    #' For each fold in the cross validation, calls Get_CVWeights to split into training and validation data
    #' and get weights, then applies Calculate_CV_Error over all candidate values for alpha. Returns array of dimensions
    #' 6 by 3 by length(parvec) containing errors. Dimensions index (1) type of error calculated
    #'          1-MSE without downweighting outliers in CV error
    #'          2-MAPE without downweighting outliers in CV error
    #'          3-MSE downweighting outliers according to BisqwtRF
    #'          4-MAPE downweighting outliers according to BisqwtRF
    #'          5-MSE downweighting outliers according to BisqwtRFL
    #'          6-MAPE downweighting outliers according to BisqwtRFL
    #' (2) applied to all cases in TRAIN2, outliers only, nonoutliers only
    #' (3) index of alpha from parvec

    #' @param  TRAIN matrix of training cases with response in last column
    #' @param  TestInd Matrix indicating which cases are in test (validation) set, TRAIN2
    #' @param  fold number of fold within cross validation
    #' @param  ndsize nodesize random forest tuning parameter for cross validation
    #' @param  OutlierInd Vector of zeros and ones indicating whether training cases came from contaminating distribution
    #' @param  parvec vector of candidate values for tuning parameter alpha
    #' @param  tol maximal change in interation for LOWESSRF weights in cross validation
    #' @param  ntreestune number of trees to use for forests involved in parameter tuning
    #' @return Returns array of dimensions
    #' 6 by 3 by length(parvec) containing errors. Dimensions index (1) type of error calculated
    #'          1-MSE without downweighting outliers in CV error
    #'          2-MAPE without downweighting outliers in CV error
    #'          3-MSE downweighting outliers according to BisqwtRF
    #'          4-MAPE downweighting outliers according to BisqwtRF
    #'          5-MSE downweighting outliers according to BisqwtRFL
    #'          6-MAPE downweighting outliers according to BisqwtRFL
    #' (2) applied to all cases in TRAIN2, outliers only, nonoutliers only
    #' (3) index of alpha from parvec
    #' @export

  Calc_fold_CV_Error <- function(TRAIN, TestInd, fold, ndsize, OutlierInd, parvec, tol=10^-4, ntreestune=ntreestune){
    CVWeightInfo <- Get_CVWeights(TRAIN=TRAIN, TestInd=TestInd, fold=fold, ndsize=ndsize, ntreestune=ntreestune)  #done for a specific fold
    CVERR <- lapply(X=1:length(parvec), FUN = Calculate_CV_Error, OOBWeights1=CVWeightInfo[[3]],
                    PredWeights1=CVWeightInfo[[4]],TRAINY=CVWeightInfo[[2]], samp=CVWeightInfo[[1]],
                    parvec=parvec, tol=tol,  BisqwtRF=CVWeightInfo[[5]],  OutlierInd=OutlierInd)
    CVERR <- abind::abind(CVERR, along=3) #convert to an array with dims 6 x 3 x length(parvec)
    return(CVERR)
  }


#' Evaluate Tuning Candidates
#'
#' This function applies Calc_fold_CV_Error over all folds of a cross validation. Use replicate if running multiple
#' repetitions of cross-validation

  #' @param  TRAIN matrix of training cases with response in last column
  #' @param  OutlierInd Vector of zeros and ones indicating whether training cases came from contaminating distribution
  #' @param  cvfolds number of folds to perform in cross validation
  #' @param  parvec vector of candidate values for tuning parameter alpha
  #' @param  ndsize nodesize random forest tuning parameter for cross validation
  #' @param  ntreestune number of trees to use for forests involved in parameter tuning
  #' @param  tol maximal change in interation for LOWESSRF weights in cross validation
  #' @return Returns array of dimensions
  #' 6 by 3 by length(parvec) by number of folds, containing errors. Dimensions index (1) type of error calculated
  #'          1-MSE without downweighting outliers in CV error
  #'          2-MAPE without downweighting outliers in CV error
  #'          3-MSE downweighting outliers according to BisqwtRF
  #'          4-MAPE downweighting outliers according to BisqwtRF
  #'          5-MSE downweighting outliers according to BisqwtRFL
  #'          6-MAPE downweighting outliers according to BisqwtRFL
  #' (2) applied to all cases in TRAIN2, outliers only, nonoutliers only
  #' (3) index of alpha from parvec
  #' (4) index of fold in cross validation
  #' @export

  Evaluate_Tuning_Candidates <- function(TRAIN, OutlierInd, cvfolds, parvec, ndsize, ntreestune=100, tol=10^-4){
  ntrain <- nrow(TRAIN)
  TestInd <- AssignFolds(ntrain, cvfolds)
  TRAINY <- TRAIN[,ncol(TRAIN)]
  CVArray <- lapply(X=1:cvfolds, Calc_fold_CV_Error, TRAIN=TRAIN, TestInd=TestInd, ndsize=ndsize, OutlierInd=OutlierInd, parvec=parvec, tol=tol, ntreestune=ntreestune )
  CVArray <- abind::abind(CVArray, along=4) #convert to an array with dims 6 x 3 x length(parvec) x number of folds
  return(CVArray)
  }


#' Predict RFLOWESS Tuning
#'
#' This is a function to be called inside TuneMultifoldCV. It makes the RFLOWESS predictions for a specific alpha
#' and records them in a matrix LPREDERR.

#' @param  OOBWeights matrix of training cases with response in last column
#' @param  PredWeights Vector of zeros and ones indicating whether training cases came from contaminating distribution
#' @param  TRAINY number of folds to perform in cross validation
#' @param  TEST Test data
#' @param  tol maximal change in interation for LOWESSRF weights in cross validation
#' @param  ind index of parameter vector to use in tuning
#' @param  parvec vector of candidate values for tuning parameter alpha
#' @param  ndsize nodesize random forest tuning parameter for cross validation
#' @param  method should Tukey or Huber weighting function be used?
#' @return Returns:
#'         LPREDERR: A length 2 vector containing MSPE and MAPE on test data using alpha parameter specified
#'         LWeights: A ntest by ntrain matrix containing LOWESSRF weights for test cases
#'         LIter: Number of iterations until convergence of LOWESSRF algorithm
#' @export

MakeLowessPred <- function(OOBWeights, PredWeights, TRAINY, TEST,  tol=tol, ind, parvec,  method="Tukey"){
  LPREDERR <- c(NA, NA)
  LPredInfo <- c(LOWESSPred(OOBWeights, PredWeights, TRAINY, tol=tol, alpha=parvec[ind], method="Tukey"))
  LPred <- LPredInfo[[1]]  #new predictions
  LWeights <- LPredInfo[[2]]  #adjusted weights
  LIter <- LPredInfo[[3]]  #iterations until convergence
  LPREDERR[1] <- mean((LPred-TEST[, ncol(TEST)])^2)  #Just do MSPE or MAPE, not weighted average on Test data
  LPREDERR[2] <- mean(abs(LPred-TEST[, ncol(TEST)]))
  return(list(LPREDERR, LWeights, LIter))
}


#' Main Parameter Tuning Function
#'
#' This function takes in information from a forest grown on a test and training set and calls Evaluate_Tuning_Candidates()
#' to assess performance for different values of alpha and then makes predictions on the test set for each alpha. It
#' records the performance on test data for each alpha chosen. We can check how the value of alpha suggested by CV performs.

#' @param  TRAIN matrix of training cases with response in last column
#' @param  TEST matrix of test cases with response in last column
#' @param  OOBWeights matrix of training cases with response in last column
#' @param  PredWeights Vector of zeros and ones indicating whether training cases came from contaminating distribution
#' @param  OutlierInd Vector of zeros and ones indicating whether training cases came from contaminating distribution
#' @param  ntrees number of trees for forest used in predictions of test data
#' @param  ntreestune number of trees for forests used to set tuning parameter alpha
#' @param  ndsize nodesize random forest tuning parameter for cross validation
#' @param  parvec vector of candidate values for tuning parameter alpha
#' @param  ntreestune number of trees to use for forests involved in parameter tuning
#' @param  cvreps number of repetitions of cross validation
#' @param  cvfolds number of folds to perform in cross validation
#' @param  tol  maximal change in interation for LOWESSRF weights in cross validation
#' @return Returns list of eight elements containing:
#'        1. LOOBERR: 6 by 3 by length (parvec) array containing OOB error
#'                 First index represents
#'                       1-MSE without downweighting outliers in CV error
#'                       2-MAPE without downweighting outliers in CV error
#'                       3-MSE downweighting outliers according to BisqwtRF
#'                       4-MAPE downweighting outliers according to BisqwtRF
#'                       5-MSE downweighting outliers according to BisqwtRFL
#'                       6-MAPE downweighting outliers according to BisqwtRFL
#'                  Second index represents whether to set alpha using  all cases in TRAIN2, outliers only, nonoutliers only
#'        2. LPREDERR: length(parvec) by 2 matrix containing MSPE in first column and MAPE in 2nd column for
#'                     predictions on test cases using each value of alpha in the parvec
#'        3. ChosenPars: 6 by 3 matrix containing the values of alpha selected using each approach with indices same as for LOOBERR
#'        4. BestPars: Vector of length 2, containing value of alpha that performed best on test data for MSPE(1) and MAPE(2)
#'        5. Diff: 6 by 3 matrix containing the difference between the value of alpha chosen by each weighting crietera and function applied
#'                 to the training data, and the one performing best on the test data
#'        6. ERR: 3 by 6 matrix containing error on test data (either MSPE in columns 1,3,5 or MAPE in columns 2,4,6) where
#'                columns denote the parameter choice determined using the 6 criteria explained in the rows for LOOBERR
#'                and rows represent using (1) All training cases, (2) Only Outliers, (3) Only nonoutliers in tuning.
#' @export
#'

TuneMultifoldCV=function(TRAIN, TEST, OOBWeights, PredWeights, OutlierInd, ntrees,ntreestune, ndsize, parvec, cvreps, cvfolds, tol){
  TRAINY <- TRAIN[,ncol(TRAIN)]
  ChosenPars <- array(NA, dim=c(4,3)) #3 is for contribution of all, outliers, and non-outliers
  BestPars <- array(NA, dim=c(2))
  Diff <- array(NA, dim=c(2,4,3) )#second dimension is for MSE or MAPE #3 is for contribution of all, outliers, and non-outliers
  LPREDERR <- array(dim=c(length(parvec), 2)) #Just one entry for MSE and one for MAPE
  ERR <- array(NA, dim=c(3,4)) #dims whether all, outliers, or nonoutliers are used, and 6 different types of downweighing
  # Get Errors for tuning by doing CV on training data. Returns errors for each parameter
  Tunerep <- replicate(cvreps, Evaluate_Tuning_Candidates(TRAIN=TRAIN, OutlierInd=OutlierInd, cvfolds=cvfolds, parvec=parvec, ndsize=ndsize, ntreestune=ntreestune))
  LOOBERR <- apply(Tunerep, c(1,2,3), mean)
  #Get RFLOWESS predictions for each test case using each tuning parameter estimate
  LPREDINFO <- sapply(X=1:length(parvec), MakeLowessPred, OOBWeights=OOBWeights,
                     PredWeights=PredWeights, TRAINY=TRAINY, TEST=TEST, tol=tol,
                     parvec=parvec, method="Tukey")
  LPREDERR <- do.call(rbind, LPREDINFO[1,]) #matrix containing MSPE and MAPE on test data using each alpha
  LWeights <- LPREDINFO[2,] #List of ntest by ntrain matrices, containing LOWESSRF weights for each value of alpha
  LIter <- unlist(LPREDINFO[3,])  #vector containing number of iterations needed using each alpha

  ChosenPars[,1] <- parvec[apply(LOOBERR[,1,], 1, which.min)] #Parameters chosen using all cases
  if(sum(is.na(LOOBERR[,,2]))==0){  #if there is contamination
    ChosenPars[,2] <- parvec[apply(LOOBERR[,2,], 1, which.min)]  #Parameters chosen using only outliers
  }
  ChosenPars[,3] <- parvec[apply(LOOBERR[,3,], 1, which.min)] #Parameters chosen using nonoutliers
  #Note sum of LOOBERR[rep,,,2] and LOOBERR[rep,,,3] is LOOBERR[rep,,,1]
  BestPars <- parvec[apply(LPREDERR, 2, which.min)]  #Best parameters for test data
  Diff[1,,] <- ChosenPars-BestPars[1]
  Diff[2,,] <- ChosenPars-BestPars[2]
  #Find prediction error corresponding to parameter choice we would have made
  #Prediction error using all cases to tune
  ERR[1,1] <- LPREDERR[which.min(LOOBERR[1,1,]),1]  #first 1 is for CVtype-unweighted/MSE, second is for using all training cases, 3rd is for MSE
  ERR[1,2] <- LPREDERR[which.min(LOOBERR[2,1,]),2] #first 2 is for CVtype-unweighted/MAPE, second is for using all training cases, 3rd is for MAPE
  ERR[1,3] <- LPREDERR[which.min(LOOBERR[3,1,]),1] #using RF downweighting
  ERR[1,4] <- LPREDERR[which.min(LOOBERR[4,1,]),2]
  #ERR[1,5] <- LPREDERR[which.min(LOOBERR[5,1,]),1] #using LOWESS downweighting
  #ERR[1,6] <- LPREDERR[which.min(LOOBERR[6,1,]),2]
  #Prediction error only outliers in TRAIN2 to tune-not sure why we'd want to do this
  if(sum(is.na(LOOBERR[,,2]))==0){
    ERR[2,1] <- LPREDERR[which.min(LOOBERR[1,2,]),1]  #first 1 is for CVtype-unweighted/MSE, second is for using all training cases, 3rd is for MSE
    ERR[2,2] <- LPREDERR[which.min(LOOBERR[2,2,]),2] #first 2 is for CVtype-unweighted/MAPE, second is for using all training cases, 3rd is for MAPE
    ERR[2,3] <- LPREDERR[which.min(LOOBERR[3,2,]),1]
    ERR[2,4] <- LPREDERR[which.min(LOOBERR[4,2,]),2]
   # ERR[2,5] <- LPREDERR[which.min(LOOBERR[5,2,]),1]
  #  ERR[2,6] <- LPREDERR[which.min(LOOBERR[6,2,]),2]
  }
  #Prediction error using only nonoutliers  in TRAIN2 to tune
  ERR[3,1] <- LPREDERR[which.min(LOOBERR[1,3,]),1]  #first 1 is for CVtype-unweighted/MSE, second is for using all training cases, 3rd is for MSE
  ERR[3,2] <- LPREDERR[which.min(LOOBERR[2,3,]),2] #first 2 is for CVtype-unweighted/MAPE, second is for using all training cases, 3rd is for MAPE
  ERR[3,3] <- LPREDERR[which.min(LOOBERR[3,3,]),1]
  ERR[3,4] <- LPREDERR[which.min(LOOBERR[4,3,]),2]
  #ERR[3,5] <- LPREDERR[which.min(LOOBERR[5,3,]),1]
  #ERR[3,6] <- LPREDERR[which.min(LOOBERR[6,3,]),2]
  #return(list(LOOBERR, LPREDERR, ChosenPars, BestPars, Diff, ERR, LWeights, LIter ))
  return(list(LOOBERR, LPREDERR, ChosenPars, BestPars, Diff, ERR))
  }
