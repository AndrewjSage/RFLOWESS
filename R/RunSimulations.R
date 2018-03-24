
#' Make Predictions
#'
#' Function to make predictions for all methods tested in simulation

#' @param DATA  object from generate_RLdata or generate_LMdata
#' @param ntrees  number of trees
#' @param ndsize nodesize
#' @param ntreestune  number of trees to use for tuning alpha
#' @param parvec vector of candidate values for tuning parameter alpha
#' @param cvreps number of repetitions to perform in cross validation
#' @param cvfolds number of folds to perform in cross validation
#' @param  tol maximal change in interation for LOWESSRF weights in cross validation
#' @return returns a list of 4 items
#'        1. Datasets (TRAIN, TEST, and Outlier Indicator)
#'        2. Matrix of 16 columns giving different predictions. Last column is true Y.
#'        3. Number of iterations
#'        4. Output from TuneMultifoldCV (a list of 8 items itself)
#' @export


Make_All_Preds <-function(DATA, ntrees, ndsize, ntreestune, parvec, cvreps, cvfolds, tol){
  TRAIN <- data.frame(DATA[[1]]) #pull out test and training set for a given rep
  TEST <- data.frame(DATA[[2]]) #pull out test and training set for a given rep
  OutlierInd <- DATA[[3]] #pull out test and training set for a given rep
  nPreds <- 16 #consider 15 different types of predictions
  Preds <- array(NA, dim=c(nrow(TEST), nPreds))  #Matrix to store predictions
  niter <- rep(NA, 7) #Create length 7 vector to store number of iterations for each type of RFLowess prediction
  TRAINY <-TRAIN[,ncol(TRAIN)]
  RF <- randomForestSRC::rfsrc(Y~., data=TRAIN, nodesize=ndsize, forest.wt="oob", membership = T)
  RFpredInfo <- predict(RF, newdata=TEST, forest.wt=TRUE, membership = T)
  RFpred <- RFpredInfo$predicted
  TrainNodes <- RF$membership
  TestNodes <- RFpredInfo$membership
  Inbag <- RF$inbag
  OOBWeights <- RF$forest.wt
  PredWeights <- RFpredInfo$forest.wt
  QRF <- quantregForest::quantregForest(x=TRAIN[,-ncol(TRAIN)], y=TRAIN[,ncol(TRAIN)], ntree=ntrees, nodesize=ndsize, keep.forest=TRUE, keep.inbag=TRUE)
  Preds[,1] <- mean(TRAIN[,ncol(TRAIN)])
  Preds[,2] <- RFpred
  Preds[,3] <- predict(QRF, newdata=TEST[,-ncol(TEST)], what=0.5)
  #Other methods
  Preds[,4] <- LiPred(OOBWeights, PredWeights, TRAINY, method="Tukey", delta=0.8,   tol=10^-6, maxiter=100)[[1]]
  Preds[,5] <- LiPred(OOBWeights, PredWeights, TRAINY, method="Huber", delta=0.005,   tol=10^-6, maxiter=100)[[1]]
  MedPreds <- Node_Tree_Agg(TrainNodes, TestNodes, Inbag, TRAINY)
  Preds[,6:8] <- as.matrix(data.frame(MedPreds[[1]],MedPreds[[2]],MedPreds[[3]]))
  Res <- TuneMultifoldCV(TRAIN, TEST, OOBWeights, PredWeights, OutlierInd, ntrees, ntreestune, ndsize, parvec, cvreps, cvfolds, tol=10^-6)
  CVERR <- Res[[1]]
  #using default tuning par of 6
  Lpred <- LOWESSPred(OOBWeights, PredWeights, TRAINY, tol=10^-6, alpha=6, method="Tukey")
  Preds[,9] <- Lpred[[1]]
  niter[1] <- Lpred[[3]]
  #using unweighted CV with MSE
  Lpred <- LOWESSPred(OOBWeights, PredWeights, TRAINY, tol=10^-6, alpha=parvec[which.min(CVERR[1,1,])], method="Tukey")
  Preds[,10] <- Lpred[[1]]
  niter[2] <- Lpred[[3]]
  #using unweighted CV with MAPE
  Lpred <- LOWESSPred(OOBWeights, PredWeights, TRAINY, tol=10^-6, alpha=parvec[which.min(CVERR[2,1,])], method="Tukey")
  Preds[,11] <- Lpred[[1]]
  niter[3] <- Lpred[[3]]
  #using weighted CV by RF resid with MSE
  Lpred <- LOWESSPred(OOBWeights, PredWeights, TRAINY, tol=10^-6, alpha=parvec[which.min(CVERR[3,1,])], method="Tukey")
  Preds[,12] <- Lpred[[1]]
  niter[4] <- Lpred[[3]]
  #using weighted CV by RF resid with MAPE
  Lpred <- LOWESSPred(OOBWeights, PredWeights, TRAINY, tol=10^-6, alpha=parvec[which.min(CVERR[4,1,])], method="Tukey")
  Preds[,13] <- Lpred[[1]]
  niter[5] <- Lpred[[3]]
  #using weighted CV by LOWESS resid with MSE
  Lpred <- LOWESSPred(OOBWeights, PredWeights, TRAINY, tol=10^-6, alpha=parvec[which.min(CVERR[5,1,])], method="Tukey")
  Preds[,14] <- Lpred[[1]]
  niter[6] <- Lpred[[3]]
  #using weighted CV by LOWESS resid with MAPE
  Lpred <- LOWESSPred(OOBWeights, PredWeights, TRAINY, tol=10^-6, alpha=parvec[which.min(CVERR[6,1,])], method="Tukey")
  Preds[,15] <- Lpred[[1]]
  niter[7] <- Lpred[[3]]
  Preds[,16] <- TEST$Y
  #return(list(DATA, Preds, niter, Res))
  return(list(Preds, niter, Res))
  }


#' Run Simulation
#'
#' Function to generate data and run full simulaton
#'
#' @param Sim Which simulation? Either "RL" for Roy Larocque (2012), or "LM" for Li, Martin (2017)
#' @param ntrain number of training cases
#' @param ntest, number of test cases
#' @param p proportion of outliers
#' @param m signal to noise parameter when Sim=="RL"
#' @param contamination Use either variance ("Var") or mean ("Mean") contamination. Only relevant for Sim =="RL".
#' @param Vartype use identity ("Id") or Toeplitz ("Toeplitz") correlation matrix. Only relevant for Sim =="LM"
#' @param DGP If Sim == "RL", which data generating process should be used? either 1 for tree-like, or 2 for non-tree
#' @param DATA  object from generate_RLdata or generate_LMdata
#' @param ntrees  number of trees
#' @param ndsize nodesize
#' @param ntreestune  number of trees to use for tuning alpha
#' @param parvec vector of candidate values for tuning parameter alpha
#' @param cvreps number of repetitions to perform in cross validation
#' @param cvfolds number of folds to perform in cross validation
#' @param  tol maximal change in interation for LOWESSRF weights in cross validation
#' @return returns a list of 4 items
#'        1. Datasets (TRAIN, TEST, and Outlier Indicator)
#'        2. Matrix of 16 columns giving different predictions. Last column is true Y.
#'        3. Number of iterations
#'        4. Output from TuneMultifoldCV (a list of 8 items itself)
#' @export


RunSimulation <- function(Sim = "RL", ntrain, ntest, p, m, contamination = "Var", Vartype="Id", DGP=2, ntrees, ndsize, ntreestune, parvec, cvreps, cvfolds, tol ){
  if (!Sim%in%c("RL", "LM")) {stop('type must be either "RL" or "LM".')}
  if (Sim == "RL"){
    DATA <- generate_RLdata(ntrain, ntest, p, m, type=contamination, DGP)
  }
  else{
    DATA <- generate_LMdata(ntrain, ntest, p, Vartype)
  }
  Preds <- Make_All_Preds(DATA=DATA, ntrees, ndsize, ntreestune, parvec, cvreps, cvfolds, tol)
  return(Preds)
}



#' Sim across m and p
#'
#' Function to repeat simulation over vectors for m and p
#'
#' @param Sim Which simulation? Either "RL" for Roy Larocque (2012), or "LM" for Li, Martin (2017)
#' @param ntrain number of training cases
#' @param ntest, number of test cases
#' @param pvec proportion of outliers
#' @param mvec value of m to use if Sim == "RL"
#' @param contamination Use either variance ("Var") or mean ("Mean") contamination. Only relevant for Sim =="RL".
#' @param Vartype use identity ("Id") or Toeplitz ("Toeplitz") correlation matrix. Only relevant for Sim =="LM"
#' @param DGP If Sim == "RL", which data generating process should be used? either 1 for tree-like, or 2 for non-tree
#' @param DATA  object from generate_RLdata or generate_LMdata
#' @param ntrees  number of trees
#' @param ndsize nodesize
#' @param ntreestune  number of trees to use for tuning alpha
#' @param parvec vector of candidate values for tuning parameter alpha
#' @param cvreps number of repetitions to perform in cross validation
#' @param cvfolds number of folds to perform in cross validation
#' @param  tol maximal change in interation for LOWESSRF weights in cross validation
#' @return returns a list of 4 items
#'        1. Datasets (TRAIN, TEST, and Outlier Indicator)
#'        2. Matrix of 16 columns giving different predictions. Last column is true Y.
#'        3. Number of iterations
#'        4. Output from TuneMultifoldCV (a list of 8 items itself)
#' @export

ApplyAcross_m_and_p <- function(Sim, ntrain, ntest, pvec, mvec, contamination="Var", Vartype="Id", DGP, ntrees, ndsize, ntreestune, parvec, cvreps, cvfolds, tol){
  Res <- sapply(mvec, simplify="array", function(m) sapply(pvec, simplify="array", function(p) RunSimulation(Sim=Sim, ntrain=ntrain, ntest=ntest,m=m, p=p, contamination=contamination, Vartype=Vartype, DGP=DGP, ntrees=ntrees, ndsize=ndsize, ntreestune=ntreestune, parvec=parvec, cvreps=cvreps, cvfolds=cvfolds, tol=tol)))
  return(Res)
}

#' Sim across m and p
#'
#' Function to repeat simulation over vector of values for p
#'
#' @param Sim Which simulation? Either "RL" for Roy Larocque (2012), or "LM" for Li, Martin (2017)
#' @param ntrain number of training cases
#' @param ntest, number of test cases
#' @param pvec proportion of outliers
#' @param m value of m to use if Sim == "RL"
#' @param contamination Use either variance ("Var") or mean ("Mean") contamination. Only relevant for Sim =="RL".
#' @param Vartype use identity ("Id") or Toeplitz ("Toeplitz") correlation matrix. Only relevant for Sim =="LM"
#' @param DGP If Sim == "RL", which data generating process should be used? either 1 for tree-like, or 2 for non-tree
#' @param DATA  object from generate_RLdata or generate_LMdata
#' @param ntrees  number of trees
#' @param ndsize nodesize
#' @param ntreestune  number of trees to use for tuning alpha
#' @param parvec vector of candidate values for tuning parameter alpha
#' @param cvreps number of repetitions to perform in cross validation
#' @param cvfolds number of folds to perform in cross validation
#' @param  tol maximal change in interation for LOWESSRF weights in cross validation
#' @return returns a list of 4 items
#'        1. Datasets (TRAIN, TEST, and Outlier Indicator)
#'        2. Matrix of 16 columns giving different predictions. Last column is true Y.
#'        3. Number of iterations
#'        4. Output from TuneMultifoldCV (a list of 8 items itself)
#' @export

ApplyAcross_p <- function(Sim, ntrain, ntest, pvec, m, contamination="Var", Vartype="Id", DGP, ntrees, ndsize, ntreestune, parvec, cvreps, cvfolds, tol){
  Res <- sapply(pvec, simplify="array", function(p) RunSimulation(Sim=Sim, ntrain=ntrain, ntest=ntest,m=m, p=p, contamination=contamination,
                                                                  Vartype=Vartype, DGP=DGP, ntrees=ntrees, ndsize=ndsize, ntreestune=ntreestune, parvec=parvec, cvreps=cvreps, cvfolds=cvfolds, tol=tol))
  return(Res)
}
