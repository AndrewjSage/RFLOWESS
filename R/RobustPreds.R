#' Tukey bisquare function
#'
#' Bisquare function applied to a vector of values, element-by-element

#' @param t & a vector of values on which to apply the Tukey Bisquare function
#' @return returns a vector the same length as t, containing values after applying Tukey bisquare function
#' @export

Tukey=function(t){
  return(max(c(1-(t)^2,0))^2)
}

#' Pseudo-Huber function
#'
#' Pseudo-Huber function, used by Li & Martin, applied to a vector of values, element-by-element

#' @param t & a vector of values on which to apply the pseudo Huber function
#' @return returns a vector the same length as t, containing values after applying pesudo Huber function
#' @export

Huber=function(t){
  return(1/sqrt(1+(t)^2))
}


#' Function to compute RFLOWESS weight multipliers for training cases
#'
#' Iteratively adjusts training case weights using residuals for out-of-bag cases.

#' @param OOBWeights   matrix of OOB training case weights
#' @param TRAINY  vector of training case response values
#' @param alpha  tuning parameter for robust weighting. Small alpha corresponds to more aggressive robust adjustment.
#' @param method  should Tukey or Huber weighting function be used?
#' @param tol  maximal tolerance for change between iterations
#' @param maxiter  maximum number of iterations
#' @return returns a list containing
#' 1- vector of multipliers for training case weights
#' 2- adjusted OOB training case predictions
#' 3- number of iterations
#' @export

Compute_deltas=function(OOBWeights, TRAINY, alpha=6, method="Tukey", tol=10^-6, maxiter=100){
  if (!method%in%c("Tukey", "Huber")) {stop('type must be either "Tukey" or "Huber".')}
  #Before going into the loop, we're just setting up the matrices and getting the original OOB
  #prediction and residuals
  d <- rep(1, length(TRAINY))
  D0 <- matrix(rep(d, length(TRAINY)), nrow=length(TRAINY), ncol=length(TRAINY), byrow=T)
  AdjWts <- D0 * OOBWeights + 10^-9  #last term adds small amount to ensure weights don't all get set to 0 when alpha is very small
  AdjWts <- AdjWts / rowSums(AdjWts)
  OOBPred <- as.matrix(AdjWts) %*% TRAINY
  resid <- TRAINY- OOBPred
  resid[is.na(resid)] <- mean(abs(resid), na.rm=T)   #if case does not come up OOB set residual to mean of residuals
  niter <- 1
  Change <- 1
  while(Change > tol & niter < maxiter){  #main loop for reweighting
    d0 <- d
    s <- median(abs(resid))
    t <- resid/(alpha*s)
    if (method=="Tukey") {
      d <- t(apply(data.frame(t), 1, Tukey)) #Tukey weights
    }else {
      d <- t(apply(data.frame(t), 1, Huber)) #Huber weights
    }
    D0 <- matrix(rep(d, length(TRAINY)), nrow=length(TRAINY), ncol=length(TRAINY), byrow=T)
    AdjWts <- D0*OOBWeights + 10^-9 # add in very small buffer so that if all are zero we don't divide by 0
    AdjWts <- AdjWts / rowSums(AdjWts)
    OOBPred0 <- OOBPred
    OOBPred <- as.matrix(AdjWts) %*% TRAINY
    resid <- TRAINY - OOBPred
    resid[is.na(resid)] <- mean(abs(resid), na.rm=T)   #if case does not come up OOB
    Change <- mean((OOBPred-OOBPred0)^2)
    niter <- niter + 1
  }
  return(list(d0, OOBPred, niter))
}


#' RFLOWESS Prediction
#'
#' Function to obtain LOWESSRF prediction using case weights for OOB predictions and test predictions
#' by calling Compute_deltas

#' @param OOBWeights  matrix of OOB training case weights
#' @param PredWeights  ntest by ntrain matrix of prediction case weights
#' @param TRAINY  vector of training case response values
#' @param method  should Tukey or Huber weighting function be used?
#' @param tol  maximal tolerance for change between iterations
#' @param maxiter  maximum number of iterations
#' @return returns a list containing
#' 1- vector of LOWESSRF predictions
#' 2- vector of adjusted case weights
#' 3- number of iterations
#' @export

LOWESSPred <- function(OOBWeights, PredWeights, TRAINY, alpha=6, method="Tukey", tol=10^-6){
  if (!method%in%c("Tukey", "Huber")) {stop('type must be either "Tukey" or "Huber".')}
  #adjust weights for training cases
  Res <- Compute_deltas(OOBWeights=OOBWeights, TRAINY=TRAINY, alpha=alpha, method=method, tol=tol)
  d <- Res[[1]]
  niter <- Res[[3]]
  #Apply adjustment to test cases
  D <- matrix(rep(d, nrow(PredWeights)), nrow=nrow(PredWeights), ncol=length(TRAINY), byrow=T)
  AdjWts <- D*PredWeights
  AdjWts <- AdjWts/rowSums(AdjWts)
  Pred <- as.matrix(AdjWts)%*%TRAINY
  return(list(Pred, AdjWts, niter))
}


#' Li, Martin Predictions
#'
#' Predictions obtained by reweighting training cases using Li & Martin's (2017) approach

#' @param OOBWeights matrix of OOB training case weights
#' @param PredWeights ntest by ntrain matrix of prediction case weights
#' @param TRAINY  vector of training case response values
#' @param method  should Tukey or Huber weighting function be used?
#' @param delta  value of tuning parameter
#' @param tol  maximal tolerance for change between iterations
#' @param maxiter  maximum number of iterations
#' @return returns a list containing
#' 1- vector of predictions
#' 2- vector of adjusted case weights
#' 3- number of iterations
#' @export

LiPred <- function(OOBWeights, PredWeights, TRAINY, method="Huber", delta=0.005,   tol=10^-6, maxiter=100){
  if (!method%in%c("Tukey", "Huber")) {stop('type must be either "Tukey" or "Huber".')}
  Weights <- PredWeights
  Pred <- PredWeights%*%scale(TRAINY) #standardize training responses
  Pred <- PredWeights%*%scale(TRAINY) #standardize training responses
  Change <- 1
  niter <- 1
  while (Change > tol & niter < maxiter){
    diffmat <- outer(as.vector(scale(TRAINY)), as.vector(Pred), "-") #matrix of differences between Y's for training cases and predicted y's for test cases
    diffmat <- diffmat/delta
    if (method == "Tukey") {
      RobWts <- t(matrix(sapply(diffmat, Tukey), nrow=nrow(diffmat), ncol=ncol(diffmat)))
    }else {
      RobWts <- t(matrix(sapply(diffmat, Huber), nrow=nrow(diffmat), ncol=ncol(diffmat)))
    }
    AdjWts <- RobWts*Weights
    AdjWts <- AdjWts/rowSums(AdjWts)
    Pred0 <- Pred
    Pred <- as.matrix(AdjWts)%*%scale(TRAINY)
    Pred[is.na(Pred)] <- Pred0[is.na(Pred)] #In case of Tukey, if all weights are 0, keep pred same
    Change <- mean((scale(Pred)-scale(Pred0))^2)
    niter <- niter+1
  }
  Pred <- Pred*sd(TRAINY) + mean(TRAINY) #Get back to original scale
  return(list(Pred, AdjWts, niter))
}

#' Median Aggregation Predictions.
#'
#' Method to get predictions from median-type aggregation methods used in paper by Roy & Larocque (2012)

#' @param TrainNodes ntrain by ntree matrix containing terminal node numbers of each training case in each tree
#' @param TestNodes ntest by ntree matrix containing terminal node numbers of each training case in each tree
#' @param Inbag ntrain by ntree matrix telling number of times training case came up in bootstrap sample used to grow each tree
#' @param TRAINY vector of training case response values
#' @return returns a list containing predictions using
#' 1- mean within terminal nodes and median across trees
#' 2- median within terminal nodes and median across trees
#' 3- median within terminal nodes and mean across trees
#' @export

Node_Tree_Agg=function(TrainNodes, TestNodes, Inbag, TRAINY){
  rownames(TrainNodes) <- 1:nrow(TrainNodes) #rename so predictions are in right order
  rownames(TestNodes) <- 1:nrow(TestNodes) #rename so predictions are in right order
  M <- reshape2::melt(TrainNodes)
  ntree <- max(M[,2])
  M$YVal <- rep(TRAINY, ntree)
  names(M) <- c("TrCase", "Tree", "Node", "Yval")
  ntree <- max(M$Tree)
  #Create a new variable that incorporates both trees and nodes. Each node gets its unique value
  #instead of being nested in tree
  M$Inbag <- reshape2::melt(Inbag)[,3]
  M <- M[M$Inbag>0,]
  M <- splitstackshape::expandRows(M, "Inbag")
  NodePreds <- data.frame(M %>%
                         group_by(Tree, Node) %>%
                         summarise (Meany = mean(Yval),Medy = median(Yval) ) )
  NodePreds$NodeKey <- NodePreds$Node+NodePreds$Tree/ntree
  TestM <- reshape2::melt(TestNodes)
  names(TestM) <- c("TeCase", "Tree", "Node")
  TestM$NodeKey <- TestM$Node+TestM$Tree/ntree
  MatchNodes <- merge(NodePreds, TestM, key="NodeKey") #Contains mean and median response for the terminal node that each test case lands in for each tree
  MedPredMat <- reshape2::acast(MatchNodes, TeCase~Tree, value.var="Medy")
  MeanPredMat <- reshape2::acast(MatchNodes, TeCase~Tree, value.var="Meany")
  Mean_Med_Predictions <- apply(MeanPredMat,1,median)       #Took mean within each node and median over trees
  Med_Med_Predictions <- apply(MedPredMat,1,median)        #Took median within each node and median over trees
  Med_Mean_Predictions <- apply(MedPredMat,1,mean)       #Took median within each node and mean over trees
  return(list(Mean_Med_Predictions, Med_Med_Predictions,Med_Mean_Predictions))
}

