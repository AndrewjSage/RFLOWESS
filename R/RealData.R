#' Assess_Folds
#'
#' Function to take partitioned data, and make predictions across all folds
#'
#'
#' @param PartitionedData  List of partitioned dataframes for training, test sets, and indicator of outlying observations
#' @param folds fold on which to assess performance
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
Assess_Folds <- function(PartitionedData, fold, ntrees, ndsize, ntreestune, parvec, cvreps, cvfolds, tol ){
DATA <- list(data.frame(PartitionedData[[1]][fold,,]), data.frame(PartitionedData[[2]][fold,,]), PartitionedData[[3]][fold,])
names(DATA[[1]])[ncol(DATA[[1]])] <- "Y"
names(DATA[[2]])[ncol(DATA[[2]])] <- "Y"
Preds <- Make_All_Preds(DATA=DATA, ntrees, ndsize, ntreestune, parvec, cvreps, cvfolds, tol)
return(Preds)
}

#' Assess_Real_Data
#'
#' Function to divide real datasets into test and training sets, possibly add contamination
#' and assess predictions
#'
#'
#' @param dataset  real dataset to work with
#' @param nfolds number of folds to partition into to assess performance
#' @param p percentage of training cases for which to add contamination (using N(0, 5*sd(Y)))
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


Assess_Real_Data <- function(dataset, nfolds, p, ntrees, ndsize, ntreestune, parvec, cvreps, cvfolds, tol ){
    FULLDATA <- Partition_Real_Data(dataset, nfolds=nfolds, p=p)
    Preds <- sapply(1:nfolds, simplify="array", Assess_Folds, PartitionedData=FULLDATA, ntrees=ntrees, ndsize=ndsize,
           ntreestune=ntreestune, parvec=parvec, cvreps=cvreps, cvfolds=cvfolds, tol=tol  )
  return(Preds)
}




