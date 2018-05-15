#' Generate Roy, Larocque data
#'
#' Function to generate data from simulations used in Roy & Larocque (2012)

#' @param ntrain  number of training cases to generate
#' @param ntest  number of test cases to generate
#' @param p  probability of coming from contaminating distribution
#' @param m  multiplier determining signal-to-noise ratio
#' @param type  should contamination be done to the variance or the mean
#' @param DGP  Which data generating process from Roy & Larocque (2012) should be used? Choices are 1 or 2.
#' @return returns a list containing the training dataset, the test dataset, and an indicator of whether training cases came from contaminating distribution.
#' @export

generate_RLdata <- function(ntrain, ntest, p,m, type="Var", DGP=1){
  if (!type%in%c("Mean", "Var")) {stop('type must be either "Var" or "Mean".')}
  if (!DGP%in%c(1,2)) {stop('DGP must either 1 or 2.')}
  X <- MASS::mvrnorm(n = ntrain, mu=c(rep(0,6)), Sigma=diag(6))
  ds <- rbinom(ntrain,1,p)  #determine which training cases are outliers
  e <- c(rep(NA, ntrain))   #initialize error vector
  e[ds==0] <- rnorm(length(e[ds==0]),0,1) #non-contaminated errors
  if (type=="Var"){
    e[ds==1]=rnorm(length(e[ds==1]),0,5) #contaminated errors for variance cont.
  } else {
    e[ds==1]=rnorm(length(e[ds==1]),5,1) #contaminated errors for mean cont.
  }
  if (DGP == 1){  #first, tree-like data generating mechanism
  Y <- m*((1*X[,1]<=0 & X[,2]<=0)+2*(X[,1]<=0 & X[,2]>0 & X[,4]<=0)+3*(X[,1]<=0 & X[,2]>0 & X[,4]>0 & X[,6]<=0)+4*(X[,1]<=0 & X[,2]>0 & X[,4]>0 & X[,6]>0)+5*(X[,1]>0 & X[,3]<=0)+6*(X[,1]>0 & X[,3]>0 & X[5]<=0) + 7*(X[,1]>0 & X[,3]>0 & X[5]>0))+e
  } else {  #second, nontree data generating mechanism
    Y=m*(X[,1]+.707*(X[,2])^2+2*(X[,3]>0)+.873*log(abs(X[,1]))*X[,3]+0.894*X[,2]*X[,4]+2*(X[,5]>0)+0.464*exp(X[,6]))+e
  }
  TRAIN <- data.frame(X, Y)         #Training Set
  Xt <- MASS::mvrnorm(n = ntest, mu=c(rep(0,6)), Sigma=diag(6))   #Generate predictors for test data
  et <- rnorm(ntest,0,1)                 #Errors for test cases-no contamination
  if (DGP == 1){  #first, tree-like data generating mechanism
    Y <- m*((1*Xt[,1]<=0 & Xt[,2]<=0)+2*(Xt[,1]<=0 & Xt[,2]>0 & Xt[,4]<=0)+3*(Xt[,1]<=0 & Xt[,2]>0 & Xt[,4]>0 & Xt[,6]<=0)+4*(Xt[,1]<=0 & Xt[,2]>0 & Xt[,4]>0 & Xt[,6]>0)+5*(Xt[,1]>0 & Xt[,3]<=0)+6*(Xt[,1]>0 & Xt[,3]>0 & Xt[5]<=0) + 7*(Xt[,1]>0 & Xt[,3]>0 & Xt[5]>0))+et
  } else {  #second, nontree data generating mechanism
    Y <- m*(Xt[,1]+.707*(Xt[,2])^2+2*(Xt[,3]>0)+.873*log(abs(Xt[,1]))*Xt[,3]+0.894*Xt[,2]*Xt[,4]+2*(Xt[,5]>0)+0.464*exp(Xt[,6]))+et
  }
  TEST <- data.frame(Xt, Y)
  return(list(TRAIN, TEST, ds))
}


#'Generate Li, Martin data
#'
#' Function to generate data from simulations used in Li & Martin (2017)

#' @param ntrain  number of training cases to generate
#' @param ntest  number of test cases to generate
#' @param p  probability of coming from contaminating distribution
#' @param Vartype  covariance matrix can be either identity matrix ("Id") or Toeplitz(0.7) matrix ("Toeplitz)
#' @return returns a list containing the training dataset, the test dataset, and an indicator of whether training cases came from contaminating distribution.
#' @export

generate_LMdata <- function(ntrain, ntest, p, Vartype="Id"){
  if (!Vartype%in%c("Id", "Toeplitz")) {stop('type must be either "Id" or "Toeplitz".')}
    if (Vartype=="Id"){
    X <- MASS::mvrnorm(n = ntrain, mu=c(rep(0,10)), Sigma=diag(10))
  } else {
    X <- MASS::mvrnorm(n = ntrain, mu=c(rep(0,10)), Sigma=toeplitz(0.7^seq.int(0, 10-1)))#Toeplitz(.7) matrix defined in Li, Martin paper
  }
  ds <- rbinom(ntrain,1,p)                  #which cases are contaminated?
  e <- rnorm(ntrain,0,1)              #non-contaminated errors
  e[ds==1] <- e[ds==1]+15*rt(sum(ds==1),2)   #non-contaminated errors
  Y <- rowSums(X^2)+e
  TRAIN <- data.frame(X, Y)                 #training set
  if (Vartype=="Id"){
    Xt <- MASS::mvrnorm(n = ntest, mu=c(rep(0,10)), Sigma=diag(10))
  } else {
    Xt=MASS::mvrnorm(n = ntrain, mu=c(rep(0,10)), Sigma=toeplitz(0.7^seq.int(0, 10-1)))#Toeplitz(.7) matrix defined in Li, Martin paper
  }
  et <- rnorm(ntest,0,1)                   #Generate test errors no contamination
  Yt <- rowSums(Xt^2)+et
  TEST <- data.frame(Xt, Yt)
  return(list(TRAIN, TEST, ds))
}


#'Partition Real Data
#'
#' Function to partition real datasets into test and training sets and possibly add contamination

#' @param dataset  real dataset to work with
#' @param nfolds number of folds to partition into to assess performance
#' @return returns a list containing the training dataset, the test dataset, and an indicator of whether training cases came from contaminating distribution.
#' @export
#'

Partition_Real_Data <- function(dataset, nfolds){
  #if (!(nrow(dataset)%%nfolds==0)) {stop('Number of folds must evenly divide number of observations')}
      orderedcases <- sample(1:nrow(dataset),replace=F)
      foldsize <- ceiling(nrow(dataset)/nfolds)
      extras <- foldsize*nfolds - length(orderedcases)
      if(extras>0){
        repeatedcases <- orderedcases[1:extras]
        orderedcases <- c(orderedcases, repeatedcases)
      }
      CVTRAINind <- array(NA, dim=c(nfolds, nrow(dataset)-foldsize))
      CVTESTind <- array(NA, dim=c(nfolds, foldsize))
      for(fold in 1:nfolds){
        CVTESTind[fold,] <- orderedcases[((fold-1)*foldsize+1):(fold*foldsize)]
        CVTRAINind[fold,] <- unique(orderedcases[!orderedcases%in%CVTESTind[fold,]])
}
    return(list(CVTRAINind, CVTESTind))
  }

