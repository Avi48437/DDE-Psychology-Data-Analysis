library(rpart)
library(pbapply)
library(parallel)

root_dir <- "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology2"

ipip_dir <- file.path(root_dir,"2. IPIP 100")
data_dir <- file.path(ipip_dir,"Analysis","Predictive_Data")
results_dir <- file.path(ipip_dir,"Analysis","Results")
scripts_dir <- file.path(ipip_dir,"Analysis","Scripts")

clf_file <- file.path(scripts_dir,"Prediction_Classifiers.R")
source(clf_file)

X <- read.csv(file.path(data_dir,"X98.csv"))
A <- read.csv(file.path(data_dir,"A98.csv"))
Y <- read.csv(file.path(data_dir,"Y.csv"))
D <- read.csv(file.path(data_dir,"D_IPIP98.csv"))

D <- D[,c("age","gender","race")]
D$gender <- factor(D$gender)
D$race <- factor(D$race)

XD_all <- cbind(D,X)
AD_all <- cbind(D,A)

prepare_outcome <- function(outcome){
  y <- Y[[outcome]]
  
  if(outcome=="educ"){
    keep <- !is.na(y)
    return(list(y=as.numeric(y[keep]),rows=which(keep),stratified=FALSE,type="ordinal"))
  }
  
  if(outcome=="marstat"){
    keep <- !is.na(y); y0 <- y[keep]; yg <- rep(NA_integer_,length(y0))
    yg[y0==1] <- 1; yg[y0 %in% c(2,3)] <- 2; yg[y0==4] <- 3; yg[y0 %in% c(5,6)] <- 4
    keep2 <- !is.na(yg)
    return(list(y=yg[keep2],rows=which(keep)[keep2],stratified=TRUE,type="nominal"))
  }
  
  if(outcome=="faminc_new"){
    keep <- !is.na(y) & y!=97; y0 <- y[keep]; yg <- rep(NA_integer_,length(y0))
    yg[y0 %in% 1:4] <- 1; yg[y0 %in% 5:8] <- 2; yg[y0 %in% 9:11] <- 3; yg[y0 %in% 12:16] <- 4
    keep2 <- !is.na(yg)
    return(list(y=yg[keep2],rows=which(keep)[keep2],stratified=TRUE,type="ordinal"))
  }
  
  if(outcome=="votereg"){
    keep <- !is.na(y) & y!=3
    return(list(y=as.numeric(y[keep]),rows=which(keep),stratified=TRUE,type="binary"))
  }
  
  if(outcome=="pew_prayer"){
    keep <- !is.na(y)
    return(list(y=as.numeric(y[keep]),rows=which(keep),stratified=TRUE,type="nominal"))
  }
  
  stop("Unknown outcome.")
}

make_split <- function(y,stratified){
  if(!stratified){
    tr <- sample(seq_along(y),floor(0.8*length(y)))
  } else {
    tr <- unlist(lapply(split(seq_along(y),y),function(ii) sample(ii,floor(0.8*length(ii)))),use.names=FALSE)
  }
  tr <- sort(tr)
  list(train=tr,test=setdiff(seq_along(y),tr))
}

outcomes <- c("educ","marstat","faminc_new","votereg","pew_prayer")
n_rep <- 100
n_cores <- 10

cl <- makeCluster(n_cores)
clusterExport(cl,"clf_file",envir=environment())
clusterEvalQ(cl,{
  library(rpart)
  source(clf_file)
  NULL
})

IPIP98_CART_Demographics <- list()

for(outcome in outcomes){
  
  cat("\nRunning:",outcome,"\n")
  
  dat <- prepare_outcome(outcome)
  y <- dat$y
  rows <- dat$rows
  
  X0 <- X[rows,,drop=FALSE]
  A0 <- A[rows,,drop=FALSE]
  XD <- XD_all[rows,,drop=FALSE]
  AD <- AD_all[rows,,drop=FALSE]
  
  splits <- lapply(1:n_rep,function(r) make_split(y,dat$stratified))
  
  clusterExport(cl,c("X0","A0","XD","AD","y","rows","splits"),envir=environment())
  
  reps <- pblapply(1:n_rep,function(r){
    
    tr <- splits[[r]]$train
    te <- splits[[r]]$test
    
    res_X <- clf_cart(X0[tr,,drop=FALSE],X0[te,,drop=FALSE],y[tr],y[te])
    res_A <- clf_cart(A0[tr,,drop=FALSE],A0[te,,drop=FALSE],y[tr],y[te])
    res_XD <- clf_cart(XD[tr,,drop=FALSE],XD[te,,drop=FALSE],y[tr],y[te])
    res_AD <- clf_cart(AD[tr,,drop=FALSE],AD[te,,drop=FALSE],y[tr],y[te])
    
    list(
      repetition=r,
      train_idx=tr,
      test_idx=te,
      train_rows=rows[tr],
      test_rows=rows[te],
      y_train=y[tr],
      y_test=y[te],
      Big5=res_X,
      DDE_A=res_A,
      Demo_Big5=res_XD,
      Demo_DDE_A=res_AD
    )
    
  },cl=cl)
  
  IPIP98_CART_Demographics[[outcome]] <- list(
    type=dat$type,
    stratified=dat$stratified,
    n=length(y),
    repetitions=reps
  )
  
  saveRDS(IPIP98_CART_Demographics,file.path(data_dir,"IPIP98_CART_Demographics_100.rds"))
  cat("Completed:",outcome,"\n")
}

stopCluster(cl)

saveRDS(IPIP98_CART_Demographics,file.path(results_dir,"IPIP98_CART_Demographics_100.rds"))
cat("\nSaved:",file.path(results_dir,"IPIP98_CART_Demographics_100.rds"),"\n")