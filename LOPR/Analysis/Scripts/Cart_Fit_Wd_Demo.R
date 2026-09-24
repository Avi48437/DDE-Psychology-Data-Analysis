library(rpart)
library(pbapply)
library(parallel)

root_dir <- "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology"

data_dir <- file.path(root_dir,"LOPR","Analysis","Predictive_Data")
results_dir <- file.path(root_dir,"LOPR","Analysis","Results")
scripts_dir <- file.path(root_dir,"LOPR","Analysis","Scripts")

clf_file <- file.path(scripts_dir,"Prediction_Classifiers.R")

source(clf_file)

X <- read.csv(file.path(data_dir,"X_rc.csv"))
A <- read.csv(file.path(data_dir,"A.csv"))
Y <- read.csv(file.path(data_dir,"Y.csv"))
D <- read.csv(file.path(data_dir,"D_LOOPR.csv"))

R <- as.matrix(D[paste0("Race_",1:6)])

D$Race <- apply(R,1,function(z){
  ii <- which(z==1)
  if(length(ii)==1) return(ii)
  if(length(ii)>1) return(7)
  NA
})

keep_race <- !is.na(D$Race)

X <- X[keep_race,,drop=FALSE]
A <- A[keep_race,,drop=FALSE]
Y <- Y[keep_race,,drop=FALSE]
D <- D[keep_race,,drop=FALSE]

D <- D[,c("Age","Sex","Race")]
D$Sex <- factor(D$Sex)
D$Race <- factor(D$Race,levels=1:7)

XD_all <- cbind(D,X)
AD_all <- cbind(D,A)

prepare_outcome <- function(outcome){
  
  y <- Y[[outcome]]
  yg <- rep(NA_integer_,length(y))
  
  if(outcome=="educ"){
    yg[y %in% 1:3] <- 1
    yg[y==4] <- 2
    yg[y==5] <- 3
    yg[y==6] <- 4
    yg[y==7] <- 5
    yg[y %in% 8:10] <- 6
    stratified <- FALSE
    type <- "ordinal"
  }
  
  if(outcome=="marstat"){
    yg[y==1] <- 1
    yg[y %in% c(3,4)] <- 2
    yg[y==2] <- 3
    yg[y==5] <- 4
    stratified <- TRUE
    type <- "nominal"
  }
  
  if(outcome=="faminc_new"){
    yg[y %in% 1:4] <- 1
    yg[y %in% 5:8] <- 2
    yg[y %in% 9:11] <- 3
    yg[y==12] <- 4
    stratified <- TRUE
    type <- "ordinal"
  }
  
  if(outcome=="pew_prayer"){
    yg <- as.numeric(y)
    stratified <- TRUE
    type <- "ordinal"
  }
  
  keep <- !is.na(yg)
  
  list(
    y=yg[keep],
    rows=which(keep),
    stratified=stratified,
    type=type
  )
}

make_split <- function(y,stratified){
  
  if(!stratified){
    tr <- sample(seq_along(y),floor(0.8*length(y)))
  } else {
    tr <- unlist(
      lapply(
        split(seq_along(y),y),
        function(ii) sample(ii,floor(0.8*length(ii)))
      ),
      use.names=FALSE
    )
  }
  
  tr <- sort(tr)
  
  list(
    train=tr,
    test=setdiff(seq_along(y),tr)
  )
}

outcomes <- c("educ","marstat","faminc_new","pew_prayer")

n_rep <- 100
n_cores <- 10

cl <- makeCluster(n_cores)

clusterExport(cl,"clf_file",envir=environment())

clusterEvalQ(cl,{
  library(rpart)
  source(clf_file)
  NULL
})

LOOPR_CART_Demographics <- list()

for(outcome in outcomes){
  
  cat("\nRunning:",outcome,"\n")
  
  dat <- prepare_outcome(outcome)
  
  y <- dat$y
  rows <- dat$rows
  
  X0 <- X[rows,,drop=FALSE]
  A0 <- A[rows,,drop=FALSE]
  XD <- XD_all[rows,,drop=FALSE]
  AD <- AD_all[rows,,drop=FALSE]
  
  splits <- lapply(
    1:n_rep,
    function(r) make_split(y,dat$stratified)
  )
  
  clusterExport(
    cl,
    c("X0","A0","XD","AD","y","rows","splits"),
    envir=environment()
  )
  
  reps <- pblapply(
    1:n_rep,
    function(r){
      
      tr <- splits[[r]]$train
      te <- splits[[r]]$test
      
      res_X <- clf_cart(
        X0[tr,,drop=FALSE],
        X0[te,,drop=FALSE],
        y[tr],
        y[te]
      )
      
      res_A <- clf_cart(
        A0[tr,,drop=FALSE],
        A0[te,,drop=FALSE],
        y[tr],
        y[te]
      )
      
      res_XD <- clf_cart(
        XD[tr,,drop=FALSE],
        XD[te,,drop=FALSE],
        y[tr],
        y[te]
      )
      
      res_AD <- clf_cart(
        AD[tr,,drop=FALSE],
        AD[te,,drop=FALSE],
        y[tr],
        y[te]
      )
      
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
      
    },
    cl=cl
  )
  
  LOOPR_CART_Demographics[[outcome]] <- list(
    type=dat$type,
    stratified=dat$stratified,
    n=length(y),
    repetitions=reps
  )
  
  saveRDS(
    LOOPR_CART_Demographics,
    file.path(results_dir,"LOOPR_CART_Demographics_100.rds")
  )
  
  cat("Completed:",outcome,"\n")
}

stopCluster(cl)

saveRDS(
  LOOPR_CART_Demographics,
  file.path(results_dir,"LOOPR_CART_Demographics_100.rds")
)

cat(
  "\nSaved:",
  file.path(results_dir,"LOOPR_CART_Demographics_100.rds"),
  "\n"
)