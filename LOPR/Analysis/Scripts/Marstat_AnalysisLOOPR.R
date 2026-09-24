## ============================================================
## LOOPR Marital Status: 100-run predictive analysis
## ============================================================

data_dir <- "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology/LOPR/Analysis/Predictive_Data"
results_dir <- "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology/LOPR/Analysis/Results"

source(file.path(scripts_dir,"Prediction_Classifiers.R"))

X <- read.csv(file.path(data_dir,"X_rc.csv"))
A <- read.csv(file.path(data_dir,"A.csv"))
Y <- read.csv(file.path(data_dir,"Y.csv"))

## ============================================================
## Marital-status grouping
##
## 1   -> Married
## 2   -> Widowed
## 3,4 -> Separated/Divorced
## 5   -> Never married
## ============================================================

mar <- Y$marstat

y <- rep(NA_character_,length(mar))

y[mar==1] <- "Married"
y[mar==2] <- "Widowed"
y[mar %in% c(3,4)] <- "Separated/Divorced"
y[mar==5] <- "Never married"

keep <- !is.na(y)

X <- X[keep,,drop=FALSE]
A <- A[keep,,drop=FALSE]

y <- factor(
  y[keep],
  levels=c(
    "Married",
    "Separated/Divorced",
    "Widowed",
    "Never married"
  )
)

table(y)

nRep <- 100
train_prop <- 0.8

## ============================================================
## Compact result storage
## ============================================================

compact_result <- function(res){
  if(is.null(res)) return(NULL)
  res$fit <- NULL
  res
}

compact_all_results <- function(x){
  lapply(x,function(z){
    list(
      Big5=compact_result(z$Big5),
      A=compact_result(z$A)
    )
  })
}

## ============================================================
## Storage
## ============================================================

LOOPR_marstat <- vector("list",nRep)

## ============================================================
## 100 repetitions
## ============================================================

for(r in seq_len(nRep)){
  
  cat("\n============================================================\n")
  cat("LOOPR marital-status repetition:",r,"/",nRep,"\n")
  cat("============================================================\n")
  
  ## ----------------------------------------------------------
  ## Stratified 80/20 split
  ## ----------------------------------------------------------
  
  train_idx <- unlist(
    lapply(split(seq_along(y),y),function(idx){
      sample(idx,floor(train_prop*length(idx)))
    }),
    use.names=FALSE
  )
  
  test_idx <- setdiff(seq_along(y),train_idx)
  
  X_train <- X[train_idx,,drop=FALSE]
  X_test <- X[test_idx,,drop=FALSE]
  
  A_train <- A[train_idx,,drop=FALSE]
  A_test <- A[test_idx,,drop=FALSE]
  
  y_train <- y[train_idx]
  y_test <- y[test_idx]
  
  ## ----------------------------------------------------------
  ## 1. Random forest
  ## ----------------------------------------------------------
  
  res_X_rf <- clf_random_forest(X_train,X_test,y_train,y_test)
  res_A_rf <- clf_random_forest(A_train,A_test,y_train,y_test)
  
  ## ----------------------------------------------------------
  ## 2. CART
  ## ----------------------------------------------------------
  
  res_X_cart <- clf_cart(X_train,X_test,y_train,y_test)
  res_A_cart <- clf_cart(A_train,A_test,y_train,y_test)
  
  ## ----------------------------------------------------------
  ## 3. DDE-structured plug-in Bayes
  ## ----------------------------------------------------------
  
  res_A_bayes <- clf_dde_structured_bayes(
    A_train,A_test,y_train,y_test
  )
  
  ## ----------------------------------------------------------
  ## 4. Multinomial logit
  ## ----------------------------------------------------------
  
  res_X_multi <- clf_multinomial_logit(
    X_train,X_test,y_train,y_test,
    formula_mode="all_pairwise"
  )
  
  res_A_multi <- clf_multinomial_logit(
    A_train,A_test,y_train,y_test,
    formula_mode="A1_A2"
  )
  
  ## ----------------------------------------------------------
  ## 5. Penalized multinomial logit
  ## ----------------------------------------------------------
  
  res_X_pmulti <- clf_penalized_multinomial(
    X_train,X_test,y_train,y_test,
    formula_mode="all_pairwise",
    alpha=0.5,
    nfolds=5,
    lambda_choice="lambda.min"
  )
  
  res_A_pmulti <- clf_penalized_multinomial(
    A_train,A_test,y_train,y_test,
    formula_mode="A1_A2",
    alpha=0.5,
    nfolds=5,
    lambda_choice="lambda.min"
  )
  
  ## ----------------------------------------------------------
  ## 6. XGBoost
  ## ----------------------------------------------------------
  
  res_X_xgb <- clf_xgboost(
    X_train,X_test,y_train,y_test,
    nrounds=150,
    max_depth=3,
    eta=0.05
  )
  
  res_A_xgb <- clf_xgboost(
    A_train,A_test,y_train,y_test,
    nrounds=150,
    max_depth=3,
    eta=0.05
  )
  
  ## ----------------------------------------------------------
  ## 7. SVM
  ## ----------------------------------------------------------
  
  res_X_svm <- clf_svm(
    X_train,X_test,y_train,y_test,
    kernel="radial"
  )
  
  res_A_svm <- clf_svm(
    A_train,A_test,y_train,y_test,
    kernel="radial"
  )
  
  ## ----------------------------------------------------------
  ## 8. Hamming-distance kNN
  ## ----------------------------------------------------------
  
  res_A_hamming <- clf_hamming_knn(
    A_train,A_test,y_train,y_test,
    k=15
  )
  
  ## ----------------------------------------------------------
  ## 9. Tree-augmented naive Bayes
  ## ----------------------------------------------------------
  
  A2_name <- grep("^A2_",names(A_train),value=TRUE)
  
  if(length(A2_name)!=1)
    stop("Expected exactly one A2 coordinate.")
  
  res_A_tan <- clf_tan_binary(
    A_train,A_test,y_train,y_test,
    root=A2_name
  )
  
  ## ----------------------------------------------------------
  ## Raw results
  ## ----------------------------------------------------------
  
  all_results <- list(
    "Random forest"=list(Big5=res_X_rf,A=res_A_rf),
    "CART"=list(Big5=res_X_cart,A=res_A_cart),
    "DDE-structured plug-in Bayes"=list(Big5=NULL,A=res_A_bayes),
    "Multinomial logit"=list(Big5=res_X_multi,A=res_A_multi),
    "Penalized multinomial logit"=list(Big5=res_X_pmulti,A=res_A_pmulti),
    "XGBoost"=list(Big5=res_X_xgb,A=res_A_xgb),
    "SVM"=list(Big5=res_X_svm,A=res_A_svm),
    "Hamming kNN"=list(Big5=NULL,A=res_A_hamming),
    "Tree-augmented naive Bayes"=list(Big5=NULL,A=res_A_tan)
  )
  
  ## ----------------------------------------------------------
  ## Store repetition
  ## ----------------------------------------------------------
  
  LOOPR_marstat[[r]] <- list(
    repetition=r,
    train_idx=train_idx,
    test_idx=test_idx,
    y_train=y_train,
    y_test=y_test,
    results=compact_all_results(all_results)
  )
  
  ## ----------------------------------------------------------
  ## Checkpoint every 5 repetitions
  ## ----------------------------------------------------------
  
  if(r %% 5==0){
    
    saveRDS(
      LOOPR_marstat,
      file.path(results_dir,"LOOPR_marstat.rds")
    )
    
    cat("Checkpoint saved through repetition",r,"\n")
  }
  
  ## ----------------------------------------------------------
  ## Cleanup
  ## ----------------------------------------------------------
  
  rm(
    res_X_rf,res_A_rf,
    res_X_cart,res_A_cart,
    res_A_bayes,
    res_X_multi,res_A_multi,
    res_X_pmulti,res_A_pmulti,
    res_X_xgb,res_A_xgb,
    res_X_svm,res_A_svm,
    res_A_hamming,res_A_tan,
    all_results
  )
  
  gc()
}

## ============================================================
## Final save
## ============================================================

saveRDS(
  LOOPR_marstat,
  file.path(results_dir,"LOOPR_marstat.rds")
)

cat("\nCompleted",length(LOOPR_marstat),"LOOPR marital-status repetitions.\n")