## ============================================================
## LOOPR Family Income: 100-run predictive analysis
## ============================================================

data_dir <- "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology/LOPR/Analysis/Predictive_Data"
results_dir <- "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology/LOPR/Analysis/Results"

source(file.path(scripts_dir,"Prediction_Classifiers.R"))

X <- read.csv(file.path(data_dir,"X_rc.csv"))
A <- read.csv(file.path(data_dir,"A.csv"))
Y <- read.csv(file.path(data_dir,"Y.csv"))

## ============================================================
## Family-income grouping
##
## 1:4  -> < $40k
## 5:8  -> $40k - < $80k
## 9:11 -> $80k - < $150k
## 12   -> >= $150k
## ============================================================

fam <- Y$faminc_new

y <- rep(NA_integer_,length(fam))
y[fam %in% 1:4] <- 1
y[fam %in% 5:8] <- 2
y[fam %in% 9:11] <- 3
y[fam==12] <- 4

keep <- !is.na(y)

X <- X[keep,,drop=FALSE]
A <- A[keep,,drop=FALSE]
y <- y[keep]

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

LOOPR_faminc <- vector("list",nRep)

## ============================================================
## 100 repetitions
## ============================================================

for(r in seq_len(nRep)){
  
  cat("\n============================================================\n")
  cat("LOOPR family-income repetition:",r,"/",nRep,"\n")
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
  ## 1. Ordinal forest
  ## ----------------------------------------------------------
  
  res_X_of <- clf_ordinal_forest(X_train,X_test,y_train,y_test)
  res_A_of <- clf_ordinal_forest(A_train,A_test,y_train,y_test)
  
  ## ----------------------------------------------------------
  ## 2. Random forest
  ## ----------------------------------------------------------
  
  res_X_rf <- clf_random_forest(X_train,X_test,y_train,y_test)
  res_A_rf <- clf_random_forest(A_train,A_test,y_train,y_test)
  
  ## ----------------------------------------------------------
  ## 3. CART
  ## ----------------------------------------------------------
  
  res_X_cart <- clf_cart(X_train,X_test,y_train,y_test)
  res_A_cart <- clf_cart(A_train,A_test,y_train,y_test)
  
  ## ----------------------------------------------------------
  ## 4. Additive ordinal logit
  ## ----------------------------------------------------------
  
  res_X_add <- clf_additive_ordinal_logit(X_train,X_test,y_train,y_test)
  res_A_add <- clf_additive_ordinal_logit(A_train,A_test,y_train,y_test)
  
  ## ----------------------------------------------------------
  ## 5. Interaction ordinal logit
  ## ----------------------------------------------------------
  
  res_X_int <- clf_interaction_ordinal_logit(
    X_train,X_test,y_train,y_test,
    interaction_mode="all_pairwise"
  )
  
  res_A_int <- clf_interaction_ordinal_logit(
    A_train,A_test,y_train,y_test,
    interaction_mode="A1_A2"
  )
  
  ## ----------------------------------------------------------
  ## 6. DDE-structured plug-in Bayes
  ## ----------------------------------------------------------
  
  res_A_bayes <- clf_dde_structured_bayes(
    A_train,A_test,y_train,y_test
  )
  
  ## ----------------------------------------------------------
  ## 7. Multinomial logit
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
  ## 8. Penalized multinomial logit
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
  ## 9. Partial proportional odds
  ## ----------------------------------------------------------
  
  res_X_ppo <- clf_partial_proportional_odds(
    X_train,X_test,y_train,y_test
  )
  
  res_A_ppo <- clf_partial_proportional_odds(
    A_train,A_test,y_train,y_test
  )
  
  ## ----------------------------------------------------------
  ## 10. XGBoost
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
  ## 11. SVM
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
  ## 12. Hamming-distance kNN
  ## ----------------------------------------------------------
  
  res_A_hamming <- clf_hamming_knn(
    A_train,A_test,y_train,y_test,
    k=15
  )
  
  ## ----------------------------------------------------------
  ## 13. Tree-augmented naive Bayes
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
    "Ordinal forest"=list(Big5=res_X_of,A=res_A_of),
    "Random forest"=list(Big5=res_X_rf,A=res_A_rf),
    "CART"=list(Big5=res_X_cart,A=res_A_cart),
    "Additive ordinal logit"=list(Big5=res_X_add,A=res_A_add),
    "Interaction ordinal logit"=list(Big5=res_X_int,A=res_A_int),
    "DDE-structured plug-in Bayes"=list(Big5=NULL,A=res_A_bayes),
    "Multinomial logit"=list(Big5=res_X_multi,A=res_A_multi),
    "Penalized multinomial logit"=list(Big5=res_X_pmulti,A=res_A_pmulti),
    "Partial proportional odds"=list(Big5=res_X_ppo,A=res_A_ppo),
    "XGBoost"=list(Big5=res_X_xgb,A=res_A_xgb),
    "SVM"=list(Big5=res_X_svm,A=res_A_svm),
    "Hamming kNN"=list(Big5=NULL,A=res_A_hamming),
    "Tree-augmented naive Bayes"=list(Big5=NULL,A=res_A_tan)
  )
  
  ## ----------------------------------------------------------
  ## Store repetition
  ## ----------------------------------------------------------
  
  LOOPR_faminc[[r]] <- list(
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
      LOOPR_faminc,
      file.path(results_dir,"LOOPR_faminc.rds")
    )
    
    cat("Checkpoint saved through repetition",r,"\n")
  }
  
  ## ----------------------------------------------------------
  ## Cleanup
  ## ----------------------------------------------------------
  
  rm(
    res_X_of,res_A_of,
    res_X_rf,res_A_rf,
    res_X_cart,res_A_cart,
    res_X_add,res_A_add,
    res_X_int,res_A_int,
    res_A_bayes,
    res_X_multi,res_A_multi,
    res_X_pmulti,res_A_pmulti,
    res_X_ppo,res_A_ppo,
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
  LOOPR_faminc,
  file.path(results_dir,"LOOPR_faminc.rds")
)

cat("\nCompleted",length(LOOPR_faminc),"LOOPR family-income repetitions.\n")