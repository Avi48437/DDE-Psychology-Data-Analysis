## ============================================================
## Prediction_Classifiers.R
##
##  1. Ordinal forest
##  2. Random forest
##  3. CART
##  4. Additive ordinal logit
##  5. Interaction ordinal logit
##  6. DDE-structured plug-in Bayes
##  7. Multinomial logit
##  8. Penalized multinomial logit
##  9. Partial/non-proportional odds
## 10. XGBoost
## 11. SVM
## 12. Hamming kNN
## 13. Tree-augmented naive Bayes
## ============================================================


## ============================================================
## Helpers
## ============================================================

.need_pkg <- function(pkg){
  if(!requireNamespace(pkg,quietly=TRUE)) stop("Package '",pkg,"' is required.")
}

.get_classes <- function(y_train,y_test=NULL){
  
  if(anyNA(y_train)) stop("y_train contains missing values.")
  
  if(is.factor(y_train)){
    classes <- levels(droplevels(y_train))
  } else {
    z <- unique(as.character(y_train))
    zn <- suppressWarnings(as.numeric(z))
    classes <- if(all(!is.na(zn))) z[order(zn)] else sort(z)
  }
  
  if(length(classes)<2) stop("At least two outcome classes are required.")
  
  if(!is.null(y_test)){
    if(anyNA(y_test)) stop("y_test contains missing values.")
    unseen <- setdiff(unique(as.character(y_test)),classes)
    if(length(unseen)>0) stop("Test classes absent from training data: ",paste(unseen,collapse=", "))
  }
  
  classes
}

.prep_xy <- function(X_train,X_test,y_train,y_test){
  
  X_train <- as.data.frame(X_train,check.names=FALSE)
  X_test <- as.data.frame(X_test,check.names=FALSE)
  
  if(nrow(X_train)!=length(y_train)) stop("Training predictor/outcome sizes do not match.")
  if(nrow(X_test)!=length(y_test)) stop("Testing predictor/outcome sizes do not match.")
  if(!identical(names(X_train),names(X_test))) stop("Training and testing predictor columns do not match.")
  if(anyNA(X_train) || anyNA(X_test)) stop("Predictor matrices contain missing values.")
  
  list(
    X_train=X_train,
    X_test=X_test,
    classes=.get_classes(y_train,y_test)
  )
}

.standardize_prob <- function(prob,classes,n){
  
  if(is.null(prob)) return(NULL)
  
  K <- length(classes)
  
  if(length(dim(prob))==3){
    if(dim(prob)[3]!=1) stop("Unexpected three-dimensional probability array.")
    prob <- prob[,,1,drop=TRUE]
  }
  
  if(length(dim(prob))>2) stop("Unsupported probability-array dimensions.")
  
  ## Binary case: some classifiers return only P(class 2)
  if(K==2){
    
    if(is.null(dim(prob)) && length(prob)==n){
      prob <- cbind(1-as.numeric(prob),as.numeric(prob))
      colnames(prob) <- classes
    }
    
    else {
      
      prob <- as.matrix(prob)
      
      if(nrow(prob)!=n && ncol(prob)==n) prob <- t(prob)
      
      if(nrow(prob)!=n)
        stop("Probability matrix has the wrong number of rows.")
      
      if(ncol(prob)==1){
        
        p <- as.numeric(prob[,1])
        
        if(!is.null(colnames(prob)) && colnames(prob)[1] %in% classes){
          
          observed_class <- colnames(prob)[1]
          other_class <- setdiff(classes,observed_class)
          
          out <- matrix(NA_real_,nrow=n,ncol=2)
          colnames(out) <- classes
          
          out[,observed_class] <- p
          out[,other_class] <- 1-p
          
          prob <- out
          
        } else {
          
          prob <- cbind(1-p,p)
          colnames(prob) <- classes
        }
      }
      
      else if(ncol(prob)==2){
        
        if(!is.null(colnames(prob)) && all(classes %in% colnames(prob))){
          prob <- prob[,classes,drop=FALSE]
        } else {
          colnames(prob) <- classes
        }
      }
      
      else {
        stop("Binary classifier returned an unexpected number of probability columns.")
      }
    }
  }
  
  ## Multiclass case
  else {
    
    if(is.null(dim(prob))){
      if(n==1 && length(prob)==K){
        prob <- matrix(prob,nrow=1)
      } else {
        stop("Could not determine probability-matrix dimensions.")
      }
    }
    
    prob <- as.matrix(prob)
    
    if(nrow(prob)!=n && ncol(prob)==n) prob <- t(prob)
    
    if(nrow(prob)!=n)
      stop("Probability matrix has the wrong number of rows.")
    
    if(ncol(prob)!=K)
      stop("Probability matrix has the wrong number of classes.")
    
    if(!is.null(colnames(prob)) && all(classes %in% colnames(prob))){
      prob <- prob[,classes,drop=FALSE]
    } else {
      colnames(prob) <- classes
    }
  }
  
  if(any(!is.finite(prob)))
    stop("Non-finite class probabilities returned.")
  
  if(any(prob < -1e-10))
    stop("Negative class probabilities returned.")
  
  prob[prob<0] <- 0
  
  rs <- rowSums(prob)
  
  if(any(rs<=0))
    stop("Invalid class-probability rows.")
  
  prob <- prob/rs
  rownames(prob) <- NULL
  
  prob
}

.pred_from_prob <- function(prob,classes){
  factor(classes[max.col(prob,ties.method="first")],levels=classes)
}

.make_result <- function(method,pred,y_train,y_test,prob=NULL,fit=NULL,extra=NULL){
  
  classes <- .get_classes(y_train,y_test)
  
  pred <- factor(as.character(pred),levels=classes)
  truth <- factor(as.character(y_test),levels=classes)
  
  if(anyNA(pred)) stop(method,": predictions contain unknown classes.")
  
  prob <- .standardize_prob(prob,classes,length(y_test))
  
  confusion <- table(
    Actual=factor(truth,levels=classes),
    Predicted=factor(pred,levels=classes)
  )
  
  confusion_rowprop <- prop.table(confusion,margin=1)
  confusion_colprop <- prop.table(confusion,margin=2)
  
  pred_num <- suppressWarnings(as.numeric(as.character(pred)))
  truth_num <- suppressWarnings(as.numeric(as.character(truth)))
  
  if(anyNA(pred_num) || anyNA(truth_num)){
    pred_num <- NULL
    truth_num <- NULL
  }
  
  list(
    method=method,
    prediction=pred,
    prediction_numeric=pred_num,
    truth=truth,
    truth_numeric=truth_num,
    probability=prob,
    confusion=confusion,
    confusion_rowprop=confusion_rowprop,
    confusion_colprop=confusion_colprop,
    classes=classes,
    train_class_counts=table(factor(as.character(y_train),levels=classes)),
    test_class_counts=table(truth),
    fit=fit,
    extra=extra
  )
}

.check_binary <- function(X){
  
  X <- as.matrix(X)
  
  if(anyNA(X)) stop("Binary predictors contain missing values.")
  
  vals <- unique(as.vector(X))
  
  if(!all(vals %in% c(0,1))) stop("This classifier requires binary 0/1 predictors.")
  
  invisible(TRUE)
}

.bt <- function(x){
  paste0("`",gsub("`","",x,fixed=TRUE),"`")
}

.resolve_interaction_mode <- function(nms,mode){
  
  if(mode!="auto") return(mode)
  
  if(any(grepl("^A1_",nms)) && any(grepl("^A2_",nms))) "A1_A2" else "all_pairwise"
}

.rhs_string <- function(nms,mode=c("additive","all_pairwise","A1_A2")){
  
  mode <- match.arg(mode)
  vars <- .bt(nms)
  
  if(mode=="additive") return(paste(vars,collapse=" + "))
  if(mode=="all_pairwise") return(paste0("(",paste(vars,collapse=" + "),")^2"))
  
  A1_names <- grep("^A1_",nms,value=TRUE)
  A2_names <- grep("^A2_",nms,value=TRUE)
  
  if(length(A1_names)==0 || length(A2_names)==0)
    stop("A1_A2 mode requires A1_ and A2_ predictor columns.")
  
  interactions <- unlist(lapply(.bt(A1_names),function(x) paste0(x,":",.bt(A2_names))),use.names=FALSE)
  
  paste(c(vars,interactions),collapse=" + ")
}

.model_matrix_pair <- function(X_train,X_test,mode="additive"){
  
  X_train <- as.data.frame(X_train,check.names=FALSE)
  X_test <- as.data.frame(X_test,check.names=FALSE)
  
  if(mode=="auto") mode <- .resolve_interaction_mode(names(X_train),"auto")
  
  f <- as.formula(paste("~",.rhs_string(names(X_train),mode)))
  
  xtr <- model.matrix(f,data=X_train)
  xte <- model.matrix(f,data=X_test)
  
  xtr <- xtr[,colnames(xtr)!="(Intercept)",drop=FALSE]
  xte <- xte[,colnames(xte)!="(Intercept)",drop=FALSE]
  
  missing_test <- setdiff(colnames(xtr),colnames(xte))
  
  if(length(missing_test)>0){
    add <- matrix(0,nrow(xte),length(missing_test),dimnames=list(NULL,missing_test))
    xte <- cbind(xte,add)
  }
  
  extra_test <- setdiff(colnames(xte),colnames(xtr))
  if(length(extra_test)>0) xte <- xte[,setdiff(colnames(xte),extra_test),drop=FALSE]
  
  xte <- xte[,colnames(xtr),drop=FALSE]
  
  storage.mode(xtr) <- "double"
  storage.mode(xte) <- "double"
  
  list(train=xtr,test=xte,formula=f,mode=mode)
}


## ============================================================
## 1. Ordinal forest
## ============================================================

clf_ordinal_forest <- function(X_train,X_test,y_train,y_test,
                               nsets=50,ntreeperdiv=30,ntreefinal=500,
                               nbest=3,num.threads=8){
  
  .need_pkg("ordinalForest")
  d <- .prep_xy(X_train,X_test,y_train,y_test)
  
  y_ord <- ordered(as.character(y_train),levels=d$classes)
  
  fit <- ordinalForest::ordfor(
    depvar="y",
    data=data.frame(y=y_ord,d$X_train,check.names=FALSE),
    perffunction="probability",
    nsets=nsets,
    ntreeperdiv=ntreeperdiv,
    ntreefinal=ntreefinal,
    nbest=nbest,
    num.threads=num.threads
  )
  
  raw <- predict(fit,newdata=d$X_test)
  
  .make_result(
    "Ordinal forest",
    raw$ypred,
    y_train,y_test,
    prob=raw$classprobs,
    fit=fit,
    extra=list(raw_prediction=raw)
  )
}


## ============================================================
## 2. Random forest
## ============================================================

clf_random_forest <- function(X_train,X_test,y_train,y_test,ntree=500,mtry=NULL){
  
  .need_pkg("randomForest")
  d <- .prep_xy(X_train,X_test,y_train,y_test)
  
  y_fac <- factor(as.character(y_train),levels=d$classes)
  
  args <- list(x=d$X_train,y=y_fac,ntree=ntree,importance=TRUE)
  if(!is.null(mtry)) args$mtry <- mtry
  
  fit <- do.call(randomForest::randomForest,args)
  
  pred <- predict(fit,d$X_test,type="response")
  prob <- predict(fit,d$X_test,type="prob")
  
  .make_result(
    "Random forest",
    pred,
    y_train,y_test,
    prob=prob,
    fit=fit,
    extra=list(importance=randomForest::importance(fit))
  )
}


## ============================================================
## 3. CART
## ============================================================

clf_cart <- function(X_train,X_test,y_train,y_test,cp=0.01,minsplit=20,maxdepth=30){
  
  .need_pkg("rpart")
  d <- .prep_xy(X_train,X_test,y_train,y_test)
  
  y_fac <- factor(as.character(y_train),levels=d$classes)
  
  fit <- rpart::rpart(
    y ~ .,
    data=data.frame(y=y_fac,d$X_train,check.names=FALSE),
    method="class",
    control=rpart::rpart.control(cp=cp,minsplit=minsplit,maxdepth=maxdepth),
    model=TRUE,x=TRUE,y=TRUE
  )
  
  pred <- predict(fit,d$X_test,type="class")
  prob <- predict(fit,d$X_test,type="prob")
  
  .make_result(
    "CART",
    pred,
    y_train,y_test,
    prob=prob,
    fit=fit,
    extra=list(cptable=fit$cptable)
  )
}


## ============================================================
## 4. Additive ordinal logit
## ============================================================

clf_additive_ordinal_logit <- function(X_train,X_test,y_train,y_test,method="logistic"){
  
  .need_pkg("MASS")
  d <- .prep_xy(X_train,X_test,y_train,y_test)
  
  y_ord <- ordered(as.character(y_train),levels=d$classes)
  
  fit <- MASS::polr(
    y ~ .,
    data=data.frame(y=y_ord,d$X_train,check.names=FALSE),
    method=method,
    Hess=FALSE,
    model=TRUE
  )
  
  pred <- predict(fit,newdata=d$X_test,type="class")
  prob <- predict(fit,newdata=d$X_test,type="probs")
  
  .make_result(
    "Additive ordinal logit",
    pred,
    y_train,y_test,
    prob=prob,
    fit=fit,
    extra=list(coefficients=coef(fit),thresholds=fit$zeta)
  )
}


## ============================================================
## 5. Interaction ordinal logit
##
## auto:
##   X      -> all pairwise X interactions
##   A      -> A1:A2 only
##   X + A  -> all main effects + A1:A2 only
## ============================================================

clf_interaction_ordinal_logit <- function(
    X_train,X_test,y_train,y_test,
    interaction_mode=c("auto","all_pairwise","A1_A2"),
    method="logistic"){
  
  .need_pkg("MASS")
  d <- .prep_xy(X_train,X_test,y_train,y_test)
  
  interaction_mode <- match.arg(interaction_mode)
  mode <- .resolve_interaction_mode(names(d$X_train),interaction_mode)
  
  y_ord <- ordered(as.character(y_train),levels=d$classes)
  f <- as.formula(paste("y ~",.rhs_string(names(d$X_train),mode)))
  
  fit <- MASS::polr(
    f,
    data=data.frame(y=y_ord,d$X_train,check.names=FALSE),
    method=method,
    Hess=FALSE,
    model=TRUE
  )
  
  pred <- predict(fit,newdata=d$X_test,type="class")
  prob <- predict(fit,newdata=d$X_test,type="probs")
  
  .make_result(
    "Interaction ordinal logit",
    pred,
    y_train,y_test,
    prob=prob,
    fit=fit,
    extra=list(
      formula=f,
      interaction_mode=mode,
      coefficients=coef(fit),
      thresholds=fit$zeta
    )
  )
}


## ============================================================
## 6. DDE-structured plug-in Bayes
##
## P(A1,A2 | Y=c)
## =
## P(A2 | Y=c) prod_j P(A1_j | A2,Y=c)
## ============================================================

clf_dde_structured_bayes <- function(A_train,A_test,y_train,y_test,alpha=1){
  
  d <- .prep_xy(A_train,A_test,y_train,y_test)
  
  A1_names <- grep("^A1_",names(d$X_train),value=TRUE)
  A2_names <- grep("^A2_",names(d$X_train),value=TRUE)
  
  if(length(A1_names)==0) stop("No A1_ coordinates found.")
  if(length(A2_names)!=1) stop("Exactly one A2 coordinate is required.")
  
  .check_binary(d$X_train[,c(A1_names,A2_names),drop=FALSE])
  .check_binary(d$X_test[,c(A1_names,A2_names),drop=FALSE])
  
  classes <- d$classes
  K <- length(classes)
  
  y_idx <- as.integer(factor(as.character(y_train),levels=classes))
  
  A1_tr <- as.matrix(d$X_train[,A1_names,drop=FALSE])
  A1_te <- as.matrix(d$X_test[,A1_names,drop=FALSE])
  A2_tr <- d$X_train[[A2_names]]
  A2_te <- d$X_test[[A2_names]]
  
  J <- ncol(A1_tr)
  
  pi_hat <- (tabulate(y_idx,nbins=K)+alpha)/(length(y_idx)+K*alpha)
  q_hat <- numeric(K)
  theta_hat <- array(NA_real_,c(K,2,J))
  
  for(cc in seq_len(K)){
    
    idx_c <- y_idx==cc
    q_hat[cc] <- (sum(A2_tr[idx_c]==1)+alpha)/(sum(idx_c)+2*alpha)
    
    for(a in 0:1){
      idx <- idx_c & A2_tr==a
      theta_hat[cc,a+1,] <- (colSums(A1_tr[idx,,drop=FALSE])+alpha)/(sum(idx)+2*alpha)
    }
  }
  
  post <- matrix(NA_real_,nrow(d$X_test),K,dimnames=list(NULL,classes))
  
  for(i in seq_len(nrow(d$X_test))){
    
    a2 <- A2_te[i]
    logp <- numeric(K)
    
    for(cc in seq_len(K)){
      logp[cc] <- log(pi_hat[cc]) +
        dbinom(a2,1,q_hat[cc],log=TRUE) +
        sum(dbinom(A1_te[i,],1,theta_hat[cc,a2+1,],log=TRUE))
    }
    
    logp <- logp-max(logp)
    post[i,] <- exp(logp)/sum(exp(logp))
  }
  
  pred <- .pred_from_prob(post,classes)
  
  .make_result(
    "DDE-structured plug-in Bayes",
    pred,
    y_train,y_test,
    prob=post,
    extra=list(
      prior=pi_hat,
      A2_probability=q_hat,
      A1_probability=theta_hat,
      A1_names=A1_names,
      A2_name=A2_names
    )
  )
}


## ============================================================
## 7. Multinomial logit
##
## auto:
##   X      -> all pairwise X interactions
##   A      -> A1:A2
##   X + A  -> all main effects + A1:A2
## ============================================================

clf_multinomial_logit <- function(
    X_train,X_test,y_train,y_test,
    formula_mode=c("auto","additive","all_pairwise","A1_A2"),
    decay=0,maxit=1000){
  
  .need_pkg("nnet")
  d <- .prep_xy(X_train,X_test,y_train,y_test)
  
  formula_mode <- match.arg(formula_mode)
  if(formula_mode=="auto") formula_mode <- .resolve_interaction_mode(names(d$X_train),"auto")
  
  y_fac <- factor(as.character(y_train),levels=d$classes)
  f <- as.formula(paste("y ~",.rhs_string(names(d$X_train),formula_mode)))
  
  fit <- nnet::multinom(
    f,
    data=data.frame(y=y_fac,d$X_train,check.names=FALSE),
    decay=decay,
    maxit=maxit,
    Hess=FALSE,
    trace=FALSE
  )
  
  prob <- predict(fit,newdata=d$X_test,type="probs")
  prob <- .standardize_prob(prob,d$classes,nrow(d$X_test))
  pred <- .pred_from_prob(prob,d$classes)
  
  .make_result(
    "Multinomial logit",
    pred,
    y_train,y_test,
    prob=prob,
    fit=fit,
    extra=list(
      formula=f,
      formula_mode=formula_mode,
      coefficients=coef(fit)
    )
  )
}


## ============================================================
## 8. Penalized multinomial logit
##
## alpha=0   ridge
## alpha=1   lasso
## alpha=.5  elastic net
## ============================================================

clf_penalized_multinomial <- function(
    X_train,X_test,y_train,y_test,
    formula_mode=c("auto","additive","all_pairwise","A1_A2"),
    alpha=0.5,nfolds=5,
    lambda_choice=c("lambda.min","lambda.1se")){
  
  .need_pkg("glmnet")
  d <- .prep_xy(X_train,X_test,y_train,y_test)
  
  formula_mode <- match.arg(formula_mode)
  lambda_choice <- match.arg(lambda_choice)
  
  if(formula_mode=="auto") formula_mode <- .resolve_interaction_mode(names(d$X_train),"auto")
  
  mm <- .model_matrix_pair(d$X_train,d$X_test,formula_mode)
  y_fac <- factor(as.character(y_train),levels=d$classes)
  
  fit <- glmnet::cv.glmnet(
    x=mm$train,
    y=y_fac,
    family="multinomial",
    alpha=alpha,
    nfolds=nfolds,
    type.measure="class"
  )
  
  prob <- predict(fit,newx=mm$test,s=lambda_choice,type="response")
  prob <- .standardize_prob(prob,d$classes,nrow(d$X_test))
  pred <- .pred_from_prob(prob,d$classes)
  
  lambda <- if(lambda_choice=="lambda.min") fit$lambda.min else fit$lambda.1se
  
  .make_result(
    "Penalized multinomial logit",
    pred,
    y_train,y_test,
    prob=prob,
    fit=fit,
    extra=list(
      formula=mm$formula,
      formula_mode=formula_mode,
      alpha=alpha,
      lambda_choice=lambda_choice,
      lambda=lambda,
      design_columns=colnames(mm$train)
    )
  )
}

## ============================================================
## 9. Partial proportional-odds model
##
## Default:
##   1. fit proportional-odds model
##   2. test each predictor for violation using nominal_test()
##   3. allow only significant predictors to be non-proportional
##
## If the selected PPO model is rank deficient or gives invalid
## probabilities, safely fall back to the proportional-odds fit.
## ============================================================

clf_partial_proportional_odds <- function(
    X_train,X_test,y_train,y_test,
    nominal_vars=NULL,
    nominal_alpha=0.05,
    p_adjust="BH",
    max_nominal=3,
    link="logit"){
  
  .need_pkg("ordinal")
  d <- .prep_xy(X_train,X_test,y_train,y_test)
  
  y_ord <- ordered(as.character(y_train),levels=d$classes)
  train_data <- data.frame(y=y_ord,d$X_train,check.names=FALSE)
  
  location_formula <- as.formula(
    paste("y ~",paste(.bt(names(d$X_train)),collapse=" + "))
  )
  
  ## Proportional-odds base model
  fit_po <- ordinal::clm(
    formula=location_formula,
    data=train_data,
    link=link
  )
  
  nominal_test_result <- NULL
  selected_nominal <- character(0)
  
  ## Automatic selection of non-proportional effects
  if(is.null(nominal_vars)){
    
    nominal_test_result <- suppressWarnings(
      ordinal::nominal_test(fit_po)
    )
    
    nt <- as.data.frame(nominal_test_result)
    pcol <- grep("^Pr\\(",names(nt),value=TRUE)
    
    if(length(pcol)>0){
      
      terms <- rownames(nt)
      pval <- nt[[pcol[1]]]
      
      keep <- terms!="<none>" & !is.na(pval) & is.finite(pval)
      
      if(any(keep)){
        
        terms <- terms[keep]
        pval <- pval[keep]
        padj <- p.adjust(pval,method=p_adjust)
        
        ord <- order(padj)
        terms <- terms[ord]
        padj <- padj[ord]
        
        selected_nominal <- terms[padj<nominal_alpha]
        
        if(length(selected_nominal)>max_nominal)
          selected_nominal <- selected_nominal[seq_len(max_nominal)]
      }
    }
    
  } else {
    
    selected_nominal <- nominal_vars
    
    bad <- setdiff(selected_nominal,names(d$X_train))
    
    if(length(bad)>0)
      stop("Unknown nominal variables: ",paste(bad,collapse=", "))
  }
  
  ## Fit PPO only when at least one effect is selected
  if(length(selected_nominal)>0){
    
    nominal_formula <- as.formula(
      paste("~",paste(.bt(selected_nominal),collapse=" + "))
    )
    
    fit_ppo <- suppressWarnings(
      ordinal::clm(
        formula=location_formula,
        nominal=nominal_formula,
        data=train_data,
        link=link
      )
    )
    
    aliased <- any(unlist(fit_ppo$aliased),na.rm=TRUE)
    
    prob_ppo <- suppressWarnings(
      predict(fit_ppo,newdata=d$X_test,type="prob")$fit
    )
    
    valid_prob <- !aliased &&
      !is.null(prob_ppo) &&
      all(is.finite(prob_ppo)) &&
      all(prob_ppo>=-1e-10) &&
      all(prob_ppo<=1+1e-10) &&
      all(rowSums(prob_ppo)>0)
    
  } else {
    
    nominal_formula <- NULL
    fit_ppo <- NULL
    prob_ppo <- NULL
    valid_prob <- FALSE
  }
  
  ## Use PPO when valid; otherwise fall back to PO
  if(valid_prob){
    
    fit <- fit_ppo
    prob <- .standardize_prob(
      prob_ppo,
      d$classes,
      nrow(d$X_test)
    )
    
    used_model <- "partial proportional odds"
    fallback <- FALSE
    
  } else {
    
    fit <- fit_po
    
    prob <- predict(
      fit_po,
      newdata=d$X_test,
      type="prob"
    )$fit
    
    prob <- .standardize_prob(
      prob,
      d$classes,
      nrow(d$X_test)
    )
    
    used_model <- "proportional odds fallback"
    fallback <- length(selected_nominal)>0
  }
  
  pred <- .pred_from_prob(prob,d$classes)
  
  .make_result(
    method="Partial proportional odds",
    pred=pred,
    y_train=y_train,
    y_test=y_test,
    prob=prob,
    fit=fit,
    extra=list(
      location_formula=location_formula,
      nominal_formula=nominal_formula,
      nominal_test=nominal_test_result,
      selected_nominal=selected_nominal,
      nominal_alpha=nominal_alpha,
      p_adjust=p_adjust,
      used_model=used_model,
      fallback=fallback
    )
  )
}
## ============================================================
## 10. XGBoost
## ============================================================

clf_xgboost <- function(
    X_train,X_test,y_train,y_test,
    nrounds=150,max_depth=3,eta=0.05,
    subsample=0.8,colsample_bytree=0.8,
    params=list(),verbose=0){
  
  .need_pkg("xgboost")
  d <- .prep_xy(X_train,X_test,y_train,y_test)
  
  mm <- .model_matrix_pair(d$X_train,d$X_test,"additive")
  
  K <- length(d$classes)
  y_idx <- as.integer(factor(as.character(y_train),levels=d$classes))-1
  
  dtrain <- xgboost::xgb.DMatrix(mm$train,label=y_idx)
  dtest <- xgboost::xgb.DMatrix(mm$test)
  
  if(K==2){
    
    default_params <- list(
      objective="binary:logistic",
      eval_metric="logloss",
      max_depth=max_depth,
      eta=eta,
      subsample=subsample,
      colsample_bytree=colsample_bytree
    )
    
  } else {
    
    default_params <- list(
      objective="multi:softprob",
      eval_metric="mlogloss",
      num_class=K,
      max_depth=max_depth,
      eta=eta,
      subsample=subsample,
      colsample_bytree=colsample_bytree
    )
  }
  
  params <- utils::modifyList(default_params,params)
  
  fit <- xgboost::xgb.train(
    params=params,
    data=dtrain,
    nrounds=nrounds,
    verbose=verbose
  )
  
  raw <- predict(fit,dtest)
  
  if(K==2){
    
    if(is.matrix(raw)) raw <- raw[,1]
    raw <- as.numeric(raw)
    prob <- cbind(1-raw,raw)
    
  } else {
    
    if(is.matrix(raw)){
      prob <- raw
    } else {
      if(length(raw)!=nrow(d$X_test)*K) stop("Unexpected XGBoost multiclass prediction dimensions.")
      prob <- matrix(as.numeric(raw),nrow=nrow(d$X_test),ncol=K,byrow=TRUE)
    }
  }
  
  prob <- .standardize_prob(prob,d$classes,nrow(d$X_test))
  pred <- .pred_from_prob(prob,d$classes)
  
  .make_result(
    "XGBoost",
    pred,
    y_train,y_test,
    prob=prob,
    fit=fit,
    extra=list(
      params=params,
      nrounds=nrounds,
      design_columns=colnames(mm$train),
      raw_prediction=raw
    )
  )
}


## ============================================================
## 11. SVM
## ============================================================

clf_svm <- function(
    X_train,X_test,y_train,y_test,
    kernel="radial",cost=1,gamma=NULL,scale=TRUE){
  
  .need_pkg("e1071")
  d <- .prep_xy(X_train,X_test,y_train,y_test)
  
  y_fac <- factor(as.character(y_train),levels=d$classes)
  
  args <- list(
    x=d$X_train,
    y=y_fac,
    type="C-classification",
    kernel=kernel,
    cost=cost,
    scale=scale,
    probability=TRUE
  )
  
  if(!is.null(gamma)) args$gamma <- gamma
  
  fit <- do.call(e1071::svm,args)
  
  raw <- predict(
    fit,
    d$X_test,
    probability=TRUE,
    decision.values=TRUE
  )
  
  prob <- attr(raw,"probabilities")
  decision_values <- attr(raw,"decision.values")
  
  prob <- .standardize_prob(prob,d$classes,nrow(d$X_test))
  
  .make_result(
    paste0("SVM (",kernel,")"),
    raw,
    y_train,y_test,
    prob=prob,
    fit=fit,
    extra=list(
      decision_values=decision_values,
      kernel=kernel,
      cost=cost,
      gamma=gamma
    )
  )
}


## ============================================================
## 12. Hamming-distance kNN
## Binary predictors only
## ============================================================

clf_hamming_knn <- function(
    X_train,X_test,y_train,y_test,
    k=15,include_ties=TRUE){
  
  d <- .prep_xy(X_train,X_test,y_train,y_test)
  
  .check_binary(d$X_train)
  .check_binary(d$X_test)
  
  Xtr <- as.matrix(d$X_train)
  Xte <- as.matrix(d$X_test)
  
  K <- length(d$classes)
  y_idx <- as.integer(factor(as.character(y_train),levels=d$classes))
  
  k <- min(k,nrow(Xtr))
  
  prob <- matrix(
    0,
    nrow=nrow(Xte),
    ncol=K,
    dimnames=list(NULL,d$classes)
  )
  
  neighbor_index <- vector("list",nrow(Xte))
  neighbor_distance <- vector("list",nrow(Xte))
  
  for(i in seq_len(nrow(Xte))){
    
    dist_i <- rowSums(sweep(Xtr,2,Xte[i,],FUN="!="))
    ord <- order(dist_i)
    
    if(include_ties){
      cutoff <- dist_i[ord[k]]
      idx <- which(dist_i<=cutoff)
    } else {
      idx <- ord[seq_len(k)]
    }
    
    counts <- tabulate(y_idx[idx],nbins=K)
    prob[i,] <- counts/sum(counts)
    
    neighbor_index[[i]] <- idx
    neighbor_distance[[i]] <- dist_i[idx]
  }
  
  pred <- .pred_from_prob(prob,d$classes)
  
  .make_result(
    "Hamming kNN",
    pred,
    y_train,y_test,
    prob=prob,
    fit=list(
      X_train=Xtr,
      y_train=y_train,
      k=k,
      include_ties=include_ties
    ),
    extra=list(
      neighbor_index=neighbor_index,
      neighbor_distance=neighbor_distance,
      k=k,
      include_ties=include_ties
    )
  )
}


## ============================================================
## TAN helpers
## ============================================================

.binary_cmi <- function(x1,x2,y_idx,K,alpha=1){
  
  class_prob <- (tabulate(y_idx,nbins=K)+alpha)/(length(y_idx)+K*alpha)
  ans <- 0
  
  for(cc in seq_len(K)){
    
    idx <- y_idx==cc
    joint <- matrix(alpha,2,2)
    
    for(a in 0:1){
      for(b in 0:1){
        joint[a+1,b+1] <- joint[a+1,b+1] + sum(x1[idx]==a & x2[idx]==b)
      }
    }
    
    p_joint <- joint/sum(joint)
    p1 <- rowSums(p_joint)
    p2 <- colSums(p_joint)
    
    ans <- ans + class_prob[cc]*sum(p_joint*log(p_joint/outer(p1,p2)))
  }
  
  ans
}

.maximum_spanning_tree <- function(W,root){
  
  p <- nrow(W)
  if(p==1) return(NA_integer_)
  
  selected <- rep(FALSE,p)
  selected[root] <- TRUE
  parent <- rep(NA_integer_,p)
  
  for(step in 2:p){
    
    best_weight <- -Inf
    best_u <- best_v <- NA_integer_
    
    for(u in which(selected)){
      for(v in which(!selected)){
        if(W[u,v]>best_weight){
          best_weight <- W[u,v]
          best_u <- u
          best_v <- v
        }
      }
    }
    
    if(is.na(best_v)) stop("Could not construct TAN spanning tree.")
    
    parent[best_v] <- best_u
    selected[best_v] <- TRUE
  }
  
  parent
}


## ============================================================
## 13. Tree-augmented naive Bayes
## Binary predictors only
## ============================================================

clf_tan_binary <- function(
    X_train,X_test,y_train,y_test,
    alpha=1,root=NULL){
  
  d <- .prep_xy(X_train,X_test,y_train,y_test)
  
  .check_binary(d$X_train)
  .check_binary(d$X_test)
  
  Xtr <- as.matrix(d$X_train)
  Xte <- as.matrix(d$X_test)
  
  p <- ncol(Xtr)
  K <- length(d$classes)
  
  y_idx <- as.integer(factor(as.character(y_train),levels=d$classes))
  
  cmi <- matrix(
    0,p,p,
    dimnames=list(colnames(Xtr),colnames(Xtr))
  )
  
  if(p>1){
    for(j in 1:(p-1)){
      for(k in (j+1):p){
        cmi[j,k] <- .binary_cmi(Xtr[,j],Xtr[,k],y_idx,K,alpha)
        cmi[k,j] <- cmi[j,k]
      }
    }
  }
  
  if(is.null(root)){
    
    root_idx <- which.max(rowSums(cmi))
    
  } else if(is.character(root)){
    
    if(length(root)!=1) stop("root must specify one predictor.")
    root_idx <- match(root,colnames(Xtr))
    if(is.na(root_idx)) stop("Requested TAN root was not found.")
    
  } else {
    
    root_idx <- as.integer(root)[1]
    if(root_idx<1 || root_idx>p) stop("Invalid TAN root index.")
  }
  
  parent <- .maximum_spanning_tree(cmi,root_idx)
  
  prior <- (tabulate(y_idx,nbins=K)+alpha)/(length(y_idx)+K*alpha)
  
  root_prob <- numeric(K)
  theta <- array(NA_real_,c(K,2,p))
  
  for(cc in seq_len(K)){
    
    idx_c <- y_idx==cc
    
    root_prob[cc] <- (
      sum(Xtr[idx_c,root_idx]==1)+alpha
    ) / (
      sum(idx_c)+2*alpha
    )
    
    for(j in seq_len(p)){
      
      if(j==root_idx) next
      
      par_j <- parent[j]
      
      for(a in 0:1){
        
        idx <- idx_c & Xtr[,par_j]==a
        
        theta[cc,a+1,j] <- (
          sum(Xtr[idx,j]==1)+alpha
        ) / (
          sum(idx)+2*alpha
        )
      }
    }
  }
  
  post <- matrix(
    NA_real_,
    nrow=nrow(Xte),
    ncol=K,
    dimnames=list(NULL,d$classes)
  )
  
  for(i in seq_len(nrow(Xte))){
    
    logp <- numeric(K)
    
    for(cc in seq_len(K)){
      
      logp[cc] <- log(prior[cc]) +
        dbinom(Xte[i,root_idx],1,root_prob[cc],log=TRUE)
      
      for(j in seq_len(p)){
        
        if(j==root_idx) next
        
        par_j <- parent[j]
        par_value <- Xte[i,par_j]
        
        logp[cc] <- logp[cc] +
          dbinom(Xte[i,j],1,theta[cc,par_value+1,j],log=TRUE)
      }
    }
    
    logp <- logp-max(logp)
    post[i,] <- exp(logp)/sum(exp(logp))
  }
  
  pred <- .pred_from_prob(post,d$classes)
  
  parent_names <- rep(NA_character_,p)
  
  for(j in seq_len(p)){
    if(!is.na(parent[j])) parent_names[j] <- colnames(Xtr)[parent[j]]
  }
  
  names(parent_names) <- colnames(Xtr)
  
  fit <- list(
    CMI=cmi,
    root=root_idx,
    root_name=colnames(Xtr)[root_idx],
    parent=parent,
    parent_names=parent_names,
    prior=prior,
    root_probability=root_prob,
    conditional_probability=theta
  )
  
  .make_result(
    "Tree-augmented naive Bayes",
    pred,
    y_train,y_test,
    prob=post,
    fit=fit,
    extra=fit
  )
}