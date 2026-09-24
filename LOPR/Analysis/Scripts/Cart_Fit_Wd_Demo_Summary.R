root_dir <- "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology"
results_dir <- file.path(root_dir,"LOPR","Analysis","Results")

R <- readRDS(file.path(results_dir,"LOOPR_CART_Demographics_100.rds"))

## ============================================================
## 1. Accuracy summaries
## ============================================================

accuracy_summary <- function(model1,model2){
  do.call(rbind,lapply(names(R),function(outcome){
    reps <- R[[outcome]]$repetitions
    
    b <- sapply(reps,function(z) mean(z[[model1]]$prediction_numeric==z[[model1]]$truth_numeric))
    a <- sapply(reps,function(z) mean(z[[model2]]$prediction_numeric==z[[model2]]$truth_numeric))
    
    data.frame(
      Outcome=outcome,
      Big5_Mean=mean(b),
      Big5_SD=sd(b),
      DDE_A_Mean=mean(a),
      DDE_A_SD=sd(a),
      Difference=mean(a)-mean(b)
    )
  }))
}

Accuracy_No_Demo <- accuracy_summary("Big5","DDE_A")
Accuracy_Demo <- accuracy_summary("Demo_Big5","Demo_DDE_A")

Accuracy_No_Demo[,-1] <- round(Accuracy_No_Demo[,-1],4)
Accuracy_Demo[,-1] <- round(Accuracy_Demo[,-1],4)

cat("\n================ WITHOUT DEMOGRAPHICS ================\n")
print(Accuracy_No_Demo,row.names=FALSE)

cat("\n================ WITH DEMOGRAPHICS ================\n")
print(Accuracy_Demo,row.names=FALSE)

## ============================================================
## 2. CART variable interpretation
## ============================================================

tree_summary <- function(reps,model){
  fits <- lapply(reps,function(z) z[[model]]$fit)
  vars <- sort(unique(unlist(lapply(fits,function(f) f$frame$var[f$frame$var!="<leaf>"]))))
  
  if(length(vars)==0)
    return(data.frame(
      Variable=character(0),
      Root_n=numeric(0),
      Root_pct=numeric(0),
      Used_n=numeric(0),
      Used_pct=numeric(0),
      Mean_splits=numeric(0)
    ))
  
  do.call(rbind,lapply(vars,function(v){
    root <- sapply(fits,function(f) f$frame$var[1]==v)
    nsplit <- sapply(fits,function(f) sum(f$frame$var==v))
    
    data.frame(
      Variable=v,
      Root_n=sum(root),
      Root_pct=100*mean(root),
      Used_n=sum(nsplit>0),
      Used_pct=100*mean(nsplit>0),
      Mean_splits=mean(nsplit)
    )
  }))
}

models <- c("Big5","DDE_A","Demo_Big5","Demo_DDE_A")

Tree_Summary <- do.call(rbind,lapply(names(R),function(outcome){
  do.call(rbind,lapply(models,function(model){
    z <- tree_summary(R[[outcome]]$repetitions,model)
    if(nrow(z)==0) return(NULL)
    data.frame(Outcome=outcome,Model=model,z)
  }))
}))

Tree_Summary$Root_pct <- round(Tree_Summary$Root_pct,1)
Tree_Summary$Used_pct <- round(Tree_Summary$Used_pct,1)
Tree_Summary$Mean_splits <- round(Tree_Summary$Mean_splits,2)

Tree_Summary <- Tree_Summary[
  order(Tree_Summary$Outcome,Tree_Summary$Model,-Tree_Summary$Root_pct,-Tree_Summary$Used_pct),
]

## ============================================================
## 3. Top CART variables
## ============================================================

for(outcome in names(R)){
  
  cat("\n====================",outcome,"====================\n")
  
  for(model in models){
    
    cat("\n",model,"\n")
    
    z <- subset(Tree_Summary,Outcome==outcome & Model==model)
    
    if(nrow(z)==0){
      cat("No splits in any of the 100 trees.\n")
    } else {
      z <- z[order(-z$Root_pct,-z$Used_pct),]
      print(head(z,10),row.names=FALSE)
    }
  }
}

models <- c("Big5","DDE_A","Demo_Big5","Demo_DDE_A")

Split_Summary <- do.call(rbind,lapply(names(R),function(outcome){
  do.call(rbind,lapply(models,function(model){
    ns <- sapply(R[[outcome]]$repetitions,function(z) sum(z[[model]]$fit$frame$var!="<leaf>"))
    
    data.frame(
      Outcome=outcome,
      Model=model,
      Mean_Splits=mean(ns),
      SD_Splits=sd(ns),
      Min_Splits=min(ns),
      Max_Splits=max(ns),
      Stump_pct=100*mean(ns==0)
    )
  }))
}))

Split_Summary[,3:7] <- round(Split_Summary[,3:7],2)
print(Split_Summary,row.names=FALSE)