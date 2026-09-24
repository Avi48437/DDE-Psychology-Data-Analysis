library(gt)

## ============================================================
## USER INPUT
## Choose one:
## "educ", "marstat", "faminc", "pew_prayer"
## ============================================================

outcome <- "educ"

## ============================================================
## Data directory
## ============================================================

data_dir <- "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology/LOPR/Analysis/Predictive_Data"

## ============================================================
## Outcome setup
## ============================================================

outcome_info <- list(
  
  educ=list(
    file="LOOPR_educ.rds",
    title="Education Prediction",
    type="ordinal"
  ),
  
  marstat=list(
    file="LOOPR_marstat.rds",
    title="Marital Status Prediction",
    type="nominal"
  ),
  
  faminc=list(
    file="LOOPR_faminc.rds",
    title="Family Income Prediction",
    type="ordinal"
  ),
  
  pew_prayer=list(
    file="LOOPR_pew_prayer.rds",
    title="Prayer Prediction",
    type="ordinal"
  )
)

if(!outcome %in% names(outcome_info))
  stop("Unknown outcome.")

info <- outcome_info[[outcome]]

sim_list <- readRDS(
  file.path(data_dir,info$file)
)

## ============================================================
## Completed repetitions
## ============================================================

completed <- sim_list[
  !vapply(sim_list,is.null,logical(1))
]

if(length(completed)==0)
  stop("No completed repetitions found.")

cat(
  "\nOutcome:",outcome,
  "\nCompleted repetitions:",length(completed),"\n\n"
)

## ============================================================
## Ordinal metrics
## ============================================================

get_ordinal_metrics <- function(res){
  
  if(is.null(res))
    return(c(
      Accuracy=NA_real_,
      MAE=NA_real_
    ))
  
  pred <- if(!is.null(res$prediction_numeric)){
    res$prediction_numeric
  } else {
    as.numeric(as.character(res$prediction))
  }
  
  truth <- if(!is.null(res$truth_numeric)){
    res$truth_numeric
  } else {
    as.numeric(as.character(res$truth))
  }
  
  c(
    Accuracy=mean(pred==truth),
    MAE=mean(abs(pred-truth))
  )
}

## ============================================================
## Nominal metrics
## ============================================================

get_nominal_metrics <- function(res){
  
  if(is.null(res))
    return(c(
      Accuracy=NA_real_,
      BalancedAccuracy=NA_real_,
      MacroF1=NA_real_
    ))
  
  cm <- res$confusion
  
  tp <- diag(cm)
  fn <- rowSums(cm)-tp
  fp <- colSums(cm)-tp
  
  recall <- tp/(tp+fn)
  precision <- tp/(tp+fp)
  f1 <- 2*precision*recall/(precision+recall)
  
  recall[!is.finite(recall)] <- 0
  f1[!is.finite(f1)] <- 0
  
  c(
    Accuracy=sum(tp)/sum(cm),
    BalancedAccuracy=mean(recall),
    MacroF1=mean(f1)
  )
}

## ============================================================
## Select metrics
## ============================================================

metric_fun <- switch(
  info$type,
  ordinal=get_ordinal_metrics,
  nominal=get_nominal_metrics
)

metric_names <- switch(
  info$type,
  ordinal=c("Accuracy","MAE"),
  nominal=c(
    "Accuracy",
    "BalancedAccuracy",
    "MacroF1"
  )
)

## ============================================================
## Summarize classifiers
## ============================================================

methods <- names(completed[[1]]$results)

summary_table <- do.call(
  rbind,
  lapply(methods,function(method){
    
    X_metrics <- do.call(
      rbind,
      lapply(completed,function(run){
        metric_fun(run$results[[method]]$Big5)
      })
    )
    
    A_metrics <- do.call(
      rbind,
      lapply(completed,function(run){
        metric_fun(run$results[[method]]$A)
      })
    )
    
    out <- data.frame(
      Classifier=method,
      NRep=length(completed),
      stringsAsFactors=FALSE
    )
    
    for(metric in metric_names){
      
      out[[paste0("Big5_",metric,"_Mean")]] <-
        mean(X_metrics[,metric],na.rm=TRUE)
      
      out[[paste0("Big5_",metric,"_SD")]] <-
        sd(X_metrics[,metric],na.rm=TRUE)
      
      out[[paste0("A_",metric,"_Mean")]] <-
        mean(A_metrics[,metric],na.rm=TRUE)
      
      out[[paste0("A_",metric,"_SD")]] <-
        sd(A_metrics[,metric],na.rm=TRUE)
    }
    
    out
  })
)

## ============================================================
## Replace NaN by NA
## ============================================================

numeric_cols <- setdiff(
  names(summary_table),
  "Classifier"
)

for(j in numeric_cols)
  summary_table[[j]][is.nan(summary_table[[j]])] <- NA_real_

## ============================================================
## DDE structured Bayes:
## do not report MAE for ordinal outcomes
## ============================================================

if(info$type=="ordinal"){
  
  bayes_idx <- summary_table$Classifier==
    "DDE-structured plug-in Bayes"
  
  summary_table$A_MAE_Mean[bayes_idx] <- NA_real_
  summary_table$A_MAE_SD[bayes_idx] <- NA_real_
}

## ============================================================
## Sort by best mean Accuracy
## ============================================================

summary_table$OrderScore <- apply(
  summary_table[,c(
    "Big5_Accuracy_Mean",
    "A_Accuracy_Mean"
  )],
  1,
  function(z){
    if(all(is.na(z))) return(NA_real_)
    max(z,na.rm=TRUE)
  }
)

summary_table <- summary_table[
  order(
    summary_table$OrderScore,
    decreasing=TRUE,
    na.last=TRUE
  ),
  ,
  drop=FALSE
]

row.names(summary_table) <- NULL

## ============================================================
## Mean (SD) formatting
## ============================================================

fmt_mean_sd <- function(mu,s){
  
  out <- rep("—",length(mu))
  
  ok <- !is.na(mu)
  
  out[ok] <- ifelse(
    is.na(s[ok]),
    sprintf("%.4f",mu[ok]),
    sprintf("%.4f (%.4f)",mu[ok],s[ok])
  )
  
  out
}

## ============================================================
## Display columns
## ============================================================

for(metric in metric_names){
  
  summary_table[[paste0("Big5_",metric)]] <-
    fmt_mean_sd(
      summary_table[[paste0("Big5_",metric,"_Mean")]],
      summary_table[[paste0("Big5_",metric,"_SD")]]
    )
  
  summary_table[[paste0("A_",metric)]] <-
    fmt_mean_sd(
      summary_table[[paste0("A_",metric,"_Mean")]],
      summary_table[[paste0("A_",metric,"_SD")]]
    )
}

## ============================================================
## GT columns
## ============================================================

display_cols <- "Classifier"

for(metric in metric_names){
  
  display_cols <- c(
    display_cols,
    paste0("Big5_",metric),
    paste0("A_",metric),
    paste0("Big5_",metric,"_Mean"),
    paste0("A_",metric,"_Mean")
  )
}

gt_data <- summary_table[,display_cols,drop=FALSE]

## ============================================================
## Build GT table
## ============================================================

result_gt <- gt_data |>
  gt() |>
  
  tab_header(
    title=md(
      paste0(
        "**LOOPR ",
        info$title,
        "**"
      )
    ),
    subtitle=paste0(
      length(completed),
      " repeated train-test analyses; entries are Mean (SD)"
    )
  )

## ============================================================
## Column labels
## ============================================================

pretty_metric <- c(
  Accuracy="Acc.",
  MAE="MAE",
  BalancedAccuracy="Bal. Acc.",
  MacroF1="Macro-F1"
)

label_args <- list(
  Classifier="Classifier"
)

for(metric in metric_names){
  
  label_args[[paste0("Big5_",metric)]] <-
    paste("Big5",pretty_metric[[metric]])
  
  label_args[[paste0("A_",metric)]] <-
    paste("DDE A",pretty_metric[[metric]])
}

result_gt <- do.call(
  cols_label,
  c(
    list(.data=result_gt),
    label_args
  )
)

## ============================================================
## Hide raw means
## ============================================================

mean_cols <- unlist(
  lapply(metric_names,function(metric){
    c(
      paste0("Big5_",metric,"_Mean"),
      paste0("A_",metric,"_Mean")
    )
  })
)

result_gt <- result_gt |>
  cols_hide(columns=all_of(mean_cols)) |>
  
  cols_align(
    align="left",
    columns=Classifier
  ) |>
  
  cols_align(
    align="center",
    columns=-Classifier
  )

## ============================================================
## Bold better Big5 vs DDE A
## ============================================================

for(metric in metric_names){
  
  X_mean <- paste0("Big5_",metric,"_Mean")
  A_mean <- paste0("A_",metric,"_Mean")
  
  X_display <- paste0("Big5_",metric)
  A_display <- paste0("A_",metric)
  
  higher_is_better <- metric!="MAE"
  
  if(higher_is_better){
    
    X_rows <- which(
      !is.na(gt_data[[X_mean]]) &
        !is.na(gt_data[[A_mean]]) &
        gt_data[[X_mean]]>=gt_data[[A_mean]]
    )
    
    A_rows <- which(
      !is.na(gt_data[[X_mean]]) &
        !is.na(gt_data[[A_mean]]) &
        gt_data[[A_mean]]>=gt_data[[X_mean]]
    )
    
  } else {
    
    X_rows <- which(
      !is.na(gt_data[[X_mean]]) &
        !is.na(gt_data[[A_mean]]) &
        gt_data[[X_mean]]<=gt_data[[A_mean]]
    )
    
    A_rows <- which(
      !is.na(gt_data[[X_mean]]) &
        !is.na(gt_data[[A_mean]]) &
        gt_data[[A_mean]]<=gt_data[[X_mean]]
    )
  }
  
  if(length(X_rows)>0){
    
    result_gt <- result_gt |>
      tab_style(
        cell_text(weight="bold"),
        cells_body(
          columns=all_of(X_display),
          rows=X_rows
        )
      )
  }
  
  if(length(A_rows)>0){
    
    result_gt <- result_gt |>
      tab_style(
        cell_text(weight="bold"),
        cells_body(
          columns=all_of(A_display),
          rows=A_rows
        )
      )
  }
}

## ============================================================
## Final formatting
## ============================================================

result_gt <- result_gt |>
  
  tab_style(
    cell_text(weight="bold"),
    cells_column_labels()
  ) |>
  
  tab_source_note(
    source_note=md(
      paste0(
        "**Bold:** better Big Five vs. DDE A mean within each paired classifier. ",
        "DDE-structured plug-in Bayes, Hamming kNN, and tree-augmented naive Bayes are A-specific."
      )
    )
  ) |>
  
  opt_row_striping() |>
  
  tab_options(
    table.font.size=11,
    heading.title.font.size=18,
    heading.subtitle.font.size=13,
    data_row.padding=6
  )

## ============================================================
## Show
## ============================================================

result_gt