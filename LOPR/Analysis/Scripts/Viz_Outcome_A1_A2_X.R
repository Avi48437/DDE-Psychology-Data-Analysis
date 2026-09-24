library(ggplot2)
library(tidyr)
library(patchwork)
library(grid)
library(scales)

plot_A1_X_by_outcome_LOOPR <- function(X,A,Y,Y_Interpretation,outcome,
                                       grouping=FALSE,group_list=NULL,
                                       Discard=FALSE,discard_codes=NULL){
  
  domains <- c("EXT","NEM","AGR","CSN","OPN")
  
  if(!all(domains %in% names(X))) stop("X must contain EXT, NEM, AGR, CSN, and OPN.")
  if(!outcome %in% names(Y)) stop("Outcome not found in Y.")
  if(nrow(X)!=nrow(A) || nrow(X)!=nrow(Y)) stop("X, A, and Y must have the same number of rows.")
  if(!all(c("Outcome","Code","Meaning") %in% names(Y_Interpretation))) stop("Y_Interpretation must contain Outcome, Code, and Meaning.")
  if(grouping && is.null(group_list)) stop("group_list must be supplied when grouping=TRUE.")
  if(Discard && is.null(discard_codes)) stop("discard_codes must be supplied when Discard=TRUE.")
  
  shorten <- function(outcome,x){
    
    maps <- list(
      
      educ=c(
        "Less than 7th grade"="<7th",
        "Junior high school (7th, 8th, or 9th grade)"="Jr high",
        "Some high school but no diploma or equivalent (10th or 11th grade)"="Some HS",
        "High school graduate (high school diploma or equivalent including GED)"="HS grad",
        "Some college (at least one year) but no degree"="Some college",
        "Associate degree in college (2-year)"="Associate",
        "Bachelor’s degree in college (4-year)"="Bachelor's",
        "Master’s degree"="Master's",
        "Doctoral degree (for example, PhD)"="Doctoral",
        "Graduate professional degree (for example, JD, MD)"="Prof degree"
      ),
      
      faminc_new=c(
        "Less than $10,000"="<10k",
        "$10,000 to $19,999"="10-20k",
        "$20,000 to $29,999"="20-30k",
        "$30,000 to $39,999"="30-40k",
        "$40,000 to $49,999"="40-50k",
        "$50,000 to $59,999"="50-60k",
        "$60,000 to $69,999"="60-70k",
        "$70,000 to $79,999"="70-80k",
        "$80,000 to $89,999"="80-90k",
        "$90,000 to $99,999"="90-100k",
        "$100,000 to $149,999"="100-150k",
        "$150,000 or more"="150k+"
      ),
      
      marstat=c(
        "Married"="Married",
        "Widowed"="Widowed",
        "Divorced"="Divorced",
        "Separated"="Separated",
        "Never married"="Never married"
      )
    )
    
    m <- maps[[outcome]]
    if(is.null(m)) return(x)
    
    ii <- x %in% names(m)
    x[ii] <- unname(m[x[ii]])
    x
  }
  
  X5 <- as.matrix(X[,domains,drop=FALSE])
  A1 <- as.matrix(A[,grep("^A1_",names(A),value=TRUE),drop=FALSE])
  A2 <- as.matrix(A[,grep("^A2_",names(A),value=TRUE),drop=FALSE])
  y <- Y[[outcome]]
  
  ## ============================================================
  ## Read category labels automatically from Y_Interpretation
  ## ============================================================
  
  lab <- Y_Interpretation[
    Y_Interpretation$Outcome==outcome &
      !is.na(Y_Interpretation$Code) &
      Y_Interpretation$Code!="",
    ,
    drop=FALSE
  ]
  
  if(nrow(lab)>0){
    
    lab$Code <- as.numeric(lab$Code)
    lab <- lab[order(lab$Code),,drop=FALSE]
    
    if(anyDuplicated(lab$Code))
      stop("Duplicated codes found in Y_Interpretation for outcome: ",outcome)
    
    labels <- setNames(as.character(lab$Meaning),as.character(lab$Code))
    
  } else {
    labels <- NULL
  }
  
  ## ============================================================
  ## Remove missing/discarded observations
  ## ============================================================
  
  keep <- !is.na(y)
  
  if(Discard)
    keep <- keep & !(y %in% discard_codes)
  
  X5 <- X5[keep,,drop=FALSE]
  A1 <- A1[keep,,drop=FALSE]
  A2 <- A2[keep,,drop=FALSE]
  y <- y[keep]
  
  if(length(y)==0)
    stop("No observations remain for outcome: ",outcome)
  
  ## ============================================================
  ## Convert codes to interpretation labels
  ## ============================================================
  
  if(is.null(labels)){
    
    original_label <- as.character(y)
    
    if(length(unique(y))>30)
      warning(outcome," has ",length(unique(y))," unique values and no coded value labels.")
    
  } else {
    
    original_label <- unname(labels[as.character(y)])
    
    missing_label <- is.na(original_label)
    original_label[missing_label] <- as.character(y[missing_label])
  }
  
  original_label <- shorten(outcome,original_label)
  
  ## ============================================================
  ## Optional grouping
  ## ============================================================
  
  if(grouping){
    
    grouped <- rep(NA_character_,length(y))
    group_names <- names(group_list)
    
    if(is.null(group_names) || any(group_names=="")){
      
      group_names <- vapply(group_list,function(g){
        
        labs <- if(is.null(labels)){
          as.character(g)
        } else {
          unname(labels[as.character(g)])
        }
        
        missing_lab <- is.na(labs)
        labs[missing_lab] <- as.character(g[missing_lab])
        
        labs <- shorten(outcome,labs)
        paste(labs,collapse="/")
        
      },character(1))
    }
    
    for(j in seq_along(group_list))
      grouped[y %in% group_list[[j]]] <- group_names[j]
    
    untouched <- is.na(grouped)
    grouped[untouched] <- original_label[untouched]
    
    group_label <- grouped
    
    if(is.numeric(y)){
      ord_codes <- order(y)
      category_order <- unique(group_label[ord_codes])
    } else {
      category_order <- unique(group_label)
    }
    
  } else {
    
    group_label <- original_label
    
    if(is.null(labels)){
      
      if(is.numeric(y)){
        code_order <- sort(unique(y))
        category_order <- as.character(code_order)
      } else {
        category_order <- unique(as.character(y))
      }
      
    } else {
      
      code_order <- as.numeric(names(labels))
      code_order <- code_order[code_order %in% unique(y)]
      
      ordered_labels <- unname(labels[as.character(code_order)])
      category_order <- shorten(outcome,ordered_labels)
      
      extra <- unique(group_label[!group_label %in% category_order])
      category_order <- c(category_order,extra)
    }
  }
  
  group <- factor(group_label,levels=category_order)
  
  keep <- !is.na(group)
  
  X5 <- X5[keep,,drop=FALSE]
  A1 <- A1[keep,,drop=FALSE]
  A2 <- A2[keep,,drop=FALSE]
  group <- droplevels(group[keep])
  group_label <- as.character(group)
  
  ## ============================================================
  ## Sort observations by outcome group
  ## ============================================================
  
  ord <- order(group)
  
  X5 <- X5[ord,,drop=FALSE]
  A1 <- A1[ord,,drop=FALSE]
  A2 <- A2[ord,,drop=FALSE]
  group <- group[ord]
  group_label <- group_label[ord]
  
  N <- nrow(A1)
  K1 <- ncol(A1)
  K2 <- ncol(A2)
  row_id <- seq_len(N)
  
  group_sizes <- as.numeric(table(group))
  boundaries <- if(length(group_sizes)>1) cumsum(group_sizes)[-length(group_sizes)]+0.5 else numeric(0)
  
  y_breaks <- if(N>=1000){
    unique(c(1,seq(500,N,500),N))
  } else {
    unique(c(1,pretty(c(1,N),5),N))
  }
  
  y_breaks <- y_breaks[y_breaks>=1 & y_breaks<=N]
  
  ## ============================================================
  ## Group means
  ## ============================================================
  
  Xavg <- matrix(NA_real_,N,5)
  A1avg <- matrix(NA_real_,N,K1)
  A2avg <- matrix(NA_real_,N,K2)
  
  for(g in levels(group)){
    
    idx <- which(group==g)
    
    Xavg[idx,] <- matrix(colMeans(X5[idx,,drop=FALSE],na.rm=TRUE),length(idx),5,byrow=TRUE)
    A1avg[idx,] <- matrix(colMeans(A1[idx,,drop=FALSE],na.rm=TRUE),length(idx),K1,byrow=TRUE)
    A2avg[idx,] <- matrix(colMeans(A2[idx,,drop=FALSE],na.rm=TRUE),length(idx),K2,byrow=TRUE)
  }
  
  make_long <- function(M,nms){
    
    z <- as.data.frame(M)
    names(z) <- nms
    z$Row <- row_id
    
    z <- pivot_longer(z,-Row,names_to="Factor",values_to="Value")
    z$Factor <- factor(z$Factor,levels=nms)
    
    z
  }
  
  Xlong <- make_long(Xavg,domains)
  A1long <- make_long(A1avg,as.character(seq_len(K1)))
  A2long <- make_long(A2avg,as.character(seq_len(K2)))
  
  X_lim <- range(Xavg,na.rm=TRUE)
  A1_lim <- range(A1avg,na.rm=TRUE)
  A2_lim <- range(A2avg,na.rm=TRUE)
  
  if(diff(X_lim)==0) X_lim <- X_lim+c(-1e-6,1e-6)
  if(diff(A1_lim)==0) A1_lim <- A1_lim+c(-1e-6,1e-6)
  if(diff(A2_lim)==0) A2_lim <- A2_lim+c(-1e-6,1e-6)
  
  X_breaks <- pretty(X_lim,5)
  X_breaks <- X_breaks[X_breaks>=X_lim[1] & X_breaks<=X_lim[2]]
  
  A1_breaks <- pretty(A1_lim,5)
  A1_breaks <- A1_breaks[A1_breaks>=A1_lim[1] & A1_breaks<=A1_lim[2]]
  
  A2_breaks <- pretty(A2_lim,5)
  A2_breaks <- A2_breaks[A2_breaks>=A2_lim[1] & A2_breaks<=A2_lim[2]]
  
  outcome_df <- data.frame(
    Row=row_id,
    Group=factor(group_label,levels=levels(group))
  )
  
  outcome_colors <- hue_pal()(nlevels(outcome_df$Group))
  names(outcome_colors) <- levels(outcome_df$Group)
  
  legend_order <- rev(levels(outcome_df$Group))
  
  heat_theme <- theme_classic(base_size=16) +
    theme(
      plot.title=element_text(size=18,face="bold",hjust=.5),
      axis.title=element_text(size=16),
      axis.text=element_text(size=13,colour="black"),
      axis.line=element_line(linewidth=.7,colour="black"),
      plot.margin=margin(5,0,5,5)
    )
  
  make_bar <- function(lo,hi,breaks,red=FALSE){
    
    z <- data.frame(x=1,y=seq(lo,hi,length.out=1000))
    z$value <- z$y
    
    p <- ggplot(z,aes(x,y,fill=value))+geom_raster()
    
    if(red){
      p <- p+scale_fill_gradient2(low="white",high="red",midpoint=mean(c(lo,hi)),limits=c(lo,hi))
    } else {
      p <- p+scale_fill_gradient(low="white",high="black",limits=c(lo,hi))
    }
    
    p+
      scale_y_continuous(breaks=breaks,position="right",expand=c(0,0))+
      scale_x_continuous(expand=c(0,0))+
      coord_cartesian(ylim=c(lo,hi),expand=FALSE)+
      guides(fill="none")+
      theme_classic(base_size=11)+
      theme(
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        axis.text.y.left=element_blank(),
        axis.ticks.y.left=element_blank(),
        axis.line.y.left=element_blank(),
        axis.text.y.right=element_text(size=11,colour="black",margin=margin(l=2)),
        plot.margin=margin(5,3,5,0)
      )
  }
  
  make_heat <- function(dat,lims,title,xlab,red=FALSE){
    
    p <- ggplot(dat,aes(Factor,Row,fill=Value))+
      geom_raster()+
      geom_hline(yintercept=boundaries,colour="grey60",linewidth=.4)
    
    if(red){
      p <- p+scale_fill_gradient2(low="white",high="red",midpoint=mean(lims),limits=lims)
    } else {
      p <- p+scale_fill_gradient(low="white",high="black",limits=lims)
    }
    
    p+
      scale_y_continuous(breaks=y_breaks,expand=c(0,0))+
      coord_cartesian(ylim=c(.5,N+.5),expand=FALSE)+
      labs(title=title,x=xlab,y="Observations")+
      guides(fill="none")+
      heat_theme
  }
  
  p1 <- make_heat(Xlong,X_lim,"Big5","Domain",TRUE)
  p2 <- make_heat(A1long,A1_lim,"A1","A1")
  p3 <- make_heat(A2long,A2_lim,"A2","A2")
  
  cb1 <- make_bar(X_lim[1],X_lim[2],X_breaks,TRUE)
  cb2 <- make_bar(A1_lim[1],A1_lim[2],A1_breaks)
  cb3 <- make_bar(A2_lim[1],A2_lim[2],A2_breaks)
  
  p4 <- ggplot(outcome_df,aes(1,Row,fill=Group))+
    geom_raster()+
    geom_hline(yintercept=boundaries,colour="white",linewidth=.6)+
    scale_fill_manual(
      values=outcome_colors,
      breaks=legend_order,
      drop=FALSE,
      guide=guide_legend(title=outcome,title.position="top")
    )+
    scale_y_continuous(expand=c(0,0))+
    scale_x_continuous(expand=c(0,0))+
    coord_cartesian(ylim=c(.5,N+.5),expand=FALSE)+
    labs(title=outcome,x=NULL,y=NULL,fill=outcome)+
    theme_classic(base_size=16)+
    theme(
      plot.title=element_text(size=18,face="bold",hjust=.5),
      axis.text=element_blank(),
      axis.ticks=element_blank(),
      axis.title=element_blank(),
      axis.line=element_line(linewidth=.7,colour="black"),
      legend.position="right",
      legend.title=element_text(size=14,face="bold"),
      legend.text=element_text(size=12),
      legend.key.height=unit(.38,"cm"),
      legend.key.width=unit(.38,"cm"),
      legend.box.spacing=unit(0,"pt"),
      legend.margin=margin(0,0,0,3),
      plot.margin=margin(5,0,5,5)
    )
  
  p1+cb1+p2+cb2+p3+cb3+p4+
    plot_layout(widths=c(1,.10,1.15,.10,.45,.10,.95))
}


## ============================================================
## Data
## ============================================================

data_dir <- "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology/LOPR/Predictive_Analysis"

X <- read.csv(file.path(data_dir,"X_rc.csv"))
A <- read.csv(file.path(data_dir,"A.csv"))
Y <- read.csv(file.path(data_dir,"Y.csv"))
Y_Interpretation <- read.csv(file.path(data_dir,"Y_Interpretation.csv"))


## ============================================================
## HOW TO CUSTOMIZE A NEW OUTCOME
## ============================================================

## First choose an outcome from the new Y.csv.
## This shows all available outcome names.
names(Y)

## Suppose the outcome you want to study is called "NEW_OUTCOME".
##
## Before plotting, inspect:
##
## 1. observed codes and frequencies
## table(Y$NEW_OUTCOME,useNA="ifany")
##
## 2. interpretation/code meanings
## Y_Interpretation[Y_Interpretation$Outcome=="NEW_OUTCOME",]
##
##
## CASE 1:
## Use all original categories exactly as stored.
##
## plot_A1_X_by_outcome_LOOPR(
##   X,A,Y,Y_Interpretation,
##   "NEW_OUTCOME"
## )
##
##
## CASE 2:
## Combine several original response codes into larger groups.
##
## Example:
##
## plot_A1_X_by_outcome_LOOPR(
##   X,A,Y,Y_Interpretation,
##   "NEW_OUTCOME",
##   grouping=TRUE,
##   group_list=list(
##     "Low"=c(1,2),
##     "Medium"=3,
##     "High"=c(4,5)
##   )
## )
##
## The numbers in group_list are the ORIGINAL codes in Y.
##
##
## CASE 3:
## Remove special/non-substantive response codes.
##
## Example: 97 = Don't know, 98 = Refused, 99 = Missing.
##
## plot_A1_X_by_outcome_LOOPR(
##   X,A,Y,Y_Interpretation,
##   "NEW_OUTCOME",
##   Discard=TRUE,
##   discard_codes=c(97,98,99)
## )
##
##
## CASE 4:
## Group categories AND discard special codes.
##
## plot_A1_X_by_outcome_LOOPR(
##   X,A,Y,Y_Interpretation,
##   "NEW_OUTCOME",
##   grouping=TRUE,
##   group_list=list(
##     "Low"=c(1,2),
##     "Medium"=3,
##     "High"=c(4,5)
##   ),
##   Discard=TRUE,
##   discard_codes=c(97,98,99)
## )
##
##
## CASE 5:
## If Y_Interpretation contains very long category names and we
## want shorter labels in the figure, add that outcome inside
## the "maps" list in shorten().
##
## Example:
##
## NEW_OUTCOME=c(
##   "Very long category name one"="Short 1",
##   "Very long category name two"="Short 2"
## )
##
## No other part of the function needs to be changed.


## ============================================================
## Existing examples
## ============================================================

plot_A1_X_by_outcome_LOOPR(
  X,A,Y,Y_Interpretation,
  "educ",
  grouping=TRUE,
  group_list=list(
    "No HS"=c(1,2,3),
    "High school"=4,
    "Some college"=5,
    "2-year"=6,
    "4-year"=7,
    "Post-grad"=c(8,9,10)
  )
)

plot_A1_X_by_outcome_LOOPR(
  X,A,Y,Y_Interpretation,
  "faminc_new",
  grouping=TRUE,
  group_list=list(
    "<40k"=c(1,2,3,4),
    "40k-80k"=c(5,6,7,8),
    "80k-150k"=c(9,10,11),
    "150k+"=c(12,13,14,15,16)
  )
)

plot_A1_X_by_outcome_LOOPR(
  X,A,Y,Y_Interpretation,
  "marstat",
  grouping=TRUE,
  group_list=list(
    "Married"=1,
    "Widowed"=2,
    "Separated/Divorced"=c(3,4),
    "Never married"=5
  )
)

plot_A1_X_by_outcome_LOOPR(
  X,A,Y,Y_Interpretation,
  "pew_prayer"
)