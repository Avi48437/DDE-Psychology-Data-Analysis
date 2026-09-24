library(ggplot2)
library(tidyr)
library(patchwork)
library(grid)
library(scales)

## ============================================================
## IPIP-100: Plot Big Five, DDE factors, and outcome groups
##
## X must contain:
## EXT, EST, AGR, CSN, OPN
##
## A must contain:
## A1_* and A2_* columns
##
## Y contains the fixed IPIP-100 outcomes listed below.
##
## Optional:
## D                = demographic data frame
## demographic_vars = demographic columns to display
##
## grouping=TRUE allows outcome categories to be combined.
## Discard=TRUE removes selected outcome codes.
## ============================================================

plot_A1_X_by_outcome <- function(X,A,Y,outcome,
                                 D=NULL,
                                 demographic_vars=NULL,
                                 grouping=FALSE,
                                 group_list=NULL,
                                 Discard=FALSE,
                                 discard_codes=NULL){
  
  ## ==========================================================
  ## Fixed IPIP-100 outcome labels
  ## ==========================================================
  
  outcome_labels <- list(
    
    educ=c(
      `1`="No HS",
      `2`="High school",
      `3`="Some college",
      `4`="2-year",
      `5`="4-year",
      `6`="Post-grad"
    ),
    
    faminc_new=c(
      `1`="<10k",
      `2`="10-20k",
      `3`="20-30k",
      `4`="30-40k",
      `5`="40-50k",
      `6`="50-60k",
      `7`="60-70k",
      `8`="70-80k",
      `9`="80-100k",
      `10`="100-120k",
      `11`="120-150k",
      `12`="150-200k",
      `13`="200-250k",
      `14`="250-350k",
      `15`="350-500k",
      `16`="500k+",
      `97`="Prefer not to say"
    ),
    
    marstat=c(
      `1`="Married",
      `2`="Separated",
      `3`="Divorced",
      `4`="Widowed",
      `5`="Never married",
      `6`="Domestic/civil"
    ),
    
    pew_prayer=c(
      `1`="Several/day",
      `2`="Once/day",
      `3`="Few/week",
      `4`="Once/week",
      `5`="Few/month",
      `6`="Seldom",
      `7`="Never",
      `8`="Don't know"
    ),
    
    votereg=c(
      `1`="Yes",
      `2`="No",
      `3`="Don't know"
    ),
    
    pid3=c(
      `1`="Democrat",
      `2`="Republican",
      `3`="Independent",
      `4`="Other",
      `5`="Not sure"
    ),
    
    pid7=c(
      `1`="Strong Democrat",
      `2`="Democrat",
      `3`="Lean Democrat",
      `4`="Independent",
      `5`="Lean Republican",
      `6`="Republican",
      `7`="Strong Republican",
      `8`="Not sure"
    ),
    
    ideo5=c(
      `1`="Very liberal",
      `2`="Liberal",
      `3`="Moderate",
      `4`="Conservative",
      `5`="Very conservative",
      `6`="Not sure"
    )
  )
  
  ## ==========================================================
  ## Check predictors and outcome
  ## ==========================================================
  
  domains <- c("EXT","EST","AGR","CSN","OPN")
  
  if(!all(domains %in% names(X)))
    stop("X must contain EXT, EST, AGR, CSN, and OPN.")
  
  if(!outcome %in% names(Y))
    stop(paste0("Outcome '",outcome,"' was not found in Y."))
  
  if(nrow(X)!=nrow(A) || nrow(X)!=nrow(Y))
    stop("X, A, and Y must have the same number of rows.")
  
  if(grouping && is.null(group_list))
    stop("grouping=TRUE requires group_list.")
  
  if(Discard && is.null(discard_codes))
    stop("Discard=TRUE requires discard_codes.")
  
  X5 <- as.matrix(X[,domains,drop=FALSE])
  A1 <- as.matrix(A[,grep("^A1_",names(A),value=TRUE),drop=FALSE])
  A2 <- as.matrix(A[,grep("^A2_",names(A),value=TRUE),drop=FALSE])
  
  y <- Y[[outcome]]
  labels <- outcome_labels[[outcome]]
  
  ## ==========================================================
  ## Optional demographics
  ##
  ## Each requested demographic variable becomes one column.
  ## Numeric variables are rescaled to [0,1].
  ## Factors/characters are first converted to numeric codes.
  ## ==========================================================
  
  use_demo <- !is.null(D) && !is.null(demographic_vars)
  
  if(use_demo){
    
    if(nrow(D)!=nrow(X))
      stop("D must have the same number of rows as X.")
    
    if(!all(demographic_vars %in% names(D)))
      stop("Some demographic_vars were not found in D.")
    
    Demo <- D[,demographic_vars,drop=FALSE]
    
    for(v in demographic_vars){
      
      z <- Demo[[v]]
      
      if(is.factor(z) || is.character(z))
        z <- as.numeric(factor(z))
      
      z <- as.numeric(z)
      
      lo <- min(z,na.rm=TRUE)
      hi <- max(z,na.rm=TRUE)
      
      if(is.finite(lo) && is.finite(hi) && hi>lo)
        z <- (z-lo)/(hi-lo)
      
      Demo[[v]] <- z
    }
    
    Demo <- as.matrix(Demo)
    
  } else {
    
    Demo <- NULL
  }
  
  ## ==========================================================
  ## Remove missing outcomes
  ## ==========================================================
  
  keep <- !is.na(y)
  
  X5 <- X5[keep,,drop=FALSE]
  A1 <- A1[keep,,drop=FALSE]
  A2 <- A2[keep,,drop=FALSE]
  
  if(use_demo)
    Demo <- Demo[keep,,drop=FALSE]
  
  y <- y[keep]
  
  ## ==========================================================
  ## Optionally discard selected outcome codes
  ##
  ## Example:
  ## Discard=TRUE, discard_codes=97
  ## removes income code 97 ("Prefer not to say").
  ## ==========================================================
  
  if(Discard){
    
    keep <- !(as.character(y) %in% as.character(discard_codes))
    
    X5 <- X5[keep,,drop=FALSE]
    A1 <- A1[keep,,drop=FALSE]
    A2 <- A2[keep,,drop=FALSE]
    
    if(use_demo)
      Demo <- Demo[keep,,drop=FALSE]
    
    y <- y[keep]
  }
  
  ## ==========================================================
  ## Original outcome-category ordering and labels
  ## ==========================================================
  
  if(is.null(labels)){
    
    category_order <- as.character(sort(unique(y)))
    label_map <- setNames(category_order,category_order)
    
  } else {
    
    category_order <- names(labels)
    label_map <- setNames(unname(labels),names(labels))
  }
  
  y_chr <- as.character(y)
  
  category_order <- category_order[
    category_order %in% unique(y_chr)
  ]
  
  ## ==========================================================
  ## Optional grouping
  ##
  ## Example:
  ##
  ## group_list=list(
  ##   "Low"=c(1,2),
  ##   "Middle"=3,
  ##   "High"=c(4,5)
  ## )
  ##
  ## Any category not supplied remains unchanged.
  ## ==========================================================
  
  if(grouping){
    
    group_codes <- lapply(group_list,as.character)
    
    all_codes <- unlist(group_codes,use.names=FALSE)
    duplicate_codes <- unique(all_codes[duplicated(all_codes)])
    
    if(length(duplicate_codes)>0)
      stop(
        paste0(
          "These outcome codes appear in more than one group: ",
          paste(duplicate_codes,collapse=", ")
        )
      )
    
    group_key <- setNames(category_order,category_order)
    group_display <- setNames(unname(label_map[category_order]),category_order)
    
    supplied_names <- names(group_list)
    
    if(is.null(supplied_names))
      supplied_names <- rep("",length(group_list))
    
    for(i in seq_along(group_codes)){
      
      codes_i <- intersect(group_codes[[i]],category_order)
      
      if(length(codes_i)==0)
        next
      
      new_key <- paste0(".GROUP_",i)
      
      if(nzchar(supplied_names[i])){
        new_label <- supplied_names[i]
      } else {
        new_label <- paste(unname(label_map[codes_i]),collapse="/")
      }
      
      group_key[codes_i] <- new_key
      group_display[new_key] <- new_label
    }
    
    group_levels <- unique(unname(group_key[category_order]))
    group_values <- unname(group_key[y_chr])
    
    group <- factor(group_values,levels=group_levels)
    
    group_labels <- vapply(
      group_levels,
      function(g){
        if(g %in% names(group_display))
          return(unname(group_display[g]))
        g
      },
      character(1)
    )
    
    names(group_labels) <- group_levels
    
  } else {
    
    group <- factor(y_chr,levels=category_order)
    
    if(is.null(labels)){
      group_labels <- setNames(category_order,category_order)
    } else {
      group_labels <- setNames(unname(labels[category_order]),category_order)
    }
  }
  
  ## ==========================================================
  ## Remove values not represented by the coding
  ## ==========================================================
  
  keep <- !is.na(group)
  
  X5 <- X5[keep,,drop=FALSE]
  A1 <- A1[keep,,drop=FALSE]
  A2 <- A2[keep,,drop=FALSE]
  
  if(use_demo)
    Demo <- Demo[keep,,drop=FALSE]
  
  group <- droplevels(group[keep])
  
  group_label <- unname(
    group_labels[as.character(group)]
  )
  
  ## ==========================================================
  ## Sort respondents by outcome group
  ## ==========================================================
  
  ord <- order(group)
  
  X5 <- X5[ord,,drop=FALSE]
  A1 <- A1[ord,,drop=FALSE]
  A2 <- A2[ord,,drop=FALSE]
  
  if(use_demo)
    Demo <- Demo[ord,,drop=FALSE]
  
  group <- group[ord]
  group_label <- group_label[ord]
  
  N <- nrow(A1)
  K1 <- ncol(A1)
  K2 <- ncol(A2)
  
  row_id <- seq_len(N)
  
  group_sizes <- as.numeric(table(group))
  
  boundaries <- if(length(group_sizes)>1){
    cumsum(group_sizes)[-length(group_sizes)]+0.5
  } else {
    numeric(0)
  }
  
  y_breaks <- if(N>=1000){
    unique(c(1,seq(500,N,500),N))
  } else {
    unique(c(1,pretty(c(1,N),5),N))
  }
  
  y_breaks <- y_breaks[y_breaks>=1 & y_breaks<=N]
  
  ## ==========================================================
  ## Replace every respondent by the mean of their outcome group
  ##
  ## This produces the block structure shown in the heatmaps.
  ## ==========================================================
  
  group_average <- function(M){
    
    out <- matrix(NA_real_,nrow(M),ncol(M))
    
    for(g in levels(group)){
      
      idx <- which(group==g)
      
      out[idx,] <- matrix(
        colMeans(M[idx,,drop=FALSE],na.rm=TRUE),
        length(idx),
        ncol(M),
        byrow=TRUE
      )
    }
    
    colnames(out) <- colnames(M)
    
    out
  }
  
  Xavg <- group_average(X5)
  A1avg <- group_average(A1)
  A2avg <- group_average(A2)
  
  if(use_demo)
    Davg <- group_average(Demo)
  
  ## ==========================================================
  ## Convert matrices to long format for ggplot
  ## ==========================================================
  
  make_long <- function(M,nms){
    
    z <- as.data.frame(M)
    names(z) <- nms
    z$Row <- row_id
    
    z <- pivot_longer(
      z,
      -Row,
      names_to="Factor",
      values_to="Value"
    )
    
    z$Factor <- factor(z$Factor,levels=nms)
    
    z
  }
  
  Xlong <- make_long(Xavg,domains)
  A1long <- make_long(A1avg,as.character(seq_len(K1)))
  A2long <- make_long(A2avg,as.character(seq_len(K2)))
  
  if(use_demo)
    Dlong <- make_long(Davg,demographic_vars)
  
  ## ==========================================================
  ## Color limits
  ## ==========================================================
  
  X_lim <- range(Xavg,na.rm=TRUE)
  A1_lim <- range(A1avg,na.rm=TRUE)
  A2_lim <- range(A2avg,na.rm=TRUE)
  
  if(diff(X_lim)==0) X_lim <- X_lim+c(-1e-6,1e-6)
  if(diff(A1_lim)==0) A1_lim <- A1_lim+c(-1e-6,1e-6)
  if(diff(A2_lim)==0) A2_lim <- A2_lim+c(-1e-6,1e-6)
  
  X_mid <- mean(X_lim)
  
  X_breaks <- pretty(X_lim,5)
  X_breaks <- X_breaks[X_breaks>=X_lim[1] & X_breaks<=X_lim[2]]
  
  A1_breaks <- pretty(A1_lim,5)
  A1_breaks <- A1_breaks[A1_breaks>=A1_lim[1] & A1_breaks<=A1_lim[2]]
  
  A2_breaks <- pretty(A2_lim,5)
  A2_breaks <- A2_breaks[A2_breaks>=A2_lim[1] & A2_breaks<=A2_lim[2]]
  
  if(use_demo){
    
    D_lim <- range(Davg,na.rm=TRUE)
    
    if(diff(D_lim)==0)
      D_lim <- D_lim+c(-1e-6,1e-6)
    
    D_breaks <- pretty(D_lim,5)
    D_breaks <- D_breaks[D_breaks>=D_lim[1] & D_breaks<=D_lim[2]]
  }
  
  ## ==========================================================
  ## Plot helpers
  ## ==========================================================
  
  heat_theme <- theme_classic(base_size=16) +
    theme(
      plot.title=element_text(size=18,face="bold",hjust=.5),
      axis.title=element_text(size=16),
      axis.text=element_text(size=13,colour="black"),
      axis.line=element_line(linewidth=.7,colour="black"),
      plot.margin=margin(5,0,5,5)
    )
  
  make_bar <- function(lo,hi,breaks,red=FALSE,mid=NULL){
    
    z <- data.frame(
      x=1,
      y=seq(lo,hi,length.out=1000)
    )
    
    z$value <- z$y
    
    p <- ggplot(z,aes(x,y,fill=value)) +
      geom_raster()
    
    if(red){
      
      if(is.null(mid))
        mid <- mean(c(lo,hi))
      
      p <- p +
        scale_fill_gradient2(
          low="white",
          high="red",
          midpoint=mid,
          limits=c(lo,hi)
        )
      
    } else {
      
      p <- p +
        scale_fill_gradient(
          low="white",
          high="black",
          limits=c(lo,hi)
        )
    }
    
    p +
      scale_y_continuous(
        breaks=breaks,
        position="right",
        expand=c(0,0)
      ) +
      scale_x_continuous(expand=c(0,0)) +
      coord_cartesian(
        ylim=c(lo,hi),
        expand=FALSE
      ) +
      guides(fill="none") +
      theme_classic(base_size=11) +
      theme(
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        axis.text.y.left=element_blank(),
        axis.ticks.y.left=element_blank(),
        axis.line.y.left=element_blank(),
        axis.text.y.right=element_text(
          size=11,
          colour="black",
          margin=margin(l=2)
        ),
        plot.margin=margin(5,3,5,0)
      )
  }
  
  make_heat <- function(dat,lims,title,xlab,red=FALSE,mid=NULL){
    
    p <- ggplot(dat,aes(Factor,Row,fill=Value)) +
      geom_raster() +
      geom_hline(
        yintercept=boundaries,
        colour="grey60",
        linewidth=.4
      )
    
    if(red){
      
      if(is.null(mid))
        mid <- mean(lims)
      
      p <- p +
        scale_fill_gradient2(
          low="white",
          high="red",
          midpoint=mid,
          limits=lims
        )
      
    } else {
      
      p <- p +
        scale_fill_gradient(
          low="white",
          high="black",
          limits=lims
        )
    }
    
    p +
      scale_y_continuous(
        breaks=y_breaks,
        expand=c(0,0)
      ) +
      coord_cartesian(
        ylim=c(.5,N+.5),
        expand=FALSE
      ) +
      labs(
        title=title,
        x=xlab,
        y="Observations"
      ) +
      guides(fill="none") +
      heat_theme
  }
  
  ## ==========================================================
  ## Big Five, A1 and A2 panels
  ## ==========================================================
  
  p1 <- make_heat(Xlong,X_lim,"Big5","Domain",TRUE,X_mid)
  p2 <- make_heat(A1long,A1_lim,"A1","A1")
  p3 <- make_heat(A2long,A2_lim,"A2","A2")
  
  cb1 <- make_bar(X_lim[1],X_lim[2],X_breaks,TRUE,X_mid)
  cb2 <- make_bar(A1_lim[1],A1_lim[2],A1_breaks)
  cb3 <- make_bar(A2_lim[1],A2_lim[2],A2_breaks)
  
  ## ==========================================================
  ## Optional demographics panel
  ## ==========================================================
  
  if(use_demo){
    
    p0 <- make_heat(
      Dlong,
      D_lim,
      "Demographics",
      NULL
    ) +
      theme(
        axis.text.x=element_text(
          angle=45,
          hjust=1
        )
      )
    
    cb0 <- make_bar(
      D_lim[1],
      D_lim[2],
      D_breaks
    )
  }
  
  ## ==========================================================
  ## Outcome panel
  ## ==========================================================
  
  outcome_df <- data.frame(
    Row=row_id,
    Group=factor(
      group_label,
      levels=unique(group_label)
    )
  )
  
  outcome_colors <- hue_pal()(
    nlevels(outcome_df$Group)
  )
  
  names(outcome_colors) <- levels(outcome_df$Group)
  
  legend_order <- rev(levels(outcome_df$Group))
  
  p4 <- ggplot(
    outcome_df,
    aes(1,Row,fill=Group)
  ) +
    geom_raster() +
    geom_hline(
      yintercept=boundaries,
      colour="white",
      linewidth=.6
    ) +
    scale_fill_manual(
      values=outcome_colors,
      breaks=legend_order,
      drop=FALSE,
      guide=guide_legend(
        title=outcome,
        title.position="top"
      )
    ) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    coord_cartesian(
      ylim=c(.5,N+.5),
      expand=FALSE
    ) +
    labs(
      title=outcome,
      x=NULL,
      y=NULL,
      fill=outcome
    ) +
    theme_classic(base_size=16) +
    theme(
      plot.title=element_text(size=18,face="bold",hjust=.5),
      axis.text=element_blank(),
      axis.ticks=element_blank(),
      axis.title=element_blank(),
      axis.line=element_line(linewidth=.7,colour="black"),
      legend.position="right",
      legend.title=element_text(size=14,face="bold"),
      legend.text=element_text(size=13),
      legend.key.height=unit(.42,"cm"),
      legend.key.width=unit(.42,"cm"),
      legend.box.spacing=unit(0,"pt"),
      legend.margin=margin(0,0,0,3),
      plot.margin=margin(5,0,5,5)
    )
  
  ## ==========================================================
  ## Final layout
  ## ==========================================================
  
  if(use_demo){
    
    demo_width <- max(.28,.18*ncol(Demo))
    
    p0 + cb0 +
      p1 + cb1 +
      p2 + cb2 +
      p3 + cb3 +
      p4 +
      plot_layout(
        widths=c(
          demo_width,.07,
          1,.07,
          1.15,.07,
          .45,.07,
          .95
        )
      )
    
  } else {
    
    p1 + cb1 +
      p2 + cb2 +
      p3 + cb3 +
      p4 +
      plot_layout(
        widths=c(
          1,.10,
          1.15,.10,
          .45,.10,
          .95
        )
      )
  }
}


## ============================================================
## USAGE
## ============================================================

## Education
plot_A1_X_by_outcome(
  X,A,Y,
  "educ"
)

## Marital status
plot_A1_X_by_outcome(
  X,A,Y,
  "marstat"
)

## Prayer
plot_A1_X_by_outcome(
  X,A,Y,
  "pew_prayer"
)

## Voting registration
plot_A1_X_by_outcome(
  X,A,Y,
  "votereg"
)

## Family income
## Code 97 = "Prefer not to say", so discard it.
plot_A1_X_by_outcome(
  X,A,Y,
  "faminc_new",
  Discard=TRUE,
  discard_codes=97
)

## Political ideology
plot_A1_X_by_outcome(
  X,A,Y,
  "ideo5"
)

## Three-category party identification
plot_A1_X_by_outcome(
  X,A,Y,
  "pid3"
)

## Seven-category party identification
plot_A1_X_by_outcome(
  X,A,Y,
  "pid7"
)


## ============================================================
## OPTIONAL GROUPING EXAMPLE
##
## Use this only when we deliberately want to merge categories.
## Original outcome codes are supplied inside group_list.
## ============================================================

# plot_A1_X_by_outcome(
#   X,A,Y,
#   "marstat",
#   grouping=TRUE,
#   group_list=list(
#     "Married"=1,
#     "Separated/Divorced"=c(2,3),
#     "Widowed"=4,
#     "Never/Domestic"=c(5,6)
#   )
# )


## ============================================================
## OPTIONAL DEMOGRAPHICS EXAMPLE
##
## D must have exactly the same respondent ordering as X, A, Y.
## Add whatever demographic columns we want to display.
## ============================================================

# plot_A1_X_by_outcome(
#   X,A,Y,
#   "educ",
#   D=D_IPIP98,
#   demographic_vars=c("age","gender","race")
# )


## ============================================================
## GROUPING + DISCARDING + DEMOGRAPHICS
## ============================================================

# plot_A1_X_by_outcome(
#   X,A,Y,
#   "faminc_new",
#   D=D_IPIP98,
#   demographic_vars=c("age","gender","race"),
#   grouping=TRUE,
#   group_list=list(
#     "<40k"=c(1,2,3,4),
#     "40k-80k"=c(5,6,7,8),
#     "80k-150k"=c(9,10,11),
#     "150k+"=c(12,13,14,15,16)
#   ),
#   Discard=TRUE,
#   discard_codes=97
# )