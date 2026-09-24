library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(grid)
library(scales)
## ============================================================
## Outcome labels
## ============================================================

outcome_labels <- list(
  
  educ=c(
    `1`="No HS",`2`="High school",`3`="Some college",
    `4`="2-year",`5`="4-year",`6`="Post-grad"
  ),
  
  faminc_new=c(
    `1`="<10k",`2`="10-20k",`3`="20-30k",`4`="30-40k",
    `5`="40-50k",`6`="50-60k",`7`="60-70k",`8`="70-80k",
    `9`="80-100k",`10`="100-120k",`11`="120-150k",
    `12`="150-200k",`13`="200-250k",`14`="250-350k",
    `15`="350-500k",`16`="500k+",`97`="Prefer not to say"
  ),
  
  marstat=c(
    `1`="Married",`2`="Separated",`3`="Divorced",
    `4`="Widowed",`5`="Never married",`6`="Domestic/civil"
  ),
  
  pew_prayer=c(
    `1`="Several/day",`2`="Once/day",`3`="Few/week",
    `4`="Once/week",`5`="Few/month",`6`="Seldom",
    `7`="Never",`8`="Don't know"
  ),
  
  votereg=c(
    `1`="Yes",`2`="No",`3`="Don't know"
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

plot_A_means <- function(outcome, drop_codes = NULL) {
  
  y <- Y[[outcome]]
  
  keep <- !is.na(y)
  
  if (!is.null(drop_codes))
    keep <- keep & !(y %in% drop_codes)
  
  A_use <- A[keep, , drop = FALSE]
  y_use <- y[keep]
  
  avg <- aggregate(A_use, by = list(Group = y_use), FUN = mean)
  
  ## Apply descriptive labels when available
  if (outcome %in% names(outcome_labels)) {
    
    lab <- outcome_labels[[outcome]]
    
    avg$GroupLabel <- unname(lab[as.character(avg$Group)])
    
  } else {
    
    avg$GroupLabel <- as.character(avg$Group)
  }
  
  avg_long <- avg |>
    select(Group, GroupLabel, everything()) |>
    pivot_longer(
      cols = starts_with("A"),
      names_to = "LatentFactor",
      values_to = "Mean"
    )
  
  ggplot(avg_long,
         aes(x = LatentFactor,
             y = factor(GroupLabel, levels = rev(unique(GroupLabel))),
             fill = Mean)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      low = "white",
      high = "black",
      limits = c(0,1)
    ) +
    labs(
      title = paste("Mean DDE latent activation by", outcome),
      x = "DDE latent factor",
      y = NULL,
      fill = "Mean"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}



plot_A1_by_outcome <- function(A,Y,outcome,drop_codes=NULL){
  
  outcome_labels <- list(
    educ=c(`1`="No HS",`2`="High school",`3`="Some college",`4`="2-year",`5`="4-year",`6`="Post-grad"),
    faminc_new=c(`1`="<10k",`2`="10-20k",`3`="20-30k",`4`="30-40k",`5`="40-50k",`6`="50-60k",
                 `7`="60-70k",`8`="70-80k",`9`="80-100k",`10`="100-120k",`11`="120-150k",
                 `12`="150-200k",`13`="200-250k",`14`="250-350k",`15`="350-500k",
                 `16`="500k+",`97`="Prefer not to say"),
    marstat=c(`1`="Married",`2`="Separated",`3`="Divorced",`4`="Widowed",`5`="Never married",`6`="Domestic/civil"),
    pew_prayer=c(`1`="Several/day",`2`="Once/day",`3`="Few/week",`4`="Once/week",`5`="Few/month",`6`="Seldom",`7`="Never",`8`="Don't know"),
    votereg=c(`1`="Yes",`2`="No",`3`="Don't know")
  )
  
  A1_names <- grep("^A1_",names(A),value=TRUE)
  A1 <- as.matrix(A[,A1_names,drop=FALSE])
  y <- Y[[outcome]]
  
  keep <- !is.na(y)
  if(!is.null(drop_codes)) keep <- keep & !(y %in% drop_codes)
  A1 <- A1[keep,,drop=FALSE]
  y <- y[keep]
  
  labels <- outcome_labels[[outcome]]
  category_order <- if(!is.null(labels)) as.numeric(names(labels)) else sort(unique(y))
  
  group <- factor(y,levels=category_order)
  keep <- !is.na(group)
  A1 <- A1[keep,,drop=FALSE]
  group <- droplevels(group[keep])
  
  if (is.null(labels)) {
    
    group_label <- as.character(group)
    
  } else {
    
    group_label <- unname(
      labels[as.character(group)]
    )
    
    ## If a label has not yet been defined,
    ## show the numeric code instead of dropping it
    missing_label <- is.na(group_label)
    
    group_label[missing_label] <-
      as.character(group[missing_label])
  }
  
  ord <- order(group)
  A1 <- A1[ord,,drop=FALSE]
  group <- group[ord]
  group_label <- group_label[ord]
  
  N <- nrow(A1)
  K <- ncol(A1)
  row_id <- 1:N
  
  group_sizes <- as.numeric(table(group))
  boundaries <- if(length(group_sizes)>1) cumsum(group_sizes)[-length(group_sizes)]+0.5 else numeric(0)
  
  y_breaks <- if(N>=1000) unique(c(1,seq(500,N,500),N)) else unique(c(1,pretty(c(1,N),5),N))
  y_breaks <- y_breaks[y_breaks>=1 & y_breaks<=N]
  
  noisy_long <- as.data.frame(A1) |>
    setNames(as.character(1:K)) |>
    mutate(Row=row_id) |>
    pivot_longer(-Row,names_to="Factor",values_to="Value")
  noisy_long$Factor <- factor(noisy_long$Factor,levels=as.character(1:K))
  
  A1_avg <- matrix(NA_real_,N,K)
  for(g in levels(group)){
    idx <- which(group==g)
    mu <- colMeans(A1[idx,,drop=FALSE],na.rm=TRUE)
    A1_avg[idx,] <- matrix(mu,length(idx),K,byrow=TRUE)
  }
  
  avg_long <- as.data.frame(A1_avg) |>
    setNames(as.character(1:K)) |>
    mutate(Row=row_id) |>
    pivot_longer(-Row,names_to="Factor",values_to="Value")
  avg_long$Factor <- factor(avg_long$Factor,levels=as.character(1:K))
  
  outcome_levels <- unique(group_label)
  outcome_df <- data.frame(Row=row_id,Group=factor(group_label,levels=outcome_levels))
  outcome_colors <- hue_pal()(length(outcome_levels))
  names(outcome_colors) <- outcome_levels
  legend_order <- rev(outcome_levels)
  
  heat_theme <- theme_classic(base_size=16) +
    theme(plot.title=element_text(size=18,face="bold",hjust=0.5),
          axis.title=element_text(size=16),
          axis.text=element_text(size=13,colour="black"),
          axis.line=element_line(linewidth=0.7,colour="black"),
          plot.margin=margin(5,0,5,5))
  
  p1 <- ggplot(noisy_long,aes(Factor,Row,fill=Value)) +
    geom_raster() +
    geom_hline(yintercept=boundaries,colour="grey75",linewidth=0.25) +
    scale_fill_gradient(low="white",high="black",limits=c(0,1)) +
    scale_y_continuous(breaks=y_breaks,expand=c(0,0)) +
    coord_cartesian(ylim=c(0.5,N+0.5),expand=FALSE) +
    labs(title=expression("Noisy "*A^{(1)}),x=expression(A[k]^{(1)}),y="Observations") +
    guides(fill="none") +
    heat_theme
  
  p2 <- ggplot(avg_long,aes(Factor,Row,fill=Value)) +
    geom_raster() +
    geom_hline(yintercept=boundaries,colour="grey70",linewidth=0.35) +
    scale_fill_gradient(low="white",high="black",limits=c(0,1)) +
    scale_y_continuous(breaks=y_breaks,expand=c(0,0)) +
    coord_cartesian(ylim=c(0.5,N+0.5),expand=FALSE) +
    labs(title=bquote(Avg(A[k]^{(1)}~"|"~.(outcome))),x=expression(A[k]^{(1)}),y="Observations") +
    guides(fill="none") +
    heat_theme
  
  make_gray_bar <- function(){
    df <- data.frame(x=1,y=seq(0,1,length.out=1000)); df$value <- df$y
    ggplot(df,aes(x,y,fill=value)) +
      geom_raster() +
      scale_fill_gradient(low="white",high="black",limits=c(0,1)) +
      scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0","0.25","0.50","0.75","1"),position="right",expand=c(0,0)) +
      coord_cartesian(ylim=c(0,1),expand=FALSE) +
      scale_x_continuous(expand=c(0,0)) +
      guides(fill="none") +
      theme_classic(base_size=12) +
      theme(axis.title=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.line.x=element_blank(),
            axis.text.y.left=element_blank(),
            axis.ticks.y.left=element_blank(),
            axis.line.y.left=element_blank(),
            axis.text.y.right=element_text(size=11,colour="black",margin=margin(l=3)),
            axis.ticks.y.right=element_line(colour="black"),
            axis.line.y.right=element_line(colour="black",linewidth=0.5),
            plot.margin=margin(5,3,5,0))
  }
  
  cb1 <- make_gray_bar()
  cb2 <- make_gray_bar()
  
  p3 <- ggplot(outcome_df,aes(1,Row,fill=Group)) +
    geom_raster() +
    geom_hline(yintercept=boundaries,colour="white",linewidth=0.5) +
    scale_fill_manual(values=outcome_colors,breaks=legend_order,drop=FALSE,
                      guide=guide_legend(title=outcome,title.position="top",label.position="right")) +
    scale_y_continuous(breaks=y_breaks,expand=c(0,0)) +
    coord_cartesian(ylim=c(0.5,N+0.5),expand=FALSE) +
    scale_x_continuous(expand=c(0,0)) +
    labs(title=outcome,x=NULL,y=NULL,fill=outcome) +
    theme_classic(base_size=16) +
    theme(plot.title=element_text(size=18,face="bold",hjust=0.5),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          axis.line=element_line(linewidth=0.7,colour="black"),
          legend.position="right",
          legend.box.spacing=unit(0,"pt"),
          legend.margin=margin(0,0,0,2),
          legend.title=element_text(size=14,face="bold"),
          legend.text=element_text(size=13),
          legend.key.height=unit(0.45,"cm"),
          legend.key.width=unit(0.45,"cm"),
          plot.margin=margin(5,0,5,5))
  
  p1 + cb1 + p2 + cb2 + p3 + plot_layout(widths=c(1,0.10,1,0.10,0.95))
}


# Following functions are added by Youlin

## ============================================================
## BFI-2 domain means for LOOPR
## ============================================================

get_loopr_domain_means <- function(X) {
  
  X <- as.data.frame(X)
  
  stopifnot(ncol(X) == 60)
  
  ## Reverse-keyed BFI-2 items
  reverse_items <- c(
    11,16,26,31,36,51,       # Extraversion
    12,17,22,37,42,47,       # Agreeableness
    3,8,23,28,48,58,         # Conscientiousness
    4,9,24,29,44,49,         # Negative Emotionality
    5,25,30,45,50,55         # Open-Mindedness
  )
  
  X_rc <- X
  
  X_rc[, reverse_items] <-
    6 - X_rc[, reverse_items]
  
  
  domain_items <- list(
    EXT = c(1,6,11,16,21,26,31,36,41,46,51,56),
    AGR = c(2,7,12,17,22,27,32,37,42,47,52,57),
    CSN = c(3,8,13,18,23,28,33,38,43,48,53,58),
    NEG = c(4,9,14,19,24,29,34,39,44,49,54,59),
    OPN = c(5,10,15,20,25,30,35,40,45,50,55,60)
  )
  
  
  X5 <- sapply(
    domain_items,
    function(idx) {
      rowMeans(
        X_rc[, idx, drop = FALSE],
        na.rm = TRUE
      )
    }
  )
  
  X5 <- as.data.frame(X5)
  
  return(X5)
}

plot_A1_X_by_outcome <- function(
    X,
    A,
    Y,
    outcome,
    drop_codes = NULL
) {
  
  ## ----------------------------------------------------------
  ## LOOPR outcome labels
  ## ----------------------------------------------------------
  
  outcome_labels <- list(
    
    education_level = c(
      `1` = "Less than high school",
      `2` = "High school",
      `3` = "Some college",
      `4` = "Associate degree",
      `5` = "Bachelor's degree",
      `6` = "Graduate degree"
    ),
    
    marriage_status = c(
      `1` = "Married",
      `2` = "Separated",
      `3` = "Divorced",
      `4` = "Widowed",
      `5` = "Never married",
      `6` = "Domestic/civil"
    )
    
    ## Add income / prayer labels here once we want
    ## exact wording from LOOPR
  )
  
  
  ## ----------------------------------------------------------
  ## Basic checks
  ## ----------------------------------------------------------
  
  stopifnot(
    nrow(X) == nrow(A),
    nrow(X) == nrow(Y),
    outcome %in% names(Y)
  )
  
  
  ## ----------------------------------------------------------
  ## Big Five means
  ## ----------------------------------------------------------
  
  X5 <- get_loopr_domain_means(X)
  
  X5 <- as.matrix(X5)
  
  
  ## ----------------------------------------------------------
  ## Extract A1 and A2
  ## ----------------------------------------------------------
  
  A1_names <- grep(
    "^A1_",
    names(A),
    value = TRUE
  )
  
  A2_names <- grep(
    "^A2_",
    names(A),
    value = TRUE
  )
  
  
  A1 <- as.matrix(
    A[, A1_names, drop = FALSE]
  )
  
  A2 <- as.matrix(
    A[, A2_names, drop = FALSE]
  )
  
  
  y <- Y[[outcome]]
  
  
  ## ----------------------------------------------------------
  ## Remove missing outcome / unwanted codes
  ## ----------------------------------------------------------
  
  keep <- !is.na(y)
  
  if (!is.null(drop_codes)) {
    keep <- keep & !(y %in% drop_codes)
  }
  
  
  X5 <- X5[keep, , drop = FALSE]
  A1 <- A1[keep, , drop = FALSE]
  A2 <- A2[keep, , drop = FALSE]
  y <- y[keep]
  
  
  ## ----------------------------------------------------------
  ## Define outcome ordering
  ## ----------------------------------------------------------
  
  labels <- outcome_labels[[outcome]]
  
  ## Always determine categories from the observed data
  category_order <- sort(unique(y))
  
  ## Remove NA just in case
  category_order <- category_order[!is.na(category_order)]
  
  group <- factor(
    y,
    levels = category_order
  )
  
  keep <- !is.na(group)
  
  
  X5 <- X5[keep, , drop = FALSE]
  A1 <- A1[keep, , drop = FALSE]
  A2 <- A2[keep, , drop = FALSE]
  
  group <- droplevels(
    group[keep]
  )
  
  
  if (is.null(labels)) {
    
    group_label <- as.character(group)
    
  } else {
    
    group_label <- unname(
      labels[as.character(group)]
    )
    
    ## If a category does not yet have a label,
    ## use its numeric code instead
    missing_label <- is.na(group_label)
    
    group_label[missing_label] <-
      as.character(group[missing_label])
  }
  
  
  ## ----------------------------------------------------------
  ## Order observations by outcome category
  ## ----------------------------------------------------------
  
  ord <- order(group)
  
  X5 <- X5[ord, , drop = FALSE]
  A1 <- A1[ord, , drop = FALSE]
  A2 <- A2[ord, , drop = FALSE]
  
  group <- group[ord]
  group_label <- group_label[ord]
  
  
  N <- nrow(A1)
  K1 <- ncol(A1)
  K2 <- ncol(A2)
  
  row_id <- seq_len(N)
  
  
  group_sizes <- as.numeric(
    table(group)
  )
  
  boundaries <- if (length(group_sizes) > 1) {
    
    cumsum(group_sizes)[
      -length(group_sizes)
    ] + 0.5
    
  } else {
    
    numeric(0)
  }
  
  
  y_breaks <- if (N >= 1000) {
    
    unique(
      c(
        1,
        seq(500, N, 500),
        N
      )
    )
    
  } else {
    
    unique(
      c(
        1,
        pretty(c(1, N), 5),
        N
      )
    )
  }
  
  
  y_breaks <- y_breaks[
    y_breaks >= 1 &
      y_breaks <= N
  ]
  
  
  ## ----------------------------------------------------------
  ## Group averages
  ## ----------------------------------------------------------
  
  Xavg <- matrix(
    NA_real_,
    N,
    ncol(X5)
  )
  
  A1avg <- matrix(
    NA_real_,
    N,
    K1
  )
  
  A2avg <- matrix(
    NA_real_,
    N,
    K2
  )
  
  
  for (g in levels(group)) {
    
    idx <- which(
      group == g
    )
    
    
    Xavg[idx, ] <- matrix(
      colMeans(
        X5[idx, , drop = FALSE],
        na.rm = TRUE
      ),
      length(idx),
      ncol(X5),
      byrow = TRUE
    )
    
    
    A1avg[idx, ] <- matrix(
      colMeans(
        A1[idx, , drop = FALSE],
        na.rm = TRUE
      ),
      length(idx),
      K1,
      byrow = TRUE
    )
    
    
    A2avg[idx, ] <- matrix(
      colMeans(
        A2[idx, , drop = FALSE],
        na.rm = TRUE
      ),
      length(idx),
      K2,
      byrow = TRUE
    )
  }
  
  
  ## ----------------------------------------------------------
  ## Convert matrices to long format
  ## ----------------------------------------------------------
  
  make_long <- function(M, nms) {
    
    z <- as.data.frame(M)
    
    names(z) <- nms
    
    z$Row <- row_id
    
    z <- tidyr::pivot_longer(
      z,
      -Row,
      names_to = "Factor",
      values_to = "Value"
    )
    
    z$Factor <- factor(
      z$Factor,
      levels = nms
    )
    
    z
  }
  
  
  domain_names <- colnames(X5)
  
  Xlong <- make_long(
    Xavg,
    domain_names
  )
  
  A1long <- make_long(
    A1avg,
    as.character(seq_len(K1))
  )
  
  A2long <- make_long(
    A2avg,
    as.character(seq_len(K2))
  )
  
  
  ## ----------------------------------------------------------
  ## Color limits
  ## ----------------------------------------------------------
  
  X_lim <- range(
    Xavg,
    na.rm = TRUE
  )
  
  A1_lim <- range(
    A1avg,
    na.rm = TRUE
  )
  
  A2_lim <- range(
    A2avg,
    na.rm = TRUE
  )
  
  
  if (diff(X_lim) == 0) {
    X_lim <- X_lim + c(-1e-6, 1e-6)
  }
  
  if (diff(A1_lim) == 0) {
    A1_lim <- A1_lim + c(-1e-6, 1e-6)
  }
  
  if (diff(A2_lim) == 0) {
    A2_lim <- A2_lim + c(-1e-6, 1e-6)
  }
  
  
  X_mid <- mean(X_lim)
  
  
  X_breaks <- pretty(
    X_lim,
    n = 5
  )
  
  X_breaks <- X_breaks[
    X_breaks >= X_lim[1] &
      X_breaks <= X_lim[2]
  ]
  
  
  A1_breaks <- pretty(
    A1_lim,
    n = 5
  )
  
  A1_breaks <- A1_breaks[
    A1_breaks >= A1_lim[1] &
      A1_breaks <= A1_lim[2]
  ]
  
  
  A2_breaks <- pretty(
    A2_lim,
    n = 5
  )
  
  A2_breaks <- A2_breaks[
    A2_breaks >= A2_lim[1] &
      A2_breaks <= A2_lim[2]
  ]
  
  
  ## ----------------------------------------------------------
  ## Outcome bar
  ## ----------------------------------------------------------
  
  outcome_df <- data.frame(
    Row = row_id,
    Group = factor(
      group_label,
      levels = unique(group_label)
    )
  )
  
  
  outcome_colors <-
    scales::hue_pal()(
      nlevels(outcome_df$Group)
    )
  
  names(outcome_colors) <-
    levels(outcome_df$Group)
  
  legend_order <-
    rev(levels(outcome_df$Group))
  
  
  ## ----------------------------------------------------------
  ## Plot theme
  ## ----------------------------------------------------------
  
  heat_theme <-
    ggplot2::theme_classic(
      base_size = 16
    ) +
    ggplot2::theme(
      plot.title =
        ggplot2::element_text(
          size = 18,
          face = "bold",
          hjust = 0.5
        ),
      axis.title =
        ggplot2::element_text(size = 16),
      axis.text =
        ggplot2::element_text(
          size = 13,
          colour = "black"
        ),
      axis.line =
        ggplot2::element_line(
          linewidth = 0.7,
          colour = "black"
        ),
      plot.margin =
        ggplot2::margin(5, 0, 5, 5)
    )
  
  
  ## ----------------------------------------------------------
  ## Color bars
  ## ----------------------------------------------------------
  
  make_bar_bwr <- function(
    lo,
    mid,
    hi,
    breaks
  ) {
    
    z <- data.frame(
      x = 1,
      y = seq(
        lo,
        hi,
        length.out = 1000
      )
    )
    
    z$value <- z$y
    
    
    ggplot2::ggplot(
      z,
      ggplot2::aes(
        x,
        y,
        fill = value
      )
    ) +
      ggplot2::geom_raster() +
      ggplot2::scale_fill_gradient2(
        low = "white",
        high = "red",
        midpoint = mid,
        limits = c(lo, hi)
      ) +
      ggplot2::scale_y_continuous(
        breaks = breaks,
        position = "right",
        expand = c(0, 0)
      ) +
      ggplot2::coord_cartesian(
        ylim = c(lo, hi),
        expand = FALSE
      ) +
      ggplot2::scale_x_continuous(
        expand = c(0, 0)
      ) +
      ggplot2::guides(fill = "none") +
      ggplot2::theme_classic(
        base_size = 11
      ) +
      ggplot2::theme(
        axis.title =
          ggplot2::element_blank(),
        axis.text.x =
          ggplot2::element_blank(),
        axis.ticks.x =
          ggplot2::element_blank(),
        axis.line.x =
          ggplot2::element_blank(),
        axis.text.y.left =
          ggplot2::element_blank(),
        axis.ticks.y.left =
          ggplot2::element_blank(),
        axis.line.y.left =
          ggplot2::element_blank(),
        axis.text.y.right =
          ggplot2::element_text(
            size = 11,
            colour = "black",
            margin =
              ggplot2::margin(l = 2)
          ),
        axis.ticks.y.right =
          ggplot2::element_line(
            colour = "black"
          ),
        axis.line.y.right =
          ggplot2::element_line(
            colour = "black"
          ),
        plot.margin =
          ggplot2::margin(
            5,
            3,
            5,
            0
          )
      )
  }
  
  
  make_bar_bw <- function(
    lo,
    hi,
    breaks
  ) {
    
    z <- data.frame(
      x = 1,
      y = seq(
        lo,
        hi,
        length.out = 1000
      )
    )
    
    z$value <- z$y
    
    
    ggplot2::ggplot(
      z,
      ggplot2::aes(
        x,
        y,
        fill = value
      )
    ) +
      ggplot2::geom_raster() +
      ggplot2::scale_fill_gradient(
        low = "white",
        high = "black",
        limits = c(lo, hi)
      ) +
      ggplot2::scale_y_continuous(
        breaks = breaks,
        position = "right",
        expand = c(0, 0)
      ) +
      ggplot2::coord_cartesian(
        ylim = c(lo, hi),
        expand = FALSE
      ) +
      ggplot2::scale_x_continuous(
        expand = c(0, 0)
      ) +
      ggplot2::guides(fill = "none") +
      ggplot2::theme_classic(
        base_size = 11
      ) +
      ggplot2::theme(
        axis.title =
          ggplot2::element_blank(),
        axis.text.x =
          ggplot2::element_blank(),
        axis.ticks.x =
          ggplot2::element_blank(),
        axis.line.x =
          ggplot2::element_blank(),
        axis.text.y.left =
          ggplot2::element_blank(),
        axis.ticks.y.left =
          ggplot2::element_blank(),
        axis.line.y.left =
          ggplot2::element_blank(),
        axis.text.y.right =
          ggplot2::element_text(
            size = 11,
            colour = "black",
            margin =
              ggplot2::margin(l = 2)
          ),
        axis.ticks.y.right =
          ggplot2::element_line(
            colour = "black"
          ),
        axis.line.y.right =
          ggplot2::element_line(
            colour = "black"
          ),
        plot.margin =
          ggplot2::margin(
            5,
            3,
            5,
            0
          )
      )
  }
  
  
  ## ----------------------------------------------------------
  ## Main heatmaps
  ## ----------------------------------------------------------
  
  p1 <- ggplot2::ggplot(
    Xlong,
    ggplot2::aes(
      Factor,
      Row,
      fill = Value
    )
  ) +
    ggplot2::geom_raster() +
    ggplot2::geom_hline(
      yintercept = boundaries,
      colour = "grey60",
      linewidth = 0.4
    ) +
    ggplot2::scale_fill_gradient2(
      low = "white",
      high = "red",
      midpoint = X_mid,
      limits = X_lim
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0.5, N + 0.5),
      expand = FALSE
    ) +
    ggplot2::labs(
      title = "Big Five",
      x = "Domain",
      y = "Observations"
    ) +
    ggplot2::guides(
      fill = "none"
    ) +
    heat_theme
  
  
  p2 <- ggplot2::ggplot(
    A1long,
    ggplot2::aes(
      Factor,
      Row,
      fill = Value
    )
  ) +
    ggplot2::geom_raster() +
    ggplot2::geom_hline(
      yintercept = boundaries,
      colour = "grey60",
      linewidth = 0.4
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "black",
      limits = A1_lim
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0.5, N + 0.5),
      expand = FALSE
    ) +
    ggplot2::labs(
      title = "A1",
      x = "A1",
      y = "Observations"
    ) +
    ggplot2::guides(
      fill = "none"
    ) +
    heat_theme
  
  
  p3 <- ggplot2::ggplot(
    A2long,
    ggplot2::aes(
      Factor,
      Row,
      fill = Value
    )
  ) +
    ggplot2::geom_raster() +
    ggplot2::geom_hline(
      yintercept = boundaries,
      colour = "grey60",
      linewidth = 0.4
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "black",
      limits = A2_lim
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0.5, N + 0.5),
      expand = FALSE
    ) +
    ggplot2::labs(
      title = "A2",
      x = "A2",
      y = "Observations"
    ) +
    ggplot2::guides(
      fill = "none"
    ) +
    heat_theme
  
  
  cb1 <- make_bar_bwr(
    X_lim[1],
    X_mid,
    X_lim[2],
    X_breaks
  )
  
  cb2 <- make_bar_bw(
    A1_lim[1],
    A1_lim[2],
    A1_breaks
  )
  
  cb3 <- make_bar_bw(
    A2_lim[1],
    A2_lim[2],
    A2_breaks
  )
  
  
  ## ----------------------------------------------------------
  ## Outcome bar
  ## ----------------------------------------------------------
  
  p4 <- ggplot2::ggplot(
    outcome_df,
    ggplot2::aes(
      1,
      Row,
      fill = Group
    )
  ) +
    ggplot2::geom_raster() +
    ggplot2::geom_hline(
      yintercept = boundaries,
      colour = "white",
      linewidth = 0.6
    ) +
    ggplot2::scale_fill_manual(
      values = outcome_colors,
      breaks = legend_order,
      drop = FALSE,
      guide =
        ggplot2::guide_legend(
          title = outcome,
          title.position = "top",
          label.position = "right"
        )
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0.5, N + 0.5),
      expand = FALSE
    ) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      title = outcome,
      x = NULL,
      y = NULL,
      fill = outcome
    ) +
    ggplot2::theme_classic(
      base_size = 16
    ) +
    ggplot2::theme(
      plot.title =
        ggplot2::element_text(
          size = 18,
          face = "bold",
          hjust = 0.5
        ),
      axis.text =
        ggplot2::element_blank(),
      axis.ticks =
        ggplot2::element_blank(),
      axis.title =
        ggplot2::element_blank(),
      axis.line =
        ggplot2::element_line(
          linewidth = 0.7,
          colour = "black"
        ),
      legend.position = "right",
      legend.box.spacing =
        grid::unit(0, "pt"),
      legend.margin =
        ggplot2::margin(
          0,
          0,
          0,
          3
        ),
      legend.title =
        ggplot2::element_text(
          size = 14,
          face = "bold"
        ),
      legend.text =
        ggplot2::element_text(
          size = 13
        ),
      legend.key.height =
        grid::unit(
          0.42,
          "cm"
        ),
      legend.key.width =
        grid::unit(
          0.42,
          "cm"
        ),
      plot.margin =
        ggplot2::margin(
          5,
          0,
          5,
          5
        )
    )
  
  
  p1 + cb1 + p2 + cb2 + p3 + cb3 + p4 +
    patchwork::plot_layout(
      widths = c(
        1,
        0.10,
        1.15,
        0.10,
        0.45,
        0.10,
        0.95
      )
    )
}


save_var_info <- function(data, survey_name, path) {
  
  var_names <- names(data)
  
  var_labels <- sapply(
    data,
    function(x) {
      lab <- attr(x, "label")
      if (is.null(lab)) NA_character_ else as.character(lab)
    }
  )
  
  var_info <- data.frame(
    index = seq_along(var_names),
    variable = var_names,
    label = unname(var_labels),
    stringsAsFactors = FALSE
  )
  
  file_name <- paste0(survey_name, "_var_info.csv")
  file_path <- file.path(path, file_name)
  
  write.csv(var_info, file_path, row.names = FALSE)
  
  message("Saved variable information to: ", file_path)
  
  return(invisible(var_info))
}

## ============================================================
## Helper: match survey rows to DDE / post-processed row order
## ============================================================

match_rows_to_dde <- function(X_survey, X_dde) {
  
  stopifnot(
    nrow(X_survey) == nrow(X_dde),
    ncol(X_survey) == ncol(X_dde)
  )
  
  make_key <- function(M) {
    apply(
      M,
      1,
      paste,
      collapse = "_"
    )
  }
  
  key_survey <- make_key(X_survey)
  key_dde <- make_key(X_dde)
  
  ## Initial match
  dde_to_survey <- match(
    key_dde,
    key_survey
  )
  
  ## Handle duplicated response profiles
  dup_keys <- unique(
    key_dde[duplicated(key_dde)]
  )
  
  for (k in dup_keys) {
    
    dde_rows <- which(
      key_dde == k
    )
    
    survey_rows <- which(
      key_survey == k
    )
    
    stopifnot(
      length(dde_rows) ==
        length(survey_rows)
    )
    
    ## Pair duplicate profiles in observed order
    dde_to_survey[dde_rows] <-
      survey_rows
  }
  
  ## Final checks
  stopifnot(
    !any(is.na(dde_to_survey)),
    length(unique(dde_to_survey)) ==
      nrow(X_dde)
  )
  
  return(dde_to_survey)
}


## ============================================================
## Helper: extract SPSS value labels
## ============================================================

extract_value_labels <- function(
    x,
    outcome
) {
  
  labs <- attr(
    x,
    "labels"
  )
  
  if (is.null(labs)) {
    return(NULL)
  }
  
  data.frame(
    outcome = outcome,
    code = as.numeric(labs),
    label = names(labs),
    stringsAsFactors = FALSE
  )
}


## ============================================================
## Helper: prepare LOOPR prediction data
## ============================================================

prepare_prediction_data <- function(
    path,
    data_dir = file.path(path, "Predictive_Analysis"),
    save_files = TRUE,
    overwrite = FALSE
) {
  
  ## ----------------------------------------------------------
  ## 1. Create output directory
  ## ----------------------------------------------------------
  
  dir.create(
    data_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  
  ## ----------------------------------------------------------
  ## 2. Load survey data
  ## ----------------------------------------------------------
  
  survey1 <- haven::read_sav(
    file.path(
      path,
      "Survey 1 - Merged - Working data file.sav"
    )
  )
  
  survey2 <- haven::read_sav(
    file.path(
      path,
      "Survey 2 - Merged - Working data file.sav"
    )
  )
  
  
  ## ----------------------------------------------------------
  ## 3. Save variable information
  ## ----------------------------------------------------------
  
  var_info1 <- save_var_info(
    survey1,
    "survey1",
    path
  )
  
  var_info2 <- save_var_info(
    survey2,
    "survey2",
    path
  )
  
  
  ## ----------------------------------------------------------
  ## 4. Extract and stack 60 BFI-2 personality items
  ## ----------------------------------------------------------
  
  personality_base <- paste0(
    "Personality_",
    1:60
  )
  
  survey1_personality <- paste0(
    "S1_",
    personality_base
  )
  
  survey2_personality <- paste0(
    "S2_",
    personality_base
  )
  
  
  stopifnot(
    all(
      survey1_personality %in%
        names(survey1)
    ),
    all(
      survey2_personality %in%
        names(survey2)
    )
  )
  
  
  survey1_analysis <- survey1[
    ,
    survey1_personality
  ]
  
  survey2_analysis <- survey2[
    ,
    survey2_personality
  ]
  
  
  ## Standardize names
  names(survey1_analysis) <-
    personality_base
  
  names(survey2_analysis) <-
    personality_base
  
  
  ## Keep survey source temporarily
  survey1_analysis$survey <- "survey1"
  survey2_analysis$survey <- "survey2"
  
  
  ## Original Survey 1 + Survey 2 order
  survey_both <- rbind(
    survey1_analysis,
    survey2_analysis
  )
  
  
  stopifnot(
    nrow(survey_both) ==
      nrow(survey1) + nrow(survey2)
  )
  
  
  ## ----------------------------------------------------------
  ## 5. Raw personality matrix in survey order
  ## ----------------------------------------------------------
  
  X_survey <- as.matrix(
    survey_both[
      ,
      personality_base
    ]
  )
  
  storage.mode(X_survey) <- "numeric"
  
  
  ## ----------------------------------------------------------
  ## 6. Load CORRECT post-processed LOOPR results
  ##    MATLAB v7.3 / HDF5
  ## ----------------------------------------------------------
  
  path_post <- file.path(path, "Analysis")
  
  mat_file <- file.path(
    path_post,
    "LOOPR_PostProcessed.mat"
  )
  
  if (!file.exists(mat_file)) {
    stop(
      paste0(
        "Cannot find:\n",
        mat_file
      )
    )
  }
  
  
  mat_h5 <- hdf5r::H5File$new(
    mat_file,
    mode = "r"
  )
  
  on.exit(
    mat_h5$close_all(),
    add = TRUE
  )
  
  
  ## Check that the higher-level LOOPR object exists
  if (!"LOOPR" %in% names(mat_h5)) {
    stop(
      paste0(
        "LOOPR was not found in LOOPR_PostProcessed.mat.\n",
        "Top-level objects are:\n",
        paste(names(mat_h5), collapse = ", ")
      )
    )
  }
  
  
  ## LOOPR is the MATLAB struct/group containing
  ## X, A1, A2, B1, B2, etc.
  LOOPR_h5 <- mat_h5[["LOOPR"]]
  
  
  ## Check required fields
  required_fields <- c(
    "X",
    "A1",
    "A2"
  )
  
  missing_fields <- setdiff(
    required_fields,
    names(LOOPR_h5)
  )
  
  if (length(missing_fields) > 0) {
    stop(
      paste0(
        "Missing field(s) inside LOOPR:\n",
        paste(missing_fields, collapse = ", "),
        "\n\nAvailable LOOPR fields are:\n",
        paste(names(LOOPR_h5), collapse = ", ")
      )
    )
  }
  
  
  ## Extract the correct post-processed variables
  X_dde <- LOOPR_h5[["X"]][, ]
  A1 <- LOOPR_h5[["A1"]][, ]
  A2 <- LOOPR_h5[["A2"]][, ]
  
  
  ## Make sure they are ordinary numeric matrices
  X_dde <- as.matrix(X_dde)
  A1 <- as.matrix(A1)
  A2 <- as.matrix(A2)
  
  storage.mode(X_dde) <- "numeric"
  storage.mode(A1) <- "numeric"
  storage.mode(A2) <- "numeric"
  
  
  cat(
    "\nLoaded LOOPR_PostProcessed.mat:\n",
    "X  :", nrow(X_dde), "x", ncol(X_dde), "\n",
    "A1 :", nrow(A1), "x", ncol(A1), "\n",
    "A2 :", nrow(A2), "x", ncol(A2), "\n\n"
  )
  
  
  ## Expected:
  ## X  = 6126 x 60
  ## A1 = 6126 x 9
  ## A2 = 6126 x 1
  
  stopifnot(
    nrow(X_survey) == nrow(X_dde),
    ncol(X_survey) == ncol(X_dde),
    nrow(A1) == nrow(X_dde),
    nrow(A2) == nrow(X_dde)
  )
  
  
  ## ----------------------------------------------------------
  ## 7. Match survey rows to LOOPR post-processed row order
  ## ----------------------------------------------------------
  
  dde_to_survey <- match_rows_to_dde(
    X_survey = X_survey,
    X_dde = X_dde
  )
  
  
  survey_aligned <- survey_both[
    dde_to_survey,
    ,
    drop = FALSE
  ]
  
  
  ## ----------------------------------------------------------
  ## 8. Construct final X
  ## ----------------------------------------------------------
  
  X_aligned <- as.matrix(
    survey_aligned[
      ,
      personality_base
    ]
  )
  
  storage.mode(X_aligned) <- "numeric"
  
  
  ## Verify exact alignment with LOOPR.X
  stopifnot(
    all(
      X_aligned == X_dde
    )
  )
  
  
  X <- as.data.frame(
    X_aligned
  )
  
  names(X) <-
    personality_base
  
  
  ## ----------------------------------------------------------
  ## 9. Construct final A
  ## ----------------------------------------------------------
  
  colnames(A1) <- paste0(
    "A1_",
    seq_len(
      ncol(A1)
    )
  )
  
  colnames(A2) <- paste0(
    "A2_",
    seq_len(
      ncol(A2)
    )
  )
  
  
  A <- data.frame(
    A1,
    A2
  )
  
  
  ## ----------------------------------------------------------
  ## 10. Construct outcomes
  ## ----------------------------------------------------------
  
  Y1 <- data.frame(
    
    education_level =
      as.numeric(
        survey1$S1_WorkEducation
      ),
    
    income =
      as.numeric(
        survey1$S1_HouseholdIncome
      ),
    
    marriage_status =
      as.numeric(
        survey1$S1_MaritalStatus
      ),
    
    frequency_of_prayer =
      as.numeric(
        survey1$S1_Spirituality_28
      )
  )
  
  
  Y2 <- data.frame(
    
    education_level =
      as.numeric(
        survey2$S2_WorkEducation
      ),
    
    income =
      as.numeric(
        survey2$S2_HouseholdIncome
      ),
    
    marriage_status =
      as.numeric(
        survey2$S2_MaritalStatus
      ),
    
    ## Frequency of prayer is not available in Survey 2
    frequency_of_prayer =
      rep(
        NA_real_,
        nrow(survey2)
      )
  )
  
  
  Y_both <- rbind(
    Y1,
    Y2
  )
  
  
  ## Reorder outcomes into the SAME LOOPR row order
  Y <- Y_both[
    dde_to_survey,
    ,
    drop = FALSE
  ]
  
  
  rownames(Y) <- NULL
  
  
  ## ----------------------------------------------------------
  ## 11. Extract outcome value labels from SPSS
  ## ----------------------------------------------------------
  
  label_list <- list(
    
    extract_value_labels(
      survey1$S1_WorkEducation,
      "education_level"
    ),
    
    extract_value_labels(
      survey1$S1_HouseholdIncome,
      "income"
    ),
    
    extract_value_labels(
      survey1$S1_MaritalStatus,
      "marriage_status"
    ),
    
    extract_value_labels(
      survey1$S1_Spirituality_28,
      "frequency_of_prayer"
    )
  )
  
  
  ## Remove NULL entries
  label_list <- Filter(
    Negate(is.null),
    label_list
  )
  
  
  if (length(label_list) > 0) {
    
    Y_labels <- do.call(
      rbind,
      label_list
    )
    
    ## Remove duplicate outcome/code entries
    Y_labels <- Y_labels[
      !duplicated(
        Y_labels[
          ,
          c(
            "outcome",
            "code"
          )
        ]
      ),
      ,
      drop = FALSE
    ]
    
    rownames(Y_labels) <- NULL
    
  } else {
    
    Y_labels <- data.frame(
      outcome = character(0),
      code = numeric(0),
      label = character(0)
    )
  }
  
  
  ## ----------------------------------------------------------
  ## 12. Final consistency checks
  ## ----------------------------------------------------------
  
  stopifnot(
    nrow(X) == nrow(A),
    nrow(X) == nrow(Y),
    nrow(X) == nrow(X_dde)
  )
  
  
  cat(
    "\nFinal prediction data:\n",
    "X:", nrow(X), "x", ncol(X), "\n",
    "A:", nrow(A), "x", ncol(A), "\n",
    "Y:", nrow(Y), "x", ncol(Y), "\n\n"
  )
  
  
  ## ----------------------------------------------------------
  ## 13. Define output files
  ## ----------------------------------------------------------
  
  X_file <- file.path(
    data_dir,
    "X.csv"
  )
  
  A_file <- file.path(
    data_dir,
    "A.csv"
  )
  
  Y_file <- file.path(
    data_dir,
    "Y.csv"
  )
  
  Y_labels_file <- file.path(
    data_dir,
    "Y_labels.csv"
  )
  
  
  ## ----------------------------------------------------------
  ## 14. Save files
  ## ----------------------------------------------------------
  
  if (save_files) {
    
    output_files <- c(
      X_file,
      A_file,
      Y_file,
      Y_labels_file
    )
    
    
    if (
      !overwrite &&
      any(
        file.exists(
          output_files
        )
      )
    ) {
      
      stop(
        paste0(
          "Prediction files already exist in:\n",
          data_dir,
          "\nSet overwrite = TRUE to replace them."
        )
      )
    }
    
    
    write.csv(
      X,
      X_file,
      row.names = FALSE
    )
    
    
    write.csv(
      A,
      A_file,
      row.names = FALSE
    )
    
    
    write.csv(
      Y,
      Y_file,
      row.names = FALSE
    )
    
    
    write.csv(
      Y_labels,
      Y_labels_file,
      row.names = FALSE
    )
    
    
    message(
      "Saved prediction data to: ",
      data_dir
    )
  }
  
  
  ## ----------------------------------------------------------
  ## 15. Return objects invisibly
  ## ----------------------------------------------------------
  
  invisible(
    list(
      X = X,
      A = A,
      Y = Y,
      Y_labels = Y_labels,
      dde_to_survey = dde_to_survey,
      survey_aligned = survey_aligned
    )
  )
}