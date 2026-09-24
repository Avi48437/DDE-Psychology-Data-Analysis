## ============================================================
## IPIP-100 downstream predictive analysis
## ============================================================

rm(list = ls())

## Folder containing X.csv, A.csv, Y.csv
# data_dir <- "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology/2. IPIP 100/Predictive_Analysis"
# data_dir <- "/Users/youlinzhang/Downloads/DDE_Fitting_Pshycology/2. IPIP 100/Predictive_Analysis"


## Load packages and helper functions
library(haven)
library(hdf5r) # to read .mat file
# stratification packages
library(rpart)
library(rpart.plot)
# library(rsample)
# packages to grow trees
library(rpart)
library(rpart.plot)
# package to read/load .mat file
library(rmatio)
source("/Users/youlinzhang/Downloads/DDE_Fitting_Pshycology/LOPR/Predictive_Analysis/Helpers.R")

## ============================================================
## Reverse-coded IPIP-100 responses
##
## Official IPIP NEO-domain scoring key:
## https://www.ipip.ori.org/newNEODomainsKey.htm
##
## X is already in the canonical order:
## EXT | EST | AGR | CSN | OPN
##
## X_rc = reverse-coded version of X.
##
## Note: the official scale is Neuroticism, whereas we use
## EST = Emotional Stability. Therefore the Neuroticism
## direction is reversed when constructing EST.
## ============================================================

# item_map <- read.csv(
#   "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology/2. IPIP 100/Analysis/IPIP100_Item_Map.csv"
# )
# item_map <- read.csv(
#   "/Users/youlinzhang/Downloads/DDE_Fitting_Pshycology/2. IPIP 100/Analysis/IPIP100_Item_Map.csv"
# )
# 
# reverse_q <- c(
#   # EXT
#   14,18,29,39,49,59,69,79,89,99,
#   # EST
#   12,17,30,40,50,60,70,80,90,100,
#   # AGR
#   2,9,22,32,42,52,62,72,82,92,
#   # CSN
#   8,20,28,38,48,58,68,78,88,98,
#   # OPN
#   4,7,24,34,44,54,64,74,84,94
# )
# 
# ## X columns are permuted, so map original q-numbers
# ## to their current canonical column positions
# 
# reverse_cols <- which(item_map$OriginalRow %in% reverse_q)
# 
# ## Sanity checks
# 
# stopifnot(ncol(X)==100)
# stopifnot(length(reverse_cols)==50)
# stopifnot(length(unique(reverse_q))==50)
# 
# print(table(
#   Domain=item_map$Domain,
#   Reverse=item_map$OriginalRow %in% reverse_q
# ))
# 
# ## Create reverse-coded copy
# 
# X_rc <- X
# X_rc[,reverse_cols] <- 6-X_rc[,reverse_cols]
# 
# 
# ## ============================================================
# plot_A1_X_by_outcome(X_rc,A,Y,"educ")
# plot_A1_X_by_outcome(X_rc,A,Y,"marstat")
# plot_A1_X_by_outcome(X_rc,A,Y,"pew_prayer")
# plot_A1_X_by_outcome(X_rc,A,Y,"votereg")
# plot_A1_X_by_outcome(X_rc,A,Y,"faminc_new",drop_codes = 97)
# plot_A1_X_by_outcome(X_rc,A,Y,"ideo5")
# plot_A1_X_by_outcome(X_rc,A,Y,"pid3")
# plot_A1_X_by_outcome(X_rc,A,Y,"pid7")
# 
# plot_A1_X_by_outcome(X,A,Y,"education_level")
## ============================================================

# Oh I see dont use that .mat file it does not contains the fit. I remeber now. It was the file save d after the tuning but for memory problem I did not saved A1 and A2 and realized later. So when I ran the DDE fitting with the best tuning setup I saw that in each iteration it reveals different structure. So I ran it with the best tuning setup and once I got the 9,1 setup I saved it as Loopr_DDE_fit.mat

## ============================================================
## LOOPR: data preparation
## ============================================================
## Load data
path <- "/Users/youlinzhang/Downloads/DDE_Fitting_Pshycology/LOPR"

data_dir <- file.path(
  path,
  "Predictive_Analysis"
)

# If X, A, Y already pre-stored
X <- read.csv(file.path(data_dir, "X.csv"))
A <- read.csv(file.path(data_dir, "A.csv"))
Y <- read.csv(file.path(data_dir, "Y.csv"))

# If not, prepare prediction data using helper function "prepare_prediction_data"
prepare_prediction_data(
  path = path,
  data_dir = data_dir,
  save_files = TRUE,
  overwrite = TRUE
)


## ============================================================
## Prediction analysis: frequency of prayer
## Baseline model using Big Five domain means
## ============================================================

## ------------------------------------------------------------
## 1. Prepare prediction data
## ------------------------------------------------------------

dat_prayer <- data.frame(
  frequency_of_prayer = as.numeric(Y$frequency_of_prayer),
  X_domain_mean
)

## Remove observations with missing outcome or predictors
dat_prayer <- dat_prayer[complete.cases(dat_prayer), ]

## Check outcome distribution
table(dat_prayer$frequency_of_prayer)
prop.table(table(dat_prayer$frequency_of_prayer))


## ------------------------------------------------------------
## 2. Stratified 80/20 train-test split
## ------------------------------------------------------------

set.seed(123)

train_idx <- unlist(
  lapply(
    split(
      seq_len(nrow(dat_prayer)),
      dat_prayer$frequency_of_prayer
    ),
    function(idx) {
      sample(
        idx,
        size = floor(0.80 * length(idx))
      )
    }
  )
)

train_prayer <- dat_prayer[train_idx, , drop = FALSE]
test_prayer  <- dat_prayer[-train_idx, , drop = FALSE]


## Check train-test distributions
table(train_prayer$frequency_of_prayer)
table(test_prayer$frequency_of_prayer)

prop.table(table(train_prayer$frequency_of_prayer))
prop.table(table(test_prayer$frequency_of_prayer))


## ------------------------------------------------------------
## 3. Grow a relatively large Big Five regression tree
## ------------------------------------------------------------

tree_prayer_big5_full <- rpart(
  frequency_of_prayer ~ .,
  data = train_prayer,
  method = "anova",
  control = rpart.control(
    cp = 0.001,
    minsplit = 30,
    minbucket = 10,
    maxdepth = 6
  )
)


## ------------------------------------------------------------
## 4. Examine tree complexity using cross-validation
## ------------------------------------------------------------

printcp(tree_prayer_big5_full)
plotcp(tree_prayer_big5_full)

cp_table_big5 <- tree_prayer_big5_full$cptable


## ------------------------------------------------------------
## 5. Minimum-CV-error tree
## ------------------------------------------------------------

row_min_big5 <- which.min(
  cp_table_big5[, "xerror"]
)

cp_min_big5 <- cp_table_big5[
  row_min_big5,
  "CP"
]

tree_prayer_big5_min <- prune(
  tree_prayer_big5_full,
  cp = cp_min_big5
)

rpart.plot(
  tree_prayer_big5_min,
  type = 2,
  extra = 101,
  fallen.leaves = TRUE,
  roundint = FALSE,
  main = "Prayer: Big Five Minimum CV Error Tree"
)


## ------------------------------------------------------------
## 6. One-standard-error tree
## ------------------------------------------------------------

xerror_min_big5 <- cp_table_big5[
  row_min_big5,
  "xerror"
]

xstd_min_big5 <- cp_table_big5[
  row_min_big5,
  "xstd"
]

threshold_1se_big5 <-
  xerror_min_big5 + xstd_min_big5

eligible_rows_big5 <- which(
  cp_table_big5[, "xerror"] <= threshold_1se_big5
)

## First eligible row = simplest tree within 1 SE
row_1se_big5 <- eligible_rows_big5[1]

cp_1se_big5 <- cp_table_big5[
  row_1se_big5,
  "CP"
]

tree_prayer_big5_1se <- prune(
  tree_prayer_big5_full,
  cp = cp_1se_big5
)

rpart.plot(
  tree_prayer_big5_1se,
  type = 2,
  extra = 101,
  fallen.leaves = TRUE,
  roundint = FALSE,
  main = "Prayer: Big Five 1-SE Tree"
)


## ------------------------------------------------------------
## 7. Evaluate minimum-CV tree on the test set
## ------------------------------------------------------------

pred_prayer_big5 <- predict(
  tree_prayer_big5_min,
  newdata = test_prayer
)

mae_prayer_big5 <- mean(
  abs(
    test_prayer$frequency_of_prayer -
      pred_prayer_big5
  )
)

rmse_prayer_big5 <- sqrt(
  mean(
    (
      test_prayer$frequency_of_prayer -
        pred_prayer_big5
    )^2
  )
)

mae_prayer_big5
rmse_prayer_big5


## ------------------------------------------------------------
## 8. Null-model comparison
## ------------------------------------------------------------

pred_prayer_null <- rep(
  mean(train_prayer$frequency_of_prayer),
  nrow(test_prayer)
)

mae_prayer_null <- mean(
  abs(
    test_prayer$frequency_of_prayer -
      pred_prayer_null
  )
)

rmse_prayer_null <- sqrt(
  mean(
    (
      test_prayer$frequency_of_prayer -
        pred_prayer_null
    )^2
  )
)


## ------------------------------------------------------------
## 9. Prediction improvement over null model
## ------------------------------------------------------------

test_R2_prayer_big5 <-
  1 -
  sum(
    (test_prayer$frequency_of_prayer -
       pred_prayer_big5)^2
  ) /
  sum(
    (test_prayer$frequency_of_prayer -
       mean(train_prayer$frequency_of_prayer))^2
  )

mae_prayer_big5
rmse_prayer_big5

mae_prayer_null
rmse_prayer_null

test_R2_prayer_big5

## ------------------------------------------------------------
## Checks
## ------------------------------------------------------------

# stopifnot(
#   ncol(X_rc) == 60,
#   ncol(X_domain_mean) == 5,
#   nrow(X_rc) == nrow(Y)
# )
# 
# dim(X_rc)
# dim(X_domain_mean)
# dim(Y)

