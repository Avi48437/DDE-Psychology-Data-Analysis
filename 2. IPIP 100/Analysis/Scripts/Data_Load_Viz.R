## ============================================================
## IPIP-100 / IPIP-98 downstream predictive analysis
##
## Conventional Big-Five scores were constructed in MATLAB
## using the official IPIP NEO-domain scoring key:
##
## https://www.ipip.ori.org/newNEODomainsKey.html
##
## X100 : 5 conventional IPIP-100 domain averages
## X98  : 5 conventional IPIP-98 domain averages
## A100 : IPIP-100 DDE representation
## A98  : IPIP-98 DDE representation
## Y    : common downstream outcomes
##
## X98 differs from X100 because BIGFIVE_44 and BIGFIVE_51
## were removed before computing the IPIP-98 OPN score.
## ============================================================
data_dir <- "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology2/2. IPIP 100/Analysis/Predictive_Data"

## Load conventional Big-Five scores
X100 <- read.csv(file.path(data_dir, "X100.csv"))
X98  <- read.csv(file.path(data_dir, "X98.csv"))

## Load DDE representations
A100 <- read.csv(file.path(data_dir, "A100.csv"))
A98  <- read.csv(file.path(data_dir, "A98.csv"))

## Common outcomes
Y <- read.csv(file.path(data_dir, "Y.csv"))

## ============================================================
## Sanity checks
## ============================================================

stopifnot(ncol(X100) == 5)
stopifnot(ncol(X98) == 5)

stopifnot(identical(names(X100), c("EXT","EST","AGR","CSN","OPN")))
stopifnot(identical(names(X98), c("EXT","EST","AGR","CSN","OPN")))

stopifnot(nrow(X100) == nrow(Y))
stopifnot(nrow(X98) == nrow(Y))
stopifnot(nrow(A100) == nrow(Y))
stopifnot(nrow(A98) == nrow(Y))

cat("X100:", dim(X100), "\n")
cat("X98 :", dim(X98), "\n")
cat("A100:", dim(A100), "\n")
cat("A98 :", dim(A98), "\n")
cat("Y   :", dim(Y), "\n")

# Visualization of the IPIP100 summary
plot_A1_X_by_outcome(X100,A100,Y,"educ")
plot_A1_X_by_outcome(X100,A100,Y,"marstat")
plot_A1_X_by_outcome(X100,A100,Y,"pew_prayer")
plot_A1_X_by_outcome(X100,A100,Y,"votereg")
plot_A1_X_by_outcome(X100,A100,Y,"faminc_new",drop_codes=97)
plot_A1_X_by_outcome(X100,A100,Y,"ideo5")
plot_A1_X_by_outcome(X100,A100,Y,"pid3")
plot_A1_X_by_outcome(X100,A100,Y,"pid7")

# Visualization of the IPIP98 summary
plot_A1_X_by_outcome(X98,A98,Y,"educ")

plot_A1_X_by_outcome(X98,A98,Y,"marstat",grouping=TRUE,
                     group_list=list("Separated/Divorced"=c(2,3),
                                     "Never married/Domestic-civil"=c(5,6)))

plot_A1_X_by_outcome(X98,A98,Y,"faminc_new",grouping=TRUE,
                     group_list=list("< $40k"=1:4,
                                     "$40k - < $80k"=5:8,
                                     "$80k - < $150k"=9:11,
                                     ">= $150k"=12:16),
                     Discard=TRUE,discard_codes=97)

plot_A1_X_by_outcome(X98,A98,Y,"votereg",Discard=TRUE,discard_codes=3)

plot_A1_X_by_outcome(X98,A98,Y,"pew_prayer",Discard=TRUE,discard_codes=8)

plot_A1_X_by_outcome(X98,A98,Y,"ideo5")
plot_A1_X_by_outcome(X98,A98,Y,"pid3")
plot_A1_X_by_outcome(X98,A98,Y,"pid7")

