## ============================================================
## IPIP-98 demographics
## ============================================================

root_dir <- "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology2"
ipip_dir <- file.path(root_dir,"2. IPIP 100")
data_dir <- file.path(ipip_dir,"Data")
predictive_dir <- file.path(ipip_dir,"Analysis","Predictive_Data")

B5_wo <- read.csv(
  file.path(data_dir,"B5_wo.csv"),
  check.names=FALSE
)

X98 <- read.csv(file.path(predictive_dir,"X98.csv"))
A98 <- read.csv(file.path(predictive_dir,"A98.csv"))
Y98 <- read.csv(file.path(predictive_dir,"Y.csv"))

## ============================================================
## Alignment checks
## ============================================================

stopifnot(
  nrow(B5_wo)==1500,
  nrow(B5_wo)==nrow(X98),
  nrow(B5_wo)==nrow(A98),
  nrow(B5_wo)==nrow(Y98)
)

## ============================================================
## Raw demographic variables
## Keep original coding
## ============================================================

demo_vars <- c(
  "caseid",
  "birthyr",
  "gender",
  "race"
)

stopifnot(all(demo_vars %in% names(B5_wo)))

D_IPIP98 <- B5_wo[,demo_vars,drop=FALSE]

## Optional age variable
D_IPIP98$age <- 2018-D_IPIP98$birthyr

D_IPIP98 <- D_IPIP98[,c(
  "caseid",
  "age",
  "birthyr",
  "gender",
  "race"
)]

## ============================================================
## Save
## ============================================================

write.csv(
  D_IPIP98,
  file.path(predictive_dir,"D_IPIP98.csv"),
  row.names=FALSE
)

cat("\nIPIP-98 demographics saved:\n")
cat(file.path(predictive_dir,"D_IPIP98.csv"),"\n")
cat("Dimensions:",nrow(D_IPIP98),"x",ncol(D_IPIP98),"\n")

summary(D_IPIP98)
