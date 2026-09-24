library(haven)

root_dir <- "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology"
loopr_dir <- file.path(root_dir,"LOPR")
data_dir <- file.path(loopr_dir,"Data")
predictive_dir <- file.path(loopr_dir,"Analysis","Predictive_Data")

master <- read.csv(file.path(data_dir,"LOOPRDataAgeGender.csv"))
S1 <- read_sav(file.path(data_dir,"Survey 1 - Merged - Working data file.sav"))
S2 <- read_sav(file.path(data_dir,"Survey 2 - Merged - Working data file.sav"))

idx1 <- match(master$CaseID,S1$CaseID)
idx2 <- match(master$CaseID,S2$CaseID)

in1 <- !is.na(idx1)
in2 <- !is.na(idx2)

stopifnot(!any(in1 & in2),all(in1 | in2))

get_var <- function(v1,v2){
  z <- rep(NA_real_,nrow(master))
  z[in1] <- as.numeric(S1[[v1]][idx1[in1]])
  z[in2] <- as.numeric(S2[[v2]][idx2[in2]])
  z
}

D_LOOPR <- data.frame(
  CaseID=master$CaseID,
  Survey=ifelse(in1,"S1","S2"),
  Age=get_var("S1_Age_1","S2_Age_1"),
  Sex=get_var("S1_Sex","S2_Sex"),
  HispanicLatino=get_var("S1_HispanicLatino","S2_HispanicLatino")
)

for(j in 1:6){
  D_LOOPR[[paste0("Race_",j)]] <- get_var(
    paste0("S1_Race_",j),
    paste0("S2_Race_",j)
  )
}

## Exact race-category meanings from SPSS
race_labels <- sapply(1:6,function(j){
  lab <- attr(S1[[paste0("S1_Race_",j)]],"label")
  if(is.null(lab) || lab=="") paste0("Race_",j) else as.character(lab)
})

print(data.frame(
  Variable=paste0("Race_",1:6),
  Meaning=race_labels
))

## Merge race indicators into one CART-friendly categorical variable
R <- as.matrix(D_LOOPR[paste0("Race_",1:6)])

D_LOOPR$Race <- apply(R,1,function(z){
  ii <- which(z==1)
  if(length(ii)==0) return(NA_character_)
  if(length(ii)>1) return("Multiple races")
  race_labels[ii]
})

## Put merged Race beside the demographic variables
D_LOOPR <- D_LOOPR[,c(
  "CaseID","Survey","Age","Sex","HispanicLatino","Race",
  paste0("Race_",1:6)
)]

write.csv(
  D_LOOPR,
  file.path(predictive_dir,"D_LOOPR.csv"),
  row.names=FALSE
)

table(D_LOOPR$Race,useNA="ifany")