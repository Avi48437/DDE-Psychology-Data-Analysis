library(haven)

root_dir <- "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology"
loopr_dir <- file.path(root_dir,"LOPR")
data_dir <- file.path(loopr_dir,"Data")
out_dir <- file.path(loopr_dir,"Analysis","Predictive_Data")

if(!dir.exists(out_dir)) dir.create(out_dir,recursive=TRUE)

## ============================================================
## Load master LOOPR data
## This determines the exact row ordering
## ============================================================

master <- read.csv(
  file.path(data_dir,"LOOPRDataAgeGender.csv"),
  check.names=FALSE
)

## ============================================================
## Load Survey 1 and Survey 2
## ============================================================

S1 <- read_sav(
  file.path(data_dir,"Survey 1 - Merged - Working data file.sav")
)

S2 <- read_sav(
  file.path(data_dir,"Survey 2 - Merged - Working data file.sav")
)

## ============================================================
## Match CaseID to master LOOPR ordering
## ============================================================

master_id <- as.character(master$CaseID)
S1_id <- as.character(S1$CaseID)
S2_id <- as.character(S2$CaseID)

if(anyDuplicated(master_id)) stop("Duplicate CaseID values in LOOPR master data.")
if(anyDuplicated(S1_id)) stop("Duplicate CaseID values in Survey 1.")
if(anyDuplicated(S2_id)) stop("Duplicate CaseID values in Survey 2.")

idx1 <- match(master_id,S1_id)
idx2 <- match(master_id,S2_id)

if(any(!is.na(idx1) & !is.na(idx2)))
  stop("Some CaseID values occur in both Survey 1 and Survey 2.")

if(any(is.na(idx1) & is.na(idx2)))
  stop("Some LOOPR CaseID values could not be matched to S1 or S2.")

survey_origin <- ifelse(!is.na(idx1),"S1","S2")

cat("\nMaster rows:",nrow(master),"\n")
cat("Survey 1 rows:",sum(survey_origin=="S1"),"\n")
cat("Survey 2 rows:",sum(survey_origin=="S2"),"\n")

cat("\nRow-order structure:\n")
print(rle(survey_origin))

## ============================================================
## Functions for exact master-order extraction
## ============================================================

combine_variable <- function(v1,v2){
  
  z <- rep(NA_real_,nrow(master))
  
  from1 <- !is.na(idx1)
  from2 <- !is.na(idx2)
  
  z[from1] <- as.numeric(S1[[v1]][idx1[from1]])
  z[from2] <- as.numeric(S2[[v2]][idx2[from2]])
  
  z
}

S1_only_variable <- function(v1){
  
  z <- rep(NA_real_,nrow(master))
  
  from1 <- !is.na(idx1)
  
  z[from1] <- as.numeric(S1[[v1]][idx1[from1]])
  
  z
}

## ============================================================
## Construct Y
##
## Josh outcomes available in LOOPR:
## educ       : S1 + S2
## faminc_new : S1 + S2
## marstat    : S1 + S2
## pew_prayer : S1 only; S2 = NA
## ============================================================

Y <- data.frame(
  educ=combine_variable(
    "S1_WorkEducation",
    "S2_WorkEducation"
  ),
  
  faminc_new=combine_variable(
    "S1_HouseholdIncome",
    "S2_HouseholdIncome"
  ),
  
  marstat=combine_variable(
    "S1_MaritalStatus",
    "S2_MaritalStatus"
  ),
  
  pew_prayer=S1_only_variable(
    "S1_Spirituality_28"
  )
)

## ============================================================
## Construct reverse-coded/scored Big-Five matrix X_rc
##
## These are already scored LOOPR domain measures.
## Do NOT reverse-code them again.
## ============================================================

X_rc <- data.frame(
  EXT=combine_variable(
    "S1_PersonalityExtraversion",
    "S2_PersonalityExtraversion"
  ),
  
  NEM=combine_variable(
    "S1_PersonalityNegativeEmotionality",
    "S2_PersonalityNegativeEmotionality"
  ),
  
  AGR=combine_variable(
    "S1_PersonalityAgreeableness",
    "S2_PersonalityAgreeableness"
  ),
  
  CSN=combine_variable(
    "S1_PersonalityConscientiousness",
    "S2_PersonalityConscientiousness"
  ),
  
  OPN=combine_variable(
    "S1_PersonalityOpenMindedness",
    "S2_PersonalityOpenMindedness"
  )
)

## ============================================================
## Final consistency checks
## ============================================================

if(nrow(Y)!=nrow(master))
  stop("Y row count does not match LOOPR master data.")

if(nrow(X_rc)!=nrow(master))
  stop("X_rc row count does not match LOOPR master data.")

cat("\nY dimensions:",nrow(Y),"x",ncol(Y),"\n")
cat("X_rc dimensions:",nrow(X_rc),"x",ncol(X_rc),"\n")

cat("\nY missing values:\n")
print(colSums(is.na(Y)))

cat("\nX_rc missing values:\n")
print(colSums(is.na(X_rc)))

cat("\nEducation:\n")
print(table(Y$educ,useNA="ifany"))

cat("\nHousehold income:\n")
print(table(Y$faminc_new,useNA="ifany"))

cat("\nMarital status:\n")
print(table(Y$marstat,useNA="ifany"))

cat("\nPrayer:\n")
print(table(Y$pew_prayer,useNA="ifany"))

## ============================================================
## Save
## ============================================================

write.csv(
  Y,
  file.path(out_dir,"Y.csv"),
  row.names=FALSE,
  na=""
)

write.csv(
  X_rc,
  file.path(out_dir,"X_rc.csv"),
  row.names=FALSE,
  na=""
)

cat("\nSaved:\n")
cat(file.path(out_dir,"Y.csv"),"\n")
cat(file.path(out_dir,"X_rc.csv"),"\n")