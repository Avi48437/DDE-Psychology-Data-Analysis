library(haven)

root_dir <- "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology"
loopr_dir <- file.path(root_dir,"LOPR")
data_dir <- file.path(loopr_dir,"Data")
analysis_dir <- file.path(loopr_dir,"Analysis","PostProcessed_Data")

if (!dir.exists(analysis_dir)) dir.create(analysis_dir,recursive=TRUE)

## ============================================================
## Load Survey 1
## ============================================================

dat <- read_sav(
  file.path(data_dir,"Survey 1 - Merged - Working data file.sav")
)

## ============================================================
## BFI-2 variables
## ============================================================

personality_vars <- paste0("S1_Personality_",1:60)

## ============================================================
## Extract SPSS question labels
## ============================================================

Question <- vapply(
  dat[personality_vars],
  function(x) {
    lab <- attr(x,"label")
    if (is.null(lab)) NA_character_ else as.character(lab)
  },
  character(1)
)

## ============================================================
## Remove common prefix
## ============================================================

Question <- sub(
  "^I am someone who\\.\\.\\. -[[:space:]]*",
  "",
  Question
)

## ============================================================
## Construct BFI-2 item lookup
## ============================================================

BFI_Item_Text <- data.frame(
  OriginalBFIItem=1:60,
  Item=paste0("BFI",1:60),
  OriginalVariable=personality_vars,
  Question=unname(Question),
  stringsAsFactors=FALSE
)

## ============================================================
## Inspect
## ============================================================

print(BFI_Item_Text)

## ============================================================
## Save
## ============================================================

output_file <- file.path(
  analysis_dir,
  "LOOPR_BFI2_Item_Text.csv"
)

write.csv(
  BFI_Item_Text,
  output_file,
  row.names=FALSE
)

cat("\nSaved to:\n",output_file,"\n")