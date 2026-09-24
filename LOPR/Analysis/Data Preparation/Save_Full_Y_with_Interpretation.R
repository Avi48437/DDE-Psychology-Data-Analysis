library(haven)

root_dir <- "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology"
loopr_dir <- file.path(root_dir,"LOPR")
data_dir <- file.path(loopr_dir,"Data")
out_dir <- file.path(loopr_dir,"Analysis","Predictive_Data")

if(!dir.exists(out_dir)) dir.create(out_dir,recursive=TRUE)

master <- read.csv(file.path(data_dir,"LOOPRDataAgeGender.csv"),check.names=FALSE)
S1 <- read_sav(file.path(data_dir,"Survey 1 - Merged - Working data file.sav"))
S2 <- read_sav(file.path(data_dir,"Survey 2 - Merged - Working data file.sav"))

## ============================================================
## Match Survey 1 and Survey 2 to exact LOOPR master order
## ============================================================

master_id <- as.character(master$CaseID)
S1_id <- as.character(S1$CaseID)
S2_id <- as.character(S2$CaseID)

if(anyDuplicated(master_id)) stop("Duplicate CaseID in LOOPR master.")
if(anyDuplicated(S1_id)) stop("Duplicate CaseID in Survey 1.")
if(anyDuplicated(S2_id)) stop("Duplicate CaseID in Survey 2.")

idx1 <- match(master_id,S1_id)
idx2 <- match(master_id,S2_id)

if(any(!is.na(idx1) & !is.na(idx2))) stop("Some CaseID values occur in both surveys.")
if(any(is.na(idx1) & is.na(idx2))) stop("Some LOOPR participants could not be matched.")

in1 <- !is.na(idx1)
in2 <- !is.na(idx2)

cat("Master N:",nrow(master),"\n")
cat("Survey 1:",sum(in1),"\n")
cat("Survey 2:",sum(in2),"\n")

## ============================================================
## Identify variables and pair S1/S2 variables by suffix
## ============================================================

v1 <- setdiff(names(S1),"CaseID")
v2 <- setdiff(names(S2),"CaseID")

suffix1 <- sub("^S1_","",v1)
suffix2 <- sub("^S2_","",v2)

map1 <- setNames(v1,suffix1)
map2 <- setNames(v2,suffix2)

all_suffix <- union(suffix1,suffix2)

## ============================================================
## Exclude variables that are not outcomes
## ============================================================

exclude_exact <- c(
  "Age_1",
  "Sex",
  "HispanicLatino",
  paste0("Race_",1:6)
)

exclude_regex <- c(
  "^Personality",
  "^BFI",
  "^CaseID$"
)

keep <- !(all_suffix %in% exclude_exact)

for(pat in exclude_regex)
  keep <- keep & !grepl(pat,all_suffix,ignore.case=TRUE)

outcomes <- sort(all_suffix[keep])

cat("\nCandidate outcomes:",length(outcomes),"\n")

## ============================================================
## Preserve names used in previous LOOPR analysis
## ============================================================

legacy_names <- c(
  WorkEducation="educ",
  HouseholdIncome="faminc_new",
  MaritalStatus="marstat",
  Spirituality_28="pew_prayer"
)

output_names <- outcomes
ii <- match(outcomes,names(legacy_names))
output_names[!is.na(ii)] <- unname(legacy_names[ii[!is.na(ii)]])

if(anyDuplicated(output_names))
  stop("Duplicate outcome names were created.")

## ============================================================
## Outcome map used internally
## ============================================================

Outcome_Map <- data.frame(
  Outcome=output_names,
  S1Variable=unname(map1[outcomes]),
  S2Variable=unname(map2[outcomes]),
  stringsAsFactors=FALSE
)

Outcome_Map$S1Variable[is.na(Outcome_Map$S1Variable)] <- ""
Outcome_Map$S2Variable[is.na(Outcome_Map$S2Variable)] <- ""

Outcome_Map$Availability <- ifelse(
  Outcome_Map$S1Variable!="" & Outcome_Map$S2Variable!="","S1 and S2",
  ifelse(Outcome_Map$S1Variable!="","S1 only","S2 only")
)

print(Outcome_Map)

## ============================================================
## Helper functions
## ============================================================

get_label <- function(x){
  z <- attr(x,"label")
  if(is.null(z)) "" else as.character(z)
}

labels_to_df <- function(x){
  z <- attr(x,"labels")
  if(is.null(z) || length(z)==0)
    return(data.frame(Code=numeric(0),Meaning=character(0)))
  
  data.frame(
    Code=as.numeric(unname(z)),
    Meaning=names(z),
    stringsAsFactors=FALSE
  )
}

combine_outcome <- function(v1="",v2=""){
  
  x1 <- if(v1!="") S1[[v1]] else NULL
  x2 <- if(v2!="") S2[[v2]] else NULL
  
  is_character <- (!is.null(x1) && is.character(x1)) ||
    (!is.null(x2) && is.character(x2))
  
  if(is_character){
    
    z <- rep(NA_character_,nrow(master))
    
    if(v1!="") z[in1] <- as.character(x1[idx1[in1]])
    if(v2!="") z[in2] <- as.character(x2[idx2[in2]])
    
  } else {
    
    z <- rep(NA_real_,nrow(master))
    
    if(v1!="") z[in1] <- as.numeric(x1[idx1[in1]])
    if(v2!="") z[in2] <- as.numeric(x2[idx2[in2]])
  }
  
  z
}

## ============================================================
## Check S1/S2 coding compatibility
## ============================================================

for(i in seq_len(nrow(Outcome_Map))){
  
  v1i <- Outcome_Map$S1Variable[i]
  v2i <- Outcome_Map$S2Variable[i]
  
  if(v1i=="" || v2i=="") next
  
  L1 <- labels_to_df(S1[[v1i]])
  L2 <- labels_to_df(S2[[v2i]])
  
  if(nrow(L1)==0 || nrow(L2)==0) next
  
  M <- merge(L1,L2,by="Code",suffixes=c("_S1","_S2"))
  
  if(any(M$Meaning_S1!=M$Meaning_S2)){
    cat("\nWARNING: Different S1/S2 value labels for",Outcome_Map$Outcome[i],"\n")
    print(M[M$Meaning_S1!=M$Meaning_S2,,drop=FALSE])
  }
}

## ============================================================
## Construct Y
## ============================================================

Y <- as.data.frame(
  lapply(seq_len(nrow(Outcome_Map)),function(i)
    combine_outcome(
      Outcome_Map$S1Variable[i],
      Outcome_Map$S2Variable[i]
    )
  ),
  check.names=FALSE
)

names(Y) <- Outcome_Map$Outcome

if(nrow(Y)!=nrow(master))
  stop("Y does not have the same number of rows as LOOPR master.")

## ============================================================
## Construct Y_Interpretation
## ============================================================

interpretation_list <- vector("list",nrow(Outcome_Map))

for(i in seq_len(nrow(Outcome_Map))){
  
  out <- Outcome_Map$Outcome[i]
  v1i <- Outcome_Map$S1Variable[i]
  v2i <- Outcome_Map$S2Variable[i]
  avail <- Outcome_Map$Availability[i]
  
  L1 <- if(v1i!="") labels_to_df(S1[[v1i]]) else data.frame()
  L2 <- if(v2i!="") labels_to_df(S2[[v2i]]) else data.frame()
  
  lab1 <- if(v1i!="") get_label(S1[[v1i]]) else ""
  lab2 <- if(v2i!="") get_label(S2[[v2i]]) else ""
  
  labs <- rbind(
    if(nrow(L1)>0) transform(L1,Source="S1") else NULL,
    if(nrow(L2)>0) transform(L2,Source="S2") else NULL
  )
  
  if(!is.null(labs) && nrow(labs)>0){
    
    codes <- sort(unique(labs$Code))
    
    rows <- lapply(codes,function(code){
      
      m1 <- labs$Meaning[labs$Code==code & labs$Source=="S1"]
      m2 <- labs$Meaning[labs$Code==code & labs$Source=="S2"]
      
      m1 <- if(length(m1)==0) "" else m1[1]
      m2 <- if(length(m2)==0) "" else m2[1]
      
      meaning <- if(m1!="" && m2!="" && m1==m2){
        m1
      } else if(m1!="" && m2!=""){
        paste0("S1: ",m1," | S2: ",m2)
      } else if(m1!=""){
        m1
      } else {
        m2
      }
      
      data.frame(
        Outcome=out,
        S1Variable=v1i,
        S2Variable=v2i,
        Availability=avail,
        Code=code,
        Meaning=meaning,
        stringsAsFactors=FALSE
      )
    })
    
    interpretation_list[[i]] <- do.call(rbind,rows)
    
  } else {
    
    meaning <- if(lab1!="" && lab2!="" && lab1==lab2){
      lab1
    } else if(lab1!="" && lab2!=""){
      paste0("S1: ",lab1," | S2: ",lab2)
    } else if(lab1!=""){
      lab1
    } else if(lab2!=""){
      lab2
    } else {
      ""
    }
    
    interpretation_list[[i]] <- data.frame(
      Outcome=out,
      S1Variable=v1i,
      S2Variable=v2i,
      Availability=avail,
      Code=NA_real_,
      Meaning=meaning,
      stringsAsFactors=FALSE
    )
  }
}

Y_Interpretation <- do.call(rbind,interpretation_list)
rownames(Y_Interpretation) <- NULL

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
  Y_Interpretation,
  file.path(out_dir,"Y_Interpretation.csv"),
  row.names=FALSE,
  na=""
)

cat("\nSaved:\n")
cat(file.path(out_dir,"Y.csv"),"\n")
cat(file.path(out_dir,"Y_Interpretation.csv"),"\n")
cat("\nY dimensions:",nrow(Y),"x",ncol(Y),"\n")