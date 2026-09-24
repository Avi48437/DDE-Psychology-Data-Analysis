# DDE Psychology Data Analysis

This repository contains the code, post-processing utilities, and analysis workflow used to fit and interpret **Deep Discrete Encoder (DDE)** models on several personality datasets.

The project is organized around a common workflow:

```text
Raw / prepared data
        ↓
DDE fitting
        ↓
Post-processing
        ↓
Canonical loading matrices and maps
        ↓
Predictive / outcome analysis
        ↓
Figures and summaries
```

Large fitted objects and raw data files are intentionally **not stored in GitHub**. They are excluded through `.gitignore` and can be distributed separately through Box.

---

## Repository Structure

```text
DDE_Fitting_Pshycology2/
│
├── Algorithms/
│   ├── get_SAEM_RL.m
│   ├── get_SAEM_RL_CSP.m
│   ├── get_SAEM_RL_CSP_D.m
│   ├── get_SAEM_RL_CSP_Gibbs.m
│   └── get_SAEM_RL_G.m
│
├── Utilities/
│   ├── fitting / simulation utilities
│   ├── loading-matrix alignment and thresholding
│   ├── plotting functions
│   ├── anchor-map functions
│   └── predictive / synthetic-data utilities
│
├── Scripts/
│   ├── B5.m
│   ├── Deeper_Sims.m
│   ├── Sim_Overspecified.m
│   └── Sim_SampleSplit.m
│
├── 1. IPIP-FFM-data/
│   ├── Data/
│   │   ├── data-final.csv                  [Box / ignored by Git]
│   │   ├── B5_Red.mat                     [Box / ignored by Git]
│   │   └── B5_Red_questions.csv           [Box / ignored by Git]
│   │
│   ├── Fitting/
│   │   ├── Scripts/
│   │   │   └── IPIPFMM_Tuning.m
│   │   └── Results/
│   │       └── IPIP_DDE_100_results6_comp.mat
│   │           [Box / ignored by Git]
│   │
│   ├── PostProcessing/
│   │   └── PostProcessingC1.m
│   │
│   ├── Analysis/
│   │   ├── PostProcessed_Data/
│   │   │   ├── IPIPFMM_CC1.mat            [Box / ignored by Git]
│   │   │   ├── IPIPFMM_CC1_Item_Map.csv
│   │   │   └── IPIPFMM_CC1_Factor_Map.csv
│   │   └── Scripts/
│   │       └── Post_Fit_Analysis.m
│   │
│   ├── Helpers/
│   │   ├── loading-matrix matching utilities
│   │   ├── key-item extraction utilities
│   │   ├── Layer-2 interpretation utilities
│   │   └── combination summaries
│   │
│   └── References/
│       └── codebook.txt
│
├── 2. IPIP 100/
│   ├── Data/
│   │   ├── B5.csv                          [Box / ignored by Git]
│   │   └── B5_wo.csv                       [Box / ignored by Git]
│   │
│   ├── Fitting/
│   │   ├── Scripts/
│   │   │   ├── DDE_Fit.m
│   │   │   ├── IPIP100_Tuning.m
│   │   │   └── IPIP_100_Permute_Rows.m
│   │   └── Results/
│   │       ├── IPIP100_DDE_fit.mat         [Box / ignored by Git]
│   │       └── B5_analysis.mat             [Box / ignored by Git]
│   │
│   ├── PostProcessing/
│   │   ├── Postprocessing_Fit.m
│   │   ├── PostProcessCPDDE98.m
│   │   └── code_map.m
│   │
│   ├── Analysis/
│   │   ├── PostProcessed_Data/
│   │   │   ├── IPIP100_PostProcessed.mat  [Box / ignored by Git]
│   │   │   ├── IPIP100_Item_Map.csv
│   │   │   ├── IPIP100_Factor_Map.csv
│   │   │   ├── IPIP98_PostProcessed.mat   [Box / ignored by Git]
│   │   │   ├── IPIP98_Item_Map.csv
│   │   │   └── IPIP98_Factor_Map.csv
│   │   ├── Data_Preparation/
│   │   │   ├── Save_Demography.R
│   │   │   └── Save_XAY.m
│   │   ├── Predictive_Data/
│   │   │   ├── X100.csv
│   │   │   ├── X98.csv
│   │   │   ├── A100.csv
│   │   │   ├── A98.csv
│   │   │   ├── rawX100.csv
│   │   │   ├── rawX98.csv
│   │   │   ├── Y.csv
│   │   │   └── D_IPIP98.csv
│   │   ├── Scripts/
│   │   │   ├── Prediction_Classifiers.R
│   │   │   ├── Educ_Analysis98.R
│   │   │   ├── Faminc_Analysis98.R
│   │   │   ├── Marstat_Analysis98.R
│   │   │   ├── Prayer_Analysis98.R
│   │   │   ├── Votereg_Analysis98.R
│   │   │   ├── Cart_Fit_Wd_Demo.R
│   │   │   ├── Cart_Fit_Wd_Demo_Summary.R
│   │   │   ├── Data_Load_Viz.R
│   │   │   ├── Viz_Outcome_A1_A2_X.R
│   │   │   ├── Summary_98.R
│   │   │   └── Post_Fit_Analysis.m
│   │   ├── Results/
│   │   │   └── *.rds                      [Box / ignored by Git]
│   │   └── Figures/
│   │       └── predictive-analysis figures
│   │
│   └── References/
│       └── WASU0003_B5_codebook.pdf
│
├── LOPR/
│   ├── Data/
│   │   ├── Survey 1 - Merged - Working data file.sav
│   │   │   [Box / ignored by Git]
│   │   ├── Survey 2 - Merged - Working data file.sav
│   │   │   [Box / ignored by Git]
│   │   └── LOOPRDataAgeGender.csv
│   │
│   ├── Fitting/
│   │   ├── Scripts/
│   │   │   └── LOPR_Tuning_and_Fit.m
│   │   └── Results/
│   │       ├── LOOPR_DDE_fit.mat
│   │       ├── LOOPR_best_DDE_tuning.mat
│   │       ├── LOOPR_DDE_tuning_checkpoint.mat
│   │       ├── LOOPR_DDE_pMSE_tuning_results.mat
│   │       └── LOOPR_DDE_pMSE_tuning_results.csv
│   │           [.mat files ignored by Git]
│   │
│   ├── PostProcessing/
│   │   ├── Item_Map.R
│   │   └── PostProcessing.m
│   │
│   ├── Analysis/
│   │   ├── PostProcessed_Data/
│   │   │   ├── LOOPR_PostProcessed.mat     [Box / ignored by Git]
│   │   │   ├── LOOPR_BFI2_Item_Text.csv
│   │   │   ├── LOOPR_Item_Map.csv
│   │   │   └── LOOPR_Factor_Map.csv
│   │   ├── Data Preparation/
│   │   │   ├── Save_A.m
│   │   │   ├── Save_Demography.R
│   │   │   ├── SaveYX_rc.R
│   │   │   └── Save_Full_Y_with_Interpretation.R
│   │   ├── Predictive_Data/
│   │   │   ├── X.csv
│   │   │   ├── X_rc.csv
│   │   │   ├── A.csv
│   │   │   ├── Y.csv
│   │   │   ├── Y_Interpretation.csv
│   │   │   └── D_LOOPR.csv
│   │   ├── Scripts/
│   │   │   ├── Prediction_Classifiers.R
│   │   │   ├── Educ_AnalysisLOOPR.R
│   │   │   ├── Faminc_AnalysisLOOPR.R
│   │   │   ├── Marstat_AnalysisLOOPR.R
│   │   │   ├── Prayer_AnalysisLOOPR.R
│   │   │   ├── Cart_Fit_Wd_Demo.R
│   │   │   ├── Cart_Fit_Wd_Demo_Summary.R
│   │   │   ├── Data_Load_Viz.R
│   │   │   ├── Summary_LOOPR.R
│   │   │   ├── Helpers.R
│   │   │   └── LOPR_Loading_Plot.m
│   │   ├── Results/
│   │   │   └── *.rds                      [Box / ignored by Git]
│   │   └── Figures/
│   │       └── predictive-analysis figures
│   │
│   └── References/
│       ├── Soto_John_2017.pdf
│       └── bfi2-form.pdf
│
├── References/
│   └── supporting personality / life-outcome papers
│
├── General_Analysis.m
├── Correlation_Plots.png
├── Loading Matrices.png
├── IPIP-98 Anchor Map.pdf
├── IPIP-FMM Anchor Map.pdf
├── LOOPR BFI-2 Anchor Map.pdf
└── .gitignore
```

---

## Main Components

### `Algorithms/`

Core MATLAB implementations of the DDE fitting procedures, including rank-likelihood and cumulative-shrinkage variants.

### `Utilities/`

Shared functions used across datasets for loading-matrix plotting, factor alignment, thresholding, anchor-map construction, simulation, synthetic-data generation, and pMSE-related calculations.

These functions are kept outside the dataset-specific folders because they are reused across analyses.

### `1. IPIP-FFM-data/`

Analysis of the 50-item IPIP Five-Factor Model questionnaire.

```text
IPIPFMM_Tuning.m
        ↓
IPIP_DDE_100_results6_comp.mat
        ↓
PostProcessingC1.m
        ↓
IPIPFMM_CC1
        ↓
Post_Fit_Analysis.m
```

The post-processing step constructs a canonical representation of the learned hierarchy, including factor ordering, item ordering, item maps, and factor maps.

### `2. IPIP 100/`

Contains the IPIP-100 DDE analysis together with the IPIP-98 comparison and predictive-analysis workflow.

This directory contains DDE fitting, post-processing, canonical mappings, predictive-feature preparation, classification / regression analyses, demographic CART analyses, and outcome-specific figures.

### `LOPR/`

Contains the DDE analysis of the LOOPR BFI-2 data.

```text
raw survey data
      ↓
DDE tuning and fitting
      ↓
post-processing
      ↓
predictive-data construction
      ↓
outcome prediction
      ↓
figures and summaries
```

### `General_Analysis.m`

Top-level comparison script used after the individual datasets have been post-processed.

It loads the canonical DDE representations from IPIP-FMM, IPIP-100 / IPIP-98, and LOOPR and compares the learned hierarchical loading structures and anchor maps.

---

## Data and Large Result Files

Large data and fitted-result files are intentionally not committed to GitHub.

The repository `.gitignore` excludes, among other things:

```text
*.mat
*.rds
*.sav
1. IPIP-FFM-data/Data/data-final.csv
1. IPIP-FFM-data/Data/B5_Red_questions.csv
2. IPIP 100/Data/B5.csv
2. IPIP 100/Data/B5_wo.csv
```

These files can be distributed separately through **Box**.

> **Box download link:** _to be added_

After downloading the external files, place them in the corresponding directories shown in the repository tree.

---

## Software

The project primarily uses:

- **MATLAB** for DDE fitting, simulation, post-processing, and loading-matrix analysis.
- **R** for predictive modeling, demographic analyses, visualization, and outcome-specific summaries.

---

## Repository Organization Principle

> **Dataset-specific code stays with the dataset; reusable DDE methodology stays at the project root.**

This separation keeps each empirical workflow reproducible while maintaining a common implementation of the underlying DDE methodology.

---

## Current Datasets

| Dataset | Personality instrument | Main role |
|---|---|---|
| IPIP-FFM | 50-item IPIP Five-Factor Model | Hierarchical DDE analysis |
| IPIP-100 / IPIP-98 | Extended IPIP personality items | Cross-dataset comparison and predictive analysis |
| LOOPR | BFI-2 | Independent personality dataset and predictive validation |

---

## Status

This repository is under active development. File organization, analysis scripts, and documentation may continue to evolve as the DDE methodology and empirical analyses are refined.
