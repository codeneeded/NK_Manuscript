# NK Manuscript — Code & Data Repository

> Complete analysis code and processed data for the publication: **A role for Fas/FasL-dependent NK-mediated cytotoxicity in infants with perinatal HIV** (*iScience*, 2025).

---

## Publication

**A role for Fas/FasL-dependent NK-mediated cytotoxicity in infants with perinatal HIV**
*iScience* (2025). https://www.cell.com/iscience/fulltext/S2589-0042(25)01898-X

---

## Overview

This repository contains all R analysis scripts, processed datasets, and output figures supporting the manuscript above. The study characterizes NK cell cytotoxic function and phenotypic subset composition in perinatally HIV-infected (HEI) and HIV-exposed uninfected (HEU) infants from the **TARA (Towards AIDS Remission Approaches)** cohort and additional cohorts, with a focus on identifying **Fas/FasL-mediated cytotoxicity** as a key mechanism by which infant NK cells kill HIV-infected targets.

Key analyses include unsupervised FlowSOM clustering of raw FCS files, linear mixed-effects modeling of NK subset predictors of specific killing, correlation analyses, and characterization of CMV serostatus and KIR profiles as modulators of NK phenotype and function.

---

## Repository Structure

```
NK_Manuscript/
│
├── R_Scripts/                          # All R analysis scripts
├── FLOWSOM/                            # FlowSOM unsupervised clustering outputs
├── FCS/                                # FCS files for FlowSOM analysis
│
├── Boxplots/                           # NK subset boxplots (HEI vs HEU, by timepoint)
├── Correlations/                       # Correlation analysis outputs
├── Corrplots/                          # Correlation plots (NK subsets vs. viral load)
├── Paired_Plots/TARA/                  # Longitudinal paired plots for the TARA cohort
├── Mixed_Effects_Model/                # Mixed-effects model coefficient plots
│
├── R_dat/                              # Intermediate R data files
├── Saved_R_Data/                       # Saved R objects (.RData) for reproducibility
│
├── NK Function MetaDATA_REVISED0724v2.xlsx   # Primary flow cytometry + metadata
├── NK Function MetaDATA_REVISED0724v2.csv    # CSV version
├── NK Function Metadata.xlsx                 # Earlier metadata version
├── NK_data_cleaned_normalized_to_CD45.csv              # Full cleaned NK data (CD45-normalized)
├── NK_Cell_data_cleaned_normalized_to_CD45_Lesley Reduced 061124.csv  # Reduced final dataset
├── Viral_Titres.csv                    # Participant HIV viral load by timepoint
├── CMV_KIR.xlsx                        # CMV serostatus and KIR data per participant
│
├── NK Manuscript Figures_Fig4.pptx     # Manuscript Figure 4
├── Statistical_Methods.docx            # Written statistical methods from the manuscript
│
├── test_1.png / test_2.png / test3.png # Development figure drafts
├── NK_Manuscript.Rproj                 # RStudio project file
├── .gitattributes / .gitignore
└── README.md
```

---

## Analyses

### 1. FlowSOM Unsupervised Clustering (`FLOWSOM/`, `FCS/`)

Raw FCS files are processed through FlowSOM self-organizing map clustering to identify NK cell phenotypic clusters in an unbiased, data-driven manner. This complements the manual gating approach and provides a high-dimensional view of NK subset structure across HEI and HEU infants.

### 2. NK Phenotype vs. Specific Killing — Mixed-Effects Models (`Mixed_Effects_Model/`, `R_Scripts/`)

Linear mixed-effects models fit per-NK-subset, with specific killing as the dependent variable:

```r
`Specific Killing` ~ `viral load` + Timepoint + gender + HIV + (1 | PID) + <NK_subset>
```

Separate model sets are run for two target cell lines:
- **HUT78/SF2** — HIV-infected T cell line (HIV-specific killing)
- **K562** — MHC-I-deficient tumor cell line (innate/missing-self killing)

A key analytical focus is **CD56dimCD16+/FasL**, which is always included in results to formally test the Fas/FasL killing hypothesis across conditions.

### 3. Longitudinal Paired Plots (`Paired_Plots/TARA/`)

Paired visualizations tracking individual participants' NK killing and subset frequencies across Entry (~1–2 months, pre-ART) and 12-month timepoints, faceted by HIV status (HEI/HEU) and treatment condition.

### 4. Boxplots — HEI vs. HEU Subset Comparisons (`Boxplots/`)

Group-level comparisons of NK subset frequencies between HEI and HEU infants at each timepoint, presented as boxplots with significance testing.

### 5. Correlation Analyses (`Correlations/`, `Corrplots/`)

Spearman or Pearson correlations between NK cell subset frequencies, functional markers, FasL expression, viral load, and specific killing outcomes. Corrplots provide matrix-level visualizations of inter-subset relationships.

### 6. CMV & KIR Analysis (`CMV_KIR.xlsx`, `R_Scripts/`)

CMV (cytomegalovirus) serostatus and KIR (killer immunoglobulin-like receptor) profiles are analyzed as covariates and modulators of NK cell phenotype and Fas/FasL-mediated killing. CMV-driven adaptive NK cell expansion is a well-established driver of elevated FasL expression, making this directly relevant to the paper's central finding.

---

## Key Data Files

| File | Description |
|---|---|
| `NK Function MetaDATA_REVISED0724v2.xlsx` | Primary dataset: flow cytometry (MFI + frequencies), specific killing, demographics, cohort, batch, timepoint |
| `NK_data_cleaned_normalized_to_CD45.csv` | Full NK panel normalized to CD45+ live lymphocytes |
| `NK_Cell_data_cleaned_normalized_to_CD45_Lesley Reduced 061124.csv` | Final reduced dataset used in manuscript analyses |
| `Viral_Titres.csv` | HIV viral load (copies/mL) per participant per timepoint |
| `CMV_KIR.xlsx` | CMV serostatus and KIR typing/expression data |
| `Statistical_Methods.docx` | Statistical methods section from the manuscript |
| `NK Manuscript Figures_Fig4.pptx` | Manuscript Figure 4 source file |

---

## Cohorts

| Cohort | Description |
|---|---|
| **TARA** | HIV-exposed infected (HEI) and uninfected (HEU) infants from Maputo, Mozambique. Primary cohort. |
| **Florah** | Additional infant cohort for comparative analysis |
| **PAVE KL** | Additional cohort for comparative analysis |

---

## NK Cell Killing Assay Design

NK cytotoxicity was measured as **specific killing (%)** at a defined effector-to-target ratio against:

| Target | Cell Line | Killing Mechanism Tested |
|---|---|---|
| **HUT78/SF2** | HIV-infected T cell line | HIV-specific cytotoxicity (Fas/FasL, NKG2D, NKp30) |
| **K562** | MHC-I-deficient tumor line | Innate missing-self cytotoxicity (NKG2D, NCRs) |

Functional NK markers measured include: FasL (CD178), IFNγ, MIP-1β (CCL4), CD107a (degranulation marker), Granzyme B, Perforin, and inhibitory receptors NKG2A and TIGIT.

---

## Dependencies

All scripts are written in **R**. Required packages:

| Package | Purpose |
|---|---|
| `FlowSOM` / `flowCore` | Unsupervised NK clustering from FCS files |
| `readxl` / `readr` | Data ingestion |
| `dplyr` / `tidyr` / `reshape2` | Data wrangling |
| `lme4` / `lmerTest` | Linear mixed-effects models |
| `broom.mixed` | Tidy model output extraction |
| `limma` / `edgeR` | Differential abundance (voom) |
| `ggplot2` / `ggrepel` / `ggpubr` | Visualization |
| `pheatmap` / `corrplot` | Heatmaps and correlation matrices |
| `extrafont` | Custom fonts for publication figures |

Install core packages via:

```r
install.packages(c("readxl", "readr", "dplyr", "tidyr", "reshape2",
                   "lme4", "lmerTest", "broom.mixed", "ggplot2",
                   "ggrepel", "ggpubr", "pheatmap", "corrplot", "extrafont"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("FlowSOM", "flowCore", "limma", "edgeR"))
```

---

## Usage

1. Clone the repository and open `NK_Manuscript.Rproj` in RStudio.
2. Update file paths within each script in `R_Scripts/` to match your local environment.
3. For FlowSOM analyses, ensure FCS files from `FCS/` are accessible.
4. Saved R objects in `Saved_R_Data/` can be loaded directly to skip reprocessing steps.

> ⚠️ **Note:** Participant-level identifiable data is not included. Raw FCS files and full metadata are available upon reasonable request from the corresponding authors.

---

## Citation

If you use this code or data, please cite:

> **A role for Fas/FasL-dependent NK-mediated cytotoxicity in infants with perinatal HIV.**
> *iScience* (2025). https://www.cell.com/iscience/fulltext/S2589-0042(25)01898-X

This work was supported by NIH/NIAID grant R01 AI127347 (S. Pahwa) and conducted using samples from the TARA cohort, Maputo, Mozambique.

---

## Contact

For questions about the code or data, please open an [issue](https://github.com/codeneeded/NK_Manuscript/issues) in this repository.
