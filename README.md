# TUBA4A – In silico, mechanistic, and genotype–phenotype analyses

This repository contains all reproducible in silico and mechanistic analyses performed for the **Johari et al. Missense variants in TUBA4A cause myo-tubulinopathies**, including:

- **FoldX** ΔΔG protein stability predictions  
- **AggreScan4D** aggregation propensity analysis  
- **Electron microscopy (EM)** morphology quantification  
- **Genotype–phenotype correlations** (domain × clinical features)  
- **Integrated mechanistic correlations** (ΔΔG, aggregation, microtubule morphology)  

All scripts are written in **R** and rely solely on open-source packages.  
The repository is fully self-contained and follows a transparent, analysis-first structure.

## 📊 Description of Datasets

### **1. FoldX datasets (`data/foldx/`)**
- `TUBA4A_FoldX_myopathy.txt`  
- `TUBA4A_FoldX_nonmyopathy.txt`  
Contain ΔΔG FoldX predictions for myopathy-associated and non-myopathy (ALS, FTD, Ataxia) variants.

### **2. Aggrescan4D datasets (`data/aggrescan4D/`)**
- `TUBA4A-aggregation4D_scores.txt`
- `variant-labels.txt`
Contain normalised aggregation propensity scores per mutant protein per amino acid position and the variant labels mapping mutants to phenotypes.

### **3. EM morphology dataset (`data/EM`)**
- `TUBA4A-EM.txt`
Per-microtubule diameter measurements for patient and control muscle biopsies.

### **4. Clinical genotype–phenotype dataset (`data/genotype_phenotype/`)**
- `TUBA4A-clinical-correlation.txt` 
Contains:
  - **MRC-scale limb strength**  
    DUL, PUL, DLL, PLL  
    (5 = normal, 1 = severe)

  - **Binary features**  
    Axial weakness, Facial weakness, Ptosis, Respiratory involvement  
    (1 = present, 0 = absent)

  - **CK elevation**  
    (0 = normal, 1 = mildly elevated, 2 = high)

  - **Domain**  
    (1 = GTPase, 2 = C-terminal, 3 = GTP-binding)

  *Severity score (mean MRC) is computed by the script, not part of input data.*

### **5. Mechanistic dataset (`data/mechanistic_correlation/`)**
- `TUBA4A-mechanistic-correlation.txt` 
Contains for each variant:
  - FoldX ΔΔG z-score  
  - Aggrescan4D aggregation score  
  - EM morphology (% normal / mild / abnormal)  
  - Domain assignment  
  - Derived morphology severity (0–2 scale)

---

## 🧪 Scripts

Each script is self-contained, thoroughly commented, and outputs to `results/`.

---

## 🚀 Running the Analyses

From the repository root:

```bash
# Clinical genotype–phenotype correlations
Rscript scripts/genotype_phenotype/20251118-correlation-geno-pheno.r

# Mechanistic correlations
Rscript scripts/mech_correlation/20251118-correlation-mechanistic.r

# Other modules (EM, FoldX, Aggrescan4D)
Rscript scripts/EM/20251118-EM.r
Rscript scripts/foldx/20251118-FoldX.r
Rscript scripts/aggrescan4D/20251118-aggresscan4D.r

Outputs appear in:
	•	results/figures/
	•	results/tables/

This repository accompanies the TUBA4A manuscript and is maintained by
Mridul Johari (2025).