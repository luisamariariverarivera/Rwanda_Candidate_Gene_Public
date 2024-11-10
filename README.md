This reposistory contains code to reproduce the analyses in Rivera et al (2024) Prenatal Exposure to the Genocide against the Tutsi in Rwanda is associated with DNA methylation at candidate genes in early adulthood: the role of trauma severity and postnatal adversity

preprint available: https://www.medrxiv.org/content/10.1101/2024.06.12.24308615v1

**Raw Data availability**
Due to restrictions in our consent form, we cannot share raw data files publicly. Requests for raw data solely for the purposes of reproduction should be made to luisa.maria.rivera@gmail.com.

**Files**
Analyses are conducted in two steps. "CandidateGene_Analysis.R" contains the main analysis code
  - We begin by loading annotated DNA methylation M-value data from our set of candidate genes and joining it with phenotypic and exposure data.
  - We run a series of linear regressions and extract significant contrast between exposure groups with and without adjustment for postnatal ACEs.
  - We correct for multiple comparisons using FDR.
  - Estimates for significant probes pre and post correction are saved.

Step 2: PROMIS-29 scores and significant probe DNAm.
- We save methylation values for the probes that survived FDR in step 2 and link them to the PROMIS29 depression and anxiety individual item endorsements.
- After setting weakly regularizing priors, we conduct a Bayesian ordinal regression to assess the relationship between PROMIS29 scores and significant probe methylation, including a random effect for item and participant and fixed effects for Mvalue and sex.
- 

**Dependencies**
library(tidyverse)
library(tidyr)
library(table1)
library(sesame)
library(sesameData)
library(data.table)
library(GenomicRanges)
library(performance)
library(dplyr)
library(factoextra)
library(flextable)
library(broom)
library(officer)
library(lme4)
library(marginaleffects)
library(tidyverse)
library(brms)
library(cubelyr)
library(bayestestR)
library(patchwork)




