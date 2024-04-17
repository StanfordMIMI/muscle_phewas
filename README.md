# Overview
Official implementation for: Abdominal CT Metrics in 17,646 Patients Reveal Associations between Myopenia, Myosteatosis, and Medical Phenotypes: A Phenome-Wide Association Study (eBioMedicine, 2024).

## Paper Website

The [paper website](https://musclephewas.github.io) includes interactive visualizations of the paper's main findings, in addition to a searchable results table.

## Paper Abstract

**Background:** Deep learning facilitates large-scale automated imaging evaluation of body composition. However, associations of body composition biomarkers with medical phenotypes have been underexplored. Phenome-wide association study (PheWAS) techniques search for medical phenotypes associated with biomarkers. A PheWAS integrating large-scale analysis of imaging biomarkers and electronic health record (EHR) data could discover previously unreported associations and validate expected associations. Here we use PheWAS methodology to determine the association of abdominal CT-based skeletal muscle metrics with medical phenotypes in a large North American cohort.

**Methods:** An automated deep learning pipeline was used to measure skeletal muscle index (SMI; biomarker of myopenia) and skeletal muscle density (SMD; biomarker of myosteatosis) from abdominal CT scans of adults between 2012-2018. A PheWAS was performed with logistic regression using patient sex and age as covariates to assess for associations between CT-derived muscle metrics and 611 common EHR-derived medical phenotypes. PheWAS P values were considered significant at a Bonferroni corrected threshold (α=0.05/1222).

**Findings:** 17,646 adults (mean age, 56 years +/- 19 [SD]; 57.5% women) were included. CT-derived SMI was significantly associated with 268 medical phenotypes; SMD with 340 medical phenotypes. Previously unreported associations with the highest magnitude of significance included higher SMI with decreased cardiac dysrhythmias (OR [95%CI], 0.59 [0.55-0.64]; P<0.0001), decreased epilepsy (OR, 0.59 [0.50-0.70]; P<0.0001), and increased elevated prostate-specific antigen (OR, 1.84 [1.47-2.31]; P<0.0001), and higher SMD with decreased decubitus ulcers (OR, 0.36 [0.31-0.42]; P<0.0001), sleep disorders (OR, 0.39 [0.32-0.47]; P<0.0001), and osteomyelitis (OR, 0.43 [0.36-0.52]; P<0.0001).

**Interpretation:** PheWAS methodology reveals previously unreported associations between CT-derived biomarkers of myopenia and myosteatosis and EHR medical phenotypes. The high-throughput PheWAS technique applied on a population scale can generate research hypotheses related to myopenia and myosteatosis and can be adapted to research possible associations of other imaging biomarkers with hundreds of EHR medical phenotypes.

## Citation

If you find this code or the findings of the paper useful, please cite:

```
@article{
zambranochaves2024musclephewas,
title = {Abdominal CT metrics in 17,646 patients reveal associations between myopenia, myosteatosis, and medical phenotypes: a phenome-wide association study},
journal = {eBioMedicine},
volume = {103},
pages = {105116},
year = {2024},
issn = {2352-3964},
doi = {https://doi.org/10.1016/j.ebiom.2024.105116},
url = {https://www.sciencedirect.com/science/article/pii/S2352396424001518},
author = {Juan M. {Zambrano Chaves} and Leon Lenchik and Isabel O. Gallegos and Louis Blankemeier and Daniel L. Rubin and Marc H. Willis and Akshay S. Chaudhari and Robert D. Boutin},
keywords = {Myopenia, Myosteatosis, Phenome wide association study, Computed tomography},
}
```
