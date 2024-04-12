library(tidyverse)
library(janitor)
library(lubridate)
library(PheWAS)
library(lmtest)


#0. get studies to include, load in segmentation data
write_dir <- '/path/to/results_dir/'
clindata_dir <- '/path/to/clinical_data/'
pop_dir <- '/path/to/population_data/'

study_data_sarcopenia <- read_csv(paste(pop_dir, 'study_data_sarcopenia_abstract.csv', sep=''))

anon_ids <- study_data_sarcopenia %>%
  select(patient_id) %>%
  pull

seg_data <- read_csv(paste(pop_dir, 'population_raw.csv', sep='')) %>%
  mutate(anon_id = str_replace(substr(filename, 1, 19), '-','_'))


#PheWAS study

#1. get diagnoses for phenotyping
dx_in_final <- function(x, pos) subset(x, `Patient Id` %in% anon_ids)

dx <- read_csv_chunked(paste(clindata_dir, '1/diagnoses.csv', sep=''), DataFrameCallback$new(dx_in_final), col_types='dcdcccccccd') %>%
  bind_rows(read_csv_chunked(paste(clindata_dir, '2/diagnoses.csv', sep=''), DataFrameCallback$new(dx_in_final), col_types='dcdcccccccd')) %>%
  bind_rows(read_csv_chunked(paste(clindata_dir, '3/diagnoses.csv', sep=''), DataFrameCallback$new(dx_in_final), col_types='dcdcccccccd'))

#1.1 window and aggregate code counts

dx_count_allcodes <- dx %>%
  clean_names %>%
  select(patient_id, age, icd9_code, icd10_code) %>%
  mutate(code_type=if_else(is.na(icd9_code) | str_detect(icd9_code, ','), "ICD10CM","ICD9CM"),
         code=if_else(is.na(icd9_code) | str_detect(icd9_code, ','), icd10_code, icd9_code))

#1.2 Get medical phenotypes (phecodes)

phecodes <- dx_count_allcodes %>% 
  left_join(select(study_data_sarcopenia, patient_id, age_at_scan), by=c("patient_id"="patient_id")) %>%
  filter(age < age_at_scan +3/12 & age>age_at_scan-1/12) %>%
  select(patient_id, code_type, code) %>%
  rename(vocabulary_id=code_type) %>%
  mapCodesToPhecodes(make.distinct=FALSE) %>%
  filter(phecode %in% (pheinfo %>% select(phecode) %>% pull)) %>%
  group_by(patient_id, phecode) %>% 
  summarize(n=n()) %>% 
  ungroup

custom_agg <- function(index){
  s <- sum(index)
  if (s >1){
    return(s)
  } else if (s == 0){
    return(s)
  } else {
    return(NA)
  }
}

study_data_phenotypes <- createPhenotypes(phecodes %>% mutate(vocabulary_id='phecode') %>% rename(code=phecode) %>% select(patient_id, vocabulary_id, code, n) %>% as.data.frame,
                                          id.sex=study_data_sarcopenia %>% mutate(sex=if_else(sex=='Female',"F","M")) %>% select(patient_id, sex) %>% as.data.frame,
                                          full.population.ids = study_data_sarcopenia %>% select(patient_id) %>% pull,
                                          aggregate.fun=custom_agg,
                                          translate=FALSE)
#2. get genotypes (SMI and SMD)
mean_sds <- study_data_sarcopenia %>% 
  select(sex, recent_height_cm, cross_sectional_area_muscle, hounsfield_unit_muscle) %>% #csa_muscle in mm^2
  mutate(recent_height_cm = as.numeric(recent_height_cm),
         l3smi = (cross_sectional_area_muscle/100)/((recent_height_cm/100)**2)) %>%
  group_by(sex) %>%
  summarize(avg_csi_muscle = mean(l3smi),
            sd_csi_muscle = sd(l3smi),
            avg_hu_muscle = mean(hounsfield_unit_muscle, na.rm = TRUE),
            sd_hu_muscle = sd(hounsfield_unit_muscle, na.rm = TRUE))
study_data_scores <- select(study_data_sarcopenia, patient_id, anon_id, sex, recent_height_cm) %>%
  left_join(seg_data, by="anon_id") %>%
  mutate(id = patient_id,
         recent_height_cm = as.numeric(recent_height_cm),
         l3smi = (cross_sectional_area_muscle/100)/((recent_height_cm/100)**2),
         smi_zscore = if_else(sex=='Female', (l3smi-mean_sds$avg_csi_muscle[1])/mean_sds$sd_csi_muscle[1],(l3smi-mean_sds$avg_csi_muscle[2])/mean_sds$sd_csi_muscle[2]),
         sm_hu_zscore = if_else(sex=='Female', (hounsfield_unit_muscle-mean_sds$avg_hu_muscle[1])/mean_sds$sd_hu_muscle[1],(hounsfield_unit_muscle-mean_sds$avg_hu_muscle[2])/mean_sds$sd_hu_muscle[2])) %>%
  select(id, smi_zscore, sm_hu_zscore)

#3. get covariates used for regression
study_data_covariates <- select(study_data_sarcopenia, patient_id, age_at_scan, sex) %>%
  mutate(id=patient_id,
         age = as.numeric(age_at_scan), 
         sex=if_else(sex=='Male', 1,0)) %>%
  select(id, age, sex) %>%
  as.data.frame

#4. PheWAS analysis
sm_phewas <-  phewas(phenotypes = study_data_phenotypes %>%rename(id=patient_id), 
                     genotypes=study_data_scores,
                     covariates=study_data_covariates,
                     additive.genotypes=FALSE,
                     significance.threshold = "bonferroni",
                     min.records=100, 
                     cores=4, 
                     MASS.confint.level =.95, 
                     return.models=TRUE)

all_results <- sm_phewas$results %>%
  select(phenotype, snp, p, OR, lower, upper, n_total, n_cases, n_controls) %>%
  left_join(PheWAS::pheinfo , by=c("phenotype"="phecode"))

studied_associations <- all_results %>%
  filter(n_cases >=100, n_controls >=100)

# 5. Evaluate nonlinear relationships
rm(dx, dx_count_allcodes, phecodes, seg_data, study_data_sarcopenia)
study_data_squared <- study_data_scores %>%
  mutate(sm_hu_zscore_sq = sm_hu_zscore**2,
         smi_zscore_sq = smi_zscore**2)

sm_phewas_xsquared_smd <- phewas(phenotypes = study_data_phenotypes %>%rename(id=patient_id), 
                                 genotypes=study_data_squared%>% select(id, sm_hu_zscore),
                                 covariates=study_data_covariates %>% left_join(select(study_data_squared, c(id, sm_hu_zscore_sq)), by="id"),
                                 additive.genotypes=FALSE,
                                 significance.threshold = "bonferroni",
                                 min.records=100,
                                 cores=4,
                                 MASS.confint.level =.95,
                                 return.models=TRUE)

significant_coeff_sq_smd <- sm_phewas_xsquared_smd$results %>% 
  filter(snp=="sm_hu_zscore_sq") %>%
  filter(p<0.05/1222) %>%
  select(phenotype) # n = 81

sm_phewas_xsquared_smi <- phewas(phenotypes = study_data_phenotypes %>%rename(id=patient_id), 
                                 genotypes=study_data_squared%>% select(id, smi_zscore),
                                 covariates=study_data_covariates %>% left_join(select(study_data_squared, c(id, smi_zscore_sq)), by="id"),
                                 additive.genotypes=FALSE,
                                 significance.threshold = "bonferroni",
                                 min.records=100,
                                 cores=4,
                                 MASS.confint.level =.95,
                                 return.models=TRUE)

significant_coeff_sq_smi <- sm_phewas_xsquared_smi$results %>% 
  filter(snp=="smi_zscore_sq") %>%
  filter(p<0.05/1222) %>%
  select(phenotype) # n = 29

# compare model types: linear vs nonlinear using likelihood ratio test
smd_lm_names <- grep("^NA", sm_phewas$models %>% names, value=TRUE) %>%
  as_tibble %>%
  filter(str_detect(value, 'sm_hu_zscore'))

smi_lm_names <- grep("^NA", sm_phewas$models %>% names, value=TRUE) %>%
  as_tibble %>%
  filter(str_detect(value, 'smi_zscore'))

nonlin_model_preferred_smd <- c()
was_not_now_is_not_significant_smd <- c()
was_not_now_is_significant_lin_only_smd <- c()
was_not_now_is_significant_sq_only_smd <- c()
was_not_now_is_significant_lin_sq_smd<- c()

was_now_is_not_significant_smd <- c()
was_now_is_significant_lin_only_smd <- c()
was_now_is_significant_sq_only_smd <- c()
was_now_is_significant_lin_sq_smd <- c()

for (i in 1:nrow(significant_coeff_sq_smd)) {
  # Accessing each row as a tibble
  phenotype_name <- significant_coeff_sq_smd[i, ]
  phenotype_name_space <- paste(phenotype_name, " ", sep="")
  
  og_model_name <- smd_lm_names %>% filter(str_detect(value, phenotype_name_space)) %>% pull
  new_model_name <- gsub("sex NA", "sex + sm_hu_zscore_sq NA", og_model_name)
  
  model_pval <- lrtest(sm_phewas$models[[og_model_name]], sm_phewas_xsquared_smd$models[[new_model_name]]) %>%
    as_tibble %>%
    slice(2) %>%
    pull('Pr(>Chisq)')
  
  old_pval <- sm_phewas$results %>%
    filter(phenotype == phenotype_name, snp=='sm_hu_zscore') %>%
    select(p) %>%
    pull
  new_pval_lin <- sm_phewas_xsquared_smd$results %>%
    filter(phenotype == phenotype_name, snp=='sm_hu_zscore') %>%
    select(p) %>%
    pull
  new_pval_sq <- sm_phewas_xsquared_smd$results %>%
    filter(phenotype == phenotype_name, snp=='sm_hu_zscore_sq') %>%
    select(p) %>%
    pull
  
  
  # Example action: print the row
  if (model_pval < 0.05/81) {
    nonlin_model_preferred_smd <- c(nonlin_model_preferred_smd, phenotype_name)
    
    if (old_pval > 0.05/1222){
      if (new_pval_lin < 0.05/1384) {
        if (new_pval_sq < 0.05/1384) {
          was_not_now_is_significant_lin_sq_smd <- c(was_not_now_is_significant_lin_sq_smd, phenotype_name)
        }
        else {
          was_not_now_is_significant_lin_only_smd <- c(was_not_now_is_significant_lin_only_smd, phenotype_name)
        }
      }
      else {
        if (new_pval_sq < 0.05/1384) {
          was_not_now_is_significant_sq_only_smd <- c(was_not_now_is_significant_sq_only_smd, phenotype_name)
        }
        else {
          was_not_now_is_not_significant_smd <- c(was_not_now_is_not_significant_smd, phenotype_name)
        }
      }
    }
    else {
      if (new_pval_lin < 0.05/1384) {
        if (new_pval_sq < 0.05/1384) {
          was_now_is_significant_lin_sq_smd <- c(was_now_is_significant_lin_sq_smd, phenotype_name)
        }
        else {
          was_now_is_significant_lin_only_smd <- c(was_now_is_significant_lin_only_smd, phenotype_name)
        }
      }
      else {
        if (new_pval_sq < 0.05/1384) {
          was_now_is_significant_sq_only_smd <- c(was_now_is_significant_sq_only_smd, phenotype_name)
        }
        else {
          was_now_is_not_significant_smd <- c(was_now_is_not_significant_smd, phenotype_name)
        }
      }
    }
  }
}

nonlin_model_preferred_smi <- c()
was_not_now_is_not_significant_smi <- c()
was_not_now_is_significant_lin_only_smi <- c()
was_not_now_is_significant_sq_only_smi <- c()
was_not_now_is_significant_lin_sq_smi <- c()

was_now_is_not_significant_smi <- c()
was_now_is_significant_lin_only_smi <- c()
was_now_is_significant_sq_only_smi <- c()
was_now_is_significant_lin_sq_smi <- c()

for (i in 1:nrow(significant_coeff_sq_smi)) {
  # Accessing each row as a tibble
  phenotype_name <- significant_coeff_sq_smi[i, ]
  phenotype_name_space <- paste(phenotype_name, " ", sep="")
  
  og_model_name <- smi_lm_names %>% filter(str_detect(value, phenotype_name_space)) %>% pull
  new_model_name <- gsub("sex NA", "sex + smi_zscore_sq NA", og_model_name)
  
  model_pval <- lrtest(sm_phewas$models[[og_model_name]], sm_phewas_xsquared_smi$models[[new_model_name]]) %>%
    as_tibble %>%
    slice(2) %>%
    pull('Pr(>Chisq)')
  
  old_pval <- sm_phewas$results %>%
    filter(phenotype == phenotype_name, snp=='smi_zscore') %>%
    select(p) %>%
    pull
  new_pval_lin <- sm_phewas_xsquared_smi$results %>%
    filter(phenotype == phenotype_name, snp=='smi_zscore') %>%
    select(p) %>%
    pull
  new_pval_sq <- sm_phewas_xsquared_smi$results %>%
    filter(phenotype == phenotype_name, snp=='smi_zscore_sq') %>%
    select(p) %>%
    pull
  
  
  # Example action: print the row
  if (model_pval < 0.05/29) {
    nonlin_model_preferred_smi <- c(nonlin_model_preferred_smi, phenotype_name)
    
    if (old_pval > 0.05/1222){
      if (new_pval_lin < 0.05/1280) {
        if (new_pval_sq < 0.05/1280) {
          was_not_now_is_significant_lin_sq_smi <- c(was_not_now_is_significant_lin_sq_smi, phenotype_name)
        }
        else {
          was_not_now_is_significant_lin_only_smi <- c(was_not_now_is_significant_lin_only_smi, phenotype_name)
        }
      }
      else {
        if (new_pval_sq < 0.05/1280) {
          was_not_now_is_significant_sq_only_smi <- c(was_not_now_is_significant_sq_only_smi, phenotype_name)
        }
        else {
          was_not_now_is_not_significant_smi <- c(was_not_now_is_not_significant_smi, phenotype_name)
        }
      }
    }
    else {
      if (new_pval_lin < 0.05/1280) {
        if (new_pval_sq < 0.05/1280) {
          was_now_is_significant_lin_sq_smi <- c(was_now_is_significant_lin_sq_smi, phenotype_name)
        }
        else {
          was_now_is_significant_lin_only_smi <- c(was_now_is_significant_lin_only_smi, phenotype_name)
        }
      }
      else {
        if (new_pval_sq < 0.05/1280) {
          was_now_is_significant_sq_only_smi <- c(was_now_is_significant_sq_only_smi, phenotype_name)
        }
        else {
          was_now_is_not_significant_smi <- c(was_now_is_not_significant_smi, phenotype_name)
        }
      }
    }
  }
}

# save final tables
get_two_signif <- function(x){
  two_signif <- formatC(x, digits=2)
  
  if(startsWith(str_sub(two_signif, start=-2), "0") | startsWith(str_sub(two_signif, start=-2), ".")){
    two_signif <- paste0(two_signif, "0")
  } 
  return(two_signif)
}
get_two_signif <- Vectorize(get_two_signif)

print_p_two_signif <- function(p) {
  
  if_else(p>=1.0,
          "1.00",
          if_else(
            p<0.0001,
            "< 0.0001",
            get_two_signif(p)
          )
  )
}

print_or_CI <- function(or, low, high){
  formatted_or <- sprintf("%.2f", or)
  formatted_low <- sprintf("%.2f", low)
  formatted_high <- sprintf("%.2f", high)
  
  paste0(formatted_or, " (", formatted_low, " - ", formatted_high, ")")
}

results_supp_table_3 <- sm_phewas_xsquared_smd$results %>%
  filter(phenotype %in% was_not_now_is_significant_lin_sq_smd | phenotype %in% was_not_now_is_significant_sq_only_smd) %>%
  left_join(PheWAS::pheinfo , by=c("phenotype"="phecode")) %>%
  mutate(bonf_pval = print_p_two_signif(p*1384),
         OR_print = round(OR, 2),
         lower_print = round(lower, 2),
         upper_print = round(upper, 2),
         OR_CI = print_or_CI(OR_print, lower_print, upper_print),
         metric = if_else(snp == "sm_hu_zscore", "SMD", "SMD * SMD")) %>%
  select(phenotype, description, group, metric, bonf_pval, OR_CI, n_total, n_cases, n_controls) %>%
  bind_rows(
    sm_phewas_xsquared_smi$results %>%
      filter(phenotype %in% was_not_now_is_significant_sq_only_smi) %>%
      left_join(PheWAS::pheinfo , by=c("phenotype"="phecode")) %>%
      mutate(bonf_pval = print_p_two_signif(p*1280),
             OR_print = round(OR, 2),
             lower_print = round(lower, 2),
             upper_print = round(upper, 2),
             OR_CI = print_or_CI(OR_print, lower_print, upper_print),
             metric = if_else(snp == "smi_zscore", "SMI", "SMI * SMI")) %>%
      select(phenotype, description, group, metric, bonf_pval, OR_CI, n_total, n_cases, n_controls)
  )

write_csv(results_supp_table_3, paste0(write_dir, 'sup_table_3_results.csv'))

results_supp_table_1 <- studied_associations %>%
  arrange(p) %>%
  mutate(bonf_pval = print_p_two_signif(p*1222),
         OR_print = round(OR, 2),
         lower_print = round(lower, 2),
         upper_print = round(upper, 2),
         OR_CI = print_or_CI(OR_print, lower_print, upper_print),
         metric = if_else(snp == "sm_hu_zscore", "SMD", "SMI")) %>%
  select(phenotype, description, group, metric, bonf_pval, OR_CI, n_total, n_cases, n_controls)

write_csv(results_supp_table_1, paste0(write_dir, 'sup_table_1_results.csv'))
