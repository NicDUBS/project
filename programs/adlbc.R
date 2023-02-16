# Name: ADLBC
#
# Label: Lab Analysis Dataset
#
# Input: adsl, lb
library(admiral)
#library(admiral.test) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)
library(xportr)
library(readxl)

# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data
lb <- haven::read_xpt(file = "sdtm/lb.xpt")

adsl <- admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values

lb <- convert_blanks_to_na(lb)

# Look-up tables ----

chemo <- lb %>%
  filter(LBTESTCD %in% c("ALB" , "ALP", "ALT",  "AST",  "BILI", "BUN",  "CA",   "CHOL", "CK",   "CL",   "CREAT","GGT",  "GLUC", "K",    "PHOS", "PROT", "SODIUM", "URATE"))

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~LBTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
  "ALB", "ALB", "Albumin (g/L)",	33,
  "ALP", "ALP", "Alkaline Phosphatase (U/L)",	22,
  "ALT", "ALT", "Alanine Aminotransferase (U/L)",	24,
  "AST", "AST", "Aspartate Aminotransferase (U/L)",	25,
  "BILI", "BILI", "Bilirubin (umol/L)",	21,
  "BUN", "BUN", "Blood Urea Nitrogen (mmol/L)",	26,
  "CA", "CA", "Calcium (mmol/L)",	30,
  "CHOL", "CHOL", "Cholesterol (mmol/L)",	34,
  "CK", "CK", "Creatine Kinase (U/L)",	35,
  "CL", "CL", "Chloride (mmol/L)",	20,
  "CREAT", "CREAT", "Creatinine (umol/L)",	27,
  "GGT", "GGT	Gamma", "Glutamyl Transferase (U/L)",	23,
  "GLUC", "GLUC", "Glucose (mmol/L)",	31,
  "K", "K", "Potassium (mmol/L)",	19,
  "PHOS", "PHOS", "Phosphate (mmol/L)",	29,
  "PROT", "PROT", "Protein (g/L)",	32,
  "SODIUM", "SODIUM", "Sodium (mmol/L)",	18,
  "URATE", "URATE", "Urate (umol/L)",	28
)

paramchg_lookup <- tibble::tribble(
  ~LBTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
  "ALB",    "_ALB",         "Albumin (g/L) change from previous visit, relative to normal range",	133,
  "ALP",    "_ALP",         "Alkaline Phosphatase (U/L) change from previous visit, relative to normal range",	122,
  "ALT",    "_ALT",         "Alanine Aminotransferase (U/L) change from previous visit, relative to normal range",	124,
  "AST",    "_AST",         "Aspartate Aminotransferase (U/L) change from previous visit, relative to normal range",	125,
  "BILI",   "_BILI",        "Bilirubin (umol/L) change from previous visit, relative to normal range",	121,
  "BUN",    "_BUN",         "Blood Urea Nitrogen (mmol/L) change from previous visit, relative to normal range",	126,
  "CA",     "_CA",          "Calcium (mmol/L) change from previous visit, relative to normal range",	130,
  "CHOL",   "_CHOL",        "Cholesterol (mmol/L) change from previous visit, relative to normal range",	134,
  "CK",     "_CK",          "Creatine Kinase (U/L) change from previous visit, relative to normal range",	135,
  "CL",     "_CL",          "Chloride (mmol/L) change from previous visit, relative to normal range",	120,
  "CREAT",  "_CREAT",       "Creatinine (umol/L) change from previous visit, relative to normal range",	127,
  "GGT",    "_GGT	Gamma",   "Glutamyl Transferase (U/L) change from previous visit, relative to normal range",	123,
  "GLUC",   "_GLUC",        "Glucose (mmol/L) change from previous visit, relative to normal range",	131,
  "K",      "_K",           "Potassium (mmol/L) change from previous visit, relative to normal range",	119,
  "PHOS",   "_PHOS",        "Phosphate (mmol/L) change from previous visit, relative to normal range",	129,
  "PROT",   "_PROT",        "Protein (g/L) change from previous visit, relative to normal range",	132,
  "SODIUM", "_SODIUM",      "Sodium (mmol/L) change from previous visit, relative to normal range",	118,
  "URATE",  "_URATE",       "Urate (umol/L) change from previous visit, relative to normal range",	128
)

# Derivations ----

# Get list of ADSL vars required for derivations
# need to add when ADSL done:
# adsl_vars <- vars(SUBJID, TRT01P,TRT01PN, TRT01A, TRT01AN, TRTSDT, TRTEDT, AGE, AGEGR1, AGEGR1N, RACE, RACEN, SEX, COMP24FL, DSRAEFL, SAFFL)
adsl_vars <- vars(SUBJID, TRT01P, TRT01A, TRTSDT, TRTEDT)

adlb1 <- chemo %>%
  # Join ADSL with LB (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = vars(STUDYID, USUBJID)
  ) %>%
  ## Calculate ADT, ADY ----
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = LBDTC
  ) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = vars(ADT))

adlb2 <- adlb1 %>%
  ## Add PARAMCD PARAM and PARAMN - from LOOK-UP table ----
  # Replace with PARAMCD lookup function
  derive_vars_merged_lookup(
    dataset_add = param_lookup,
    new_vars = vars(PARAMCD, PARAM, PARAMN),
    by_vars = vars(LBTESTCD),
    check_type = "none",
    print_not_mapped = FALSE
  ) %>%
  ## Calculate PARCAT1 AVAL AVALC ANRLO ANRHI ----
  mutate(
    PARCAT1 = 'CHEM',
    AVAL = LBSTRESN,
    AVALC = LBSTRESC,
    ANRLO = LBSTNRLO,
    ANRHI = LBSTNRHI,
    A1LO = LBSTNRLO,
    A1HI = LBSTNRHI
  )


## Get Visit Info ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#visits)
adlb3 <- adlb2 %>%
  # Derive Timing
  mutate(
    AVISIT = case_when(
      str_detect(VISIT, "SCREEN") ~ "Baseline",
      !is.na(VISIT) ~ str_to_title(VISIT),
      TRUE ~ NA_character_
    ),
    AVISITN = case_when(
      AVISIT == "Baseline" ~ 0,
      !is.na(VISITNUM) ~ VISITNUM
    )
  )



## Calculate ANRIND : requires the reference ranges ANRLO, ANRHI ----
adlb4 <- adlb3 %>%
  derive_var_anrind()

## Derive baseline flags ----
adlb5 <- adlb4 %>%
  # Calculate BASETYPE
  mutate(
    BASETYPE = "LAST"
  ) %>%
  # Calculate ABLFL
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(STUDYID, USUBJID, BASETYPE, PARAMCD),
      order = vars(ADT, VISITNUM, LBSEQ),
      new_var = ABLFL,
      mode = "last"
    ),
    filter = (!is.na(AVAL) & ADT <= TRTSDT & !is.na(BASETYPE))
  )

## Derive baseline information ----
adlb6 <- adlb5 %>%
  # Calculate BASE
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVAL,
    new_var = BASE
  ) %>%
  # Calculate BASEC
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVALC,
    new_var = BASEC
  ) %>%
  # Calculate BNRIND
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = ANRIND,
    new_var = BNRIND
  ) %>%
  # Calculate CHG
  derive_var_chg()
  #%>%
  # Calculate PCHG
  #derive_var_pchg()


## Calculate lab grading ----

# Assign ATOXDSCL and ATOXDSCH to hold lab grading terms
# ATOXDSCL and ATOXDSCH hold terms defined by NCI-CTCAEv4.
# See (https://pharmaverse.github.io/admiral/articles/lab_grading.html#implement_ctcv4)
#grade_lookup <- tibble::tribble(
#  ~PARAMCD, ~ATOXDSCL, ~ATOXDSCH,
#  "ALB", "Hypoalbuminemia", NA_character_,
#  "ALKPH", NA_character_, "Alkaline phosphatase increased",
#  "ALT", NA_character_, "Alanine aminotransferase increased",
#  "AST", NA_character_, "Aspartate aminotransferase increased",
#  "BILI", NA_character_, "Blood bilirubin increased",
#  "CA", "Hypocalcemia", "Hypercalcemia",
#  "CHOLES", NA_character_, "Cholesterol high",
#  "CK", NA_character_, "CPK increased",
#  "CREAT", NA_character_, "Creatinine increased",
#  "GGT", NA_character_, "GGT increased",
#  "GLUC", "Hypoglycemia", "Hyperglycemia",
#  "HGB", "Anemia", "Hemoglobin increased",
#  "POTAS", "Hypokalemia", "Hyperkalemia",
#  "LYMPH", "CD4 lymphocytes decreased", NA_character_,
#  "PHOS", "Hypophosphatemia", NA_character_,
#  "PLAT", "Platelet count decreased", NA_character_,
#  "SODIUM", "Hyponatremia", "Hypernatremia",
#  "WBC", "White blood cell decreased", "Leukocytosis",
#)
#
# Assign grade criteria
# metadata atoxgr_criteria_ctcv4 used to implement NCI-CTCAEv4
# user could change to atoxgr_criteria_ctcv5 to implement NCI-CTCAEv5
# Note: Hyperglycemia and Hypophosphatemia not defined in NCI-CTCAEv5 so
# user would need to amend look-up table grade_lookup
# See (https://pharmaverse.github.io/admiral/articles/lab_grading.html#implement_ctcv5)
#grade_crit <- atoxgr_criteria_ctcv4


# Add ATOXDSCL and ATOXDSCH
#adlb <- adlb %>%
#  derive_vars_merged(
#    dataset_add = grade_lookup,
#    by_vars = vars(PARAMCD)
#  ) %>%
#  # Derive toxicity grade for low values ATOXGRL
#
#  derive_var_atoxgr_dir(
#    meta_criteria = grade_crit,
#    new_var = ATOXGRL,
#    tox_description_var = ATOXDSCL,
#    criteria_direction = "L",
#    get_unit_expr = extract_unit(PARAM)
#  ) %>%
#  # Derive toxicity grade for low values ATOXGRH
#  # default metadata atoxgr_criteria_ctcv4 used
#  derive_var_atoxgr_dir(
#    meta_criteria = grade_crit,
#    new_var = ATOXGRH,
#    tox_description_var = ATOXDSCH,
#    criteria_direction = "H",
#    get_unit_expr = extract_unit(PARAM)
#  ) %>%
#  # (Optional) derive overall grade ATOXGR (combining ATOXGRL and ATOXGRH)
#  derive_var_atoxgr() %>%
#  # Derive baseline toxicity grade for low values BTOXGRL
#  derive_var_base(
#    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
#    source_var = ATOXGRL,
#    new_var = BTOXGRL
#  ) %>%
#  # Derive baseline toxicity grade for high values BTOXGRH
#  derive_var_base(
#    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
#    source_var = ATOXGRH,
#    new_var = BTOXGRH
#  ) %>%
#  # Derive baseline toxicity grade for for overall grade BTOXGR
#  derive_var_base(
#    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
#    source_var = ATOXGR,
#    new_var = BTOXGR
#  )
#

## Calculate R2BASE, R2ANRLO and R2ANRHI ----
adlb7 <- adlb6 %>%
  derive_var_analysis_ratio(
    numer_var = AVAL,
    denom_var = A1LO
  ) %>%
  derive_var_analysis_ratio(
    numer_var = AVAL,
    denom_var = A1HI
  )

## SHIFT derivation ----
#adlb <- adlb %>%
#  # Derive shift from baseline for analysis indicator
#  derive_var_shift(
#    new_var = SHIFT1,
#    from_var = BNRIND,
#    to_var = ANRIND
#  ) %>%
#  # Derive shift from baseline for overall grade
#  restrict_derivation(
#    derivation = derive_var_shift,
#    args = params(
#      new_var = SHIFT2,
#      from_var = BTOXGR,
#      to_var = ATOXGR
#    ),
#    filter = !is.na(ATOXDSCL) | !is.na(ATOXDSCH)
#  )
#
## Flag variables (ANL01FL, LVOTFL) ----
# ANL01FL: Flag last result within an AVISIT for post-baseline records
# LVOTFL: Flag last valid on-treatment record
#adlb <- adlb %>%
#  restrict_derivation(
#    derivation = derive_var_extreme_flag,
#    args = params(
#      by_vars = vars(USUBJID, PARAMCD, AVISIT),
#      order = vars(ADT, AVAL),
#      new_var = ANL01FL,
#      mode = "last"
#    ),
#    filter = !is.na(AVISITN) & ONTRTFL == "Y"
#  ) %>%
#  restrict_derivation(
#    derivation = derive_var_extreme_flag,
#    args = params(
#      by_vars = vars(USUBJID, PARAMCD),
#      order = vars(ADT, AVAL),
#      new_var = LVOTFL,
#      mode = "last"
#    ),
#    filter = ONTRTFL == "Y"
#  )

## Get treatment information ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#treatment_bds)
adlbc <- adlb7 %>%
  # Assign TRTA, TRTP
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  )

## Get extreme values ----
#adlb <- adlb %>%
#  # get MINIMUM value
#  derive_extreme_records(
#    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
#    order = vars(AVAL, ADT, AVISITN),
#    mode = "first",
#    # "AVISITN < 9997" to evaluate only real visits
#    filter = (!is.na(AVAL) & ONTRTFL == "Y" & AVISITN < 9997),
#    set_values_to = vars(
#      AVISITN = 9997,
#      AVISIT = "POST-BASELINE MINIMUM",
#      DTYPE = "MINIMUM"
#    )
#  ) %>%
#  # get MAXIMUM value
#  derive_extreme_records(
#    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
#    order = vars(desc(AVAL), ADT, AVISITN),
#    mode = "first",
#    # "AVISITN < 9997" to evaluate only real visits
#    filter = (!is.na(AVAL) & ONTRTFL == "Y" & AVISITN < 9997),
#    set_values_to = vars(
#      AVISITN = 9998,
#      AVISIT = "POST-BASELINE MAXIMUM",
#      DTYPE = "MAXIMUM"
#    )
#  ) %>%
#  # get LOV value
#  derive_extreme_records(
#    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
#    order = vars(ADT, AVISITN),
#    mode = "last",
#    # "AVISITN < 9997" to evaluate only real visits
#    filter = (ONTRTFL == "Y" & AVISITN < 9997),
#    set_values_to = vars(
#      AVISITN = 9999,
#      AVISIT = "POST-BASELINE LAST",
#      DTYPE = "LOV"
#    )
#  )
#
### Get ASEQ ----
#adlb <- adlb %>%
#  # Calculate ASEQ
#  derive_var_obs_number(
#    new_var = ASEQ,
#    by_vars = vars(STUDYID, USUBJID),
#    order = vars(PARAMCD, ADT, AVISITN, VISITNUM),
#    check_type = "error"
#  )

# Add all ADSL variables
#adlb <- adlb %>%
#  derive_vars_merged(
#    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
#    by_vars = vars(STUDYID, USUBJID)
#  )

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...

# Save output ----

#dir <- tempdir() # Change to whichever directory you want to save the dataset in
#saveRDS(adlb, file = file.path(dir, "adlb.rds"), compress = "bzip2")
#saveRDS(adlbc, file = "./adam/adlbc.rds", compress = "bzip2")

var_spec <- readxl::read_xlsx(path="./metadata/specs.xlsx", sheet = "Variables") %>%
  dplyr::rename(type = "Data Type") %>%  rlang::set_names(tolower)


adlbcall <- adlbc %>%
  xportr::xportr_type(var_spec, "ADLBC", "message") %>%
  xportr::xportr_length(var_spec, "ADLBC", "message") %>%
  xportr::xportr_label(var_spec, "ADLBC", "message") %>%
  xportr::xportr_order(var_spec, "ADLBC", "message") %>%
  xportr::xportr_format(var_spec, "ADLBC", "message")

adlbcall %>%
  xportr::xportr_write(path="./adam/adlbc.xpt", label = "Subject-Level Analysis Dataset")

