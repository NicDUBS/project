# Name: ADSL
#
# Label: Subject Level Analysis Dataset
#

# Input: dm, ex, ds
library(admiral)
library(admiral.test) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)

# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

dm <- haven::read_xpt(file = "sdtm/dm.xpt")
ds <- haven::read_xpt(file = "sdtm/ds.xpt")
ex <- haven::read_xpt(file = "sdtm/ex.xpt")
ae <- haven::read_xpt(file = "sdtm/ae.xpt")
lb <- haven::read_xpt(file = "sdtm/lb.xpt")
qs <- haven::read_xpt(file = "sdtm/qs.xpt")
sv <- haven::read_xpt(file = "sdtm/sv.xpt")

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values

dm <- admiral::convert_blanks_to_na(dm)
ds <- admiral::convert_blanks_to_na(ds)
ex <- admiral::convert_blanks_to_na(ex)
ae <- admiral::convert_blanks_to_na(ae)
lb <- admiral::convert_blanks_to_na(lb)
qs <- admiral::convert_blanks_to_na(qs)
sv <- admiral::convert_blanks_to_na(sv)

# User defined functions ----

# Here are some examples of how you can create your own functions that
#  operates on vectors, which can be used in `mutate`.

# Grouping

format_agegr1 <- function(x) {
  case_when(
    x < 65 ~ 1,
    between(x, 18, 80) ~ 2,
    x > 80 ~ 3,
    TRUE ~ NA_integer_
  )
}

format_bmiblgr1 <- function(x) {
  case_when(
    x < 25 ~ '<25',
    x >= 25 & x <30 ~ '25-<30',
    x >= 30 ~ '>=30',
    TRUE ~ NA_character_
  )
}

format_sitegrp1 <- function(x){
  countTable <- table(dm$SITEID)
}

# DCSREAS:	Reason for Discontinuation from Study
# Grouping of DCDECOD values to support summarizing study completion status and reason for discontinuation
format_DCSREAS <- function(x) {
  y<-case_when(
    x == "COMPLETED" ~ "Completed",
    x == "ADVERSE EVENT" ~ "Adverse Event",
    x == "DEATH" ~ "Death",
    x == "SCREEN FAILURE" ~ "I/E Not Met",
    x == "LACK OF EFFICACY" ~ "Lack of Efficacy",
    x == "LOST TO FOLLOW-UP" ~ "Lost to Follow-up",
    x == "PHYSICIAN DECISION" ~ "Physician Decision",
    x == "PROTOCOL VIOLATION" ~ "Protocol Violation",
    x == "STUDY TERMINATED BY SPONSOR" ~ "Sponsor Decision",
    x == "WITHDRAWAL BY SUBJECT" ~ "Withdrew Consent",
    TRUE ~ NA_integer_
  )
  if(is.na(y))
    {
    warning(paste("WARNING: *USER* Please add :", x))
  }
  return(y)
}


# EOSSTT: End of Study Status
# COMPLETED if ADSL.DCDECOD='COMPLETED'. DISCONTINUED if ADSL.DCDECOD not equal to COMPLETED.

# EOSSTT mapping
format_eoxxstt <- function(x) {
  case_when(
    x %in% c("COMPLETED") ~ "COMPLETED",
    !is.na(x) ~ "DISCONTINUED",
  )
}



# Codelist ----
## Race code list ----
race_lookup <- tibble::tribble(
  ~RACE, ~RACEN,
  "AMERICAN INDIAN OR ALASKA NATIVE", 1,
  "ASIAN", 2,
  "BLACK OR AFRICAN AMERICAN", 3,
  "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER", 5,
  "WHITE", 6
)

## Age group code list ----
agegr1_lookup <- tibble::tribble(
  ~AGEGR1, ~AGEGR1N,
  "<65", 1,
  "65-80", 2,
  ">80", 3
)

## ARM code list ----
arm_lookup <- tibble::tribble(
  ~ARM, ~ARMN,
  "Placebo", 1,
  "Xanomeline Low Dose", 2,
  "Xanomeline High Dose", 3
)

# Derivations ----

## copy the data from dM ----
adsl01 <- dm

## add the Pooled Site Group 1 ----

sitegr1_01 <-dm %>%
  dplyr::count(SITEID, ARM) %>%
  dplyr::filter(n<3) %>%
  mutate(SITEGR1=900) %>%
  select(SITEID, SITEGR1) %>%
  distinct()

adsl02_01 <- derive_vars_merged(
  adsl01,
  dataset_add = sitegr1_01,
  by_vars = vars(SITEID)
)
adsl02 <- adsl02_01 %>% mutate(SITEGR1 = ifelse(is.na(adsl02$SITEGR1),adsl02$SITEID, adsl02$SITEGR1))

  adsl02 <- derive_vars_merged(
    adsl01,
    dataset_add = sitegr1_01,
    by_vars = vars(SITEID)
  )

  ## derive treatment variables (TRT01P, TRT01A) ----
  # See also the "Visit and Period Variables" vignette
  # (https://pharmaverse.github.io/admiral/articles/visits_periods.html#treatment_adsl)
  adsl03 <- adsl02 %>% mutate(TRT01P = ARM, TRT01A = ACTARM)

  adsl04 <- derive_vars_merged(
    adsl03,
    dataset_add = dplyr::rename(arm_lookup,TRT01P=ARM, TRT01PN=ARMN),
    by_vars = vars(TRT01P)
  )
  adsl05 <- derive_vars_merged(
    adsl04,
    dataset_add = dplyr::rename(arm_lookup,TRT01A=ARM, TRT01AN=ARMN),
    by_vars = vars(TRT01A)
  )

  # impute start and end time of exposure to first and last respectively, do not impute date
  sv %>% filter(VISITNUM==3) %>%
    select(USUBJID, SVSTDTC) %>%
    admiral::derive_vars_dt(
      dtc = SVSTDTC,
      new_vars_prefix = "TRTS"
    ) %>% select(USUBJID, TRTSDT)


  ## derive treatment start date (TRTSDTM) ----
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
      (EXDOSE == 0 &
        str_detect(EXTRT, "PLACEBO"))) &
      !is.na(EXSTDTM),
    new_vars = vars(TRTSDTM = EXSTDTM, TRTSTMF = EXSTTMF),
    order = vars(EXSTDTM, EXSEQ),
    mode = "first",
    by_vars = vars(STUDYID, USUBJID)
  ) %>%
  ## derive treatment end date (TRTEDTM) ----
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
      (EXDOSE == 0 &
        str_detect(EXTRT, "PLACEBO"))) & !is.na(EXENDTM),
    new_vars = vars(TRTEDTM = EXENDTM, TRTETMF = EXENTMF),
    order = vars(EXENDTM, EXSEQ),
    mode = "last",
    by_vars = vars(STUDYID, USUBJID)
  ) %>%
  ## Derive treatment end/start date TRTSDT/TRTEDT ----
  derive_vars_dtm_to_dt(source_vars = vars(TRTSDTM, TRTEDTM)) %>%
  ## derive treatment duration (TRTDURD) ----
  derive_var_trtdurd()

## Disposition dates, status ----
# convert character date to numeric date without imputation
ds_ext <- derive_vars_dt(
  ds,
  dtc = DSSTDTC,
  new_vars_prefix = "DSST"
)

# Screen fail date
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(SCRFDT = DSSTDT),
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD == "SCREEN FAILURE"
  ) %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(EOSDT = DSSTDT),
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD != "SCREEN FAILURE"
  ) %>%
  # EOS status
  derive_var_disposition_status(
    dataset_ds = ds_ext,
    new_var = EOSSTT,
    status_var = DSDECOD,
    format_new_var = format_eoxxstt,
    filter_ds = DSCAT == "DISPOSITION EVENT"
  ) %>%
  # Last retrieval date
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(FRVDT = DSSTDT),
    filter_add = DSCAT == "OTHER EVENT" & DSDECOD == "FINAL RETRIEVAL VISIT"
  ) %>%
  # Derive Randomization Date
  derive_vars_merged(
    dataset_add = ds_ext,
    filter_add = DSDECOD == "RANDOMIZED",
    by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(RANDDT = DSSTDT)
  ) %>%
  # Death date - impute partial date to first day/month
  derive_vars_dt(
    new_vars_prefix = "DTH",
    dtc = DTHDTC,
    highest_imputation = "M",
    date_imputation = "first"
  ) %>%
  # Relative Day of Death
  derive_vars_duration(
    new_var = DTHADY,
    start_date = TRTSDT,
    end_date = DTHDT
  ) %>%
  # Elapsed Days from Last Dose to Death
  derive_vars_duration(
    new_var = LDDTHELD,
    start_date = TRTEDT,
    end_date = DTHDT,
    add_one = FALSE
  )

## Last known alive date ----
ae_start_date <- date_source(
  dataset_name = "ae",
  date = AESTDT
)
ae_end_date <- date_source(
  dataset_name = "ae",
  date = AEENDT
)
lb_date <- date_source(
  dataset_name = "lb",
  date = LBDT,
  filter = !is.na(LBDT)
)
trt_end_date <- date_source(
  dataset_name = "adsl",
  date = TRTEDT
)

# impute AE start and end date to first
ae_ext <- ae %>%
  derive_vars_dt(
    dtc = AESTDTC,
    new_vars_prefix = "AEST",
    highest_imputation = "M"
  ) %>%
  derive_vars_dt(
    dtc = AEENDTC,
    new_vars_prefix = "AEEN",
    highest_imputation = "M"
  )

# impute LB date to first
lb_ext <- derive_vars_dt(
  lb,
  dtc = LBDTC,
  new_vars_prefix = "LB",
  highest_imputation = "M"
)

adsl <- adsl %>%
  derive_var_extreme_dt(
    new_var = LSTALVDT,
    ae_start_date, ae_end_date, lb_date, trt_end_date,
    source_datasets = list(ae = ae_ext, lb = lb_ext, adsl = adsl),
    mode = "last"
  ) %>%
  derive_var_merged_exist_flag(
    dataset_add = ex,
    by_vars = vars(STUDYID, USUBJID),
    new_var = SAFFL,
    condition = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO")))
  ) %>%
  ## Groupings and others variables (RACEGR1, AGEGR1, REGION1,LDDTHGR1,DTH30FL, DTHA30FL, DTHB30FL  ) ----
  mutate(
    RACEGR1 = format_racegr1(RACE),
    AGEGR1 = format_agegr1(AGE),
    REGION1 = format_region1(COUNTRY),
    LDDTHGR1 = format_lddthgr1(LDDTHELD),
    DTH30FL = if_else(LDDTHGR1 == "<= 30", "Y", NA_character_),
    DTHA30FL = if_else(LDDTHGR1 == "> 30", "Y", NA_character_),
    DTHB30FL = if_else(DTHDT <= TRTSDT + 30, "Y", NA_character_),
    DOMAIN = NULL
  )

## derived MMSETOT:	MMSE Total ----
# sum of QS.QSORRES values for the subject when QSCAT = MINI-MENTAL STATE
mental <- qs %>%
  filter(qs$QSCAT=="MINI-MENTAL STATE")

summental <- mental %>%
  group_by(USUBJID) %>%
  summarise(mmsetot=sum(QSSTRESN))



# DCDECOD:	Standardized Disposition Term
# DS.DSDECOD where DSCAT='DISPOSITION EVENT'
DCDECOD <- ds %>%
  filter(DSCAT=='DISPOSITION EVENT') %>%
  rename(DCDECOD=DSDECOD) %>%
  select (USUBJID, DCDECOD)

# RFENDT: Date of Discontinuation/Completion
# RFENDTC converted to SAS date
# a terliner RFENDTC <- ymd(RFENDTC)

# Save output ----

dir <- tempdir() # Change to whichever directory you want to save the dataset in
saveRDS(adsl, file = file.path(dir, "adsl.rds"), compress = "bzip2")

