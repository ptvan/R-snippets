############################
# Pharmaverse clinical data
############################
library(admiral)
library(pharmaversesdtm)
library(pharmaverseadam)
library(metatools)
library(metacore)
library(xportr)
library(lubridate)
library(stringr)
library(dplyr)

# extract the provided example ADSL
adsl <- admiral::admiral_adsl

# alternately, generate some ADSL from the example SDTM
dm <- pharmaversesdtm::dm %>%
      convert_blanks_to_na() %>%
      select(-DOMAIN)

# merge the ADSL to EX using `STUDYID` and `USUBJID`
# for a subset of variables

adsl_vars <- exprs(TRTSDT, TRTSDTM, TRTEDT, TRTEDTM)
adex <- derive_vars_merged(
  ex,
  dataset_add = adsl,
  new_vars = adsl_vars,
  by_vars = get_admiral_option("subject_keys")
)

# derive datetime and analysis day:
adex <- derive_vars_dt(adex, new_vars_prefix = "AST", dtc = EXSTDTC)
adex <- derive_vars_dt(adex, new_vars_prefix = "AEN", dtc = EXENDTC)

adex <- derive_vars_dtm(
  adex,
  dtc = EXSTDTC,
  highest_imputation = "M",
  new_vars_prefix = "AST"
)

adex <- derive_vars_dtm(
  adex,
  dtc = EXENDTC,
  highest_imputation = "M",
  date_imputation = "last",
  new_vars_prefix = "AEN"
)

###########################################################
# Observational Medical Outcomes Partnership (OMOP) data
###########################################################
library(DBI)
library(duckdb)
library(duckplyr)
library(CDMConnector)
library(dplyr)

# connect to a Synthea-generated CDM (from OHDSI's Eunomia project)
# These datasets are stored in $EUNOMIA_DATA_FOLDER, as DuckDB databases,
# specified in $HOME/.Renviron
db <- dbConnect(duckdb::duckdb(), 
                dbdir = eunomiaDir(datasetName = "GiBleed"))
cdm <- cdmFromCon(db, cdmSchema = "main", writeSchema = "main")

# reminder that this is synthetic data
cdm$cdm_source

# list raw DuckDB tables
listTables(db)

# glimpse patient-level demographics table
cdm$person %>% glimpse()

# `$` is for internal use only, use pull() to get values from tables
cdm$person %>% pull(gender_concept_id) %>% unique()

# summarize patients' gender by birth year by merging `person` with `concept`
cdm$person %>%
  group_by(year_of_birth, gender_concept_id) %>%
  summarize(n = n(), .groups= "drop") %>%
  inner_join(cdm$concept, by=c("gender_concept_id" = "concept_id")) %>%
  select(year_of_birth, concept_name, n) 
  
# the dplyr call above is just wrapping SQL
gender_by_year_SQL <- cdm$person %>%
  group_by(year_of_birth, gender_concept_id) %>%
  summarize(n = n(), .groups= "drop") %>%
  inner_join(cdm$concept, by=c("gender_concept_id" = "concept_id")) %>%
  select(year_of_birth, concept_name, n) 
gender_by_year_SQL %>% show_query()