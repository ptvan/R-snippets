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

# generate some ADSL from the example SDTM
dm <- pharmaversesdtm::dm %>%
      convert_blanks_to_na() %>%
      select(-DOMAIN)

###########################################################
# Observational Medical Outcomes Partnership (OMOP) data
###########################################################
library(DBI)
library(duckdb)
library(duckplyr)
library(CDMConnector)
# connect to a Synthea-generated CDM (via OHDSI Eunomia project)
# These datasets are stored in $EUNOMIA_DATA_FOLDER, specified in $HOME/.Renviron
# as DuckDB databases
db <- dbConnect(duckdb::duckdb(), 
                dbdir = eunomiaDir(datasetName = "GiBleed"))
cdm <- cdmFromCon(db, cdmSchema = "main", writeSchema = "main")
