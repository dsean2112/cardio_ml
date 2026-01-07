# Helper functions for Cohort analysis

#' @description
#' Input is a vector of ICD codes, output is a vector of record IDs
#' 
icd2recid <- function(icds) {
  # Input: icd codes, .txt file. Output: dataframe of record id
  library(arrow)
  library(dplyr)
  # Input: list of icds. 
  icds <- paste0(icds, collapse = "|") 
  # Output: dataframe of mrns
  files <- list.files(
    '/mmfs1/projects/cardio_darbar_chi/common/data/Diagnosis/',
    recursive = TRUE,
    full.names = TRUE
  )
  
  recids <- c()
  for (i in c(1:length(files))) {
    batch <- arrow::open_dataset(files[i]) |>
      dplyr::filter(grepl(x = icd_code, pattern = icds)) |>
      collect()
    
    recids <- c(recids, batch$record_id)
    
    print(paste0('Finished ', i, '. ', nrow(batch)))
    rm(batch)
  }
  recids <- unique(recids)
  return(data.frame(record_id = recids))
}


# icdOrMrn2demo_table -----------------------------------------------------
#' Input: dataframe with either MRN (name 'mrn') or Record ID (name record_id)
icdOrMrn2demo_table <- function(table) {
  
  # Load and edit master demographic table:
  library(dplyr)
  demo <- read.csv('/mmfs1/projects/cardio_darbar_chi/common/data/demographics-0.csv')
  demo <- demo %>% select(-first_name, -last_name, -ssn, -age)
  demo$mrn <- as.numeric(demo$mrn)
  
  demo$gender[demo$gender == 0] <- "NA"
  demo$gender[demo$gender == 1] <- "m"
  demo$gender[demo$gender == 2] <- "f"
  
  demo$race[demo$race == 1] <- "AI_AN" # american indian / alaska native
  demo$race[demo$race == 2] <- "Asian"
  demo$race[demo$race == 3] <- "NH_PI" # native hawaiian / pacific islander
  demo$race[demo$race == 4] <- "Black"
  demo$race[demo$race == 5] <- "White"
  demo$race[demo$race == 6] <- "Other"
  demo$race[demo$race == 7] <- "Decline"
  demo$race[demo$race == 8] <- "NA"
  
  demo$ethnicity[demo$ethnicity == 1] <- "Not_HL" # hispanic/latino
  demo$ethnicity[demo$ethnicity == 2] <- "HL"
  demo$ethnicity[demo$ethnicity == 3] <- "decline"
  demo$ethnicity[demo$ethnicity == 4] <- "NA"
  
  demo$insurance_type[demo$insurance_type == 1] <- "private"
  demo$insurance_type[demo$insurance_type == 2] <- "medicaid"
  demo$insurance_type[demo$insurance_type == 3] <- "medicare"
  demo$insurance_type[demo$insurance_type == 4] <- "selfpay"
  demo$insurance_type[demo$insurance_type == 5] <- "other"
  demo$insurance_type[demo$insurance_type == 6] <- "NA"
  
  # Filter master demographic table based on mrn or record_id:
  # --- Case 1: table has MRN ---
  if ("mrn" %in% names(table)) {
    # Ensure mrn is integer
    table$mrn <- as.integer(table$mrn)
    
    # If demo has record_id, remove it
    if ("record_id" %in% names(demo)) {
      demo <- demo %>% select(-mrn)
    }
    
    # Ensure demo$mrn exists and is integer
    if ("mrn" %in% names(demo)) {
      demo$mrn <- as.integer(demo$mrn)
    }
    
    # Merge on MRN
    out <- merge(table, demo, by = "mrn", all.x = TRUE)
    return(out)
  }
  
  # --- Case 2: table has RECORD_ID ---
  if ("record_id" %in% names(table)) {
    # Ensure record_id is integer
    table$record_id <- as.integer(table$record_id)
    
    # Ensure demo$record_id exists and is integer
    if ("record_id" %in% names(demo)) {
      demo$record_id <- as.integer(demo$record_id)
    }
    
    # Merge on record_id
    out <- merge(table, demo, by = "record_id", all.x = TRUE)
    return(out)
  }
  
  # --- If neither ID column exists ---
  stop("Neither 'mrn' nor 'record_id' found in table.")
}


# icdRecid2date -----------------------------------------------------------
icdRecid2date <- function(table,icds,name) {
  # Input: table of recids; ICD code/s; name of new column
  # Output: table with column for date of initial dx
  # Note: can handle multiple ICDs for single dx
  
  icds <- paste0(icds, collapse = "|") 
  
  files <- list.files(
    '/mmfs1/projects/cardio_darbar_chi/common/data/Diagnosis/',
    recursive = TRUE,
    full.names = TRUE
  )
  
  all_years <- c()
  for (i in c(1:length(files))) {
    batch <- arrow::open_dataset(files[i]) |>
      dplyr::filter(grepl(x = icd_code, pattern = icds) & 
                      record_id %in% table$recid) |>
      select(record_id, date) |>
      collect()
    
    # recids <- c(recids, batch$record_id)
    
    print(paste0('Finished ', i, '. ', nrow(batch)))
    all_years <- rbind(all_years,batch)
    rm(batch)
  }
  
  min_dates <- all_years %>%
    group_by(record_id) %>%
    summarize(date = min(date, na.rm = TRUE))  # Get min date per record_id
  
  min_dates <- min_dates %>% rename(recid = record_id) # rename
  min_dates <- min_dates %>% rename(!!name := date)
  
  final_table <- table %>%
    left_join(min_dates, by = "recid")  # Perform a left join on 'rec_id'
  
  return(final_table)
}