mrn_copy_wfdb <- function(mrnFile, folderName) {
  # Copy ECGs of a given mrn list to a specified folder
  
  # MRN file: csv, list of mrns
  # folderName: output folder
  
  # Set paths
  home <- fs::path_expand('/shared/projects/cardio_darbar/common/')
  main <- fs::path('data','cohorts', 'f-wave')
  wfdb <- fs::path(home, 'data','uic', 'wfdb')
  
  library(vroom)
  library(fs)
  library(dplyr)
  
  # Load mrn file
  mrnData = fs::path(home, main, mrnFile) |>
    vroom::vroom(delim=',')
  
  mrnData = mrnData[,2]
  
  # Create output folder
  outputFolder <- fs::path(home, main, folderName)
  if (!fs::dir_exists(outputFolder)) {
    fs::dir_create(outputFolder)
  }
  
  # Load log file with wfdb file locations
  logFile <- fs::path(home,'data','cohorts','dm_afib','logs','wfdb_new', ext = 'log')
  logData <- vroom::vroom(file = logFile, delim = ',')
    # system2(command = 'wc',
    #         args = paste('-l', logFile),
    #         stdout = TRUE) |>
    # readr::parse_number() |>
    # {\(.x) .x - 1 }() # Subtract off header
  
  colnames(logData) <- c('mrn','MUSE_ID','WFDB_PATH')
  colnames(mrnData) <- c('mrn')
    
  batchData <- inner_join(mrnData, logData, by = "mrn")
  
  # Copy files to output folder:
  for (i in 1:nrow(batchData)) {
    fn <- as.character(batchData[i, 2])
    
    fp <-
      fs::path(home, batchData[i, 3]) |>
      fs::path_dir()
    
    files <- list.files(
      path = fp,
      pattern = as.character(batchData[i, 2]),
      full.names = TRUE
    )
    fs::file_copy(path = files,
                  new_path = outputFolder,
                  overwrite = TRUE)
    
    if (length(files) >= 1) {
      cat("\tCopying:", fn, "\n")
    }
    
  }
}

mrn_wfdb_table <- function(mrnFile) {
  # Create a table of ECG file names, locations and associated muse file names
  # from list of mrns
  
  # MRN file: csv, list of mrns
  
  home <- fs::path_expand('/shared/projects/cardio_darbar/common/')
  main <- fs::path('data','cohorts', 'f-wave')
  wfdb <- fs::path(home, 'data','uic', 'wfdb')
  
  library(vroom)
  library(fs)
  library(dplyr)
  
  # Load mrn file
  mrnData = fs::path(home, main, mrnFile) |>
    vroom::vroom(delim=',')
  mrnData = mrnData[,2]
  
  # Load log file containing all WFDB files:
  wfdb_logFile <- fs::path(home,'data','cohorts','dm_afib','logs','wfdb_new', ext = 'log')
  wfdb_logData <- vroom::vroom(file = wfdb_logFile, delim = ',')
  
  
  # Coordinate column names and join
  colnames(wfdb_logData) <- c('mrn','MUSE_ID','WFDB_PATH')
  colnames(mrnData) <- c('mrn')
  
  batchData <- inner_join(mrnData, wfdb_logData, by = "mrn")
  
  # Load MUSE file log
  muse_logFile <- fs::path(home,'data','cohorts','dm_afib','logs','muse_new', ext = 'log')
  muse_logData <- vroom::vroom(file = muse_logFile)
  
  colnames(muse_logData) <- c('MUSE_ID','MUSE_PATH')
  batchData <- left_join(batchData, muse_logData, by='MUSE_ID')
  
  # Remove duplicate muse files:
  batchData <- batchData[!duplicated(batchData$MUSE_ID),]
  
  return(batchData)
}

muse_dx_line <- function(muse_ids, muse_path = '/mmfs1/projects/cardio_darbar_chi/common/data/muse/muse.log') {
  # Input: muse files list
  # output: diagnosis line
  
  library(XML)
  library(xml2)
  if (class(muse_path) == 'character') {
    muse_log <- read.csv(muse_path)
  } else if (class(muse_path) == 'data.frame') { # for greater efficiency, allow user to set muse_path to be the actual muse log df
    muse_log <- muse_path
  }
  
  
  muse_log <- muse_log %>% dplyr::filter(FILE_NAME %in% muse_ids)
  # FINISH CODING: 
  #   Load muse_path csv: read.csv
  #   For each muse id
  #     Locate muse file using muse_log
  #     Complete loop completed below:
  
  base_path <- '/mmfs1/projects/cardio_darbar_chi/common'
  xml <- read_xml(fs::path(base_path, muse_log$PATH, ext = 'xml'))
  diagnosis <- xml_text(xml_find_all(xml, ".//DiagnosisStatement"))
  
  diagnosis <- xml_text(xml_find_all(xml, ".//DiagnosisStatement"))
  diagnosis <- gsub("ENDSLINE|USERINSERT", "", diagnosis) # Remove extraneous text
  diagnosis <- paste(diagnosis, collapse = " ") # Collapse into a single string
  
  return(diagnosis)
}
  
icd2mrn <- function(icds_path) {
  # OUTDATED, use icd2recid
  # library(arrow)
  # library(dplyr)
  # library(vroom)
  # 
  # diagnosis_path <- '/mmfs1/projects/cardio_darbar_chi/common/data/Diagnosis/'
  # demographics_path <- '/mmfs1/projects/cardio_darbar_chi/common/data/demographics-0.csv'
  # 
  # # Input: path to .txt file with list of mrns
  # icds <- readLines(icds_path)
  # icds <- paste0("^", icds, collapse = "|")
  # 
  # 
  # # Output: dataframe of mrns
  # files <- list.files(diagnosis_path, 
  #                     recursive = TRUE, 
  #                     full.names = TRUE)
  # 
  # recids <- c()
  # for (i in c(1:length(files))) {
  #   batch <- arrow::open_dataset(files[i]) |>
  #     dplyr::filter(grepl(x = ICD_CODE, pattern = icds)) |>
  #     collect()
  #   
  #   recids <- c(recids,batch$RECORD_ID)
  #   recids <- unique(recids)
  #   
  #   print(paste('Finished',i,'.',nrow(batch)))
  #   rm(batch)
  # }
  # recids <- unique(recids)
  # 
  # demographics <- vroom(demographics_path, col_select = c("RECORD_ID", "MRN"))
  # output <- demographics |> filter(RECORD_ID %in% recids)
  # return(output)
}


#' @description
#' Input: vector of ICD codes (ICD10 only, ICD9s are converted). Optional: start and end dates.
#' Output: dataframe of record ids
#' 
icd2recid <- function(icds, start_date = NULL, end_date = NULL) {
  library(arrow)
  library(dplyr)
  library(lubridate)
  
  # Collapse ICDs into regex
  icd_pattern <- paste0(icds, collapse = "|")
  
  # List all diagnosis files
  files <- list.files(
    '/mmfs1/projects/cardio_darbar_chi/common/data/Diagnosis/',
    recursive = TRUE,
    full.names = TRUE
  )
  
  # Extract year from filenames like "Diagnosis_2010.parquet"
  file_years <- files |>
    basename() |>
    sub("Diagnosis_", "", x = _) |>
    sub(".parquet", "", x = _) |>
    as.integer()
  
  # If dates provided, filter files by year
  if (!is.null(start_date)) {
    start_date <- as.Date(start_date)
    start_year <- year(start_date)
  } else {
    start_year <- min(file_years)
  }
  
  if (!is.null(end_date)) {
    end_date   <- as.Date(end_date)
    end_year   <- year(end_date)
  } else {
    end_year <- max(file_years)
  }
  
  keep_idx <- file_years >= start_year & file_years <= end_year
  files <- files[keep_idx]
  file_years <- file_years[keep_idx]
  
  recids <- c()
  
  for (i in seq_along(files)) {
    yr <- file_years[i]
    ds <- arrow::open_dataset(files[i])
    
    # Base filter: ICD codes
    ds <- ds |> filter(grepl(icd_code, pattern = icd_pattern))
    
    # If date filtering needed for boundary years
    if (!is.null(start_date)) {
      if (yr == start_year) {
        ds <- ds |> filter(date >= start_date)
      }
    }
    if (!is.null(end_date)) {
      if (yr == end_year) {
        ds <- ds |> filter(date <= end_date)
      }
    }
    
    batch <- ds |> collect()
    
    recids <- c(recids, batch$record_id)
    
    message("Finished ", i, ". Rows: ", nrow(batch))
  }
  
  recids <- unique(recids)
  return(data.frame(record_id = recids))
}

icdRecid2date <- function(table,icds,name) {
  # Input: table of recids; ICD code/s; name of new column
  # Output: table with column for date of initial dx
  # Note: can handle multiple ICDs for single dx
  
  icds <- paste0(icds, collapse = "|") 
  
  if (any(names(table) %in% 'recid')) {
    recid_colname <- 'recid'
  } else {recid_colname <- 'record_id'}
  
  files <- list.files(
    '/mmfs1/projects/cardio_darbar_chi/common/data/Diagnosis/',
    recursive = TRUE,
    full.names = TRUE
  )
  
  all_years <- c()
  for (i in c(1:length(files))) {
    batch <- arrow::open_dataset(files[i]) |>
      dplyr::filter(grepl(x = icd_code, pattern = icds) & 
                      record_id %in% table[[recid_colname]]) |>
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
  
  if (recid_colname == 'recid') {
    min_dates <- min_dates %>% rename(recid = record_id) # rename
  }
  
  min_dates <- min_dates %>% rename(!!name := date)
  
  final_table <- table %>%
    left_join(min_dates, by = recid_colname)  # Perform a left join on 'rec_id'
  
  return(final_table)
}

recid2demo <- function(ids, id_method = 'record_id', demo_path = '/home/dseaney2/cardio_darbar_chi_link/common/data/demographics-0.csv') {
  library(arrow)
  library(dplyr)
  demo = read.csv(demo_path)
  # demo = demo[,c('record_id','mrn','gender','race','ethnicity')]
  # demo = demo[,c('RECORD_ID','MRN','GENDER','RACE','ETHNICITY')]
  demo$mrn = as.numeric(demo$mrn)
  
  # If input is a vector, transform to a dataframe
  if (!class(ids) == 'data.frame') {
    table = data.frame(id = ids)
    colnames(table) = id_method
  } else {
    table <- ids
  }
  
  table <- left_join(table,demo, by = id_method)
  
  races <- c('AI_AN','Asian','NH_PI','Black','White','Other','Declined','NA')
  genders <- c('Unknown','Male','Female')
  ethnicities <- c('Not_Hispanic','Hispanic','Declined','NA')
  
  table$race <- races[table$race]
  table$ethnicity <- ethnicities[table$ethnicity]
  table$gender <- genders[table$gender + 1]
  
  race_summary <- array(NA,length(races))
  for (i in 1:length(races)) {
    race_summary[i] <- round(sum(table$race == races[i]) / nrow(table),3)
  }
  
  gender_summary <- array(NA,length(genders))
  for (i in 1:length(genders)) {
    gender_summary[i] <- round(sum(table$gender == genders[i]) / nrow(table),3)
  }
  
  ethnicity_summary <- array(NA,length(ethnicities))
  for (i in 1:length(ethnicities)) {
    ethnicity_summary[i] <- round(sum(table$ethnicity == ethnicities[i]) / nrow(table),3)
  }
  
  summary <- data.frame(matrix(c(nrow(table),gender_summary,race_summary,ethnicity_summary), nrow = 1))
  colnames(summary) <- c('total',genders,races,ethnicities)
  
  summary$total = nrow(table)
  
  output = list(table,summary)
  names(output) = c('table','summary')
  
  return(output)
}