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
  muse_log <- read.csv(muse_path)
  
  # FINISH CODING: 
  #   Load muse_path csv: read.csv
  #   For each muse id
  #     Locate muse file using muse_log
  #     Complete loop completed below:
  
  base_path <- '/mmfs1/projects/cardio_darbar_chi/common'
  xml <- read_xml(fs::path(base_path, muse_path, ext = 'xml'))
  diagnosis <- xml_text(xml_find_all(xml, ".//DiagnosisStatement"))
  
  return(diagnosis)
}
  
icd2mrn <- function(icds_path) {
  library(arrow)
  library(dplyr)
  library(vroom)
  
  diagnosis_path <- '/mmfs1/projects/cardio_darbar_chi/common/data/backup_clinical_data/Diagnosis/'
  demographics_path <- '/mmfs1/projects/cardio_darbar_chi/common/data/demographics-0.csv'
  
  # Input: path to .txt file with list of mrns
  icds <- readLines(icds_path)
  icds <- paste0("^", icds, collapse = "|")
  
  
  # Output: dataframe of mrns
  files <- list.files(diagnosis_path, 
                      recursive = TRUE, 
                      full.names = TRUE)
  
  recids <- c()
  for (i in c(1:length(files))) {
    batch <- arrow::open_dataset(files[i]) |>
      dplyr::filter(grepl(x = ICD_CODE, pattern = icds)) |>
      collect()
    
    recids <- c(recids,batch$RECORD_ID)
    recids <- unique(recids)
    
    print(paste('Finished',i,'.',nrow(batch)))
    rm(batch)
  }
  recids <- unique(recids)
  
  demographics <- vroom(demographics_path, col_select = c("RECORD_ID", "MRN"))
  output <- demographics |> filter(RECORD_ID %in% recids)
  return(output)
}

recid2demo <- function(ids, id_method = 'RECORD_ID', demo_path = '/home/dseaney2/cardio_darbar_chi_link/common/data/demographics-0.csv') {
  library(arrow)
  library(dplyr)
  demo = read.csv(demo_path)
  # demo = demo[,c('record_id','mrn','gender','race','ethnicity')]
  demo = demo[,c('RECORD_ID','MRN','GENDER','RACE','ETHNICITY')]
  demo$MRN = as.numeric(demo$MRN)
  
  table = data.frame(id = ids)
  colnames(table) = id_method
  
  table <- left_join(table,demo, by = id_method)
  
  races <- c('AI_AN','Asian','NH_PI','Black','White','Other','Declined','NA')
  genders <- c('Unknown','Male','Female')
  ethnicities <- c('Not_Hispanic','Hispanic','Declined','NA')
  
  table$RACE <- races[table$RACE]
  table$ETHNICITY <- ethnicities[table$ETHNICITY]
  table$GENDER <- genders[table$GENDER + 1]
  
  race_summary <- array(NA,length(races))
  for (i in 1:length(races)) {
    race_summary[i] <- round(sum(table$RACE == races[i]) / nrow(table),3)
  }
  
  gender_summary <- array(NA,length(genders))
  for (i in 1:length(genders)) {
    gender_summary[i] <- round(sum(table$GENDER == genders[i]) / nrow(table),3)
  }
  
  ethnicity_summary <- array(NA,length(ethnicities))
  for (i in 1:length(ethnicities)) {
    ethnicity_summary[i] <- round(sum(table$ETHNICITY == ethnicities[i]) / nrow(table),3)
  }
  
  summary <- data.frame(matrix(c(nrow(table),gender_summary,race_summary,ethnicity_summary), nrow = 1))
  colnames(summary) <- c('total',genders,races,ethnicities)
  
  summary$total = nrow(table)
  
  output = list(table,summary)
  names(output) = c('table','summary')
  
  return(output)
}