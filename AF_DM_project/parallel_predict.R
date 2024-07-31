# Prep parallel -----------------------------------------------------------

library(vroom)
library(fs)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)
library(EGM)
library(tools)
options(wfdb_path = 'wsl /usr/local/bin') 

# For each batch:
# Load batch
# For each lead:
# filter
# spectro
# predict
# save probabilties

outputfile_path <- 'C:/Users/darre/OneDrive/Documents/UICOM Research/AF_DM/predictions/af_predx/'
outputfile_list <-paste0('raw_probabilities',1:48,'.RData')

coordinator <- c(array(1,6),
                 array(2,6),
                 array(3,6),
                 array(4,6),
                 array(5,6),
                 array(6,6),
                 array(7,6),
                 array(8,6))

signal_batch_path <- 'C:/Users/darre/OneDrive/Documents/UICOM Research/AF_DM/signal/af_predx/'
signal_batch_names <- paste0('batch',1:48,'.RData')

model_path <- "C:/Users/darre/OneDrive/Documents/UICOM Research/ECG Segmentation/saved_model/"
model_names <- c("spectro_bilstm_i","spectro_bilstm_ii","spectro_bilstm_iii",
                 "spectro_bilstm_avr","spectro_bilstm_avl","spectro_bilstm_avf",
                 "spectro_bilstm_v1","spectro_bilstm_v2","spectro_bilstm_v3",
                 "spectro_bilstm_v4","spectro_bilstm_v5","spectro_bilstm_v6")

setwd("C:/Users/darre/OneDrive/Documents/UICOM Research/AF_DM/")


# Run parallel ------------------------------------------------------------

# request_cores <- 8
# cluster <- makeCluster(request_cores)
# registerDoParallel(cluster)
# ncores <- foreach::getDoParWorkers()
# 
# foreach(i = 1:ncores, .packages = c('EGM','keras')) %dopar% {
#   # Prep:
#   options(wfdb_path = 'wsl /usr/local/bin')
#   source('C:/Users/darre/OneDrive/Documents/UICOM Research/ECG Segmentation/ML_functions.R')
#   
#   
#   # Log file
#   logfile_name <- paste0('log_', i, '.log')
#   
#   # Set which files the core will iterate thru:
#   file_groups <- which(coordinator == i)
#   
#   for (batch in file_groups) {
#     write.table(paste0(
#       "Beginning batch",batch),logfile_name,
#       append = TRUE,col.names = FALSE)
#     
#     signal_file <- paste0(signal_batch_path,signal_batch_names[batch])
#     load(signal_file) # (with name wfdb_array)
#     
#     predictions_array <- array(0, c(dim(wfdb_array), 4))
#     
#     for (lead in 1:12) {
#       input <- wfdb_array[, , lead]
#       input_filtered <- filter_samples(input)
#       input_spectro <- build_spectrogram(input_filtered)
#       
#       model <- load_model_tf(paste0(model_path, model_names[lead]))
#       predictions_array[, , lead, ] <- model %>% predict(input_spectro)
#       
#       write.table(
#         paste0("Finished batch ",batch,' lead ',lead),
#         logfile_name,append = TRUE,col.names = FALSE)
#       
#     }
#     
#     output_file <- paste0(outputfile_path,outputfile_list[batch])
#     save(predictions_array, file = output_file) # 500
#     write.table(
#       paste0("Saved batch",batch),
#       logfile_name,append = TRUE,col.names = FALSE)
#     
#   }
#     
# }
# 
# stopCluster(cluster)


# Data clean up -----------------------------------------------------------

# Issues:
messup <- data.frame(batch = integer(),sample = integer(),time = integer(),lead = integer(),value = integer())
for (batch in 1:48) {
  load(batch_names[batch])
  raw <- which(is.na(predictions_array), arr.ind = T)
  raw <- raw[,c(1,3,4)]
  # df[(nrow(messup)+1):]
  raw <- raw[which(!duplicated(raw)),]
  raw <- cbind(array(batch,nrow(raw)),raw)
  messup <- rbind(messup,raw)
}

messup_clean <- messup[!duplicated(messup[,c(1:3)]),c(1:3)]
colnames(messup_clean) <- c('batch','sample','lead')

model_path <- "C:/Users/darre/OneDrive/Documents/UICOM Research/ECG Segmentation/saved_model/"
model_names <- c("spectro_bilstm_i","spectro_bilstm_ii","spectro_bilstm_iii",
                 "spectro_bilstm_avr","spectro_bilstm_avl","spectro_bilstm_avf",
                 "spectro_bilstm_v1","spectro_bilstm_v2","spectro_bilstm_v3",
                 "spectro_bilstm_v4","spectro_bilstm_v5","spectro_bilstm_v6")

outputfile_path <- 'C:/Users/darre/OneDrive/Documents/UICOM Research/AF_DM/raw_predictions/af_predx/'
outputfile_list <-paste0('raw_probabilities',1:48,'.RData')



# for (lead in 1:12) {
#   model <- load_model_tf(paste0(model_path, model_names[lead]))
#   
#   unique_batches <- messup_clean[messup_clean$lead == lead,]
#   for (i in 1:nrow(unique_batches)) {
#     
#     signal_file <- paste0(signal_batch_path,signal_batch_names[unique_batches$batch[i]])
#     load(signal_file) # (with name wfdb_array)
#     
#     output_file <- paste0(outputfile_path,outputfile_list[unique_batches$batch[i]])
#     load(output_file)
#     
#     #     save(predictions_array, file = output_file)
#     
#     input <- wfdb_array[unique_batches$sample[i], , lead]
#     input_filtered <- filter_samples(input)
#     input_spectro <- build_spectrogram(input_filtered)
#     
#     if (sum(is.na(predictions_array[unique_batches$sample[i], , lead, ])) > 0) {
#       predictions_array[unique_batches$sample[i], , lead, ] <- model %>% predict(input_spectro)
#       save(predictions_array, file = output_file)
#     } else {
#       print(paste('Aborting: not writing prediction for batch',unique_batches$batch[i],'sample',unique_batches$sample[i],'lead',lead,'. No NAN detected.'))
#     }
#     
#     print(paste('Finished lead',lead))
#     
#   }
# }



# for (batch in unique(messup_clean$batch)) {
#   
#   # batch <- 47
#   unique_batches <- messup_clean[messup_clean$batch == batch,]
#   unique_batches
#   
#   signal_file <- paste0(signal_batch_path,signal_batch_names[batch])
#   load(signal_file) # (with name wfdb_array)
#   
#   output_file <- paste0(outputfile_path,outputfile_list[batch])
#   load(output_file)
#   
#   # i <- 1
#   # plot_func(y = wfdb_array[unique_batches$sample[i], , unique_batches$lead[i]])
#   
#   for (i in 1:nrow(unique_batches)) {
#     model <- load_model_tf(paste0(model_path, model_names[unique_batches$lead[i]]))
#     
#     input <- wfdb_array[unique_batches$sample[i], , unique_batches$lead[i]]
#     input_filtered <- filter_samples(input)
#     input_spectro <- build_spectrogram(input_filtered)
#     
#     if (sum(is.na(predictions_array[unique_batches$sample[i], , unique_batches$lead[i], ])) > 0) {
#       predictions_array[unique_batches$sample[i], , unique_batches$lead[i], ] <- model %>% predict(input_spectro)
#     } else {
#       print(paste('Aborting: not writing prediction for batch',unique_batches$batch[i],'sample',unique_batches$sample[i],'lead',lead,'. No NAN detected.'))
#     }
#     
#     
#   }
#   save(predictions_array, file = output_file)
#   print(paste('Finished batch',batch))
# }
  

(messedup$batch - 1)*500 + messedup$sample
# Write master_lead -------------------------------------------------------

# Single sample example
# load('raw_predictions/af_predx/raw_probabilities1.RData')
# sample1 <- predictions_array[1,,,]
# sample1_reorg <- aperm(sample1, c(2,1,3))
# sample1_combined <- colSums(sample1_reorg)
# 
# master_lead1 <- apply(sample1_combined,1,which.max)
# 
batch_names <- paste0('raw_predictions/af_predx/raw_probabilities',1:48,'.RData')
master_lead_names <- paste0('master_lead/af_predx/master_lead',1:48,'.RData')



for (batch in 1:48) {
  load(batch_names[batch])
  # name of file: 'predictions_array'
  
  master_lead <- array(0, c(nrow(predictions_array), 5000))
  
  # start <- 500*(batch-1) + 1
  # stop <- start + nrow(predictions_array) - 1
  
  for (sample in 1:nrow(predictions_array)) {
    # Reorganize matrix dimmensions, add the 'lead' dimmension together, find max by column:
    master_lead[sample,] <- apply(colSums(aperm(predictions_array[sample,,,],c(2,1,3)),na.rm = TRUE),1,which.max) - 1
  }
  
  save(master_lead, file = master_lead_names[batch])
  
}



single_prediction <- apply(predictions_array[1,,1,],1,which.max) - 1

load('signal/af_predx/batch48.RData')

