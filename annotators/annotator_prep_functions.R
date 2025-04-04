# Build training and testing  sets---------------------------------------
prep_ludb <- function(lead,
                      annotator_style = 2,
                      split = 0.7) {
  # Load LUDB set
  load('ludb_set.RData')
  
  
  # Select lead**
  lead <- 1 # number 1 thru 12. Note: in wfdb$signal, first column is
  # Select annotation style**
  
  # Assign the annotation style for ML input:
  use_func <- paste0('ann_wfdb2continuous', annotator_style)
  #         ann_wfdb2continuous1: 1 0 0 0 1 0 0 0 2 0 0 2 ...
  #         ann_wfdb2continuous2: 1 1 1 1 1 0 0 0 2 2 2 2 ...
  #         ann_wfdb2continuous3: 1 2 2 2 3 0 0 0 4 5 5 6 ...
  
  # Split into train/test
  split <- 0.7 # % of samples going to training set
  
  leads <- c("i","ii","iii","avr","avl","avf","v1","v2","v3","v4","v5","v6")
  
  length <- 5000
  sampling_frequency <- 500
  annotation_function <- get(use_func) # Get the function object
  
  # Split sets:
  samples <- length(ludb_set)
  sample_size <- ceiling(split * samples)
  training_samples <- sample(samples, sample_size)
  testing_samples <- 1:samples
  testing_samples <- testing_samples[-training_samples]
  
  training_set <- ludb_set[c(training_samples)]
  testing_set <- ludb_set[c(testing_samples)]
  
  training_signal      <- array(NA, c(length(training_set), length))
  training_annotations <- array(NA, c(length(training_set), length))
  testing_signal       <- array(NA, c(length(testing_set), length))
  testing_annotations  <- array(NA, c(length(testing_set), length))
  
  # Prep for calculating indices for time range between 2 and 8 seconds in training set
  start_index <- 2 * sampling_frequency + 1
  end_index <- 8 * sampling_frequency
  interval_length <- 4 * sampling_frequency
  
  
  # Training Set
  for (sample in 1:length(training_set)) {
    # **Take a random 4 second interval of signal**:
    
    # Randomly select the starting point for the 4 second interval
    random_start <- sample(start_index:(end_index - interval_length), 1)
    # Define the 2-second interval
    random_interval <- random_start:(random_start + interval_length - 1)
    
    signal <- training_set[[sample]]$signal[[leads[lead]]]
    ann <- annotation_function(training_set[[sample]]$annotation[[leads[lead]]])
    
    signal[-random_interval] <- 0
    ann[-random_interval] <- 0
    
    
    
    # **Now, randomly shift the non zero signal/ann values**:
    
    # Determine a random starting point for the shifted interval
    vector_length <- length(signal)
    shifted_start <- sample(1:(vector_length - length(random_interval)), 1)
    
    # Generate the new shifted indices
    shifted_indices <- shifted_start:(shifted_start + length(random_interval) - 1)
    
    # Create new vectors initialized with zeros
    shifted_signal <- rep(0, vector_length)
    shifted_ann <- rep(0, vector_length)
    
    # Shift the non-zero values to the new interval
    shifted_signal[shifted_indices] <- signal[random_interval]
    shifted_ann[shifted_indices] <- ann[random_interval]
    
    # Add sample to matrix
    training_signal[sample, ] <- shifted_signal # signal
    training_annotations[sample, ] <- shifted_ann # ann
  }
  
  # Testing Set (no need to select a random interval)
  for (sample in 1:length(testing_set)) {
    signal <- testing_set[[sample]]$signal[[leads[lead]]]
    ann <- annotation_function(testing_set[[sample]]$annotation[[leads[lead]]])
    
    
    # Set values of first and last 1000 indices to 0 (could also do 1st and last annotated indices)
    signal[c(1:1000, 4001:5000)] <- 0
    ann[c(1:1000, 4001:5000)] <- 0
    
    # Add sample to matrix
    testing_signal[sample, ] <- signal
    testing_annotations[sample, ] <- ann
  }
  
  return(list(training_signal,training_annotations,
              testing_signal,testing_annotations,
              training_samples, testing_samples))
}

# ann_wfdb2continuous1 ----------------------------------------------------
#' @description format [1 0 0 0 1 0 0 0 0 2 0 0 0 2]
ann_wfdb2continuous1 <- function(object) {
  if (any(class(object) == 'egm')) {
    length <- nrow(object$signal)
    wfdb_ann <- object$annotation
  } else if (any(class(object) == 'annotation_table')) {
    wfdb_ann <- object
  }
  
  output <- rep(0,length)
  
  
  # P wave --------------------------
  # Find P indices, then verify onset and offset
  p_ind <- which(wfdb_ann$type == 'p')
  p_ind <- p_ind[wfdb_ann$type[p_ind - 1] == '(' & wfdb_ann$type[p_ind + 1] == ')']
  
  # For every p_ind, create array of integers like: p_onset : p_offset
  p_continuous <- unlist(lapply(p_ind, function(ind) {
    c(wfdb_ann$sample[(ind - 1)],wfdb_ann$sample[(ind + 1)])
  }))
  
  output[p_continuous] <- 1
  
  
  # QRS -------------------------- 
  # Find QRS indices, then verify onset and offset
  qrs_ind <- which(wfdb_ann$type == 'N')
  qrs_ind <- qrs_ind[wfdb_ann$type[qrs_ind - 1] == '(' & wfdb_ann$type[qrs_ind + 1] == ')']
  
  # For every qrs_ind, create array of integers like: qrs_onset : qrs_offset
  qrs_continuous <- unlist(lapply(qrs_ind, function(ind) {
    c(wfdb_ann$sample[(ind - 1)],wfdb_ann$sample[(ind + 1)])
  }))
  
  output[qrs_continuous] <- 2
  
  
  # T wave --------------------------  
  # Find T indices, then verify onset and offset
  t_ind <- which(wfdb_ann$type == 't')
  t_ind <- t_ind[wfdb_ann$type[t_ind - 1] == '(' & wfdb_ann$type[t_ind + 1] == ')']
  
  # For every t_ind, create array of integers like: t_onset : t_offset
  t_continuous <- unlist(lapply(t_ind, function(ind) {
    c(wfdb_ann$sample[(ind - 1)],wfdb_ann$sample[(ind + 1)])
  }))
  
  output[t_continuous] <- 3
  
  
  return(output)
  
}

# ann_wfdb2continuous2 ----------------------------------------------------
#' @description format [1 1 1 1 0 0 0 0 2 2 2]
#' @param object Either 'egm' or 'annotation_table'. If annotation table, should specify 'length'
#' 
#' @param length Length of corresponding ECG signal. Not necessary if 'object' is of class 'egm'. Assumed to be 5000.
ann_wfdb2continuous2 <- function(object, length=5000) {
  
  if (any(class(object) == 'egm')) {
    length <- nrow(object$signal)
    wfdb_ann <- object$annotation
  } else if (any(class(object) == 'annotation_table')) {
    wfdb_ann <- object
  }
  
  output <- rep(0,length)
  
  
  # Find P indices, then verify onset and offset
  p_ind <- which(wfdb_ann$type == 'p')
  p_ind <- p_ind[wfdb_ann$type[p_ind - 1] == '(' & wfdb_ann$type[p_ind + 1] == ')']
  
  # For every p_ind, create array of integers like: p_onset : p_offset
  p_continuous <- unlist(lapply(p_ind, function(ind) {
    wfdb_ann$sample[(ind - 1)]:wfdb_ann$sample[(ind + 1)]
  }))
  
  output[p_continuous] <- 1
  
  
  # Find QRS indices, then verify onset and offset
  qrs_ind <- which(wfdb_ann$type == 'N')
  qrs_ind <- qrs_ind[wfdb_ann$type[qrs_ind - 1] == '(' & wfdb_ann$type[qrs_ind + 1] == ')']
  
  # For every qrs_ind, create array of integers like: qrs_onset : qrs_offset
  qrs_continuous <- unlist(lapply(qrs_ind, function(ind) {
    wfdb_ann$sample[(ind - 1)]:wfdb_ann$sample[(ind + 1)]
  }))
  
  output[qrs_continuous] <- 2
  
  
  # Find T indices, then verify onset and offset
  t_ind <- which(wfdb_ann$type == 't')
  t_ind <- t_ind[wfdb_ann$type[t_ind - 1] == '(' & wfdb_ann$type[t_ind + 1] == ')']
  
  # For every t_ind, create array of integers like: t_onset : t_offset
  t_continuous <- unlist(lapply(t_ind, function(ind) {
    wfdb_ann$sample[(ind - 1)]:wfdb_ann$sample[(ind + 1)]
  }))
  
  output[t_continuous] <- 3
  
  
  return(output)
}


# ann_wfdb2continuous3 ----------------------------------------------------
#' @description format [1 2 2 3 0 0 0 0 4 5 5 6]
ann_wfdb2continuous3 <- function(object) {
  if (any(class(object) == 'egm')) {
    length <- nrow(object$signal)
    wfdb_ann <- object$annotation
  } else if (any(class(object) == 'annotation_table')) {
    wfdb_ann <- object
  }
  
  output <- rep(0,length)
  
  
  # Find P indices, then verify onset and offset
  p_ind <- which(wfdb_ann$type == 'p')
  p_ind <- p_ind[wfdb_ann$type[p_ind - 1] == '(' & wfdb_ann$type[p_ind + 1] == ')']
  
  # For every p_ind, create array of integers like: p_onset : p_offset
  p_continuous <- unlist(lapply(p_ind, function(ind) {
    wfdb_ann$sample[(ind - 1)]:wfdb_ann$sample[(ind + 1)]
  }))
  # Note P onset and offset values
  p_onset <- wfdb_ann$sample[p_ind-1]
  p_offset <- wfdb_ann$sample[p_ind+1]
  
  output[p_continuous] <- 2
  output[p_onset] <- 1
  output[p_offset] <- 3
  
  
  # Find QRS indices, then verify onset and offset
  qrs_ind <- which(wfdb_ann$type == 'N')
  qrs_ind <- qrs_ind[wfdb_ann$type[qrs_ind - 1] == '(' & wfdb_ann$type[qrs_ind + 1] == ')']
  
  # For every qrs_ind, create array of integers like: qrs_onset : qrs_offset
  qrs_continuous <- unlist(lapply(qrs_ind, function(ind) {
    wfdb_ann$sample[(ind - 1)]:wfdb_ann$sample[(ind + 1)]
  }))
  # Note QRS onset and offset values
  qrs_onset <- wfdb_ann$sample[qrs_ind-1]
  qrs_offset <- wfdb_ann$sample[qrs_ind+1]
  
  output[qrs_continuous] <- 5
  output[qrs_onset] <- 4
  output[qrs_offset] <- 6
  
  
  # Find T indices, then verify onset and offset
  t_ind <- which(wfdb_ann$type == 't')
  t_ind <- t_ind[wfdb_ann$type[t_ind - 1] == '(' & wfdb_ann$type[t_ind + 1] == ')']
  
  # For every t_ind, create array of integers like: t_onset : t_offset
  t_continuous <- unlist(lapply(t_ind, function(ind) {
    wfdb_ann$sample[(ind - 1)]:wfdb_ann$sample[(ind + 1)]
  }))
  # Note T onset and offset values
  t_onset <- wfdb_ann$sample[t_ind-1]
  t_offset <- wfdb_ann$sample[t_ind+1]
  
  output[t_continuous] <- 8
  output[t_onset] <- 7
  output[t_offset] <- 9
  
  return(output)
  
}