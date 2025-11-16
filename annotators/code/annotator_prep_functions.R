# Build training and testing  sets---------------------------------------
#' @descrption Input: custom list of LUDB annotations and signal. Splits into training/testing. 
#' Takes random continuous segment between 2 and 8 seconds, zeros the rest. Randomly shifts the segment across the vector
#' Then normalizes the signal for each sample.
#' 
#' @param dilate_range For training set: 2 element vector. Defines the range the samples will be randomly. Recommend: c(0.03,0.05) (ie 3 to 5%)
#' 
#' @param max_noise For training set: maximum standard deviation of the added Gaussian noise. ECGs are scaled from 0 to 1. Recommend: between 0.01 to 0.07
#' 
#' @param rounds number of times the training dataset will be duplicated. Note: a value of 2 will double the length of the training set. A value of 1 will not change the length. Recommend: 3 or 4
#' 
#' @param filter if the signal should be run thru a bandpass filter to remove baseline wander, high frequency noise. Recommend: FALSE
#' 
#' @param number_of_derivs if/how many derivatives should be included. If set to 1: input would be the original signal AND the first derviative (matrix form)
#' 
#' @param split percentage of samples to go in the training set (remaining go into testing set)
#' 
#' @param data data input. default 'ludb' sources the LUDB dataset. Alternative is to set input as the dataset you would like to be split
prep_ludb <- function(lead,
                      annotator_style = 2,
                      split = 0.7,
                      dilate_range = c(0.03,0.05),
                      max_noise = 0.05,
                      rounds = 1,
                      filter,
                      number_of_derivs,
                      mask_value = -1,
                      normalize = TRUE,
                      data = 'ludb') {
  library(stats)
  library(signal)
  
  # Load input data
  if (any(data == 'ludb')) {
    load('../ludb_set.RData')
    input_data <- ludb_set
    leads <- c("i","ii","iii","avr","avl","avf","v1","v2","v3","v4","v5","v6")
  } else {
    input_data <- data
    leads <- c('I','II','III','AVR','AVL','AVF','V1','V2','V3','V4','V5','V6')
  }

  
  # Assign the annotation style for ML input:
  use_func <- paste0('ann_wfdb2continuous', annotator_style)
  
  #         ann_wfdb2continuous1: 1 0 0 0 1 0 0 0 2 0 0 2 ...
  #         ann_wfdb2continuous2: 1 1 1 1 1 0 0 0 2 2 2 2 ...
  #         ann_wfdb2continuous3: 1 2 2 2 3 0 0 0 4 5 5 6 ...
  

  
  length <- 5000
  sampling_frequency <- 500
  annotation_function <- get(use_func) # Get the function object
  
  # Split sets:
  samples <- length(input_data)
  sample_size <- ceiling(split * samples)
  training_samples <- sample(samples, sample_size)
  testing_samples <- 1:samples
  testing_samples <- testing_samples[-training_samples]
  
  training_set <- input_data[c(training_samples)]
  testing_set <- input_data[c(testing_samples)]
  
  training_signal      <- array(NA, c(length(training_set)*rounds, length, number_of_derivs+1))
  training_annotations <- array(NA, c(length(training_set)*rounds, length))
  testing_signal       <- array(NA, c(length(testing_set), length, number_of_derivs+1))
  testing_annotations  <- array(NA, c(length(testing_set), length))
  
  # Prep for calculating indices for time range between 2 and 8 seconds in training set
  start_index <- 2 * sampling_frequency + 1
  end_index <- 8 * sampling_frequency
  interval_length <- 4 * sampling_frequency
  
  
  # Training Set
  for (round in 1:rounds) {
    for (sample in 1:length(training_set)) {
      # **Take a random 4 second interval of signal**:
      
      # Randomly select the starting point for the 4 second interval
      random_start <- sample(start_index:(end_index - interval_length), 1)
      # Define the 4-second interval
      random_interval <- random_start:(random_start + interval_length - 1)
      
      signal <- training_set[[sample]]$signal[[leads[lead]]]
      ann <- annotation_function(training_set[[sample]]$annotation[[leads[lead]]])
      
      if (filter) {
        signal <- ecg_filter(signal = signal)
      }
      
      # Set indices outside of interval to 0
      # signal[-random_interval] <- 0
      # ann[-random_interval] <- 0
      
      # Random interval signal/ann vectors
      signal_iso <- signal[random_interval]
      ann_iso <- ann[random_interval]
      
      
      # Compress signal:
      # Generate a random compression/dilation factor
      
      if (round != 1) {
        compression_factor <- 1 + sample(c(-1, 1), size = 1) * runif(1, min = dilate_range[1], max = dilate_range[2]) # Random factor between 1 +/- (3 to 5)%
        new_length <- round(length(signal_iso) * compression_factor) # Define new length for the signal based on compression/dilation factor
        # Resample the ECG signal to the new length
        compressed_signal <- resample(signal_iso, p = new_length, q = length(signal_iso))
        # Cannot use same interpolation to do ann vector, which contains *integers* (not decimals). ie, must ensure termini of QRS, T don't become 1s or 2s if rounding
        new_indices <- seq(from = 1,
                           to = length(ann_iso),
                           length.out = new_length) # Create new indices from 1 to length(ann_iso) with the new length
        compressed_ann <- ann_iso[round(new_indices)] # For each new index, use the nearest original annotation
      } else {
        # if on the 1st round, do not dilate/compress, or add noise
        compressed_signal <- signal_iso
        compressed_ann <- ann_iso
        new_length <- length(compressed_signal)
      }
      
      
      
      
      # Generate Gaussian random noise
      # Compress signal:
      # Generate a random compression/dilation factor
      if (round != 1) {
        # if on the 1st round, do not dilate/compress, or add noise
        sigma <- runif(1, min = 0, max = max_noise) # vary the noise SD from 0 to max_noise parameter
        range <- range(compressed_signal)
        noise <- rnorm(length(compressed_signal),
                       mean = 0,
                       sd = sigma * abs(range[1] - range[2])) # create random noise
        compressed_signal <- compressed_signal + noise # Add the noise to the original ECG signal
      } 
      
      # Normalize signal
      if (normalize) {
        compressed_signal <- (compressed_signal - min(compressed_signal)) / (max(compressed_signal) - min(compressed_signal)) * 100
      }
      
      # Now, randomly place the random signal/ann interval:
      # Determine a random starting point for the shifted interval
      vector_length <- length(signal)
      shifted_start <- sample(1:(vector_length - new_length), 1)
      
      # Generate the new shifted indices
      shifted_indices <- shifted_start:(shifted_start + new_length - 1)
      
      # Create new vectors initialized with zeros
      shifted_signal <- array(mask_value, c(vector_length,number_of_derivs+1))
      shifted_ann <- rep(0, vector_length)
      
      # Add the derivatives as needed
      if (number_of_derivs > 0) {
        compressed_signal <- add_derivs(signal = compressed_signal,
                                        number_of_derivs = number_of_derivs,
                                        mask_value = mask_value)
      }
      
      # Shift the non-zero values to the new interval
      shifted_signal[shifted_indices,] <- compressed_signal
      shifted_ann[shifted_indices] <- compressed_ann
 
      
      # Add sample to matrix
      index <- sample + (round-1)*length(training_set)
      training_signal[index, , ] <- shifted_signal # signal
      training_annotations[index, ] <- shifted_ann # ann
    }
  }
  
  # Testing Set (no need to select a random interval)
  for (sample in 1:length(testing_set)) {
    signal <- testing_set[[sample]]$signal[[leads[lead]]]
    ann <- annotation_function(testing_set[[sample]]$annotation[[leads[lead]]])
    
    # Filter signal:
    signal <- ecg_filter(signal = signal)
    
    # Normalize signal
    if (normalize) {
      signal <- (signal - min(signal)) / (max(signal) - min(signal)) * 100
    }
    
    if (number_of_derivs > 0) {
      signal <- add_derivs(signal = signal,
                                      number_of_derivs = number_of_derivs,
                                      mask_value = mask_value)
    }
    
    # Set values of first and last 1000 indices to 0 (could also do 1st and last annotated indices)
    signal[c(1:1000, 4001:5000),] <- mask_value
    ann[c(1:1000, 4001:5000)] <- 0
    
    # Add sample to matrix
    testing_signal[sample, , ] <- signal
    testing_annotations[sample, ] <- ann
  }
  
  # Prepare output
  output <- list(
    training_signal = training_signal,
    training_annotations = training_annotations,
    testing_signal = testing_signal,
    testing_annotations = testing_annotations,
    training_samples = training_samples,
    testing_samples = testing_samples,
    mask_value = mask_value,
    normalize = normalize
  )
  
  return(output)
  
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
  # } else if (any(class(object) == 'annotation_table')) {
  #   wfdb_ann <- object
  } else {
    wfdb_ann <- object
  }
  
  if (nrow(object) == 0) {
    output <- rep(NA,length)
    return(output)
    warning('Annotation file is empty')
    break
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


# ann wfdb/compact --------------------------------------------------------
ann_wfdb2compact <- function(ann_wfdb) {
  library(dplyr)
  class(ann_wfdb) <- "data.frame"
  ann_wfdb <- ann_wfdb %>% mutate(idx = row_number())
  
  # Extract waves
  ann_compact <- ann_wfdb %>%
    mutate(
      type_lead = lead(type),
      type_lag = lag(type),
      sample_lead = lead(sample),
      sample_lag = lag(sample)
    ) %>%
    filter(type %in% c("p", "N", "t"), type_lag == "(", type_lead == ")") %>%
    transmute(
      type,
      onset = sample_lag,
      peak = sample,
      offset = sample_lead
    )
  return(ann_compact)
}

ann_compact2continuous <- function(ann_compact) {
  # Assuming your dataframe is called wave_table
  ann_continuous <- rep(0, 5000)
  
  for (i in seq_len(nrow(ann_compact))) {
    start <- ann_compact$onset[i]
    end <- ann_compact$offset[i]
    
    value <- switch(ann_compact$type[i],
                    "p" = 1,
                    "N" = 2,
                    "t" = 3,
                    "V" = 4)
    
    ann_continuous[start:end] <- value
  }
  return(ann_continuous)
}

# Filter ------------------------------------------------------------------
ecg_filter <- function(signal, frequency = 500, low = 0.5, high = 40) {
  library(signal)
  if (is.vector(signal)) {
    signal <- matrix(signal, nrow = 1)
  }
  
  # Apply the filter row-wise
  t(apply(signal, 1, function(row) {
    butterFilter <- butter(4, c(low / (frequency / 2), high / (frequency / 2)), type = "pass")
    
    # Eliminate phase distortion
    filtfilt(butterFilter, row)
    
  }))
}


# Difference/Derivative functions -----------------------------------------
add_derivs <- function(signal, number_of_derivs = 2,mask_value=0) {
  if (number_of_derivs == 0) {
    return(signal)
    break
  }
  # If input is of size 1 x 5000, change to vector:
  if (length(dim(signal)) > 1) {
    signal <- c(signal)
    
  }
  
  # Do first derivative manually
  diff_signal <- cbind(signal, c(mask_value, diff(signal)))
  
  if (number_of_derivs > 1) {
    for (i in 2:number_of_derivs) {
      last_deriv <- c(mask_value, diff(diff_signal[, ncol(diff_signal)]))
      diff_signal <- cbind(diff_signal, last_deriv)
    }
  }
  
  return(diff_signal)
}


# Confusion matrix analysis ----------------------------------------------------------------
confusion_analysis <- function(predictions = predictions_integer,
                               actual = testing_annotations) {
  # Quick analysis of predicted values
  
  library(caret)
  levels <- sort(unique(as.vector(actual))) # number of value types must be the same when factoring
  
  confusion_matrix <- confusionMatrix(
    data      = factor(predictions, levels = levels),
    reference = factor(actual, levels = levels)
  )
  
  confusion_table <- t(t(confusion_matrix$table) / colSums(confusion_matrix$table)) * 100
  
  print(round(confusion_matrix$byClass,3))
  print(round(confusion_table, 1))
  # return(confusion_matrix)
  
  output <- list(round(confusion_matrix$byClass,2), round(confusion_table,1))
  names(output) <- c('byClass','table')
  
  # output <- list(confusion_matrix,confusion_table)
  names(output) <- c('matrix','table')
  return(output)
  
}


# LUDB wave analysis --------------------------------------------------------

# Count number of P, QRS, T waves in each sample, cross compare
# Output: which samples / number of samples which do not have the right number of waves

# For each wave, compare onset / offset times
check_ann_prog_ludb <- function(predicted, ludb) {
  # Check annotation progression against LUDB (ie against a ground truth)
  
  # Input: predicted and LUDB annotations. Either matrix form or LUDB
  # Psuedocode:
  #   move into wfdb format as needed
  #   **use LUDB R-R intervals as ground truth reference markers**
  #   define a tolerance window- different for QRS vs p/t...
  #     Could do a tolerance test for both onset offset, 
  #     or could fine the midpoint between them, then compare length?
  
  #   verify predicted QRS are within tolerance 
  #     also verify number of QRS are the same 
  #     if extra QRS outside the tolerance, note. 
  #     if multiple QRS within the tolerance, note
  
  #   Check p and t waves
  #     for each R-R: 
  #       verify correct number and p/t-waves. If extra/lacking, note
  #       verify they're within tolerance. If outside tolerance, note
  #       should do both tasks for indicies prior to 1st R, after last R too
  
  # If annotation is in matrix format (not WFDB), convert to wfdb
  if (class(predicted) == 'numeric') {
    predicted <- ann_continuous2wfdb(predicted)
  }
  
  if (class(ludb) == 'numeric') {
    ludb <- ann_continuous2wfdb(ludb)
  }
  
  # To compare ground truths, remove predicted annotations prior to first LUDB annotation, after last LUDB annotation
  
  
  frequency <- 500
  # QRS tolerance: avg duration = 120 ms. --> tolerance = 60 ms
  qrs_tolerance <- 60 / 1000 * frequency
  p_tolerance <- 60 / 1000 * frequency
  t_tolerance <- 60 / 1000 * frequency
  
  
  tolerance <- c(qrs_tolerance,p_tolerance,t_tolerance)
  waves <- c('N','p','t')
  # Initialize points counter (for each wave type)
  points <- array(0,3)
  duplicate_points <- array(0,3)
  
  # For each wave type (QRS/P/T)
  for (wave in 1:length(waves)) {
    # Define a time stamp for a wave via the mean of the onset and offset
    pred_wave <- unlist(lapply(which(predicted$type == waves[wave]), function(idx) {
      round(mean(predicted$sample[idx + 1], predicted$sample[idx + 1]))
    }))
    
    ludb_wave <- unlist(lapply(which(ludb$type == waves[wave]), function(idx) {
      round(mean(ludb$sample[idx + 1], ludb$sample[idx + 1]))
    }))
    
    # Initialize a vector to track paired ground truth indices
    paired_ludb <- rep(FALSE, length(ludb_wave)) # FALSE means unpaired initially

    
    # Loop through each predicted QRS index
    for (i in seq_along(pred_wave)) {
      pred_time <- pred_wave[i]
      
      # Find all ground truth indices within the tolerance range
      valid_pairs <- which(abs(ludb_wave - pred_time) <= qrs_tolerance)
      
      if (length(valid_pairs) == 0) {
        # No valid pairing found, mark as unpaired
        points[wave] <- points[wave] + 1
        
        # If wave is duplicated, add to separate count:
        # Find index in original wfdb dataframe:
        row <- which(predicted$type == waves[wave])[i]
        # Safeguard for undefined indices
        valid_row_minus_3 <- (row - 3 >= 1) 
        valid_row_plus_3 <- (row + 3 <= nrow(predicted))
        # Modified condition with checks for valid indices
        if ((valid_row_minus_3 && predicted$type[row - 3] == waves[wave]) ||
            (valid_row_plus_3 && predicted$type[row + 3] == waves[wave])) {
          duplicate_points[wave] <- duplicate_points[wave] + 1
        }
        
        
      } else {
        # Check if any of the valid pairs are already paired
        paired_indices <- which(paired_ludb[valid_pairs])
        unpaired_indices <- setdiff(valid_pairs, paired_indices)
        
        if (length(unpaired_indices) > 0) {
          # Pair to the first unpaired ground truth event
          paired_ludb[unpaired_indices[1]] <- TRUE
        } else {
          # If all valid pairs are already paired, add +1 point
          points[wave] <- points[wave] + 1
          
          # If wave is duplicated, add to separate count:
          # Find index in original wfdb dataframe:
          row <- which(predicted$type == waves[wave])[i]
          # Safeguard for undefined indices
          valid_row_minus_2 <- (row - 2 >= 1) 
          valid_row_plus_2 <- (row + 2 <= nrow(predicted))
          # Modified condition with checks for valid indices
          if ((valid_row_minus_2 && predicted$type[row - 2] == waves[wave]) ||
              (valid_row_plus_2 && predicted$type[row + 2] == waves[wave])) {
            duplicate_points[wave] <- duplicate_points[wave] + 1
          }
          
        }
      }
    }
    
    # Unneeded- paired_ludb
    # # Count additional pairings to the same ground truth
    # points[wave] <- points[wave] + sum(paired_ludb > 1)
    
  }
  
    
  points_df <- data.frame(
    qrs_unpair = points[1] - duplicate_points[1], 
    p_unpair = points[2] - duplicate_points[2], 
    t_unpair = points[3] - duplicate_points[3],
    qrs_dupl = duplicate_points[1], 
    p_dupl = duplicate_points[2], 
    t_dupl = duplicate_points[3]
  )

  return(points_df)
  
}

check_ann_prog_ludb_helper <- function(testing_annotations = testing_annotations, 
                                       predictions_integer = predictions_integer) {
  table <- c()
  for (sample in 1:nrow(testing_annotations)) {
    table <- rbind(table,check_ann_prog_ludb(ludb = testing_annotations[sample,], 
                                               predicted = predictions_integer[sample,])) 
  }
  
  samples_w_missed_wave <- data.frame(
    p = round((sum(table$p_unpair) + sum(table$p_dupl)) / nrow(table), 3),
    qrs = round((sum(table$qrs_unpair) + sum(table$qrs_dupl)) / nrow(table), 3),
    t = round((sum(table$t_unpair) + sum(table$t_dupl)) / nrow(table), 3)
  )
  
  missed_waves_per_sample <- data.frame(
    p = round(length(which(table$p_unpair != 0 | table$p_dupl != 0)) / nrow(table), 3),
    qrs = round(length(which(table$qrs_unpair != 0 |table$qrs_dupl != 0)) / nrow(table), 3),
    t = round(length(which(table$t_unpair != 0 | table$t_dupl != 0)) / nrow(table), 3)
  )
  
  
  output <- list(samples_w_missed_wave, missed_waves_per_sample)
  names(output) <- c('samples_w_missed_wave','missed_waves_per_sample')
  print(output)
  return(output)
}

ann_continuous2wfdb <- function(annotations, Fs = 500) {
  # Not up to date, largely defunct for now. Used to convert index to index
  # annotation (value for each time point) to wave onset/peak/offset format. 
  # WFDB format
  
  if (sum(is.na(annotations)) == length(annotations)) {
    annotation_table <- data.frame(time = c(), 
                                   sample = c(), 
                                   type = c(), 
                                   subtype = c(), 
                                   channel = c(), 
                                   number = c()
                                   )
    return(annotation_table)
    break
  }
  
  if (any(annotations == 1)) {
    p_waves <- which(annotations == 1)
    p_change <- ((p_waves[-1] - p_waves[1:(length(p_waves) - 1)]))
    new_pwave <- c(1, ((which(p_change != 1) + 1)))
    new_pwave <- c(new_pwave, length(p_waves) + 1) # account for last value
    
    p_on <- c()
    p_off <- c()
    p_max <- c()
    
    for (i in 1:(length(new_pwave) - 1)) {
      p_on <- c(p_on, p_waves[new_pwave[i]])
      p_off <- c(p_off, p_waves[new_pwave[i + 1] - 1])
      p_max <- c(p_max, round(mean(c(p_on[i], p_off[i])))) # temporary, will need to find peak point
    }
  }
    
  if (any(annotations == 2)) {
    qrs_waves <- which(annotations == 2)
    qrs_change <- ((qrs_waves[-1] - qrs_waves[1:(length(qrs_waves) - 1)]))
    new_qrswave <- c(1, ((which(qrs_change != 1) + 1)))
    new_qrswave <- c(new_qrswave, length(qrs_waves) + 1) # account for last value
    
    qrs_on <- c()
    qrs_off <- c()
    qrs_max <- c()
    
    for (i in 1:(length(new_qrswave) - 1)) {
      qrs_on <- c(qrs_on, qrs_waves[new_qrswave[i]])
      qrs_off <- c(qrs_off, qrs_waves[new_qrswave[i + 1] - 1])
      qrs_max <- c(qrs_max, round(mean(c(
        qrs_on[i], qrs_off[i]
      )))) # temporary, will need to find peak point
    }
  }
  
  if (any(annotations == 3)) {
    t_waves <- which(annotations == 3)
    t_change <- ((t_waves[-1] - t_waves[1:(length(t_waves) - 1)]))
    new_twave <- c(1, ((which(t_change != 1) + 1)))
    new_twave <- c(new_twave, length(t_waves) + 1) # account for last value
    
    t_on <- c()
    t_off <- c()
    t_max <- c()
    
    for (i in 1:(length(new_twave) - 1)) {
      t_on <- c(t_on, t_waves[new_twave[i]])
      t_off <- c(t_off, t_waves[new_twave[i + 1] - 1])
      t_max <- c(t_max, round(mean(c(t_on[i], t_off[i])))) # temporary, will need to find peak point
    }
  }
  
  t <- (1 : length(annotations) ) / Fs
  
  
  # Build Table
  type <- c()

  sample <- c()
  if (exists('p_on')) {
    type <- c(type,array(c("(","p",")"),length(p_on)*3))
    for (i in 1:length(p_on)) {
      sample <- c(sample, p_on[i], p_max[i], p_off[i])
    }
  }
  if (exists('qrs_on')) {
    type <- c(type,array(c("(","N",")"),length(qrs_on)*3))
    for (i in 1:length(qrs_on)) {
      sample <- c(sample, qrs_on[i], qrs_max[i], qrs_off[i])
    }
  }
  
  if (exists('t_on')) {
    type <- c(type,array(c("(","t",")"),length(t_on)*3))
    for (i in 1:length(t_on)) {
      sample <- c(sample, t_on[i], t_max[i], t_off[i])
    }
  }
  
  time <- t[sample]
  subtype <- array(0,length(sample))
  channel <- array(0,length(sample))
  number <- array(0,length(sample))
  
  annotation_table <- data.frame(
    time = time,
    sample = sample,
    type = type,
    subtype = subtype,
    channel = channel,
    number = number)
  
  class(annotation_table) <- c('data.frame','annotation_table')
  
  if (nrow(annotation_table > 0)) {
    annotation_table <- annotation_table[order(annotation_table$sample), ]
    row.names(annotation_table) <- NULL
  }
  
  return(annotation_table)
}



# UIH relative analysis ---------------------------------------------------
# ie- functions when we don't have an answer key

# Check if waves go: P > QRS > T > P...
check_ann_prog <- function(annotation) {
  ann <- ann_continuous2wfdb(annotation)
  
  # Remove first and last 0.5s
  ann <- ann[(ann$sample <= 4750) & (ann$sample >= 250),]
  
  progression <- ann$type[ann$type %in% c('p','N','t')]
  
  # Create logical checks
  check_p <- which(progression == "p")  # Indices of "p"
  check_N <- which(progression == "N")  # Indices of "N"
  check_t <- which(progression == "t")  # Indices of "t"
  
  # Ensure "p" is followed by "N"
  valid_pN <- all(progression[check_p[check_p < length(progression)] + 1] == "N")
  
  # Ensure "N" is followed by "t"
  valid_Nt <- all(progression[check_N[check_N < length(progression)] + 1] == "t")
  
  # Ensure "t" is followed by "p"
  valid_tp <- all(progression[check_t[check_t < length(progression)] + 1] == "p")
  
  # Final check
  if (valid_pN & valid_Nt & valid_tp & length(progression > 1)) {
    # Check if statements, and check the annotation has more than 1 wave
    
    # print("The progression follows the expected pattern.")
    output <- 1
  } else {
    # print("The progression does NOT follow the expected pattern.")
    output <- 0
  }
  return(output)
}

# If there are two waves of the same type (ie pwave) which are within 10 indicies,
# and those indicies are labelled as 0 (non-waves), convert the 0s to pwave values
fill_wave_gaps <- function(annotation = predictions_integer, max_gap = 10) {
  # Find nonzero indices
  nonzero_indices <- which(annotation != 0)
  
  # Iterate through the nonzero indices
  for (i in seq_along(nonzero_indices)) {
    if (i == length(nonzero_indices)) break  # Stop if at last index
    
    current_val <- annotation[nonzero_indices[i]]
    next_idx <- nonzero_indices[i + 1]
    
    # Check if the next wave is of the same value and gap is ≤ max_gap
    if (annotation[next_idx] == current_val && (next_idx - nonzero_indices[i] <= max_gap)) {
      annotation[(nonzero_indices[i] + 1):(next_idx - 1)] <- current_val
    }
  }
  
  return(annotation)
}

remove_short_waves <- function(annotation, max_wave_length=10) {
  # annotation: vector of integers, not in wfdb format
  
  # Run-length encoding
  rle_obj <- rle(annotation)
  
  # Identify runs that are P/QRS/T (1,2,3) AND shorter than threshold
  too_short <- rle_obj$values %in% c(1, 2, 3) & rle_obj$lengths < max_wave_length
  
  # Replace those runs with 0
  rle_obj$values[too_short] <- 0
  
  # Reconstruct the vector
  cleaned <- inverse.rle(rle_obj)
  return(cleaned)
}


check_ann_prog_RPeaks <- function(signal,
                                  predictions_integer, 
                                  pr_tolerance = 200,
                                  rt_tolerance = 170,
                                  terminal_boundary = 250
) {
  # pr_tolerance <- 0.2 * 500 * 2 # 0.2 sec, fudge factor of 100%
  # rt_tolerance <- 0.170 * 500 * 2 # ST: 0.120 sec, RS: ~0.05 sec, fudge factor of 100%
  
  Rpeaks_lower_bound <- terminal_boundary
  Rpeaks_upper_bound <- 5000 - terminal_boundary
  
  # Convert predictions_integer sample to wfdb format
  ann <- ann_continuous2wfdb(predictions_integer)
  sig <- signal
  
  ## Identify and match QRS: ##
  
  Rpeaks <- EGM::detect_QRS(ecg_filter(sig), frequency = 500)
  # Remove Rpeaks if they are too close to the beginning/end
  Rpeaks <- Rpeaks[Rpeaks > Rpeaks_lower_bound &
                     Rpeaks < Rpeaks_upper_bound]
  
  # Find wave indicies/onsets/offsets
  QRS_wave_indices <- which(ann$type == "N")
  QRS_wave_onsets <- ann$sample[QRS_wave_indices - 1]
  QRS_wave_offsets <- ann$sample[QRS_wave_indices + 1]
  p_indices <- which(ann$type == "p")
  p_wave_onsets <- ann$sample[p_indices - 1]
  t_indices <- which(ann$type == "t")
  t_wave_onsets <- ann$sample[t_indices - 1]
  
  # If QRS onset is <250, offset > 4750, remove
  if (any(QRS_wave_onsets < Rpeaks_lower_bound)) {
    logic <- QRS_wave_onsets < Rpeaks_lower_bound
    QRS_wave_onsets <- QRS_wave_onsets[!logic]
    QRS_wave_offsets <- QRS_wave_offsets[!logic]
  }
  if (any(QRS_wave_offsets > Rpeaks_upper_bound)) {
    logic <- QRS_wave_offsets > Rpeaks_upper_bound
    QRS_wave_onsets <- QRS_wave_onsets[!logic]
    QRS_wave_offsets <- QRS_wave_offsets[!logic]
  }
  
  matched_Rpeaks <- list()
  unmatched_Rpeaks <- Rpeaks  # Start with all Rpeaks, remove matched ones later
  unmatched_QRS_waves <- data.frame(onset = integer(), offset = integer())
  
  # Find the matched/unmatched QRS/Rpeaks
  for (i in seq_along(QRS_wave_onsets)) {
    # Find Rpeaks within onset/offset range
    contained <- Rpeaks[Rpeaks >= QRS_wave_onsets[i] &
                          Rpeaks <= QRS_wave_offsets[i]]
    
    if (length(contained) > 0) {
      matched_Rpeaks[[i]] <- contained
      unmatched_Rpeaks <- setdiff(unmatched_Rpeaks, contained)  # Remove matched Rpeaks
    } else {
      unmatched_QRS_waves <- rbind(
        unmatched_QRS_waves,
        data.frame(onset = QRS_wave_onsets[i], offset = QRS_wave_offsets[i])
      )
    }
  }
  
  # matched_waves <- data.frame(Rpeaks = unlist(matched_Rpeaks))
  matched_waves <- data.frame(Rpeaks = unlist(Rpeaks))
  matched_waves$p <- NA
  matched_waves$t <- NA
  
  # Find missed pwaves, twaves
  
  for (i in seq_along(matched_waves$Rpeaks)) {
    peak <- matched_waves$Rpeaks[i]
    
    # Check for P-wave onset within pr_tolerance
    p_match <- p_wave_onsets[p_wave_onsets >= (peak - pr_tolerance) &
                               p_wave_onsets < peak]
    
    # If there are multiple matches, pick the closest one
    if (length(p_match >= 1)) {
      matched_waves$p[i] <- max(p_match)
    }
    
    # Check for T-wave onset within rt_tolerance
    t_match <- t_wave_onsets[t_wave_onsets <= (peak + rt_tolerance) &
                               t_wave_onsets > peak]
    
    # If there are multiple matches, pick the closest one
    if (length(t_match >= 1)) {
      matched_waves$t[i] <- min(t_match)
    }
  }
  
  # Find duplicated pwaves, twaves
  duplicate_pwaves <- c()
  duplicate_twaves <- c()
  
  for (i in seq_along(matched_waves$Rpeaks)) {
    if (i > 1) {
      # Skip first peak
      prev_peak <- matched_waves$Rpeaks[i - 1]
      curr_peak <- matched_waves$Rpeaks[i]
      pwaves <- p_wave_onsets[p_wave_onsets > prev_peak &
                                p_wave_onsets < curr_peak]
      
      if (length(pwaves) > 1) {
        unmatched_p <- setdiff(pwaves, matched_waves$p[i])
        duplicate_pwaves <- c(duplicate_pwaves, unmatched_p)
      }
    }
    
    if (i < nrow(matched_waves)) {
      # Skip last peak
      curr_peak <- matched_waves$Rpeaks[i]
      next_peak <- matched_waves$Rpeaks[i + 1]
      twaves <- t_wave_onsets[t_wave_onsets > curr_peak &
                                t_wave_onsets < next_peak]
      
      if (length(twaves) > 1) {
        unmatched_t <- setdiff(twaves, matched_waves$t[i])
        duplicate_twaves <- c(duplicate_twaves, unmatched_t)
      }
    }
  }
  
  output <- data.frame(
    missed_QRS = length(unmatched_Rpeaks) / nrow(matched_waves),
    missed_P = sum(is.na(matched_waves$p)) / nrow(matched_waves),
    missed_T = sum(is.na(matched_waves$t)) / nrow(matched_waves),
    duplicate_QRS = nrow(unmatched_QRS_waves) / nrow(matched_waves),
    duplicate_P = length(duplicate_pwaves) / nrow(matched_waves),
    duplicate_T = length(duplicate_twaves) / nrow(matched_waves)
  )
  
  return(output)
  # print(unmatched_Rpeaks)
  # print(unmatched_QRS_waves)
  # print(duplicate_pwaves)
  # print(duplicate_twaves)
  # matched_waves
  
  # Pseudocode:  
  # Need to adjust for waves which occur too close to beginning/end: 
  #   eliminate Rpeaks if <250, >4750
  #   there should be one labelled p-wave onset prior to  QRS (earliest: index ~150)
  #     ignore waves prior to this p-wave
  #   there should be one labelled t-wave onset after the QRS
  #     ignore waves after this t-wave
  
  # For "normal" cases not at the beginning/end:
  #   Check if there is p-wave before, t-wave after
  #   Note if there are duplicates:
  #     more than 1 p-wave before the prior QRS
  #     more than 1 t-wave after  the QRS
  #   Note if there are unmatched:
  #     no p-wave prior to QRS
  #     no t-wave after QRS
  
  # Toggle based on Rpeak- so it still looks for a p-wave even if QRS is missed, for example
  
  # Outputs:
  #   % missed P, QRS, T
  #   % duplicated P, QRS, T
  
}
# Custom plotting function ------------------------------------------------
plot_func <- function(y, 
                      color = 0,
                      linewidth=0.5, 
                      pointsize = 0.5, 
                      ylim = NULL, 
                      plotly = 'yes', 
                      x) {
  # Custom plotting function, Can toggle between ggplot (plotly = 'no') vs
  # plotly 
  library(ggplot2)
  library(plotly)

  # If annotations are wfdb format, change to vector format
  if (any(class(color) == "annotation_table")) {
    color <- ann_wfdb2continuous2(color, length = length(y))
  }
  
  
  y <- c(y)
  
  color[color == 1] <- 'p'
  color[color == 2] <- 'N'
  color[color == 3] <- 't'
  color[color == 4] <- 'V'
  
  x <- 1:length(y)
  
  frame <- data.frame(Time = x, Signal = y)
  plot <-
    ggplot(frame, aes(Time, Signal, color = color)) +
    geom_path(linewidth = linewidth, aes(group = 1)) + geom_point(size = pointsize) +
    scale_x_continuous(breaks = seq(0, 10, 1)) + 
    theme(legend.position = "none") + 
    theme(legend.title = element_blank()) +
    theme_bw() + 
    coord_cartesian(ylim = ylim)
  
  if (plotly == 'yes') {
    plot <- ggplotly(plot)
  }
  
  # original linewidth = 0.25
  return(plot)
}


plot_func2 <- function(y, 
                       color = 0,
                       linewidth = 0.5, 
                       pointsize = 0.5, 
                       ylim = NULL, 
                       plotly = 'yes',
                       legend = TRUE,
                       x_units = 'index') {
  #' @param x_units: can either be 'sec' for seconds (time is from 0 to 10 sec on standard ECG). OR index, for index by index
  library(ggplot2)
  library(plotly)
  
  # If annotations are wfdb format, change to vector format
  if (any(class(color) == "annotation_table")) {
    color <- ann_wfdb2continuous2(object = color, length = length(y))
  }
  
  # If signal is an array (ie dim 1 x 5000), reduce to vector 
  y <- c(y) 
  
  # Map numeric codes to labels
  color[color == 0] <- '0'
  color[color == 1] <- 'P'
  color[color == 2] <- 'QRS'
  color[color == 3] <- 'T'
  color[color == 4] <- 'V'
  
  if (x_units == 'index') {
    x <-  seq_along(y)
  } else if (x == 'sec') {
    x <- seq_along(y) / 500
  }
  
  frame <- data.frame(Time = x, Signal = y, Wave = color)
  
  # Create a grouping variable that increments whenever Wave changes
  frame$segment <- cumsum(c(TRUE, diff(as.numeric(as.factor(frame$Wave))) != 0))
  
  plot <- ggplot(frame, aes(Time, Signal, color = Wave, group = segment)) +
    geom_path(linewidth = linewidth) +
    geom_point(size = pointsize) +
    scale_color_manual(
      values = c("0" = "darkgray", "P" = "#1A759F", "QRS" = "#28B000", "T" = "#BC4B51", "V" = "red")
    ) +
    theme_bw() +
    coord_cartesian(ylim = ylim) +
    theme(
      legend.title = element_blank(),
      legend.position = if (legend) "right" else "none"
    )
  
  if (plotly == 'yes') {
    plot <- ggplotly(plot)
  }
  
  return(plot)
}

# bilstm_cnn model --------------------------------------------------------
build_bilstm_cnn <- function(bilstm_layers, 
                             number_of_derivs,
                             mask_value,
                             kernel,
                             filters,
                             activation = 'softmax',
                             num_classes = num_classes,
                             input_length = 5000) {
  library(keras)
  library(tensorflow)
  
  # Model parameters
  bilstm_layers <- bilstm_layers         # LSTM units
  num_channels <- number_of_derivs+1  # single lead input
  mask_value <- out$mask_value        # if you wish to mask zero values in the input
  kernel <- kernel                    # size of kernel for each convolutional layer
  filters <- filters                  # number of filters for each convolutional layer
  activation <- activation            # final activation for multi-class predictions
  
  
  # Define the model input
  inputs <- layer_input(shape = c(input_length, num_channels), dtype = 'float32') |>
    layer_masking(mask_value = mask_value) # Masking input values equal to 0
  
  # -- Convolutional Feature Extraction Block --
  # This block applies 1D convolutions to extract local features from the ECG slice.
  cnn_block <- inputs %>%
    # First convolution: using a larger kernel can help capture the broader wave morphology.
    layer_conv_1d(
      filters = filters[1],
      kernel_size = kernel[1],
      activation = 'relu',
      padding = "same"
    ) %>%
    # Max pooling to downsample and focus on the most prominent features.
    layer_max_pooling_1d(pool_size = 2) %>%
    # Second convolution: deeper features at a possibly finer resolution.
    layer_conv_1d(
      filters = filters[2],
      kernel_size = kernel[2],
      activation = 'relu',
      padding = "same"
    ) %>%
    layer_max_pooling_1d(pool_size = 2)
  
  # -- Recurrent Block for Temporal Modeling --
  # The output of the convolutional block is still sequential.
  # The bidirectional LSTM helps capture signals from both past and future contexts.
  lstm_block <- cnn_block %>%
    bidirectional(layer_lstm(
      units = bilstm_layers,
      return_sequences = TRUE,  # important to output a sequence for per-time-step predictions
      activation = 'tanh'
    ))
  
  # -- Upsampling Layer --
  upsampled_block <- lstm_block %>%
    layer_upsampling_1d(size = 4)  # Upsample by a factor of 2 to restore time steps
  
  # -- Output Layer --
  # A dense layer produces classification probabilities at every time step.
  outputs <- upsampled_block %>%
    layer_dense(units = num_classes, activation = activation, name = "predictions")
  
  # Define, compile, and summarize the model
  model <- keras_model(inputs = inputs, outputs = outputs, name = "CNN_LSTM_ECG_Model")
  
  model %>% compile(
    optimizer = optimizer_adam(),
    loss = "sparse_categorical_crossentropy",
    metrics = "accuracy"
  )
  
  return(model)
}



# bilstm_cnn_branched model -----------------------------------------------


build_bilstm_cnn_branched <- function(bilstm_layers,
                                      number_of_derivs,
                                      mask_value,
                                      kernels = c(7, 31, 81),
                                      branch_filters = c(32, 32, 32),
                                      proj_filters = 128,
                                      activation = 'softmax',
                                      num_classes,
                                      input_length = 5000) {
  library(keras)
  library(tensorflow)
  
  num_channels <- number_of_derivs + 1
  
  inputs <- layer_input(shape = c(input_length, num_channels), dtype = 'float32') %>%
    layer_masking(mask_value = mask_value)
  
  # Branch A: small kernel (local edges and peaks)
  branch_a <- inputs %>%
    layer_conv_1d(filters = branch_filters[1],
                  kernel_size = kernels[1],
                  padding = 'same',
                  activation = 'relu')
  
  # Branch B: medium kernel (morphology)
  branch_b <- inputs %>%
    layer_conv_1d(filters = branch_filters[2],
                  kernel_size = kernels[2],
                  padding = 'same',
                  activation = 'relu')
  
  # Branch C: large kernel (whole-wave context)
  branch_c <- inputs %>%
    layer_conv_1d(filters = branch_filters[3],
                  kernel_size = kernels[3],
                  padding = 'same',
                  activation = 'relu')
  
  # Optionally include a pooled branch for even coarser context
  branch_pool <- inputs %>%
    layer_max_pooling_1d(pool_size = 4, padding = 'same') %>%
    layer_conv_1d(filters = proj_filters / 4,
                  kernel_size = 11,
                  padding = 'same',
                  activation = 'relu') %>%
    layer_upsampling_1d(size = 4)
  
  # Concatenate multi-scale features
  multi_scale <- layer_concatenate(list(branch_a, branch_b, branch_c, branch_pool), axis = -1)
  
  # Projection conv to mix channels and control dimensionality
  projected <- multi_scale %>%
    layer_conv_1d(filters = proj_filters,
                  kernel_size = 1,
                  padding = 'same',
                  activation = 'relu') %>%
    layer_batch_normalization()
  
  # Optional light downsampling to reduce sequence length for LSTM compute
  # Keep stride = 1 for per-sample output alignment; use pooling only if you upsample later
  # Here we keep time resolution and let the BiLSTM operate on full-resolution features
  lstm_block <- projected %>%
    bidirectional(layer_lstm(units = bilstm_layers,
                             return_sequences = TRUE,
                             activation = 'tanh'))
  
  # If pooling was used earlier and time dimension was reduced, upsample back here
  # For this architecture we preserve time dimension and directly predict per time step
  outputs <- lstm_block %>%
    layer_dense(units = num_classes, activation = activation, name = 'predictions')
  
  model <- keras_model(inputs = inputs, outputs = outputs, name = 'multiscale_cnn_bilstm_ecg')
  
  model %>% compile(
    optimizer = optimizer_adam(),
    loss = 'sparse_categorical_crossentropy',
    metrics = c('accuracy')
  )
  
  return(model)
}

# Find Isoelectric -------------------------------------------------------------
isoelec_find <- function(signal,annotations) {
  # Finds the mean value of the T-P intervals within the given sample lead
  # Output: single mV value
  # take mean vs. median of T-P intervals?
  
  if (sum(annotations == 1) > 0 & sum(annotations == 3) > 0) {
    pwaves <- make_wave_table(annotations, wave_value = 1)
    twaves <- make_wave_table(annotations, wave_value = 3)
  } else if (sum(annotations == 3) > 0 & sum(annotations == 2)) {
    # If there are no p-waves, but there are t-waves, find T-QRS interval:
    twaves <- make_wave_table(annotations, wave_value = 3)
    # print(paste("No pwaves found. Using T-R interval for isoelec point"))
    pwaves <- make_wave_table(annotations, wave_value = 2)
    pwaves$wave_type <- "p"
  } else {
    # If no t-waves, find median of ECG:
    warning("No twaves found. Using median value for isoelec point")
    return(median(signal))
    break
  }
  
  combined <- rbind(pwaves,twaves)
  combined <- combined[order(combined$wave_on),]
  
  
  isoelectric_line <- array(0,length(unique(combined$sample)))
  for (i in 1:length(unique(combined$sample))) {
    # for each unique sample
    sample_table <- combined[combined$sample == i,]
    p_ind <- which(sample_table$wave_type == "p")
    
    # skip first p-wave if there's no preceeding t-wave:
    if (p_ind[1] == 1) {
      start <- 2
      
      # If there's only one P-wave, which is the first wave in the ECG:
      if (length(p_ind) < 2) {
        isoelectric_line[i] <- median(signal)
        break
      }
    } else{
      start <- 1
    }
    
    # build array of p-t interval points
    pt_ind <- c()
    for (j in start:length(p_ind)) {
      if (sample_table$wave_type[p_ind[j] - 1] == "t") {
        new_pt_ind <- sample_table$wave_off[p_ind[j] - 1] : sample_table$wave_on[p_ind[j]]
        pt_ind <- c(pt_ind,new_pt_ind)
      }
    }
    
    # In the rare case that there are p and t waves, but not in sequential order,
    # find median of signal:
    if (length(pt_ind) == 0) {
      isoelectric_line[i] <- median(signal)
    } else {
      isoelectric_line[i] <- mean(signal[pt_ind])
    }
  }
  
  
  
  return(isoelectric_line)
  
  # for each p wave, find each which is preceeded by t wave
  # of those, find T-P interval 
  
  # find which p waves are preceeded by t waves
  
  
  # make wavetable of just p and t waves
  # sort into time order
  # find indices which switch from t to p wave
  # find average/median of the 0 waves contained within those indices^
  
  
}



# Wave Table --------------------------------------------------------------
make_wave_table <- function(annotations,  wave_value) {
  # Gives table of specified on/offset times for any wave- P/QRS/T 
  # Multi-sample input handling using "sample" column
  
  # Could include time values, in addition to indices
  if (is.vector(annotations) == TRUE) {
    annotations <- matrix(annotations, nrow = 1)
  }
  
  wave_classes <- c(0,"p","N","t")
  
  wave_on <- c()
  wave_off <- c()
  sample <- c()
  wave_type <- c()
  
  for (i in 1:dim(annotations)[[1]]) {
    wave_cluster <-  which(annotations[i, ] == wave_value)
    if (length(wave_cluster) == 0) {
      wave_type = c(wave_type, NA)
      wave_on <- c(wave_on, NA)
      wave_off <- c(wave_off, NA)
      sample <- c(sample, i)
      break
    }
    
    change <- (wave_cluster[-1] - wave_cluster[1:(length(wave_cluster) - 1)])
    
    wave_on <- c(wave_on, wave_cluster[1])
    wave_on <- c(wave_on, wave_cluster[which(change != 1) + 1])
    
    wave_off <- c(wave_off, wave_cluster[which(change != 1)])
    wave_off <- c(wave_off, wave_cluster[length(wave_cluster)])
    
    sample <- c(sample, array(i, sum(change != 1) + 1))
    
    wave_type <- c(wave_type,array(wave_classes[wave_value+1],(length(wave_on) - length(wave_type))))
    
  }
  
  
  # wave_type <- array(wave_classes[wave_value+1],length(wave_on))
  
  wave_table <- data.frame(wave_type = wave_type, wave_on = wave_on, wave_off = wave_off, sample = sample)
  
  return(wave_table)
}

# Generate random ECGs ----------------------------------------------------
generate_ecgs <- function(size=10,dx_pattern='sinus|afib') {
  library(fs)
  library(EGM)
  library(stringr)
  library(XML)
  library(xml2)
  library(dplyr)
  options(wfdb_path = '/mmfs1/home/dseaney2/wfdb/bin')
  # source('../../ECG Segmentation/code/cluster_functions.R') # unused?
  
  muse_path = '/mmfs1/projects/cardio_darbar_chi/common/data/muse/muse.log'
  muse_log <- read.csv(muse_path)
  
  wfdb_path = '/mmfs1/projects/cardio_darbar_chi/common/data/wfdb/wfdb.log'
  wfdb_log <- read.csv(wfdb_path) 
  
  base_path <-  "/mmfs1/projects/cardio_darbar_chi/common/"
  
  temp_folder_name <- 'temp_ecg_generator'
  dir.create(temp_folder_name)
  uih_ecgs <- list(NA,c(size))
  
  counter <- 1
  row <- 1
  while (row <= size) {
    idx <- sample(1:nrow(wfdb_log), 1, replace = FALSE)
    
    
    dir <- fs::path(base_path,dirname(wfdb_log$PATH[idx]))
    file <- wfdb_log$FILE_NAME[idx]
    
    muse_idx <- which((muse_log$FILE_NAME == wfdb_log$FILE_NAME[idx]))
    xml <- read_xml(fs::path(base_path, dirname(muse_log$PATH[muse_idx]), muse_log$FILE_NAME[muse_idx],ext = 'xml'))
    diagnosis <- xml_text(xml_find_all(xml, ".//DiagnosisStatement"))
    
    # Check if sinus rhythm
    if (any(grepl(diagnosis,pattern=dx_pattern))) {
      test <- read_wfdb(record = file, record_dir = dir)
      sig <- test$signal$I
      # Check if ECG is 500 Hz sampling
      if (length(sig) == 5000) {
        # Need to move ecgpuwave file to same folder as header and data file
        file.copy(from = fs::path(dir,file,ext='hea'), to = fs::path(temp_folder_name))
        file.copy(from = fs::path(dir,file,ext='dat'), to = fs::path(temp_folder_name))
        
        # Copy ecgpuwave ECG
        ecgpuwave_dir <- sub("wfdb", "ecgpuwave", dir)
        file.copy(from = fs::path(ecgpuwave_dir,file,ext='ecgpuwave'), to = fs::path(temp_folder_name))
        
        # Read ECG 
        wfdb <- read_wfdb(record = file, record_dir = temp_folder_name, annotator = 'ecgpuwave')
        uih_ecgs[row] <- list(wfdb)
        
        row <- row+1
        file.remove(list.files(temp_folder_name, full.names=TRUE))
      }
    }
    
    counter <- counter+1
    
    if (!row %% 10) {
      print(paste0('Added: ',row,', Iterated thru: ', counter))
    }
  }
  unlink("temp")
  return(uih_ecgs)
}



# Predict -----------------------------------------------------------------
#' @param input: for input_class = 'wfdb': variable is of class 'list', where one index is list of **wfdb** format
#'               for input_class = 'array' (or unlabeled): variable is an array, where columns are for the sample number, and rows are for each time step
#' 
#' @param lead_number: integer, where they follow the order of 'leads', see below. If set to lead_number = 'all', all 12 leads will be predicted
#' 
#' @param model_number: model number from the model_log.RData file, where the number represents a row in the file.
#'  best models for each lead: 
#'    c(861, 856, 851,  846,  836,  841,  826, 821, 866, 871, 876, 881)
#'    c('I','II','III','AVR','AVL','AVF','V1','V2','V3','V4','V5','V6')
#'    If model_number is set to 'best', the function will use the requested lead to pick the best model for that lead
#'    
#' @param do_predictionInteger_threshold: method used to convert raw ML output probabilities to integers representing wave values
#'    If true, a predictionInteger_threshold is required (default of 0.5)
#'    If false, an absolute method will be used- take the greatest probability wave (or no wave) from each time step. 
#'        Roughly equivalent to a threshold method with a threshold value of 0.5
#' 
#' @param do_fill_wave_gaps: if a single wave (ie P wave) has a gap of less than wave_gap_threshold, fill the gap with P wave markers
#' @param wave_gap_threhold: see above
#' 
#' @param do_remove_short_waves: if there is a short wave (ie P wave) that is less than short_wave_threshold (ie 10 indices), remove the wave
#' @param short_wave_threshold: see above
#' 
#' @return ML predictions for all requested leads. This can be in wfdb format (recommended), or matrix format ([number_of_samples x time_steps]). If wfdb format, you may simply loop the function over all 12 leads, and it will add automatically

predict_ecgs <- function(input,
                         lead_number='all',
                         model_number='best',
                         model_log_path='../models/model_log.RData',
                         model_folder_path='../models/',
                         input_class='wfdb',
                         filter=TRUE,
                         do_predictionInteger_threshold = TRUE,
                         predictionInteger_threshold=0.5,
                         do_fill_wave_gaps=TRUE,
                         wave_gap_threshold=20,
                         do_remove_short_waves=TRUE,
                         short_wave_threshold=10) {
  library(keras)
  
  load(model_log_path)

  lead_name_list <- c('I','II','III','AVR','AVL','AVF','V1','V2','V3','V4','V5','V6')
  
  if (any(lead_number == 'all')) {
    lead_number <- 1:12
  }
  
    
  if (input_class == 'wfdb') {
    # If single ECG input, transform to list format to allow for compatability
    if (any(class(input) == 'egm')) {
      single_ecg_input <- TRUE
      input <- list(input)
    }
    
    # retain signal values for output and remove pre-existing annotations, if in wfdb format
    output <- lapply(input, function(x) { # remove any pre-existing annotation
      x$annotation <- NULL
      x
    })
  } else {
    output <- array(NA,c(length(input),length(input[[1]]$signal[[1]]),length(lead_number)))
  }
    
  counter <- 1
  for (lead in lead_number) {
    lead_name <- lead_name_list[lead]
    
    # If model_number is defined as 'best', pick the best performing model for the specified lead:
    if (model_number == 'best') {
      best_models <-  c(861, 856, 851, 846, 836, 841, 826, 821, 866, 871, 876, 881)
      model_number <- best_models[lead]
    }
    
    # Change ECG input from list to array
    
    input_signal <- do.call(rbind, lapply(1:length(input), function(idx)
      input[[idx]]$signal[[lead_name]]))
    
    # Filter
    if (filter) {
      for (i in 1:nrow(input_signal)) {
        input_signal[i,] <- ecg_filter(input_signal[i, ])
      }
    }
    
    # Normalize from 0 to 100
    if (model_log$normalize[model_number]) {
      for (i in 1:nrow(input_signal)) {
        input_signal[i, ] <- (input_signal[i, ] - min(input_signal[i, ])) / (max(input_signal[i, ]) - min(input_signal[i, ])) * 100
      }
    }
    
    # Add derivatives if needed
    number_of_derivs <- model_log$derivs[model_number]
    if (number_of_derivs > 0) {
      input_old <- input_signal
      input_signal <- array(NA, c(dim(input_old), number_of_derivs + 1))
      for (i in 1:nrow(input_signal)) {
        input_signal[i, , ] <- add_derivs(signal = input_old[i, ], number_of_derivs = number_of_derivs)
      }
    }
    
    # Predict
    model <- load_model_tf(paste0(model_folder_path, model_log$name[model_number], '.h5')) # keras
    
    predictions <- predict_ecg_raw(
      model = model,
      input_signal = input_signal,
      window_size = 5000,
      overlap_length = 500,
      ignore = 250
    )
    
    # predictions <- model %>% predict(input_signal) # keras
    
    # Convert probabilities to integer values - either use threshold, or greatest probability
    if (do_predictionInteger_threshold) {
      # Treshold method:
      predictions_integer <- predictions2integer_threshold(predictions, 
                                                           threshold = predictionInteger_threshold)
    } else {
      # Absolute (greatest probability) method:
      predictions_integer <- array(0, c(nrow(predictions), ncol(predictions)))
      for (i in 1:nrow(predictions)) {
        predictions_integer[i, ] <- max.col(predictions[i, , ])
      }
      # convert from dimension value 1,2,3,4 to 0,1,2,3
      predictions_integer <- predictions_integer - 1
    }
    
    
    # Fill gaps as needed
    if (do_fill_wave_gaps) {
      for (i in 1:nrow(predictions_integer)) {
        predictions_integer[i, ] <- fill_wave_gaps(predictions_integer[i, ], wave_gap_threshold)
      }
    }
    
    # Remove short waves as needed
    if (do_remove_short_waves) {
      for (i in 1:nrow(predictions_integer)) {
        predictions_integer[i, ] <- remove_short_waves(predictions_integer[i, ], max_wave_length = short_wave_threshold)
      }
    }
    
    # Remove waves less than threshold
    
    # Transform to wfdb format, if input is wfdb format
    if (input_class == 'wfdb') {
      for (i in seq_len(nrow(predictions_integer))) {
        
        # Ensure annotation slot exists
        if (is.null(output[[i]]$annotation)) {
          output[[i]]$annotation <- list()
        }
        # Create dataframe for this lead
        output[[i]]$annotation[[lead_name]] <- ann_continuous2wfdb(predictions_integer[i, ])
      }
    } else {
      output[,,counter] <- predictions_integer
    }
    
    print(paste('Finished lead',lead_name_list[lead]))
    counter <- counter+1
  }
  
  # If input was a single ECG, convert output format from a list of WFDBs to a single wfdb
  if (exists('single_ecg_input')) {
    if (single_ecg_input) {
      output <- output[[1]]
    }
  }
  
  return(output)
}

# predict_ecgs_raw --------------------------------------------------------
# Function to predict ECGs, specifically built to handle signals > 5000 (sliding window)

#' @description function which takes in a matrix of ECG signal, and outputs raw prediction. Able to handle signals > 5000 in length via windowing techniques
#' 
#' @param model tensor flow model from keras
#' @param input_signal ECG signal in array format
#' @param window_size length of each individual prediction window 
#' @param overlap_length length that is shared between each prediction window
#' @param ignore the length at both beginning and end of the prediction window that is left unlabeled (ie due to boundary issues). Note: this is ignored for the beginning of the first window, and end of last
predict_ecg_raw <- function(model, input_signal, 
                            window_size = 5000, 
                            overlap_length = 500, 
                            ignore = 250) {
  
  N <- dim(input_signal)[[2]]
  n_samples <- dim(input_signal)[[1]]
  n_classes <- 4
  step_size <- window_size - overlap_length
  
  # Initial starts
  starts <- seq(1, N - window_size + 1, by = step_size)
  
  # Force last window to align with end
  last_start <- N - window_size + 1
  if (tail(starts, 1) != last_start) {
    starts <- c(starts, last_start)
  }
  
  # Accumulator and counter
  acc <- array(0, c(n_samples,N,n_classes))
  counts <- rep(0, N)
  
  for (i in seq_along(starts)) {
    start <- starts[i]
    end <- start + window_size - 1
    
    # Define chunk for prediction
    if (length(dim(input_signal)) == 2) {
      # 2D case: [samples, time]
      chunk <- input_signal[, start:end, drop = FALSE]
    } else if (length(dim(input_signal)) == 3) {
      # 3D case: [samples, time, derivs]
      chunk <- input_signal[, start:end, , drop = FALSE]
    }
    
    # Run prediction
    pred <- predict(object = model, x = chunk)  # output shape: n_samples x 5000 x 4
    
    # Determine indices to keep
    keep_start <- if (i == 1) 1 else (ignore + 1)
    keep_end   <- if (i == length(starts)) window_size else (window_size - ignore)
    
    keep_idx <- keep_start:keep_end
    global_idx <- start + keep_idx - 1
    
    # Add to accumulator
    acc[,global_idx, ] <- acc[,global_idx, ] + pred[,keep_idx, ]
    counts[global_idx] <- counts[global_idx] + 1
  }
  
  # Average overlaps
  final_pred <- acc / counts
  
  return(final_pred)
}

# prediction_integer threshold --------------------------------------------
predictions2integer_threshold <- function(predictions, threshold = 0.5) {
  # Dimensions
  n_samples   <- dim(predictions)[1]
  n_timesteps <- dim(predictions)[2]
  n_classes   <- dim(predictions)[3]  # should be 4
  
  # Initialize output with 0s (or 1s if you want "no wave" as default)
  predictions_integer <- matrix(0, nrow = n_samples, ncol = n_timesteps)
  
  # Loop over samples
  for (s in 1:n_samples) {
    # Extract probabilities for this sample: timesteps × classes
    sample_probs <- predictions[s,,]   # dims: n_timesteps × n_classes
    
    # Focus only on classes 2:4 (P, QRS, T)
    wave_probs <- sample_probs[, 2:4, drop = FALSE]
    
    # For each timestep, check which classes exceed threshold
    above_thresh <- wave_probs > threshold
    
    # Count how many classes exceed threshold per timestep
    n_above <- rowSums(above_thresh)
    
    # Case 1: exactly one class above threshold → assign that class 
    one_class_idx <- which(n_above == 1)
    if (length(one_class_idx) > 0) {
      class_ids <- apply(above_thresh[one_class_idx, , drop = FALSE], 1, which.max)
      predictions_integer[s, one_class_idx] <- class_ids
    }
    
    # Case 2: multiple classes above threshold → assign 5 (mixed wave)
    multi_class_idx <- which(n_above > 1)
    if (length(multi_class_idx) > 0) {
      predictions_integer[s, multi_class_idx] <- 4
    }
    
    # Case 3: no class above threshold → leave as 0 
  }
  
  return(predictions_integer)
}



# check annotations -------------------------------------------------------
#' @description Checks wave progression of ECG annotations. Uses a Pan-Tompkins method from EGM to check for QRS, then adjacent windows for P and T waves. Checks for missed or duplicate waves
#'  NOTE: input can either be a wfdb format ECG and specified lead (ie 'I') OR individual signal and annotation files for a single lead (as well as fs rate)
#'    eg: validate_ecg_waves(wfdb = wfdb, lead = 'I') OR validate_ecg_waves(signal = wfdb$signal$I, ann_wfdb = wfdb$annotation$I)
#' @param wfdb wfdb object with sig, ann and header. If defined, 'lead' should also be defined
#' @param lead lead name. ie c('I'). To be used with 'wfdb' style input
#' 
#' @param signal signal of the ECG. Simple vector. To be used with 'ann_wfdb', as well as fs
#' @param ann_wfdb wfdb format annotation of the ECG. To be used with 'signal', as well as fs
#' @param fs sampling rate of ECG. Only needed when using signal / ann_wfdb inputs. Not needed for wfdb / lead input
#' 
#' @param terminal_boundary the flanking indices you would like removed from the annotation check. 
#' 
#' @param pr_window defined window to search for the P-waves, distance from the Rpeak. ie c(80,  350). NOTE: by default, these are dynamically calculated based on RR interval
#' @param rt_window defined window to search for the T-waves, distance from the Rpeak. ie c(120, 500). NOTE: by default, these are dynamically calculated based on RR interval
#' 
#' @return Returns a dataframe with columns for QRS, P and T wave flags and flag indices, and check if in range of windows. The flag columns specify missing or duplicate (or NA) 
#'    flag: NA (normal wave), duplicate or missing
#'    p/t_flag_idx: index of the duplicate waves, OR the index of where the P/T wave may be located
#'    p/t_in_range: if the wave is contained within the defined window


validate_ecg_waves <- function(wfdb = NULL, 
                               lead = NULL, 
                               signal = NULL, 
                               ann_wfdb = NULL, 
                               fs = 500,
                               terminal_boundary = 250,
                               pr_window = NULL,   # ms before R
                               rt_window = NULL) { # ms after R
  
  if (exists('wfdb')) {
    signal <- wfdb$signal[[lead]]
    ann_wfdb <- wfdb$annotation[[lead]]
    fs <-  attributes(wfdb$header)$record_line$frequency
  }
  
  sig_length <- length(signal) # update this in future
  
  # Collapse annotations
  ann_compact <- ann_wfdb2compact(ann_wfdb)
  
  # Detect R peaks using detect_QRS, remove those at termini
  rpeaks <- EGM::detect_QRS(signal, frequency = fs)
  rpeaks <- rpeaks[rpeaks < (length(signal) - terminal_boundary)]
  rpeaks <- rpeaks[rpeaks > terminal_boundary]
  
  if (!length(rpeaks)) {
    empty <- data.frame(rpeak = numeric(0),
                        p_flag = character(0),
                        p_flag_idx = character(0),
                        t_flag = character(0),
                        t_flag_idx = character(0),
                        qrs_flag = character(0),
                        qrs_flag_idx = character(0),
                        stringsAsFactors = FALSE)
    return(empty)
  }
  
  calc_pr_window_ms <- function(rr_ms) {
    hr <- 60000 / rr_ms
    hr <- min(max(hr, 40), 160)
    pr_center <- 200 - 0.8 * (hr - 60)
    pr_center <- max(min(pr_center, 220), 110)
    lower <- max(pr_center - ifelse(hr < 70, 70, 55), 60)
    upper <- min(pr_center + ifelse(hr < 70, 80, 65), 260)
    if (upper <= lower) upper <- lower + 30
    c(lower, upper)
  }
  
  # Similarly, calculate RR interval
  # Probably take an average of Frederica and Bazzett for now
  calc_rt_window_ms <- function(rr_ms) {
    rr_s <- rr_ms / 1000
    qt_baz <- 440 * sqrt(rr_s)
    qt_frd <- 430 * rr_s^(1/3)
    qt_center <- (qt_baz + qt_frd) / 2
    qt_center <- max(min(qt_center, 500), 300)
    qt_low <- max(qt_center - 90, 240)
    qt_high <- min(qt_center + 90, 520)
    rt_low <- max(qt_low - 60, 120)
    rt_high <- min(qt_high - 20, 520)
    if (rt_high <= rt_low) rt_high <- rt_low + 40
    c(rt_low, rt_high)
  }
  
  rr_intervals_ms <- diff(rpeaks) * 1000 / fs
  rr_valid <- rr_intervals_ms[rr_intervals_ms >= 300 & rr_intervals_ms <= 2000]
  rr_reference <- if (length(rr_valid)) stats::median(rr_valid) else 1000
  
  # Make sure RR intervals are valid and reasonable
  sanitize_rr <- function(rr_vec) {
    rr_vec[!is.finite(rr_vec) | rr_vec < 300 | rr_vec > 2000] <- rr_reference
    rr_vec
  }
  
  rr_prev <- sanitize_rr(c(rr_reference, rr_intervals_ms))
  rr_next <- sanitize_rr(c(rr_intervals_ms, rr_reference))
  rr_local <- sanitize_rr((rr_prev + rr_next) / 2)
  
  # --- QRS CHECK SECTION ---
  # Identify QRS from annotations, remove those at termini
  qrs_ann <- ann_compact |> dplyr::filter(type == "N")
  qrs_ann <- qrs_ann |> dplyr::filter(peak < (length(signal) - terminal_boundary) & peak > terminal_boundary)
  
  # For each rpeak, check if it falls inside any QRS annotation
  rpeak_matches <- sapply(rpeaks, function(r) {
    any(r >= qrs_ann$onset & r <= qrs_ann$offset)
  })
  
  # For each QRS annotation, check if it has a matching rpeak
  qrs_matches <- sapply(1:nrow(qrs_ann), function(i) {
    any(rpeaks >= qrs_ann$onset[i] & rpeaks <= qrs_ann$offset[i])
  })
  
  # Missed = rpeak without QRS
  missed_rpeaks <- rpeaks[!rpeak_matches]
  
  # Duplicate = QRS without rpeak
  duplicate_qrs <- qrs_ann[!qrs_matches, ]
  
  # Collapse duplicate QRS indices into one string
  duplicate_idx <- if (nrow(duplicate_qrs) > 0) {
    paste(round(rowMeans(duplicate_qrs[, c("onset", "offset")])), collapse = ";")
  } else NA
  
  qrs_flags <- data.frame(
    rpeak = rpeaks,
    qrs_flag = ifelse(rpeaks %in% missed_rpeaks, "missed", NA),
    qrs_flag_idx = ifelse(rpeaks %in% missed_rpeaks, as.character(rpeaks), NA),
    stringsAsFactors = FALSE
  )
  
  # Add one row for all duplicate QRS (collapsed)
  if (!is.na(duplicate_idx)) {
    dup_row <- data.frame(
      rpeak = NA,
      qrs_flag = "duplicate",
      qrs_flag_idx = duplicate_idx,
      stringsAsFactors = FALSE
    )
    qrs_flags <- dplyr::bind_rows(qrs_flags, dup_row)
  }
  # --- END QRS CHECK SECTION ---
  
  # Initialize results as a list of rows
  results <- lapply(seq_along(rpeaks), function(idx) {
    r <- rpeaks[idx]
    prior_r <- if (idx == 1) {terminal_boundary} else {rpeaks[idx-1]} # define prior R peak. If undefined, set as terminal boundary
    next_r <- if (idx == length(rpeaks)) {sig_length - terminal_boundary} else {rpeaks[idx+1]} # define next R peak. If undefined, set as terminal boundary
    
    pr_samples <- if (is.null(pr_window)) {
      win <- round(calc_pr_window_ms(rr_prev[idx]) * fs / 1000)
      if (win[2] <= win[1]) win[2] <- win[1] + 1
      win
    } else {
      win <- round(pr_window * fs / 1000)
      if (length(win) != 2) stop("pr_window must be length 2 if provided.")
      if (win[2] <= win[1]) win[2] <- win[1] + 1
      win
    }
    
    rt_samples <- if (is.null(rt_window)) {
      win <- round(calc_rt_window_ms(rr_local[idx]) * fs / 1000)
      if (win[2] <= win[1]) win[2] <- win[1] + 1
      win
    } else {
      win <- round(rt_window * fs / 1000)
      if (length(win) != 2) stop("rt_window must be length 2 if provided.")
      if (win[2] <= win[1]) win[2] <- win[1] + 1
      win
    }
    
    # Define search windows
    pr_range <- (r - pr_samples[2]):(r - pr_samples[1])
    pr_range <- pr_range[pr_range >= 1]
    rt_range <- (r + rt_samples[1]):(r + rt_samples[2])
    rt_range <- rt_range[rt_range >= 1 & rt_range <= length(signal)]
    
    # Restrict to annotations within RR intervals
    
    p_candidates <- ann_compact |> dplyr::filter(type == "p", onset %in% prior_r:r)
    t_candidates <- ann_compact |> dplyr::filter(type == "t", offset %in% r:next_r)
    
    # Flags
    p_flag <- NA
    p_flag_idx <- NA
    p_in_range <- NA
    if (nrow(p_candidates) == 0) {
      p_flag <- "missing"
      p_flag_idx <- if (length(pr_range)) paste(round(mean(pr_range))) else NA
    } else if (nrow(p_candidates) > 1) {
      p_flag <- "duplicate"
      p_flag_idx <- paste(p_candidates$peak, collapse = ";")
    }
    if(any(p_candidates$peak %in% pr_range)) {p_in_range <- TRUE} else {p_in_range <- FALSE}
    
    
    t_flag <- NA
    t_flag_idx <- NA
    if (nrow(t_candidates) == 0) {
      t_flag <- "missing"
      t_flag_idx <- if (length(rt_range)) paste(round(mean(rt_range))) else NA
    } else if (nrow(t_candidates) > 1) {
      t_flag <- "duplicate"
      t_flag_idx <- paste(t_candidates$peak, collapse = ";")
    }
    if(any(t_candidates$peak %in% rt_range)) {t_in_range <- TRUE} else {t_in_range <- FALSE}
    
    return(data.frame(rpeak = r,
                      p_flag = p_flag,
                      p_flag_idx = p_flag_idx,
                      p_in_range = p_in_range,
                      t_flag = t_flag,
                      t_flag_idx = t_flag_idx,
                      t_in_range = t_in_range,
                      stringsAsFactors = FALSE))
  })
  
  results_df <- dplyr::bind_rows(results)
  
  # Merge QRS flags into results
  # final <- results_df |>
  #   left_join(qrs_flags, by = "rpeak")
  
  final <- 
    results_df |>
    dplyr::left_join(qrs_flags |> 
                       dplyr::filter(!is.na(rpeak)), by = "rpeak")
  
  # Add the duplicate summary row back in
  dup_row <- qrs_flags |> dplyr::filter(is.na(rpeak))
  final <- dplyr::bind_rows(final, dup_row)
  
  
  return(final)
}
# Remove waves less than ___ ms 
# Check for QRS --> check for P, T waves. Use detect_QRS
# signal confidence? 
# ^ similar to confidence: compare all 12 leads, look for outliers and flag. 
# score for each wave, each time point. if < 0.25 for a wave, remove.
# Remove waves less than ___ ms 
# Check for QRS --> check for P, T waves. Use detect_QRS
# signal confidence? 
# ^ similar to confidence: compare all 12 leads, look for outliers and flag. 
# score for each wave, each time point. if < 0.25 for a wave, remove.