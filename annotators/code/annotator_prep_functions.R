# Build training and testing  sets---------------------------------------
#' @descrption Input: custom list of LUDB annotations and signal. Splits into training/testing. 
#' Takes random continuous segment between 2 and 8 seconds, zeros the rest. Randomly shifts the segment across the vector
#' Then normalizes the signal for each sample.
#' 
#' @param dilate_range For training set: 2 element vector. Defines the range the samples will be randomly 
#' 
#' @param max_noise For training set: maximum standard deviation of the added Gaussian noise. ECGs are scaled from 0 to 1. 
#' 
#' @param rounds number of times the training dataset will be duplicated. Note: a value of 2 will double the length of the training set. A value of 1 will not change the length
prep_ludb <- function(lead,
                      annotator_style = 2,
                      split = 0.7,
                      dilate_range = c(0.03,0.05),
                      max_noise = 0.05,
                      rounds = 1) {
  # Load LUDB set
  load('../ludb_set.RData')
  
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
  
  training_signal      <- array(NA, c(length(training_set)*rounds, length))
  training_annotations <- array(NA, c(length(training_set)*rounds, length))
  testing_signal       <- array(NA, c(length(testing_set), length))
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
      
      # DO NOT Filter signal:
      # signal <- ecg_filter(signal = signal)
      
      # Normalize signal
      signal <- (signal - min(signal)) / (max(signal) - min(signal))
      
      # Set indices outside of interval to 0
      # signal[-random_interval] <- 0
      # ann[-random_interval] <- 0
      
      # Random interval signal/ann vectors
      signal_iso <- signal[random_interval]
      ann_iso <- ann[random_interval]
      
      
      # Compress signal:
      # Generate a random compression/dilation factor
      
      if (round != 1) {
        compression_factor <- 1 + sample(c(-1, 1), size = 1) * runif(1, min = min, max = max) # Random factor between 1 +/- (3 to 5)%
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
        noise <- rnorm(length(compressed_signal),
                       mean = 0,
                       sd = sigma) # create random noise
        compressed_signal <- compressed_signal + noise # Add the noise to the original ECG signal
      } 
      
      
      # Now, randomly place the random signal/ann interval:
      # Determine a random starting point for the shifted interval
      vector_length <- length(signal)
      shifted_start <- sample(1:(vector_length - new_length), 1)
      
      # Generate the new shifted indices
      shifted_indices <- shifted_start:(shifted_start + new_length - 1)
      
      # Create new vectors initialized with zeros
      shifted_signal <- rep(0, vector_length)
      shifted_ann <- rep(0, vector_length)
      
      # Shift the non-zero values to the new interval
      shifted_signal[shifted_indices] <- compressed_signal
      shifted_ann[shifted_indices] <- compressed_ann
      
      # Add sample to matrix
      index <- sample + (round-1)*length(training_set)
      training_signal[index, ] <- shifted_signal # signal
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
    (signal - min(signal)) / (max(signal) - min(signal))
    
    # Set values of first and last 1000 indices to 0 (could also do 1st and last annotated indices)
    signal[c(1:1000, 4001:5000)] <- 0
    ann[c(1:1000, 4001:5000)] <- 0
    
    # Add sample to matrix
    testing_signal[sample, ] <- signal
    testing_annotations[sample, ] <- ann
  }
  
  # Prepare output
  output <- list(
    training_signal = training_signal,
    training_annotations = training_annotations,
    testing_signal = testing_signal,
    testing_annotations = testing_annotations,
    training_samples = training_samples,
    testing_samples = testing_samples
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


# Further analysis --------------------------------------------------------

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

check_ann_prog <- function(predicted, ludb) {
  # Check annotation progression
  #   ie- P > QRS > T
  
  # NOTE: this is for use on non-LUDB ECGs. ie those without a ground truth annotation
  #   
  
  # Input: single ECG annotation - vector or wfdb format
  
  # Pseudocode:
  #   filter ECGs
  #   could use EGM::qrs_detect() for ~ground truth. Could use ecgpuwave 
  
  # If annotation is in matrix format (not WFDB), convert to wfdb
  annotations <- ann_continuous2wfdb(annotations)
  
  wave_progression <- annotations$type[!annotations$type %in% c('(',')')]
  
  p <- which(wave_progression == 'p')
  qrs <- which(wave_progression == 'N')
  t <- which(wave_progression == 't')
  
  (p+1) %in% qrs
  (qrs+1) %in% t
  (t+1) %in% p
  
}

# Could only focus on annotations after the first R peak, before last Rpeak (use egm::find_Rpeaks)

ann_continuous2wfdb <- function(annotations, Fs = 500) {
  # Not up to date, largely defunct for now. Used to convert index to index
  # annotation (value for each time point) to wave onset/peak/offset format. 
  # WFDB format
  
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
  
  t <- (0 : length(annotations) - 1 ) / Fs
  
  
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
  
  annotation_table <- annotation_table[order(annotation_table$sample),]
  row.names(annotation_table) <- NULL
  
  return(annotation_table)
}


# Custom plotting function ------------------------------------------------
plot_func <- function(y, 
                      color = 0,
                      linewidth=0.5, 
                      pointsize = 1.5, 
                      ylim = NULL, 
                      plotly = 'yes', 
                      x) {
  # Custom plotting function, Can toggle between ggplot (plotly = 'no') vs
  # plotly 
  library(ggplot2)
  library(plotly)
  color <- c(color)
  y <- c(y)
  
  color[color == 1] <- 'p'
  color[color == 2] <- 'N'
  color[color == 3] <- 't'
  
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

