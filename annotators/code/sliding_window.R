# interpolate -------------------------------------------------------------
# Function to resample ECG signal from 360 Hz to 500 Hz using cubic spline interpolation
resample_ecg <- function(signal, old_fs, new_fs) {
  # Original time points
  t_old <- seq(0, (length(signal) - 1) / old_fs, by = 1 / old_fs)
  
  # New time points
  t_new <- seq(0, max(t_old), by = 1 / new_fs)
  
  # Perform cubic spline interpolation
  signal_new <- spline(x = t_old, y = signal, xout = t_new)$y
  
  return(signal_new)
}


# sliding window ----------------------------------------------------------
# Test function
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



# alpha sandbox -----------------------------------------------------------------

# Sys.info()['user']
options(wfdb_path = 'wsl /usr/local/bin')
source('annotator_prep_functions.R')
# Example:
# ecg_resampled <- resample_ecg(signal, old_fs = 360, new_fs = 500)
# length(ecg_resampled)  # should be ~10000


# Example usage:
# ecg_resampled <- resample_ecg(signal, old_fs = 360, new_fs = 500)

test = EGM::read_wfdb(record = 100, 
                      record_dir = 'C:/Users/darre/OneDrive/Documents/UICOM Research/ECG Segmentation/raw_datasets/mit-bih', 
                      annotator = 'atr' )

sig <- test$signal$V5[501:7701]

test2 <- test
test2$signal <-  test$signal[501:7701,]
V5 <- resample_ecg(signal = test2$signal$V5,old_fs = 360, new_fs = 500)
test2$signal <- data.frame(sample = 1:length(V5),
                   V5 = V5)
output <- predict_ecgs(input = test2, lead_number = 11)

# STEPS (within each lead loop)
# 1. interpolate as needed (use wfdb attributes, or manual specification?)
# 2. check length, if > 5000: 
#       divide into 5000 chunks, overlapping by specified amount (min 500)
#       remove terminal 250
#       if some windows overlap, average the raw probabilities
#       the last window should overlap by a factor to make it have length 5000 and have the last index match the global last index
# 3. 

input_signal <- array(data = sig2, dim = c(1,length(sig2)))
lead <- 11

# Load model
best_models <-  c(861, 856, 851, 846, 836, 841, 826, 821, 866, 871, 876, 881)
model_number <- best_models[lead]
model <- keras::load_model_tf(paste0(model_folder_path, model_log$name[model_number], '.h5'))

# Filter
for (i in 1:nrow(input_signal)) {
  input_signal[i,] <- ecg_filter(input_signal[i, ])
}

# Remove distorted signal due to filtering issues
input_signal <- input_signal[,1001:9550, drop = FALSE]

# Normalize signal
for (i in 1:nrow(input_signal)) {
  input_signal[i, ] <- (input_signal[i, ] - min(input_signal[i, ])) / (max(input_signal[i, ]) - min(input_signal[i, ])) * 100
}

# Add derivatives as needed
number_of_derivs <- model_log$derivs[model_number]
if (number_of_derivs > 0) {
  input_old <- input_signal
  input_signal <- array(NA, c(dim(input_old), number_of_derivs + 1))
  for (i in 1:nrow(input_signal)) {
    input_signal[i, , ] <- add_derivs(signal = input_old[i, ], number_of_derivs = number_of_derivs)
  }
}

# Output raw prediction probabilities
predictions <- predict_ecg_raw(model = model, input_signal = input_signal)

# Convert raw prediction probabilities to integer values
predictionInteger_threshold <- 0.5
predictions_integer <- predictions2integer_threshold(predictions, 
                                                     threshold = predictionInteger_threshold)

# Plot
plot_func2(input_signal[,,1],c(predictions_integer))


# beta testing/sandbox -----------------------------------------------------------------
# Notes:
# 1. ensure predict_ecg() function handles a single WFDB sample correctly (status: completed?)
# 2. need to improve filtering function on longer samples (distortion at the beginning and ends --> throws off annotators)


# Sys.info()['user']
options(wfdb_path = 'wsl /usr/local/bin')
source('annotator_prep_functions.R')

# Load ECG
test = EGM::read_wfdb(record = 100, 
                      record_dir = 'C:/Users/darre/OneDrive/Documents/UICOM Research/ECG Segmentation/raw_datasets/mit-bih', 
                      annotator = 'atr' )

# Take interval
test$signal <-  test$signal[1001:9001,]
V5 <- resample_ecg(signal = test$signal$V5,old_fs = 360, new_fs = 500)
test$signal <- data.frame(sample = 1:length(V5),
                           V5 = V5)
output <- predict_ecgs(input = test, lead_number = 11)
plot_func2(ecg_filter(output$signal$V5),output$annotation$V5)
