# Filter Samples ----------------------------------------------------------
filter_samples <- function(signal, frequency = 500, low = 0.5, high = 40) {
  library(dplR)
  
  if (is.vector(signal) == TRUE) {
    filtered <- pass.filt(signal, c(f/high, f/low), "pass", method = c("Butterworth"), n = 4, Rp = 0.1)
  } else {
    filtered <- array(0,dim(signal))
    for (i in 1:dim(signal[[1]])) {
      filtered[i,] <- pass.filt(signal[i,], c(f/high, f/low), "pass", method = c("Butterworth"), n = 4, Rp = 0.1)
    }
  }
}

# Build Spectrogram -------------------------------------------------------
build_spectrogram <- function(input, Fs = 500, window_size = 128, high_freq_cutoff = 45) {
  # Input is a *single* vector. Currently does not handle
  
  # Currently only works with 500 Hz
  
  library(signal)
  
  overlap <- window_size - 1
  padding <- window_size / 2
  
  # output <- array(0, )
  
  if (is.vector(input) == TRUE) {
    # Single vector handling:
    
    # Create spectrogram:
    complex <- specgram(
      c(array(0, padding), input, array(0, padding)),
      n = window_size,
      Fs = Fs,
      overlap = overlap
    )
    
    # Eliminate high frequency noise:
    index = length(complex[complex$f < high_freq_cutoff])
    
    complex$S <- complex$S[1:index, ]
    complex$f <- complex$f[complex$f < high_freq_cutoff]
    
    # Normalize values:
    normalized  <- (complex$S - mean(complex$S)) / sd(complex$S) #separately average real and imaginary?
    
    # Separate complex number into real and imaginary parts, then combine:
    real <- Re(normalized)
    imaginary <- Im(normalized)
    
    # Align dimmensions for ML input:
    dim <- dim(t(rbind(real, imaginary)))
    
    output <- array(0, dim = c(1,dim))
    output[1,,] <- t(rbind(real, imaginary))
  }
  
  else {
    # Array handling:
    
    # Determine size of output array:
    complex <- specgram(
      c(array(0, padding), input[1, ], array(0, padding)),
      n = window_size,
      Fs = Fs,
      overlap = overlap
    )
    
    index = length(complex[complex$f < high_freq_cutoff])
    output <- array(0, c(dim(input), index * 2))
    
    for (i in 1:dim(input)[[1]]) {
      
      complex <- specgram(
        c(array(0, padding), input[i, ], array(0, padding)),
        n = window_size,
        Fs = Fs,
        overlap = overlap
      )
      
      # Eliminate high frequency noise:
      index = length(complex[complex$f < high_freq_cutoff])
      
      complex$S <- complex$S[1:index, ]
      complex$f <- complex$f[complex$f < high_freq_cutoff]
      
      # Normalize values:
      normalized  <- (complex$S - mean(complex$S)) / sd(complex$S) #separately average real and imaginary?
      
      # Separate complex number into real and imaginary parts, then combine:
      real <- Re(normalized)
      imaginary <- Im(normalized)
      
      output[i, , ] <- t(rbind(real, imaginary))
    }
  }
  
  return(output)
}

# Predict Samples ---------------------------------------------------------
predict_samples <- function(signal, model_name, model_path) {
  # bug: must make compatible when submitting a single sample. ie array of length x must be converted to matrix of 
  # size: 1 by x
  # OR single spectrogram of size x by 24 to size 1 by x by 24
  
  # if (class(signal) == "integer" | class(signal == "numeric"))
  
  library(keras)
  
  # if (NROW(sample_data) == 1 & dim(sample_data[[1]] != 1)) {
  #   sample_data <- array(sample_data,c(1,))
  # }
  
  # rm("model")
  
  model <- load_model_tf(paste0(model_path, model_name))
  
  predictions_raw <- model %>% predict(signal)
  
  # Convert probabilities to integer values:
  dim <- dim(predictions_raw)
  predictions_integer <- array(0, dim[1:2])
  
  for (i in 1:dim[1]) {
    predictions_integer[i, ] <- max.col(predictions_raw[i, , ])
  }
  
  annotations <- predictions_integer - 1
  
  return(annotations)
}

# R Peak Isolation --------------------------------------------------------
R_peak_isolation <- function(signal, annotations, R_value = 2) {
  # Determine R peak values
  
  #input signal: raw signal, not spectrogram
  
  # Use 10-sec duration model ( predict_samples() ) to determine annotations
  
  if (is.vector(signal) == TRUE) {
    signal <- matrix(signal, nrow = 1)
  }
  
  Rpeaks_full <- array(0,dim(signal)[[1]])
  
  for (sample in 1:dim(signal)[[1]]) {
    # Find each continuous QRS segment:
    # QRS_predict: which time points are QRS
    # QRS_cluster_index: indices of QRS_predict which jump to the next QRS interval
    QRS_predict <-  which(annotations[sample, ] == R_value)
    change <- (QRS_predict[-1] - QRS_predict[1:(length(QRS_predict) - 1)])
    
    QRS_cluster_index <- c(1)
    QRS_cluster_index <- c(QRS_cluster_index, (which(change != 1) + 1))
    QRS_cluster_index <- c(QRS_cluster_index, length(QRS_predict)) # wrong: need to add in last index
    
    # Find R peak value within each QRS segment:
    Rpeaks <- array(0, length(QRS_cluster_index) - 1)
    for (i in 1:(length(QRS_cluster_index) - 1)) {
      sample_range <- signal[sample, QRS_predict[QRS_cluster_index[i]:(QRS_cluster_index[i + 1] - 1)]]
      Rpeaks[i] <- which.max(sample_range) + QRS_predict[QRS_cluster_index[i]] - 1
    }
    Rpeaks_full[sample] <- list(Rpeaks)
  }
  return(Rpeaks_full)
}

# RR Table ----------------------------------------------------------------
RR_table_make <- function(Rpeaks) {
  # Helper function for formatting
  
  R_on <- c()
  R_off <- c()
  Rpeaks_sample <- c()
  # test_rows_sample <- c()
  
  for (i in 1:length(Rpeaks)) {
    Rpeaks_single <- unlist(Rpeaks[i])
    
    R_on <- c(R_on, Rpeaks_single[1:(length(Rpeaks_single) - 1)])
    R_off <- c(R_off, Rpeaks_single[-1])
    Rpeaks_sample <- c(Rpeaks_sample, array(i, length(Rpeaks_single) - 1))
    
  }
  
  RR_info <- data.frame(Onset = R_on,
                        Offset = R_off,
                        Sample = Rpeaks_sample)
  
  return(RR_info)
}

# RR_extraction -----------------------------------------------------------
splice_signal <- function(RR_info,signal) {
  
  # in if statement: use RR_info samples. If multiple sample numbers, then...
  # Using R peak values, splice samples into individual R to R interval samples
  
  if (length(dim(signal)) == 2) {
    # Signal is not a spectrogram:
    signal_spliced <- array(0, c(dim(RR_info)[[1]], 1024))
    
    for (i in 1:dim(RR_info)[[1]]) {
      length <- length(RR_info$Onset[i]:RR_info$Offset[i])
      
      if (length > 1024) {
        stop(paste("R-R interval within Sample",RR_info$Sample[i],"is greater than 1024 limit"))
        # stop()
      }
      
      signal_spliced[i, 1:length] <- signal[RR_info$Sample[i], RR_info$Onset[i]:RR_info$Offset[i]]
    }
    
  } else if (length(dim(signal)) == 3) {
    signal_spliced <- array(0, c(dim(RR_info)[[1]], 1024, dim(signal)[[3]]))
    
    for (i in 1:dim(RR_info)[[1]]) {
      length <- length(RR_info$Onset[i]:RR_info$Offset[i])
      
      if (length > 1024) {
        stop(paste("R-R interval within Sample",RR_info$Sample[i],"is greater than 1024 limit"))
        # stop()
      }
      
      signal_spliced[i, 1:length, ] <- signal[RR_info$Sample[i], RR_info$Onset[i]:RR_info$Offset[i], ]
    }
    
    # else if dim(signal) == NULL
  }
  return(signal_spliced)
}

# Single_waveform_prediction -----------------------------------------------------------
RR_waveform_prediction <- function(input,model_name,model_path) {
  library(caret)
  library(keras)
  
  model <- load_model_tf(paste0(model_path, model_name))
  
  predictions <- model %>% predict(input)
  
  a <- dim(predictions)
  predictions_integer <- array(0, a[1:2])
  
  # use max.col to determine most probable value  for each time point
  for (i in 1:a) {
    predictions_integer[i, ] <- max.col(predictions[i, , ])
  }
  
  #convert from dimension value 1,2,3,4 to 0,1,2,3
  annotations <- predictions_integer - 1
  
  return(annotations)
  
}

# Stitch RR Annotations into their respective samples ---------------------
RR_stitch_ann <- function(annotations, RR_info) {

  annotations_stitched <- array(0,c(max(RR_info$Sample),5000))
  
  for (i in 1:max(RR_info$Sample)) {
    for (j in which(RR_info$Sample == i)) {
      range <- RR_info$Onset[j]:RR_info$Offset[j]
      annotations_stitched[i,range] <- annotations[j,1:length(range)]
    }
  }
  return(annotations_stitched)
}

# Crude 250 Hz upscale
upscale_250 <- function(signal) {
  
  if (is.vector(signal) == TRUE) {
    signal_upsample <- array(0, length(signal))
    
    for (i in 1:2499) {
      signal_upsample[(i * 2 - 1):(i * 2)] <- c(signal[i], mean(signal[i:(i + 1)]))
    }
    signal_upsample[4999:5000] <- signal[2500]
  } else {
    
    signal_upsample <- array(0, dim(signal))
    
    for (row in 1:dim(signal)[[1]]) {
      for (i in 1:2499) {
        signal_upsample[row, (i * 2 - 1):(i * 2)] <- c(signal[row, i], mean(signal[row, i:(i + 1)]))
      }
      signal_upsample[row, 4999:5000] <- signal[row, 2500]
    }
  }
  return(signal_upsample)
}



# WFDB Annotation Format --------------------------------------------------
annotation_table_creator <- function(annotations, Fs = 500) {
  
  p_waves <- which(annotations == 1)
  p_change <- ((p_waves[-1] - p_waves[1:(length(p_waves) - 1)]))
  new_pwave <- c(1, ((which(p_change!= 1) + 1)))
  new_pwave <- c(new_pwave, length(p_waves) + 1) # account for last value
  
  p_on <- c()
  p_off <- c()
  p_max <- c()
  
  for (i in 1:(length(new_pwave) - 1) ) {
    p_on <- c(p_on, p_waves[new_pwave[i]])
    p_off <- c(p_off, p_waves[new_pwave[i + 1] - 1])
    p_max <- c(p_max, round(mean(c(p_on[i], p_off[i])))) # temporary, will need to find peak point
  }
  
  
  qrs_waves <- which(annotations == 2)
  qrs_change <- ((qrs_waves[-1] - qrs_waves[1:(length(qrs_waves) - 1)]))
  new_qrswave <- c(1, ((which(qrs_change!= 1) + 1)))
  new_qrswave <- c(new_qrswave, length(qrs_waves) + 1) # account for last value
  
  qrs_on <- c()
  qrs_off <- c()
  qrs_max <- c()
  
  for (i in 1:(length(new_qrswave) - 1) ) {
    qrs_on <- c(qrs_on, qrs_waves[new_qrswave[i]])
    qrs_off <- c(qrs_off, qrs_waves[new_qrswave[i + 1] - 1])
    qrs_max <- c(qrs_max, round(mean(c(qrs_on[i], qrs_off[i])))) # temporary, will need to find peak point
  }
  
  t_waves <- which(annotations == 3)
  t_change <- ((t_waves[-1] - t_waves[1:(length(t_waves) - 1)]))
  new_twave <- c(1, ((which(t_change!= 1) + 1)))
  new_twave <- c(new_twave, length(t_waves) + 1) # account for last value
  
  t_on <- c()
  t_off <- c()
  t_max <- c()
  
  for (i in 1:(length(new_twave) - 1) ) {
    t_on <- c(t_on, t_waves[new_twave[i]])
    t_off <- c(t_off, t_waves[new_twave[i + 1] - 1])
    t_max <- c(t_max, round(mean(c(t_on[i], t_off[i])))) # temporary, will need to find peak point
  }
  
  t <- (0 : length(annotations) - 1 ) / Fs
  

  # Build Table
  type <- array(c("(","p",")"),length(p_on)*3)
  type <- c(type,array(c("(","N",")"),length(qrs_on)*3))
  type <- c(type,array(c("(","t",")"),length(t_on)*3))
  
  sample <- c()
  for (i in 1:length(p_on)) {
    sample <- c(sample, p_on[i], p_max[i], p_off[i])
  }
  for (i in 1:length(qrs_on)) {
    sample <- c(sample, qrs_on[i], qrs_max[i], qrs_off[i])
  }
  for (i in 1:length(t_on)) {
    sample <- c(sample, t_on[i], t_max[i], t_off[i])
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
}

# Custom Plot -------------------------------------------------------------

plot_func <- function(x, y, color = 0,linewidth=0.5) {
  library(ggplot2)
  library(plotly)
  color <- c(color)
  
  color[color == 1] <- 'p'
  color[color == 2] <- 'N'
  color[color == 3] <- 't'
  
  x <- 1:length(y)
  
  frame <- data.frame(Time = x, Signal = y)
  plot <-
    ggplotly(
      ggplot(frame, aes(Time, Signal, color = color)) +
        geom_path(linewidth = linewidth, aes(group = 1)) + geom_point() +
        scale_x_continuous(breaks = seq(0, 10, 1)) + theme(legend.position =
                                                             "none") + theme(legend.title = element_blank())
    )
  
  # original linewidth = 0.25
  return(plot)
}


