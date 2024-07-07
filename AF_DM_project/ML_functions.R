# Filter Samples ----------------------------------------------------------
filter_samples <- function(signal, frequency = 500, low = 0.5, high = 40) {
  # Butterworth bandpass filter
  library(dplR)
  
  if (is.vector(signal) == TRUE) {
    signal <- matrix(signal, nrow = 1)
  }
  
  filtered <- array(0, dim(signal))
  for (i in 1:dim(signal)[[1]]) {
    filtered[i, ] <- pass.filt(
      signal[i, ],
      c(f / high, f / low),
      "pass",
      method = c("Butterworth"),
      n = 4,
      Rp = 0.1
    )
  }
return(filtered)
}

# Build Spectrogram -------------------------------------------------------
build_spectrogram <- function(input, Fs = 500, window_size = 128, high_freq_cutoff = 45) {
  # Builds a spectrogram, outputs real and imaginary values separately
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
  # Prediction using model of choice
  
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
peak_isolation <- function(signal, annotations, wave_value = 2) {
  # Determine peak values of a wave of choice. Within a single continuous wave, 
  # Finds points which have a chance in slope value. Of those, find the value
  # the furthest from the isoelectric line
  
  # input signal: raw signal, not spectrogram
  # Use 10-sec duration model ( predict_samples() ) to determine annotations
  
  if (is.vector(signal) == TRUE) {
    signal <- matrix(signal, nrow = 1)
  }
  
  if (is.vector(annotations) == TRUE) {
    annotations <- matrix(annotations, nrow = 1)
  }
  
  peaks_full <- array(0,dim(signal)[[1]])
  
  
  for (sample in 1:dim(signal)[[1]]) {
    
    # midline <- median(signal)
    midline <- isoelec_find(signal[sample,], annotations[sample,])
    # Find each continuous QRS segment:
    
    # wave_cluster: which time points are QRS
    # wave_bounds: indices of QRS_predict which jump to the next QRS interval
    
    wave_cluster <-  which(annotations[sample, ] == wave_value)
    change <- (wave_cluster[-1] - wave_cluster[1:(length(wave_cluster) - 1)])
    
    wave_bounds <- c(1)
    wave_bounds <- c(wave_bounds, (which(change != 1) + 1))
    wave_bounds <- c(wave_bounds, length(wave_cluster))
    
    # Find R peak value within each QRS segment:
    peaks <- array(0, length(wave_bounds) - 1)
    for (i in 1:(length(wave_bounds) - 1)) {
      sample_range <- signal[sample, wave_cluster[wave_bounds[i]:(wave_bounds[i + 1] - 1)]]
      
      if (length(sample_range) < 40) {
        #window is 1/3 of range, odd number
        window <- ceiling(2*round(length(sample_range) / 6) + 1)
      } else {
      window <- 11
      }
      
      sample_range_mean <- rollmean_custom(sample_range,window)
      
      # find which indices where the slope changes value
      candidates <- which(diff(sign(diff(sample_range_mean))) != 0)+ wave_cluster[wave_bounds[i]]
      
      # If there are no candidate points:
      if (length(candidates) ==0) {
        peaks[i] <- NA
        text <- paste("No candidate points chosen on sample", sample, "wave number",i)
        warning(text)
      } else {
      peaks[i] <- candidates[which.max(abs(signal[candidates] - c(midline)))]
      }
      
      # Previous method using max value:
      # peak_max <- which.max(sample_range) + wave_cluster[wave_bounds[i]] - 1
      # peak_min <- which.min(sample_range) + wave_cluster[wave_bounds[i]] - 1
      # 
      # if ( abs(signal[peak_max] - midline) > abs(signal[peak_min] - midline)) {
      #   peaks[i] <- peak_max
      # } else {
      #   peaks[i] <- peak_min
      # }
    
    }
    peaks_full[sample] <- list(peaks)
  }
  return(peaks_full)
}



# Isoelectric Line --------------------------------------------------------
isoelec_find <- function(signal,annotations) {
  # Finds the mean value of the T-P intervals within the given sample lead
  # take mean vs. median of T-P intervals?
  
  pwaves <- make_wave_table(annotations, wave_value = 1)
  twaves <- make_wave_table(annotations, wave_value = 3)

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
    isoelectric_line[i] <- mean(signal[pt_ind])
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



# RR Table ----------------------------------------------------------------
RR_table_make <- function(Rpeaks) {
    # Helper function for formatting
    # Creates table of R-R intervals, dividing a multi-beat sample
  
  if (class(Rpeaks) != "list") {
    Rpeaks <- list(Rpeaks)
  }
  
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
  # Using the R-R helper function, signal is divided into R-R interval chunks
  # For R-R interval annotation models (as opposed to 10 sec annotation models)
  
  
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
        paste("R-R interval within Sample",RR_info$Sample[i],"is greater than 1024 limit")

      } else {
      signal_spliced[i, 1:length, ] <- signal[RR_info$Sample[i], RR_info$Onset[i]:RR_info$Offset[i], ]
      }
    }
    
    # else if dim(signal) == NULL
  }
  return(signal_spliced)
}

# Single_waveform_prediction -----------------------------------------------------------
RR_waveform_prediction <- function(input,model_name,model_path) {
  # R-R annotation prediction (as opposed to 10 sec annotation prediction)
  library(caret) # probably don't need. Probably just need keras
  library(keras) 
  
  # will want to avoid calling libraries inside functions for toher users 
  
  # default label for confidence, add parameter to be __% sure that it is a p/qrs/t wave.
  # or add in a function for confidence interval of each time point
  #Could combine these^ across all 12 leads... 
  
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
  # After R-R model annotations, this functions combines the R-R intervals back
  # into their full combined length. Able to handle multiple full-length 
  # samples

  annotations_stitched <- array(0,c(max(RR_info$Sample),5000))
  
  for (i in 1:max(RR_info$Sample)) {
    for (j in which(RR_info$Sample == i)) {
      range <- RR_info$Onset[j]:RR_info$Offset[j]
      if (length(range) > 1024) {
        paste("Sample",i,"contains a segment longer than 1024")
        annotations_stitched[i,range] <- array(0,length(range))
      } else{
      annotations_stitched[i,range] <- annotations[j,1:length(range)]
      }
    }
  }
  return(annotations_stitched)
}


# Wave Table --------------------------------------------------------------
make_wave_table <- function(annotations,  wave_value) {
  # Gives table of specified on/offset times for any wave- P/QRS/T 
  # Multi-sample input handling using "sample" column
  
  # Could include time values, in addition to indeces
  if (is.vector(annotations) == TRUE) {
    annotations <- matrix(annotations, nrow = 1)
  }
  
  wave_on <- c()
  wave_off <- c()
  sample <- c()
  
  for (i in 1:dim(annotations)[[1]]) {
    wave_cluster <-  which(annotations[i, ] == wave_value)
    
    change <- (wave_cluster[-1] - wave_cluster[1:(length(wave_cluster) - 1)])
    
    wave_on <- c(wave_on, wave_cluster[1])
    wave_on <- c(wave_on, wave_cluster[which(change != 1) + 1])
    
    wave_off <- c(wave_off, wave_cluster[which(change != 1)])
    wave_off <- c(wave_off, wave_cluster[length(wave_cluster)])
    
    sample <- c(sample, array(i, sum(change != 1) + 1))
    
  }
  
  wave_classes <- c(0,"p","N","t")
  wave_type <- array(wave_classes[wave_value+1],length(wave_on))
  
  wave_table <- data.frame(wave_type = wave_type, wave_on = wave_on, wave_off = wave_off, sample = sample)
  
  return(wave_table)
}


# Crude 250 Hz upscale -------------------------------------------------------------------------
upscale_250 <- function(signal) {
  # Used only for model testing. Simple upsampling from 250 Hz to 500 Hz for 
  # MIT-BIH dataset
  
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
  # Not up to date, largely defunct for now. Used to convert index to index
  # annotation (value for each time point) to wave onset/peak/offset format. 
  # WFDB format
  
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


# Custom rollmean ---------------------------------------------------------
rollmean_custom <- function(signal,window) {
  # Custom rolling average function to avoid clipping terminal values.
library(zoo)
signal_mean <- rollmean(x = signal, k = window, fill = NA)

end <- length(signal_mean)

if (window > 1) {
for (i in 1: ((window - 1)/2)) {
  signal_mean[i] <- mean(signal[1: (i + (window - 1)/2) ])
  signal_mean[end - i + 1] <- mean(signal[(end - i + 1 - (window - 1)/2) : end ])
}
}

return(signal_mean)

}

# Custom 12 lead plot -----------------------------------------------------
plot_func12 <- function(y, color = 0) {
  # Input 12 lead data for  1 sample
  library(ggplot2)
  library(ggpubr)

  i   <- plot_func(y = y[,1], color = color[,1], plotly = "no")
  ii  <- plot_func(y = y[,2], color = color[,2], plotly = "no")
  iii <- plot_func(y = y[,3], color = color[,3], plotly = "no")
  avr <- plot_func(y = y[,4], color = color[,4], plotly = "no")
  avl <- plot_func(y = y[,5], color = color[,5], plotly = "no")
  avf <- plot_func(y = y[,6], color = color[,6], plotly = "no")
  v1  <- plot_func(y = y[,7], color = color[,7], plotly = "no")
  v2  <- plot_func(y = y[,8], color = color[,8], plotly = "no")
  v3  <- plot_func(y = y[,9], color = color[,9], plotly = "no")
  v4  <- plot_func(y = y[,10], color = color[,10], plotly = "no")
  v5  <- plot_func(y = y[,11], color = color[,11], plotly = "no")
  v6  <- plot_func(y = y[,12], color = color[,12], plotly = "no")
  
  
  
  # subplot(i,ii,iii,avr,avl,avf,v1,v2,v3,v4,v5,v6, nrows = 4)
  ggarrange(i,ii,iii,avr,avl,avf,v1,v2,v3,v4,v5,v6, ncol = 4, nrow = 3)
}

# Custom Plot -------------------------------------------------------------
plot_func <- function(x, y, color = 0,linewidth=0.5, plotly = 'yes') {
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
        geom_path(linewidth = linewidth, aes(group = 1)) + geom_point() +
        scale_x_continuous(breaks = seq(0, 10, 1)) + theme(legend.position =
                                                             "none") + theme(legend.title = element_blank())
  
  if (plotly == 'yes') {
  plot <- ggplotly(plot)
  }
  
  # original linewidth = 0.25
  return(plot)
}


