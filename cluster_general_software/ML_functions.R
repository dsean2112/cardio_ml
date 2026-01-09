# Filter Samples ----------------------------------------------------------
# filter_samples <- function(signal, frequency = 500, low = 0.5, high = 40) {
#   # Butterworth bandpass filter
#   library(dplR)
#   
#   # If vector, set to array
#   if (is.vector(signal) == TRUE) {
#     signal <- matrix(signal, nrow = 1)
#   }
#   
#   # Filter each line
#   filtered <- array(0, dim(signal))
#   for (i in 1:dim(signal)[[1]]) {
#     filtered[i, ] <- pass.filt(
#       signal[i, ],
#       c(frequency / high, frequency / low),
#       "pass",
#       method = c("Butterworth"),
#       n = 4,
#       Rp = 0.1
#     )
#   }
# return(filtered)
# }

# NEW function:
filter_samples <- function(signal, frequency = 500, low = 0.5, high = 40) {
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

# Build Spectrogram -------------------------------------------------------
build_spectrogram <- function(input, Fs = 500, window_size = 128, high_freq_cutoff = 45) {
  # Builds a spectrogram, outputs real and imaginary values separately
  # Input is a *single* vector. Currently does not handle
  
  # Currently only works with 500 Hz
  
  library(signal)
  
  # Coordinate overlap, padding and window size
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
  # Output: single integer for each time point
  
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
simple_Rpeak_isolation <- function(signal, annotations) {
  wave_value= 2

  if (is.vector(signal) == TRUE) {
    signal <- matrix(signal, nrow = 1)
  }
  
  if (is.vector(annotations) == TRUE) {
    annotations <- matrix(annotations, nrow = 1)
  }
  
  peaks_full <- array(0,dim(signal)[[1]])
  
  
  for (sample in 1:dim(signal)[[1]]) {
    
 
    # Find each continuous QRS segment:
    
    # wave_cluster: which time points are QRS
    # wave_bounds: indices of QRS_predict which jump to the next QRS interval
    
    wave_cluster <-  which(annotations[sample, ] == wave_value)
    
    if (length(wave_cluster) <= 1) {
      waves <- c('p','N','t')
      print(paste('No',waves[wave_value],'wave detected'))
      break
    }
    
    change <- (wave_cluster[-1] - wave_cluster[1:(length(wave_cluster) - 1)])
    
    wave_bounds <- c(1)
    wave_bounds <- c(wave_bounds, (which(change != 1) + 1))
    wave_bounds <- c(wave_bounds, length(wave_cluster))
    
    midline <- median(signal[sample])
    # Find R peak value within each QRS segment:
    peaks <- array(0, length(wave_bounds) - 1)
    for (i in 1:(length(wave_bounds) - 1)) {
      sample_range <- signal[sample, wave_cluster[wave_bounds[i]:(wave_bounds[i + 1] - 1)]]
      

      # find which indices where the slope changes value
      candidates <- which(diff(sign(diff(sample_range))) != 0)+ wave_cluster[wave_bounds[i]]
      
      # If there are no candidate points:
      if (length(candidates) ==0) {
        peaks[i] <- NA
        text <- paste("No candidate points chosen on sample", sample, "wave number",i)
        # warning(text)
      } else {
        peaks[i] <- candidates[which.max(signal[candidates] - c(midline))]
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
  

peak_isolation <- function(signal, annotations, wave_value = 2, midline_calc = "isoelec") {
  # Determine peak values of a wave of choice. Within a single continuous wave, 
  # Finds points which have a chance in slope value. Of those, find the value
  # the furthest from the isoelectric line
  
  # input signal: raw signal, not spectrogram
  # Use 10-sec duration model ( predict_samples() ) to determine annotations
  
  # Midline_calc option: a midline is needed to determine peaks. For samples 
  # without dependable P/T annotations for T-P interval isoelectric point 
  # calculation, set midline_calc to FALSE. In this case, the median of the 
  # entire ECG will be used instead. This is meant mainly for R peak calculation
  
  if (is.vector(signal) == TRUE) {
    signal <- matrix(signal, nrow = 1)
  }
  
  if (is.vector(annotations) == TRUE) {
    annotations <- matrix(annotations, nrow = 1)
  }
  
  peaks_full <- array(0,dim(signal)[[1]])
  
  
  for (sample in 1:dim(signal)[[1]]) {
    
    # midline <- median(signal)
    if (midline_calc == "isoelec") {
      midline <- isoelec_find(signal[sample, ], annotations[sample, ])
    } else {
      midline <- median(signal[sample,])
    }
    
    if (!any(annotations == wave_value)) {
      peaks_full[sample] <- NA
      next
    }
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
      
      # Set window for rolling average
      if (length(sample_range) < 40) {
        #window is 1/3 of range, odd number
        window <- ceiling(2*round(length(sample_range) / 6) + 1)
      } else {
      window <- 11
      }
      
      # Rolling average
      sample_range_mean <- rollmean_custom(sample_range,window)
      
      # find which indices where the slope changes value
      candidates <- which(diff(sign(diff(sample_range_mean))) != 0)+ wave_cluster[wave_bounds[i]]
      
      # If there are no candidate points:
      if (length(candidates) ==0) {
        peaks[i] <- NA
        text <- paste("No candidate points chosen on sample", sample, "wave number",i)
        # warning(text)
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
  
  # will want to avoid calling libraries inside functions for other users 
  
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


# Crude 250 Hz upscale -------------------------------------------------------------------------
upscale_250 <- function(signal) {
  # Used only for model testing. Simple upsampling from 250 Hz to 500 Hz for 
  # MIT-BIH dataset
  
  if (is.vector(signal) == TRUE | length(dim(signal)) == 1) {
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

  i   <- plot_func(y = y[,1], color = color[,1], plotly = "no") + theme_void() + theme(legend.position = "none")
  ii  <- plot_func(y = y[,2], color = color[,2], plotly = "no") + theme_void() + theme(legend.position = "none")
  iii <- plot_func(y = y[,3], color = color[,3], plotly = "no") + theme_void() + theme(legend.position = "none")
  avr <- plot_func(y = y[,4], color = color[,4], plotly = "no") + theme_void() + theme(legend.position = "none")
  avl <- plot_func(y = y[,5], color = color[,5], plotly = "no") + theme_void() + theme(legend.position = "none")
  avf <- plot_func(y = y[,6], color = color[,6], plotly = "no") + theme_void() + theme(legend.position = "none")
  v1  <- plot_func(y = y[,7], color = color[,7], plotly = "no") + theme_void() + theme(legend.position = "none")
  v2  <- plot_func(y = y[,8], color = color[,8], plotly = "no") + theme_void() + theme(legend.position = "none")
  v3  <- plot_func(y = y[,9], color = color[,9], plotly = "no") + theme_void() + theme(legend.position = "none")
  v4  <- plot_func(y = y[,10], color = color[,10], plotly = "no") + theme_void() + theme(legend.position = "none")
  v5  <- plot_func(y = y[,11], color = color[,11], plotly = "no") + theme_void() + theme(legend.position = "none")
  v6  <- plot_func(y = y[,12], color = color[,12], plotly = "no") + theme_void() + theme(legend.position = "none")
  
  # ****Revised version:
  # plot_list <- lapply(seq_len(nrow(y)), function(i) {
  #   plot_func(y[i, ],plotly = FALSE) +  + theme_void() + theme(legend.position = "none") # Pass each row as a vector to plot_func
  # })
  # ggarrange(plotlist = plot_list)
  
  plot_names <- c('i','ii','iii','avr','avl','avf','v1','v2','v3','v4','v5','v6')
  
  
  
  # subplot(i,ii,iii,avr,avl,avf,v1,v2,v3,v4,v5,v6, nrows = 4)
  ggarrange(i,ii,iii,avr,avl,avf,v1,v2,v3,v4,v5,v6, ncol = 4, nrow = 3)
}

# Custom Plot -------------------------------------------------------------
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


# Confidence Plot ---------------------------------------------------------
plot_confidence <- function(y,color=0) {
  # library(ggplot2)
  library(plotly)
  library(RColorBrewer)
  
  color <- c(color)
  color <- round(color,2)
  
  y <- c(y)
  
  x <- 1:length(y)
  
  frame <- data.frame(Time = x, Signal = y, Confidence = color)
  
  fig <- plot_ly(frame, x = ~Time)
  fig <- fig %>% add_trace(data = frame, y = ~Signal, type = 'scatter', 
                           mode = 'lines', color = I('black')) %>% 
    add_trace(y = ~Signal, color = ~Confidence, size = 10,
              type = 'scatter', mode = 'markers',
              text = ~paste("Conf:",Confidence))
    
  
  return(fig)
}


# WFDB to Continuous Annotation -------------------------------------------
ann_wfdb_to_complete <- function(wfdb) {
  # Transform wfdb annotation format into a continuous annotation format for
  # Plotting, ML, etc. 
  # Continuous: 1 annotation per time step
  
  continuous_ann <- array(0,5000)
  
  open_parenth_index <- which(wfdb$type == '(')
  skip_index <- which( wfdb$type[open_parenth_index + 1] == ')' )
  if (length(skip_index > 0)) {
    print(paste0("No wave at ()", "position ",skip_index))
    
    # for (n in 1:length(skip_index)) {
    #   missed_waves[nrow(missed_waves) + 1, ] = c(AFib_index[rec], leads[ann], skip_index[n])
    # }
    
  }
  #build solutions array:
  # dimmensions <- dim(test$signal)
  # markers <- matrix(0, dimmensions[1])
  
  # for each pwave, determine on/off points, and fill the corresponding time range with '1' in markers array
  pwaves <- which(wfdb$type == 'p')
  if (length(pwaves) > 0) {
    for (i in 1:length(pwaves)) {
      
      # Verify proper ( p ) format:
      
      logical <- c(wfdb$type[pwaves[i]-1] == '(', 
                   wfdb$type[pwaves[i]+1] == ')',
                   pwaves[i] != 1)
      
      if (!(sum(logical) !=3 | is.na(sum(logical)) == TRUE)) {
        bounds <- wfdb$sample[c(pwaves[i] - 1, pwaves[i] + 1)]
        continuous_ann[(bounds[1] + 1):(bounds[2] + 1)] <- 1
      } else {
        print(paste('Improper wave format on Position ', pwaves[i]))
        # missed_waves[nrow(missed_waves) + 1,] = c(AFib_index[rec],leads[ann],pwaves[i])
      }
      
      # markers[(bounds[1] + 1):(bounds[2] + 1), ] <-
      #   1 #'p' #remove "+1"??
      
    }
  } else{
    # print(paste("No P-Wave on Sample", rec))
  }
  
  # for each qrs, fill range with '2'
  qrswaves <- which(wfdb$type == 'N')
  if (length(qrswaves) > 0) {
    for (i in 1:length(qrswaves)) {
      
      
      logical <- c(wfdb$type[qrswaves[i]-1] == '(', 
                   wfdb$type[qrswaves[i]+1] == ')',
                   qrswaves[i] != 1)
      
      if (!(sum(logical) !=3 | is.na(sum(logical)) == TRUE)) {
        bounds <- wfdb$sample[c(qrswaves[i] - 1, qrswaves[i] + 1)]
        continuous_ann[(bounds[1] + 1):(bounds[2] + 1)] <- 2
      } else {
        print(paste('Improper wave format on Position ',qrswaves[i]))
        # missed_waves[nrow(missed_waves) + 1,] = c(AFib_index[rec],leads[ann],qrswaves[i])
      }
      
      # markers[(bounds[1] + 1):(bounds[2] + 1), ] <- 2 #'N'
      continuous_ann[(bounds[1] + 1):(bounds[2] + 1)] <- 2
    }
  } else{
    # print(paste("No QRS-Complex on Sample", rec))
  }
  
  
  # for each twave, fill the range with '3'
  twaves <- which(wfdb$type == 't')
  if (length(twaves) > 0) {
    for (i in 1:length(twaves)) {
      
      logical <- c(wfdb$type[twaves[i]-1] == '(', 
                   wfdb$type[twaves[i]+1] == ')',
                   twaves[i] != 1)
      
      if (!(sum(logical) !=3 | is.na(sum(logical)) == TRUE)) {
        bounds <- wfdb$sample[c(twaves[i] - 1, twaves[i] + 1)]
        continuous_ann[(bounds[1] + 1):(bounds[2] + 1)] <- 3
      } else {
        print(paste('Improper wave format on Position ', twaves[i]))
        # missed_waves[nrow(missed_waves) + 1,] = c(AFib_index[rec],leads[ann],twaves[i])
      }
      
    }
  } else{
    # print(paste("No T-Wave on Sample", rec))
  }
  
  return(continuous_ann)
  
}
 