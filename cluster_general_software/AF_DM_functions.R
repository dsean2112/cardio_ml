# HRV ---------------------------------------------------------------------
HRV_RMSSD <- function(Rpeaks, Hz = 500, annotations = 0) {
  # Find RMSSD for 10 sec ECG
  # Input should be index values, not time values
  
  # Annotations are used to check if terminal QRS waves should be removed
  
  if (length(unlist(Rpeaks)) <= 1) {
    print(paste('No Rpeaks in input'))
    return(NA)
    stop()
  }
  
  # Multi-sample input:
  if (class(Rpeaks) == "list") {
    RMSSD <- array(0,length(Rpeaks))
    for (i in (1:length(Rpeaks))) {
      Rpeaks_indv <- unlist(Rpeaks[i])
      
      if (is.na(Rpeaks_indv[1]) == TRUE) {
        Rpeaks_indv <- Rpeaks_indv[-1]
      }
      
      if (is.na(Rpeaks_indv[length(Rpeaks_indv)] == TRUE)) {
        Rpeaks_indv <- Rpeaks_indv[-(length(Rpeaks_indv))]
      }
      
      
      if (length(Rpeaks_indv) < 3) {
        print(paste('Not enough Rpeaks'))
        return(NA)
        stop()
      }
      
      Rpeaks_indv <- Rpeaks_indv / Hz
      interval <- (Rpeaks_indv[-1] - Rpeaks_indv[1: (length(Rpeaks_indv) - 1)]) * 1000 # convert to ms
      RMSSD[i] <- sqrt( 1/(length(interval)- 1) * sum((interval[1: (length(interval) - 1)] - interval[-1])^2))
    }
    
  # Single-sample input:  
  # } else if (is.vector(Rpeaks) == TRUE) {
  } else if (class(Rpeaks) != 'list') {

    if (is.na(Rpeaks[1]) == TRUE) {
      Rpeaks <- Rpeaks[-1]
    }
    
    if (is.na(Rpeaks[length(Rpeaks)] == TRUE)) {
      Rpeaks <- Rpeaks[-(length(Rpeaks))]
    }
    
    if (length(Rpeaks) < 3) {
      print(paste('Not enough Rpeaks'))
      return(NA)
      stop()
    }
    
    Rpeaks = Rpeaks / Hz
    interval <- (Rpeaks[-1] - Rpeaks[1:(length(Rpeaks) - 1)]) * 1000 # convert to ms
    RMSSD <- sqrt(1 / (length(interval) - 1) * sum((interval[1:(length(interval) - 1)] - interval[-1]) ^2))
  }
  
  return(RMSSD)
 
}

# Wave Axis ---------------------------------------------------------------
find_wave_axis <- function(signal12, ann12, wave_value) {
  # Input: filtered signal**
  # Output: 3D spatial vectors for each beat
  # Find wave axis for each beat. Input: 12 lead signal and annotations.
  
  # Rpeaks <- peak_isolation(signal12[,1],ann12[,1])
  
  # waves <- c()
  
  
  if (length(dim(ann12)) < 2) {
    ann12 <- array(rep(ann12, 12), c(5000,12))
  }
  
  # If there are very few values with the specified wave value, stop the function:
  if (sum(ann12 == wave_value, na.rm=TRUE) < 24) {
    output <- list(array(c(NA,NA,NA),c(1,3)),array(c(NA,NA,NA),c(1,3)))
    names(output) <- c('every_wave','mean_wave')
    warning('No waves of value',wave_value,'were found. Exiting.')
    return(output)
    break
    
  }
  
  # If any leads used for the Kors transformation (I-II,V1-V6) are NA, break
  if (any(is.na(ann12[,c(1:2,7:12)]))) {
    output <- list(array(c(NA,NA,NA),c(1,3)),array(c(NA,NA,NA),c(1,3)))
    names(output) <- c('every_wave','mean_wave')
    warning('No waves of value ',wave_value,' were found. Exiting.')
    return(output)
    break
    
  }
  
  # Lists of AUC and wave onset values for all 12 leads
  AUC_full <- array(0,12)
  wave_on_full <- array(0,12)
  
  # Remove waves within specified value of beginning/end of the ECG:
  terminal_window <- 75
  
  for (i in 1:12) {
    wave_table <- find_wave_AUC(signal = signal12[,i], annotations = ann12[,i], wave_value = 1)
    
    remove_rows <- which(wave_table$wave_on < terminal_window | wave_table$wave_off > (5000 - terminal_window))
    if (length(remove_rows) > 0) {
    wave_table <- wave_table[-remove_rows,]
    }
    
    AUC_full[i] <- list(round(wave_table$AUC,2))
    wave_on_full[i] <- list(wave_table$wave_on)
    
  }
  
  # window: aligning the wave across each lead over __ frames- ie waves don't
  # occur at identical times. 
  # Because we can't assume each p-wave is the same index (one lead could miss 
  # a wave), the leads are cross-chekced relative to lead 2 (most accurate lead)
  window <- 75
  
  # For each wave (using lead 2 as reference):
  #   For each lead:
  #     Find AUC for each wave, then find 3D spatial vector
  wave_3D_vector <- array(0,c(length(unlist(AUC_full[2])), 3))
  AUC_12lead <- array(NA,c(length(unlist(AUC_full[2])),12))
  
  # If lead 2 has no waves of interest, exit:
  if (length(unlist(wave_on_full[2])) == 0 | sum(is.na(unlist(wave_on_full[2]))) > 0) {
    output <- list(array(c(NA,NA,NA),c(1,3)),array(c(NA,NA,NA),c(1,3)))
    names(output) <- c('every_wave','mean_wave')
    warning('ECG has no p-waves.')
    return(output)
    break
  }
  
  for (i in 1:length(unlist(AUC_full[2]))) {
    # Use most accurate lead- lead 2
    
    
    index <- unlist(wave_on_full[2])[[i]] # Find wave onset index
    
    # Find AUC for each lead which is present within 50 indices of lead 2 onset
    
    for (j in 1:12) {
      # Isolate lead of interest
      wave_on_indv <- unlist(wave_on_full[j])
      AUC_indv <- unlist(AUC_full[j])
      
      # Find which index within the lead corresponds to lead 2 index:
      wave_on_index <- which(wave_on_indv < index + window & wave_on_indv > index - window)
      
      # If there are multiple waves within range, pick the closest:
      if (length(wave_on_index) > 1) {
      wave_on_index <- wave_on_index[which.min(abs(wave_on_indv[wave_on_index] - index))]
      }
      
      # Print warning message if there is no corresponding wave in the lead:
      if (length(wave_on_index) == 0) {
        text <- paste("Warning: Lead #",j,"Beat #",i,"has no wave!")
        # warning(text)
      } else {
      # Array of AUC for each lead
      AUC_12lead[i,j] <- AUC_indv[wave_on_index]
      }
    }
  }
  
  # Calculate 3D vector for each beat:
  #   set rows 3 to 6 equal to zero- to remove NA values from unnecessary leads
  AUC_12lead_adjusted <- AUC_12lead
  AUC_12lead_adjusted[, c(3:6)] <- 0
  
  for (i in 1:nrow(AUC_12lead_adjusted)) {
    wave_3D_vector[i, ] <- kors(AUC_12lead_adjusted[i, ])
  }
  
  AUC_mean <- array(colMeans(AUC_12lead,na.rm = TRUE),c(1,12))
  AUC_mean[, c(3:6)] <- 0
  mean_vector <- round(kors(AUC_mean),2)
  
  
  output <- list(wave_3D_vector,mean_vector)
  names(output) <- c('every_wave','mean_wave')
  
  return(output)
  # For every p-wave in lead II:
  #   Find AUC of p wave in other leads
  
  
  # Verify wave_on is within ~50 indices
  
  
  # compare beat to beat pwaves AUC
  
  # As a check, can compare AUC to QR, RS intervals
}

# Wave AUC ----------------------------------------------------------------

find_wave_AUC <- function(signal, annotations, wave_value) {
  # NOTE: use filtered signal as input
  # Wave_value: 1,2 or 3 for P, QRS or T.
  
  # Finds the average AUC of a given wave type, for a single sample
  # baseline for the AUC is the median of the T-P interval
  
  # Currently, AUC does not use time values for x, but this can be adjusted.
  
  # Find isoelectric point: 
  isoelec <- isoelec_find(signal, annotations)
  library(DescTools)
  
  if (sum(annotations == wave_value) == 0) {
    warning('No waves of specified value in ECG.')
    wave_table <- make_wave_table(annotations, wave_value)
    wave_table$AUC = NA
    
  } else {
    # Make table to keep track of onset/offset of wave of interest:
    wave_table <- make_wave_table(annotations, wave_value)
    
    # # Remove first and last wave:
    # wave_table <- wave_table[2:(nrow(wave_table)-1),]
    
    # For each wave, find AUC and add to wave_table dataframe:
    integral_indv <- array(0, length(wave_table$wave_type))
    for (i in 1:length(wave_table$wave_type)) {
      y <- signal[wave_table$wave_on[i]:wave_table$wave_off[i]] - c(isoelec)
      x <- 1:length(y)
      integral_indv[i] <- AUC(x = x, y = y)
    }
    wave_table$AUC <- integral_indv
  }
  return(wave_table)
  
  # simple: find variance over the 10 sec interval, average across all leads?
  # for loop over all 12 leads?
  
  # Could coordinate across leads using R values. Between 0 to R1, R1 to R2...
  # R_last to end
  # It seems reasonable that R peak occurs after pwave for any combo of leads
  # Error message if there are multiple p waves
}
# P Amplitude / Duration --------------------------------------------------
wave_character <- function(signal, annotations, wave_value, Hz = 500) {
  # Uses make_wave_table. For each wave type specified (p/QRS/t), find onset,
  # offset, duration and amplitude of the wave. Output in dataframe format
  
  wave_table <- make_wave_table(annotations, wave_value)
  
  # If table returns NA, exit
  if (sum(is.na(wave_table) > 0)) {
    wave_table$duration <- NA
    wave_table$amplitude <- NA
    return(wave_table)
    exit('Sample had no waves of specified type')
  }
  
  
  wave_ind <- which(annotations != 0)
  # Remove first wave, if no other types of waves precede it:
  if (wave_ind[1] == wave_table$wave_on[1]) {
    if (nrow(wave_table) > 1) {
      wave_table <- wave_table[-1, ]
    } else {
      # In addition to above condition, if there is only 1 wave, exit function:
      wave_table$duration <- NA
      wave_table$amplitude <- NA
      return(wave_table)
      exit('Only one wave found at terminus of annotations.')
    }
  }
  
  # Remove last wave, if no other types of waves follow it:
  if (wave_ind[length(wave_ind)] == wave_table$wave_off[nrow(wave_table)]) {
    if (nrow(wave_table) > 1) {
      wave_table <- wave_table[-nrow(wave_table), ]
    } else {
      # In addition to above condition, if there is only 1 wave, exit function:
      wave_table$duration <- NA
      wave_table$amplitude <- NA
      return(wave_table)
      exit('Only one wave found at terminus of annotations.')
    }
  }
  
  duration <- (wave_table$wave_off - wave_table$wave_on) / Hz
  
  amplitude <- array(NA, nrow(wave_table))
  for (i in 1:nrow(wave_table)) {
  max <- max(signal[wave_table$wave_on[i]:wave_table$wave_off[i]])
  min <- min(signal[wave_table$wave_on[i]:wave_table$wave_off[i]])
  
  amplitude[i] <- max - min
  
  }
  
  wave_table$duration <- duration
  wave_table$amplitude <- amplitude
  
  return(wave_table)
}
  





# GEH ---------------------------------------------------------------------
geh <- function(XYZ_M, origin_point, GEH_Ronset, GEH_Rpeak, GEH_Roffset, GEH_Tpeak, GEH_Toffset, fs = 500, amp_r = 1) {
  # Translated from MATLAB using Dr. Tereshchenko's code:
  # https://github.com/Tereshchenkolab/Global-Electrical-Heterogeneity/blob/master/GEH_analysis_git.m
  # **Difficult to check for accuracy. Unsure if results seem correct**
  
  fs_d <- fs
  # amp_r: amplitude resolution in microVolts...
  
  XYZ_median      <- XYZ_M * amp_r # median
  R_VM            <- GEH_Rpeak # r peak ~300
  q_points_VM     <- GEH_Ronset # QRS onset.
  s_points_VM     <- GEH_Roffset # QRS offset
  tp_points_VM    <- GEH_Tpeak; # t peak?
  te_points_VM    <- GEH_Toffset; # t offset?
  OriginPoint_idx <- origin_point # origin

  if (sum(is.na(c(R_VM,q_points_VM,s_points_VM,tp_points_VM,te_points_VM,OriginPoint_idx))) > 0) {
    warning('One or more components are NA. Breaking.')
    return(array(NA,20))
  }
  
  # If t_peak is further than 600 index, return NA:
  if (tp_points_VM > 600) {
    warning('T peak is of value ',tp_points_VM,'. Returning NA.')
    return(NA)
  }
  
  if (te_points_VM > 600) {
    te_points_VM <- 600
    warning('T offset is of value ',te_points_VM,'. Reset to 600.')
  }
  
  # Var Calculation ---------------------------------------------------------
  
  # % calculation of Vector Magnitudes (Euclidian norm) using the origin point
  
  VecMag <- array(0,length(XYZ_median[,1]))
  
  for (ii in 1:length(XYZ_median[,1])) {
    VecMag[ii] <- norm(c(XYZ_median[OriginPoint_idx,1] - XYZ_median[ii,1], 
                         XYZ_median[OriginPoint_idx,2] - XYZ_median[ii,2], 
                         XYZ_median[OriginPoint_idx,3] - XYZ_median[ii,3]),
                       type = "2");
  }
  
  # XYZ must be an array
  
  # % find R peak in median XYZ beat
  Rx_val <-  max(XYZ_median[1:500,1])
  Ry_val <-  max(XYZ_median[1:500,2])
  Rz_val <-  max(XYZ_median[1:500,3])
  
  Rx <-  which.max(XYZ_median[1:500,1])
  Ry <-  which.max(XYZ_median[1:500,2])
  Rz <-  which.max(XYZ_median[1:500,3])
  
  # % define R peak and T peak as R axis and T axis
  Raxis = XYZ_median[R_VM,]
  Taxis = XYZ_median[tp_points_VM[1],]
  
  # % ========================== Calculate AUC on Vector Magnitude ===========================
  # AUC from QRS thru T wave end
  library(DescTools)
  x <- 1:length(q_points_VM[1]:te_points_VM[1])
  
  spac_incr <-  1000/fs #% spacing increment for trapz calculation
  AUC_VM_QT <- 0
  AUC_VM_QT <- AUC(x = x, y = abs(VecMag[q_points_VM[1]:te_points_VM[1]])) * spac_incr
  
  
  # % ======================= GEH Variable Calculation =============================
  # %  origin point values
  CP <- XYZ_median[OriginPoint_idx,]
  # % Y axis vector
  Ynew <- c(0,1,0)
  
  
  
  # % QRS and T integration: for Wilson SVG calculation
  # AUC from Q to T end on each lead
  x <- 1:length(q_points_VM[1]:te_points_VM[1])
  SumVGx <- AUC(x = x, y = XYZ_median[(q_points_VM[1] : te_points_VM[1]),1])*spac_incr
  SumVGy <- AUC(x = x, y = XYZ_median[(q_points_VM[1] : te_points_VM[1]),2])*spac_incr
  SumVGz <- AUC(x = x, y = XYZ_median[(q_points_VM[1] : te_points_VM[1]),3])*spac_incr
  
  
  # % QRS and T integration for area vectors
  # AUC of QRS, S to T_end for each spatial lead
  
  qs <- 1:length(q_points_VM[1]:s_points_VM[1])
  st <- 1:length(s_points_VM[1]:te_points_VM[1])
  
  meanVxQ <- AUC(x=qs, y=XYZ_median[q_points_VM[1]:s_points_VM[1],1])*spac_incr
  meanVxT <- AUC(x=st, y=XYZ_median[s_points_VM[1]:te_points_VM[1],1])*spac_incr
  
  meanVyQ <- AUC(x=qs, y=XYZ_median[q_points_VM[1]:s_points_VM[1],2])*spac_incr
  meanVyT <- AUC(x=st, y=XYZ_median[s_points_VM[1]:te_points_VM[1],2])*spac_incr
  
  meanVzQ <- AUC(x=qs, y=XYZ_median[q_points_VM[1]:s_points_VM[1],3])*spac_incr
  meanVzT <- AUC(x=st, y=XYZ_median[s_points_VM[1]:te_points_VM[1],3])*spac_incr
  
  
  
  # % QRS area and T area vectors based on integrals
  # Find average AUC for QRS, S to t_off
  MEAN_QRSO <- c(meanVxQ, meanVyQ, meanVzQ)
  MEAN_TO <- c(meanVxT, meanVyT, meanVzT)
  
  
  # %% QT interval in sec
  timeM <- ((1:length(VecMag))/fs)*1000
  QT_interval <- timeM[te_points_VM[1]] - timeM[q_points_VM[1]]
  
  # % peak vectors QRS and T amplitude
  # Raxis: peak QRS values. Taxis: peak T values
  QRS_amp <- sqrt(Raxis[1]^2 + Raxis[2]^2 + Raxis[3]^2)
  T_amp <- sqrt(Taxis[1]^2 + Taxis[2]^2 + Taxis[3]^2)
  
  # % peak SVG vector and mean SVG vector calculation as vector sum of QRS and T vectors 
  
  SVG_axis <- rowSums(array(c(Taxis,Raxis),c(3,2)))
  SVG_MO <- rowSums(array(c(MEAN_TO,MEAN_QRSO),c(3,2)))
  
  # % Origin Point-P, Q-S and S-T vector calculation
  qs3 = XYZ_median[q_points_VM[1]:s_points_VM[1],3]
  qs1 = XYZ_median[q_points_VM[1]:s_points_VM[1],1]
  qs2 = XYZ_median[q_points_VM[1]:s_points_VM[1],2]
  st3 = XYZ_median[s_points_VM[1]:te_points_VM[1],3]
  st1 = XYZ_median[s_points_VM[1]:te_points_VM[1],1]
  st2 = XYZ_median[s_points_VM[1]:te_points_VM[1],2]
  
  
  
  
  # Angle Calcs -------------------------------------------------------------
  
  library(DescTools)
  
  # % peak QRS-T angle
  QRSTang <- round(RadToDeg(acos( Dot(Raxis,Taxis) / (sqrt(Raxis[1]^2 + Raxis[2]^2 + Raxis[3]^2) * sqrt(Taxis[1]^2 + Taxis[2]^2 + Taxis[3]^2)))),1)
  
  # % mean QRS-T angle
  QRSTang_M <- round(RadToDeg(acos(Dot(MEAN_QRSO,MEAN_TO) / (sqrt(MEAN_QRSO[1]^2 + MEAN_QRSO[2]^2 + MEAN_QRSO[3]^2) * sqrt(MEAN_TO[1]^2 + MEAN_TO[2]^2 + MEAN_TO[3]^2)))),1)
  
  # % Azimuth of QRS: peak, area 
  AZ_OQ <- round( (RadToDeg(acos(Raxis[1]/sqrt(Raxis[1]^2 + Raxis[3]^2)))) * (((Raxis[3]<0)*-1) + ((Raxis[3]>0)*1)) ,1)
  AZ_OQ_M <- round( (RadToDeg(acos( MEAN_QRSO[1]/sqrt(MEAN_QRSO[1]^2 + MEAN_QRSO[3]^2)))) * (((MEAN_QRSO[3]<0)*-1) + ((MEAN_QRSO[3]>0)*1)) ,1)
  
  # % Azimuth of T: peak, area 
  AZ_OT <- round( (RadToDeg(acos(Taxis[1]/sqrt(Taxis[1]^2 + Taxis[3]^2)))) * (((Taxis[3]<0)*-1) + ((Taxis[3]>0)*1)),1)
  AZ_OT_M <- round( (RadToDeg(acos(MEAN_TO[1]/sqrt(MEAN_TO[1]^2 + MEAN_TO[3]^2)))) * (((MEAN_TO[3]<0)*-1) + ((MEAN_TO[3]>0)*1)) ,1)
  
  # % AzCimuth of SVG: peak, area 
  AZ_SVG <- round( (RadToDeg(acos(SVG_axis[1]/sqrt(SVG_axis[1]^2 + SVG_axis[3]^2)))) * (((SVG_axis[3]<0)*-1) + ((SVG_axis[3]>0)*1)),1)
  AZ_SVG_M <- round( (RadToDeg(acos(SVG_MO[1]/sqrt(SVG_MO[1]^2 + SVG_MO[3]^2)))) * (((SVG_MO[3]<0)*-1)+((SVG_MO[3]>0)*1)),1)
  
  # % Elevation of QRS: peak, area 
  EL_OQ <- round(RadToDeg(acos(Dot(Raxis,Ynew)/(sqrt(Raxis[1]^2+Raxis[2]^2+Raxis[3]^2)*sqrt(Ynew[1]^2+Ynew[2]^2+Ynew[3]^2)))),1)
  EL_OQ_M <- round(RadToDeg(acos(Dot(MEAN_QRSO,Ynew)/(sqrt(MEAN_QRSO[1]^2 + MEAN_QRSO[2]^2 + MEAN_QRSO[3]^2)*sqrt(Ynew[1]^2+Ynew[2]^2+Ynew[3]^2)))),1)
  
  # % Elevation of T: peak, area 
  EL_OT <- round(RadToDeg(acos(Dot(Taxis,Ynew)/(sqrt(Taxis[1]^2+Taxis[2]^2+Taxis[3]^2)*sqrt(Ynew[1]^2+Ynew[2]^2+Ynew[3]^2)))),1)
  EL_OT_M <- round(RadToDeg(acos(Dot(MEAN_TO,Ynew)/(sqrt(MEAN_TO[1]^2+MEAN_TO[2]^2+MEAN_TO[3]^2)*sqrt(Ynew[1]^2+Ynew[2]^2+Ynew[3]^2)))),1)
  
  # % Elevation of SVG: peak, area 
  EL_SVG <- round(RadToDeg(acos(Dot(SVG_axis,Ynew)/(sqrt(SVG_axis[1]^2+SVG_axis[2]^2+SVG_axis[3]^2)*sqrt(Ynew[1]^2+Ynew[2]^2+Ynew[3]^2)))),1)
  EL_SVG_M <- round(RadToDeg(acos(Dot(SVG_MO,Ynew)/(sqrt(SVG_MO[1]^2+SVG_MO[2]^2+SVG_MO[3]^2)*sqrt(Ynew[1]^2+Ynew[2]^2+Ynew[3]^2)))),1)
  
  
  # % ========================== Magnitudes Calculation ============================
  # % Magnitude of QRS: peak, area
  QRS_Mag <- round( sqrt(Raxis[1]^2+Raxis[2]^2+Raxis[3]^2) )
  QRS_Mag_M <- round( sqrt(MEAN_QRSO[1]^2+MEAN_QRSO[2]^2+MEAN_QRSO[3]^2) )
  
  # % Magnitude of T: peak, area
  T_Mag <- round( sqrt(Taxis[1]^2+Taxis[2]^2+Taxis[3]^2) )
  T_Mag_M <- round( sqrt(MEAN_TO[1]^2+MEAN_TO[2]^2+MEAN_TO[3]^2) )
  
  # % Magnitude of SVG: peak
  SVG_Mag <- round( sqrt(SVG_axis[1]^2+SVG_axis[2]^2+SVG_axis[3]^2) )
  
  # % Magnitude of WVG: Wilson's Ventricular Gradient
  
  WVG <- round( sqrt((SumVGx^2) + (SumVGy^2) + (SumVGz^2)) )
  
  
  # GEH_df = data.frame(QRSTang = QRSTang, QRSTang_M = QRSTang_M,
  #            AZ_OQ = AZ_OQ, AZ_OQ_M = AZ_OQ_M,
  #            AZ_OT = AZ_OT, AZ_OT_M = AZ_OT_M,
  #            AZ_SVG = AZ_SVG, AZ_SVG_M = AZ_SVG_M,
  #            EL_OQ = EL_OQ, EL_OQ_M = EL_OQ_M,
  #            EL_OT = EL_OT, EL_OT_M = EL_OT_M,
  #            EL_SVG = EL_SVG, EL_SVG_M = EL_SVG_M,
  #            QRS_Mag = QRS_Mag, QRS_Mag_M = QRS_Mag_M,
  #            T_Mag = T_Mag, T_Mag_M = T_Mag_M, 
  #            SVG_Mag = SVG_Mag, WVG = WVG)
  
  GEH_vector <- c(QRSTang, QRSTang_M, AZ_OQ, AZ_OQ_M, AZ_OT, AZ_OT_M,
                  AZ_SVG, AZ_SVG_M, EL_OQ, EL_OQ_M, EL_OT, EL_OT_M, 
                  EL_SVG, EL_SVG_M, QRS_Mag, QRS_Mag_M, 
                  T_Mag, T_Mag_M, SVG_Mag, WVG)
  return(GEH_vector)
  
}

# GEH QRS --------------------------------------------------------------------
find_geh_QRS_intervals <- function(leads12, ann12, skip_leads) {
  # Find average QR, RS intervals across all leads 
  # Can specify specific leads rather than all 12
  # Input prep for GEH
  
  if (length(dim(ann12)) < 2) {
    ann12 <- array(rep(ann12, 12), c(5000,12))
  }
  
  QR_interval <- array(NA,dim(ann12)[[2]])
  RS_interval <- array(NA,dim(ann12)[[2]])
  
  for (i in 1:dim(ann12)[[2]]) {
    # If there are no QRS, skip to next sample
    Rwaves <- make_wave_table(ann12[,i], wave_value = 2)
    if (nrow(Rwaves) == 1) {
      if (is.na(Rwaves$wave_type)) {
        QR_interval[i] <- NA
        RS_interval[i] <- NA
        next
      }
    }
    
    Rwaves$peak <- unlist(peak_isolation(leads12[, i], ann12[, i], wave_value = 2))
    
    # If there are no QRS, switch to next in loop:
    if (sum((Rwaves$wave_type == 'N'), na.rm = TRUE) < 1) {
      QR_interval[i] <- list(NA)
      RS_interval[i] <- list(NA)
      next
    }
    
    # if (Rwaves$wave_on[1] == 1) { # commenting out this code, as annotators should be good enough at signal termini with QRSs
    #   start <- 2
    # } else{
    #   start <- 1
    # }
    start <- 1
    
    QR_indv <- c()
    RS_indv <- c()
    
    for (j in start:length(Rwaves[, 1])) {
      QR_indv <- c(QR_indv, Rwaves$peak[j] - Rwaves$wave_on[j])
      RS_indv <- c(RS_indv, Rwaves$wave_off[j] - Rwaves$peak[j])
    }
    
    QR_interval[i] <- list(QR_indv)
    RS_interval[i] <- list(RS_indv)
  }
  
  if (exists('skip_leads')) {
    QR_interval[skip_leads] <- lapply(QR_interval[skip_leads], function(x) rep(NA))
    RS_interval[skip_leads] <- lapply(RS_interval[skip_leads], function(x) rep(NA))
  }
  
  # QR_median <- round(median(unlist(QR_interval), na.rm = TRUE))
  # RS_median <- round(median(unlist(RS_interval), na.rm = TRUE))
  
  # Weight each lead evenly to account for any rogue annotations
  QR_median <- median(sapply(QR_interval, median, na.rm = TRUE),na.rm=TRUE) 
  RS_median <- median(sapply(RS_interval, median, na.rm = TRUE),na.rm=TRUE) 
  
  return(data.frame(QR = QR_median, RS = RS_median))
  
}


# GEH Toff ----------------------------------------------------------------
find_geh_RToff_interval <- function(leads12, ann12, skip_leads) {
  # Find average R - T-offset interval across all 12 leads. 
  # Can specify specific leads rather than all 12
  # GEH input prep
  
  # If there is only one lead of annotion
  if (length(dim(ann12)) < 2) {
    ann12 <- array(rep(ann12, 12), c(5000,12))
  }
  
  Rtoff_interval <- array(0,dim(ann12)[[2]])
  # single lead
  for (i in 1:dim(ann12)[[2]]) {
    # If there are no QRS, skip to next sample
    Rwaves <- make_wave_table(ann12[,i], wave_value = 2)
    if (nrow(Rwaves) == 1) {
      if (is.na(Rwaves$wave_type)) {
        Rtoff_interval[i] <- NA
        next
      }
    }
    Rwaves$peak <- unlist(peak_isolation(leads12[,i], ann12[,i], wave_value = 2))
    Twaves <- make_wave_table(ann12[,i], wave_value = 3)
    Twaves$peak <- unlist(peak_isolation(leads12[,i], ann12[,i], wave_value = 3))
    
    combined <- rbind(Rwaves, Twaves)
    combined <- combined[order(combined$wave_on),]
    
    # skip first t-wave if there's no preceding R-wave:
    
    t_ind <- which(combined$wave_type == "t")
    
    # If there are no QRS, skip:
    if (sum((combined$wave_type == 'N'), na.rm = TRUE) < 1) {
      Rtoff_interval[i] <- list(NA)
      next
    }
    
    # If there are no t-waves, skip:
    if (length(t_ind) == 0) {
      Rtoff_interval[i] = NA
      next
    }
    
    # If there are no QRS prior to the first t-wave, skip
    if (t_ind[1] == 1) {
      start <- 2
      # If there are no other t-waves with a QRS prior, skip to next loop:
      if (length(t_ind) == 1) {
        Rtoff_interval[i] = list(NA)
        next
      }
      
    } else{
      start <- 1
    }
    
    RToff_indv <- c()
    for (j in start:length(t_ind)) {
      if (combined$wave_type[t_ind[j] - 1] == "N") {
        interval <- combined$wave_off[t_ind[j]] - combined$peak[t_ind[j] - 1]
        RToff_indv <- c(RToff_indv, interval)
      }
    }
    
    Rtoff_interval[i] <- list(RToff_indv)
    
  }
  
  # If user supplies leads to skip over, replace those lead values with NA. Not currently used
  # if (exists('skip_leads')) {
  #   Rtoff_interval[skip_leads] <- lapply(Rtoff_interval[skip_leads], function(x) rep(NA))
  # }
  
  # Weight each lead evenly to account for any rogue annotations
  overall_median <- median(sapply(Rtoff_interval, median, na.rm = TRUE),na.rm=TRUE) # weight each lead evenly to account for any rogue annotations
  return(round(overall_median))
}

# GEH Tpeak --------------------------------------------------------------
find_geh_RTpeak_interval <- function(leads12, ann12, skip_leads) {
  # Find R - T-peak
  # handles one sample
  # skip_leads: when finding the RT interval, skip the leads given. Vector format
  
  # If there is only one lead of annotion
  if (length(dim(ann12)) < 2) {
    ann12 <- array(rep(ann12, 12), c(5000,12))
  }
  
  # For each lead
  Rt_interval <- array(NA,dim(ann12)[[2]])
  for (i in 1:dim(ann12)[[2]]) {
    # Find R and T peaks
    
    # If there are no QRS, skip to next sample
    Rwaves <- make_wave_table(ann12[,i], wave_value = 2)
    if (nrow(Rwaves) == 1) {
      if (is.na(Rwaves$wave_type)) {
        Rt_interval[i] <- NA
        next
      }
    }
    Rwaves$index <- unlist(peak_isolation(leads12[,i], ann12[,i], wave_value = 2))
    Twaves <- make_wave_table(ann12[,i], wave_value = 3)
    Twaves$index <- unlist(peak_isolation(leads12[,i], ann12[,i], wave_value = 3))
    
    combined <- rbind(Rwaves, Twaves)
    combined <- combined[order(combined$wave_on), ]
    
    t_ind <- which(combined$wave_type == "t")
    
    # If there are no t-waves, skip:
    if (length(t_ind) == 0) {
      Rt_interval[i] = NA
      next
    }
    
    # skip first t-wave if there's no preceding QRS-wave:
    if (t_ind[1] == 1) {
      start <- 2
      # If there are no other t-waves with a QRS prior, skip to next loop:
      if (length(t_ind) == 1) {
        Rt_interval[i] = NA
        next
      }
    } else{
      start <- 1
    }
    
    
    Rt_indv <- c()
    for (j in start:length(t_ind)) {
      if (combined$wave_type[t_ind[j] - 1] == "N") {
        interval <- combined$index[t_ind[j]] - combined$index[t_ind[j] - 1]
        Rt_indv <- c(Rt_indv, interval)
      }
    }
    Rt_interval[i] <- list(Rt_indv)
    
  }
  
  if (exists('skip_leads')) {
    Rt_interval[skip_leads] <- lapply(Rt_interval[skip_leads], function(x) rep(NA))
  }
  
  # Weight each lead evenly to account for any rogue annotations
  Rt_interval_median <- median(sapply(Rt_interval, median, na.rm = TRUE),na.rm=TRUE) 
  
  return(round(Rt_interval_median))
  
}

# GEH RPeak ----------------------------------------------------------------
find_geh_Rpeak <- function(XYZ_M, origin) {
  # Find R peak value across all 3D spatial leads (XYZ)
  # Input prep for GEH, added to the above GEH functions to find Q, S, T_peak, 
  # T_off indices
  
  # Finds max in the middle third of the 600-long recording. Windowing method
  # is done to prevent error in tachycardic patients, which have multiple
  # QRSs in one recording
  
  if (is.na(origin)) {
    return(NA)
    stop('Origin point is NA.')
  }
  
  Rpeaks <- array(0,3)
  for (i in 1:3){
    Rpeaks[i] <- which.max(abs(XYZ_M[201:400,i] - XYZ_M[origin,i])) + 200
  }
  median <- round(median(Rpeaks))
  return(median)
}

# Median Beat Function ----------------------------------------------------
find_median_beat <- function(XYZ, Rpeaks, fs = 500) {
  # Inputs: 10 second XYZ, Rpeaks for a single lead
  # Output: one single average beat for all 3 spatial leads, centered on the 
  # Rpeak
  
  # Definitions:
  # Rpeaks: vector of Rpeaks on a reference single lead of a 12 lead ECG
  # Rx: vector of Rpeaks from the X dimension of a 3D spatial lead
  # I_dv: vector of max slopes of the Rpeaks, from the X dimension of a 3D spatial lead 
  
  XYZ_O <- XYZ
  
  if (class(Rpeaks) == "list") {
    Rpeaks <- unlist(Rpeaks)
  }
  Rpeaks <- na.omit(Rpeaks)
  
  if (length(Rpeaks) == 0) {
    warning('No Rpeaks. Returning NA')
    return(NA)
  }
  
  # Find Rpeaks within 3D spatial lead, using original Rpeaks as reference - double checks peak values
    # Rpeaks: Rpeak on a reference single lead of a 12 lead ECG
    # Rx: Rpeak from the X dimension of a 3D spatial lead
  window <- 0.1 # seconds
  window_ind <- round(window*fs)
  
  Rx <- array(0,length(Rpeaks))
  for (i in 1:length(Rpeaks)) {
    
    # Ranges for search
    start <- (Rpeaks[i] - window_ind)
    end <- (Rpeaks[i] + window_ind)
    
    start[start<1] <- 1
    end[end>length(XYZ_O[,1])] <- length(XYZ_O[,1])

    # Find max value of the first dimension of the 3D spatial lead
    Rx[i] <- which.max(abs(XYZ_O[start:end,1])) + start - 1
  }
  
  
  Total_samples <- dim(XYZ_O)[[1]] # verify dimension used for length
  Beat_length <- as.integer(round(fs*1.2)) # added 20 index safety factor
  
  # If 3D spatial Rpeaks are too close to beginning, skip it:
  if (sum(Rx < ((Beat_length/2) + 20)) > 0) {
    too_early_peaks <- which(Rx < ((Beat_length/2) + 20))
    Rx <- Rx[-too_early_peaks]
  }

  # If 3D spatial peaks are too close to end, skip it:
  if (sum(Rx > (Total_samples - (Beat_length/2) + 20)) > 0) {
    too_late_peaks <- which(Rx > (Total_samples - (Beat_length/2) + 20))
    Rx <- Rx[-too_late_peaks]
  }
  
  # If too few peaks in the X dimmension (Rx), exit:
  if (length(Rx) < 3) {
    warning('Too few usable RPeaks. Returning NA')
    return(NA)
  }
  
  I_dv=array(0,length(Rx)); # I_dv: max slope change centered in each beat
  
  # For the included Rpeaks, find maximum slope change (I_dv) within the window:
  
  #   Set window size based on HR (12-lead Rpeaks) to search for max slope:
  if (length(Rpeaks) < 20) {
    window <- 100
  } else {
    window <- 75
  }
  
  # If search range for the last max slope is out of range, don't use it
  if (Rx[length(Rx)] + window > 5000) {
    Rx = Rx[-length(Rx)]
    I_dv=array(0,(length(Rx)-1))
  }
  
  # Find I_dv, the max slope centered around the Rpeak, in the X dimension of the 3D spatial lead
  for (ii in 1:length(Rx)) {
    I_dv_temp <- which.max(abs(diff(XYZ_O[ (Rx[ii] - (window - 1)) : (Rx[ii] + window), 1])))
    # [value, index]
    # [~,I_dv_temp]=max(abs(diff(XYZ_O(     Rx(ii)-199:Rx(ii)+200,1    ))));
    I_dv[ii]=Rx[ii] + I_dv_temp - (window - 1)
  }
  
  # If centering point indices are less than half the XYZ_M size, eliminate it:
  if (sum(I_dv < 300) > 0 || sum(I_dv > 4700) > 0) {
    I_dv <- I_dv[-which(I_dv < 300 | I_dv > 4700)]
  }
  
  if (length(I_dv) < 3) {
    warning('Too few usable centering points after filtering. Returning NA')
    return(NA)
  }
  
  # FOr each 3D spatial lead, create median beats.
    # ie averaging each beat, centering on the maximum dV/dt in the X dimmension (ie I_dv)
    # Beat_length: the length of the median beat
    # total_samples: the total number of Rpeaks to be averaged
  
  total_samples <- length(I_dv)
  
  # First create matrices for each dimension, one column for each Rpeak:
  Beats_x_T <- array(0, c(Beat_length, total_samples)) # vector for centered beat ranges
  Beats_y_T <- array(0, c(Beat_length, total_samples))
  Beats_z_T <- array(0, c(Beat_length, total_samples))
  
  for (ii in 1:length(I_dv)) {
    Beats_x_T[,ii] <- XYZ_O[(I_dv[ii] - (Beat_length/2 - 1)) : (I_dv[ii] + (Beat_length/2)), 1]
    Beats_y_T[,ii] <- XYZ_O[(I_dv[ii] - (Beat_length/2 - 1)) : (I_dv[ii] + (Beat_length/2)), 2]
    Beats_z_T[,ii] <- XYZ_O[(I_dv[ii] - (Beat_length/2 - 1)) : (I_dv[ii] + (Beat_length/2)), 3]
  }
  
  # Then average these matrices to create an average vector for each lead (X, Y and Z).
    # The individual lead vectors are combined to create a single matrix combining all 3 leads
  XYZ_M_T <- array(0, c(Beat_length, 3))
  for (i in 1 : Beat_length) {
    XYZ_M_T[i,1] <- median(Beats_x_T[i,])
    XYZ_M_T[i,2] <- median(Beats_y_T[i,])
    XYZ_M_T[i,3] <- median(Beats_z_T[i,]) 
  }
  
  # Calculate the Vector Magnitude (Euclidean norm) - currently unused 
  VecMag_T <- sqrt(rowSums(XYZ_M_T^2))
  
  timeM=((1:length(VecMag_T))/fs)*1000 # - currently unused 
  
  return(XYZ_M_T)
}

# Origin Function ---------------------------------------------------------
find_origin <- function(XYZ_M, Rpeaks, fs = 500) {
  # Finds the isoelectric point across all 3D spatial leads
  # Input: median XYZ, Rpeaks
  # Translated from: https://github.com/Tereshchenkolab/Origin/blob/master/origin_point_github.m#L41C1-L140C8
  
  # If XYZ_M = NA (ie due to error), break:
  if (sum(is.na(XYZ_M)) == length(XYZ_M)) {
    return(NA)
    break
  }
  
  if (class(Rpeaks) == "list") {
    Rpeaks <- unlist(Rpeaks)
  }
  
  Rpeaks <- na.omit(Rpeaks)
  
  if (length(Rpeaks) < 3) {
    return(NA)
    break
  }
  
  Rx <- Rpeaks
  
  # Rx is the peak of the median beat** --> ~300

  
  # %% Calculate various windows analysis
  w <- floor( (mean(diff(Rx))*.8) - (372.69*(fs/1000)))
  
  # Add in minimum size for window:
  if (w < 2) {
    w = 2
  }
  
  # %% Is the RR interval < 600ms (0.6*fs)
  if (mean(diff(Rx)) < (600*(fs/1000)) ) { # (A.Label_beat~='S') ??? 
    cent <- (600-260)*(fs/1000)
  } else {
    cent <- (600-320)*(fs/1000)
  }
  
  
  w1_start <- round(cent-w)
  w1_end <- round(cent+w)
  
  w2_start <- round(cent-(160*(fs/1000)))
  w2_end <- round(cent+(160*(fs/1000)))
  
  if (w1_start < w2_start) {
    w1_start <- w2_start
  }
  
  if (w1_end > w2_end) {
    w1_end <- w2_end
  }
  
  # %% cut xyz median beat based on the window
  xyz_w1 = XYZ_M[w1_start:w1_end,]
  xyz_w2 = XYZ_M[w2_start:w2_end,]
  
  # if (is.vector(xyz_w1)) {
  #   xyz_w1 <- xyz_w2
  #   warning('xyz_w1 is a vector... investigate further')
  # }
  
  # %% calculate origin point with clustering algorithm ori_clustering
  output <- ori_clustering(xyz_median = xyz_w1)
  cluster_pt <- unlist(output[1])
  
  if (is.na(cluster_pt)) { # is "NA" or is empty?
    output <- ori_clustering(xyz_w2)
    cluster_pt <- unlist(output[1])
    Ori_pt <- cluster_pt + w2_start
  } else {
    Ori_pt <- cluster_pt + w1_start
  }
  
  # % calculate the first origin point and the modified vector magnitude
  VecMag_M1 <- array(0, length(XYZ_M[,1]))
  
  for (ii in 1:length(XYZ_M[, 1])) {
    VecMag_M1[ii] <- norm(c(XYZ_M[Ori_pt[1], 1] - XYZ_M[ii, 1], XYZ_M[Ori_pt[1], 2] - XYZ_M[ii, 2], XYZ_M[Ori_pt[1], 3] - XYZ_M[ii, 3]), type = "2")
  }
  
  
  # % calculate the sum of the absolute gradient on XYZ within w1
  output <- abs_grad_fun(xyz_median = xyz_w1)
  abs_grad_ori <- unlist(output[1])
  
  # % calcaulate the second origin point
  Ori_pt2 <- abs_grad_ori + w1_start
  
  # % calculate the second origin point and the modified vector magnitude
  VecMag_M2 <- array(NA, length(XYZ_M[,1]))
  
  
  for (ii in 1:length(XYZ_M[, 1])) {
    VecMag_M2[ii] <- norm(c(XYZ_M[Ori_pt2[1], 1] - XYZ_M[ii, 1], XYZ_M[Ori_pt2[1], 2] - XYZ_M[ii, 2], XYZ_M[Ori_pt2[1], 3] - XYZ_M[ii, 3]), type = "2")
  }
  
  # % calculate the area under the curve of the given angle
  library(DescTools)
  x <- 1:length(w1_start:w1_end)
  
  # Area_m1 <- AUC(x = x, y = VecMag_M1[w1_start:w1_end])
  # Area_m2 <- AUC(x = x, y = VecMag_M2[w1_start:w1_end])
  
  Area_m1 <- DescTools::AUC(x = x, y = VecMag_M1[w1_start:w1_end], method = "trapezoid")
  Area_m2 <- DescTools::AUC(x = x, y = VecMag_M2[w1_start:w1_end], method = "trapezoid")
  
  # % modify the XYZ median beat with the two origon point calculated
  XYZ_M1 <- t(array(c(XYZ_M[Ori_pt[1],]), c(3,dim(XYZ_M)[[1]]))) - XYZ_M
  XYZ_M2 <- t(array(c(XYZ_M[Ori_pt2[1],]), c(3,dim(XYZ_M)[[1]]))) - XYZ_M
  
  
  # % calculate the sum of the absolute gradient on modified XYZ median beats
  data_grad_M1 <- data_processing(XYZ_M1,'gradsum',5)
  data_grad_M2 <- data_processing(XYZ_M2,'gradsum',5)
  
  # % calcualte the average of the sum of the absolute gradient about the
  # % origin points with 3 points on either side
  Sum_M1 <- mean(data_grad_M1[ (Ori_pt-3) : (Ori_pt+3)]);
  Sum_M2 <- mean(data_grad_M2[ (Ori_pt2-3) : (Ori_pt2+3)]);
  
  # % decide which origin point to use based on the area under the curve and
  # % the average gradient about the origin points
  if ((Area_m1 < Area_m2) && (Sum_M1 < Sum_M2)) {
    dat_ori <- Ori_pt
    VecMag_ML <- VecMag_M1
  } else {
    dat_ori <- Ori_pt2
    VecMag_ML <- VecMag_M2
  }
  
  origin_point <- dat_ori
  
  return(origin_point)
}

# Origin Supporting Functions -------------------------------------------------------------
abs_grad_fun <- function(xyz_median) {
  
  if (is.vector(xyz_median)) {
    xyz_median <- array(xyz_median,c(1,length(xyz_median)))
  }
  
  if (nrow(xyz_median) < 3) {
    return(NA)
    break
  }
  
  diff_w1 <- data_processing(data = xyz_median,process = 'gradsum',process_window = 2)
  iter_num <- 0
  min_idx <- which.min(diff_w1)
  
  # Loop thru the windows until the condition is met, then stop
  for (ii in 1:(length(diff_w1)-1)) {
    if (abs(diff_w1[length(diff_w1)-ii] - diff_w1[min_idx]) < 0.1) {
      iter_num <- length(diff_w1) - ii
      break
    }
  }
  
  ori_pt_idx <- iter_num
  ori_val <- xyz_median[iter_num,]
  
  output <- list(ori_pt_idx, ori_val)
  return(output)
}

ori_clustering <- function(xyz_median) {
  if (is.vector(xyz_median)) {
    xyz_median <- array(xyz_median,c(1,length(xyz_median)))
  }
  
  processed_data  <-  data_processing(data = xyz_median,
                                      process = 'norm_mvvar',
                                      process_window = 10) # change to 10?
  
  bindwidth <- 10
  
  # MATLAB-style bin edges
  edges <- seq(min(processed_data), max(processed_data) + bindwidth, by = bindwidth)
  
  # Use cut() with left-closed, right-open intervals
  bin_label <- cut(
    processed_data,
    breaks = edges,
    include.lowest = TRUE,
    right = FALSE,   # MATLAB uses left-closed, right-open
    labels = FALSE
  )
  
  # Count elements per bin
  N <- tabulate(bin_label, nbins = length(edges) - 1)
  
  xyz_norm  <-  array(NA,c(dim(xyz_median)[[1]], 1))
  
  for (ii in 1:length(xyz_median[,1])) {
    xyz_norm[ii] <- norm(c(xyz_median[ii, 1], xyz_median[ii, 2], xyz_median[ii, 3]), type = "2")
  }
  
  idx <- which(N == max(N)) # which bin has the most points
  candidate_idx <-  which(bin_label == idx[1]) # find points contained in the bin
  idx_diff <- diff(candidate_idx) # Are these points consecutive?
  diff_idx <- c(0, which(idx_diff >= 4)) #If consecutive, length(diff_idx) = 1
  
  if (length(diff_idx) == 1) {
    ori_pt_cluster <- candidate_idx
    
    ori_pt_idx <- round(median(candidate_idx))
  } else {
    cluster_ele_count <- array(0, c(length(diff_idx), 1))
    
    for (ii in 1:length(diff_idx)) {
      if (ii == length(diff_idx)) {
        cluster_ele_count[ii] <- sum(idx_diff[ (diff_idx[ii] + 1) : length(idx_diff)])
      } else {
        cluster_ele_count[ii] <- sum(idx_diff[ (diff_idx[ii] + 1) : (diff_idx[ii + 1] - 1) ])
      }
    }
    I <- order(cluster_ele_count, decreasing = TRUE)
    
    if (I[1] == max(I)) {
      cluster_group1 <- candidate_idx[ (diff_idx[I[1]] + 1) : length(candidate_idx)]
    } else {
      cluster_group1 <- candidate_idx[ (diff_idx[I[1]] + 1) : diff_idx[I[1] + 1] - 1]
    }
    
    if (I[2] == max(I)) {
      cluster_group2 <- candidate_idx[ (diff_idx[I[2]] + 1) : length(candidate_idx)]
    } else {
      cluster_group2 <- candidate_idx[ (diff_idx[I[2]] + 1) : diff_idx[I[2] + 1] - 1]
    }
    
    cluster_avg_slope <- array(0, c(2))
    data_gradient_sum = data_processing(data = xyz_median, process = 'gradsum',process_window = 5)
    cluster_avg_slope[1] = mean(data_gradient_sum[cluster_group1]);
    cluster_avg_slope[2] = mean(data_gradient_sum[cluster_group2]);
    
    bestgroup <- which.min(cluster_avg_slope)
    
    if (bestgroup == 1) {
      ori_pt_cluster <-  cluster_group1
      ori_pt_idx <- round(median(cluster_group1))
    } else {
      ori_pt_cluster <- cluster_group2
      ori_pt_idx <-  round(median(cluster_group2))
    }
    
    
  }
  
  output <- list(ori_pt_idx,ori_pt_cluster)
  return(output)
}

data_processing <- function(data, process, process_window) {
  # This function outputs two different methods of findings areas of change (ie spikes) vs consistency within all 3 spatial vectors
  # norm_mvvar: finds the Euclidean norm (sqrt(x^2+y^2+z^2)), then the variance of a given window of the euclidean norm
    # better at finding slow drift of signal vs. constant
  # gradsum: calculates the absolute value of gradient across all 3 dimensions
    # better at finding quick spikes vs. constant
  library(zoo)
  
  if (is.vector(data)) {
    data <- array(data,c(1,length(data)))
  }
  
  if (process == "norm_mvvar") {
    
    data_norm <- matrix(data = NA, dim(data)[[1]])
    for (ii in 1:length(data[, 1])) {
      data_norm[ii] <- norm(c(data[ii, 1], data[ii, 2], data[ii, 3]), type = "2")
      
    }
    if (process_window > nrow(data)) {
      process_window <- nrow(data)
    }
    
    if (process_window > 1) {
      processed_data <- zoo::rollapply(
        data_norm,
        width = process_window,
        FUN = var,
        align = "center",
        partial = TRUE
      )
      
      end <- length(processed_data)
      for (i in 1:((process_window - 1) / 2)) {
        processed_data[i] <- var(data_norm[1:(i + (process_window - 1) / 2)])
        processed_data[end - i + 1] <- var(data_norm[(end - i + 1 - (process_window - 1) /
                                                        2):end])
      }
    } else {
      processed_data <- data_norm
    }
    
    
    # gradient:
  } else if (process == "gradsum") {
    # **need to take averages while keeping them in their columns. Then take abs, gradient, and add up
    # pre_processed_data <- abs(gradient(rollmean(data[, 1], process_window, fill = NA))) + abs(gradient(rollmean(data[, 2], process_window, fill = NA))) + abs(gradient(rollmean(data[, 3], process_window, fill = NA)))
    
    if (process_window > nrow(data)) {
      process_window <- nrow(data)
    }
    
    if (nrow(data) < 3) {
      return(NA)
      break
    }
    
    pre_processed_data_x <- zoo::rollmean(data[,1], process_window, fill = NA, align = "center")
    pre_processed_data_y <- zoo::rollmean(data[,2], process_window, fill = NA, align = "center")
    pre_processed_data_z <- zoo::rollmean(data[,3], process_window, fill = NA, align = "center")
    
    end <- length(pre_processed_data_x)
    for (i in 1: ((process_window - 1)/2)) {
      pre_processed_data_x[i] <- mean(data[1: (i + (process_window - 1)/2), 1])
      pre_processed_data_x[end - i + 1] <- mean(data[(end - i + 1 - (process_window - 1)/2) : end, 1])
      
      pre_processed_data_y[i] <- mean(data[1: (i + (process_window - 1)/2), 2])
      pre_processed_data_y[end - i + 1] <- mean(data[(end - i + 1 - (process_window - 1)/2) : end, 2])
      
      pre_processed_data_z[i] <- mean(data[1: (i + (process_window - 1)/2), 3])
      pre_processed_data_z[end - i + 1] <- mean(data[(end - i + 1 - (process_window - 1)/2) : end, 3])
    }
    
    processed_data <- abs(gradient(pre_processed_data_x)) + abs(gradient(pre_processed_data_y)) + abs(gradient(pre_processed_data_z))
  }
  
  return(processed_data)
}

gradient <- function(data) {
  grad <- array(0,length(data))
  end <- length(data)
  
  grad[1] <- data[2] - data[1]
  grad[end] <- data[end] - data[end-1]
  
  grad[2:(end-1)] <- (data[3:end] - data[1:(end-2)]) / 2
  # for (i in 2:(length(data)-1)) {
  #   grad[i] <- (data[i+1] - data[i - 1]) / 2
  # }
  return(grad)
}


# Kors function -----------------------------------------------------------
kors <- function(leads12) {
  # Transforms 12 lead into 3 perpendicular leads
korsMatrix <- matrix(c(0.38, -0.07, 0, 0, 0, 0, -0.13, 0.05, -0.01, 0.14, 0.06,0.54,  
                      -0.07, 0.93, 0, 0, 0, 0, 0.06, -0.02, -0.05, 0.06, -0.17, 0.13,  
                      0.11, -0.23, 0, 0, 0, 0, -0.43, -0.06,-0.14,-0.20,-0.11,0.31), 
                    c(3,12),byrow = TRUE)

transform_k <- leads12 %*% t(korsMatrix)

return(transform_k)
}

# Function to find average distance from Rpeak to QRS onset, offset; T peak, onset
# average across all beats, in all 12 leads
# Align to **~300 frame** / Rpeak of median beat / max value? 
# input: samples, annotations, which leads to average

# when finding peak values, will need to account for negative deflection leads (ie QRS of avf(?))
# to do so: 
# input: bandpass filtered signal of wave indicies of interest
# find isoelectric point of each sample lead / each beat(?), subtract from wave index values
# find max of abs value of isoelectric bp_signal