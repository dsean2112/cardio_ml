# HRV ---------------------------------------------------------------------
HRV_RMSSD <- function(Rpeaks, Hz = 500) {
  # Find RMSSD for 10 sec ECG
  # Input should be index values, not time values
  
  # Multi-sample input:
  if (class(Rpeaks) == "list") {
    RMSSD <- array(0,length(Rpeaks))
    for (i in (1:length(Rpeaks))) {
      Rpeaks_indv <- unlist(Rpeaks[i])
      Rpeaks_indv <- Rpeaks_indv / Hz
      interval <- Rpeaks_indv[-1] - Rpeaks_indv[1: (length(Rpeaks_indv) - 1)]
      RMSSD[i] <- sqrt( 1/(length(interval)- 1) * sum((interval[1: (length(interval) - 1)] - interval[-1])^2))
    }
    
  # Single-sample input:  
  } else if (is.vector(Rpeaks = TRUE)) {
    Rpeaks = Rpeaks / Hz
    interval <- Rpeaks[-1] - Rpeaks[1:(length(Rpeaks) - 1)]
    RMSSD <- sqrt(1 / (length(interval) - 1) * sum((interval[1:(length(interval) - 1)] - interval[-1]) ^2))
  }
  
  return(RMSSD)
 
}




# P Amplitude / Duration --------------------------------------------------
wave_character <- function(signal, annotations, wave_value, Hz = 500) {
  # Uses make_wave_table. For each wave type specified (p/QRS/t), find onset,
  # offset, duration and amplitude of the wave. Output in dataframe format
  
  wave_table <- make_wave_table(annotations, wave_value)

  duration <- (wave_table$wave_off - wave_table$wave_on) / Hz
  
  amplitude <- array(0, dim(wave_table)[[1]])
  for (i in 1:dim(wave_table)[[1]]) {
  max <- max(signal[wave_table$sample[i], wave_table$wave_on[i]:wave_table$wave_off[i]])
  min <- min(signal[wave_table$sample[i], wave_table$wave_on[i]:wave_table$wave_off[i]])
  
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
  
  # R_VM            <- 298 # r peak ~300
  # q_points_VM     <- 287 # QRS onset.
  # s_points_VM     <- 308 # QRS offset
  # tp_points_VM    = 460; # t peak?
  # te_points_VM    = 507; # t offset?
  
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
  library(DescTools)
  x <- 1:length(q_points_VM[1]:te_points_VM[1])
  
  spac_incr <-  1000/fs #% spacing increment for trapz calculation
  AUC_VM_QT <- 0
  AUC_VM_QT <- AUC(x = x, y = abs(VecMag[q_points_VM[1]:te_points_VM[1]])) * spac_incr
  
  
  # % ======================= GEH Variable Calculation =============================
  # %  origin point
  CP <- XYZ_median[OriginPoint_idx,]
  # % Y axis vector
  Ynew <- c(0,1,0)
  
  
  
  # % QRS and T integration: for Wilson SVG calculation
  x <- 1:length(q_points_VM[1]:te_points_VM[1])
  SumVGx <- AUC(x = x, y = XYZ_median[(q_points_VM[1] : te_points_VM[1]),1])*spac_incr
  SumVGy <- AUC(x = x, y = XYZ_median[(q_points_VM[1] : te_points_VM[1]),2])*spac_incr
  SumVGz <- AUC(x = x, y = XYZ_median[(q_points_VM[1] : te_points_VM[1]),3])*spac_incr
  
  
  # % QRS and T integration for area vectors
  
  qs <- 1:length(q_points_VM[1]:s_points_VM[1])
  st <- 1:length(s_points_VM[1]:te_points_VM[1])
  
  meanVxQ <- AUC(x=qs, y=XYZ_median[q_points_VM[1]:s_points_VM[1],1])*spac_incr
  meanVxT <- AUC(x=st, y=XYZ_median[s_points_VM[1]:te_points_VM[1],1])*spac_incr
  meanVyQ <- AUC(x=qs, y=XYZ_median[q_points_VM[1]:s_points_VM[1],2])*spac_incr
  meanVyT <- AUC(x=st, y=XYZ_median[s_points_VM[1]:te_points_VM[1],2])*spac_incr
  meanVzQ <- AUC(x=qs, y=XYZ_median[q_points_VM[1]:s_points_VM[1],3])*spac_incr
  meanVzT <- AUC(x=st, y=XYZ_median[s_points_VM[1]:te_points_VM[1],3])*spac_incr
  
  
  
  # % QRS area and T area vectors based on integrals
  MEAN_QRSO <- c(meanVxQ, meanVyQ, meanVzQ)
  MEAN_TO <- c(meanVxT, meanVyT, meanVzT)
  
  
  # %% QT interval
  timeM <- ((1:length(VecMag))/fs)*1000
  QT_interval <- timeM[te_points_VM[1]] - timeM[q_points_VM[1]]
  
  ##LEAVE OFF:
  # % peak vectors QRS and T amplitude
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
  QRSTang <- RadToDeg(acos( Dot(Raxis,Taxis) / (sqrt(Raxis[1]^2 + Raxis[2]^2 + Raxis[3]^2) * sqrt(Taxis[1]^2 + Taxis[2]^2 + Taxis[3]^2))))
  
  # % mean QRS-T angle
  QRSTang_M <- RadToDeg(acos(Dot(MEAN_QRSO,MEAN_TO) / (sqrt(MEAN_QRSO[1]^2 + MEAN_QRSO[2]^2 + MEAN_QRSO[3]^2) * sqrt(MEAN_TO[1]^2 + MEAN_TO[2]^2 + MEAN_TO[3]^2))))
  
  # % Azimuth of QRS: peak, area 
  AZ_OQ <- (RadToDeg(acos(Raxis[1]/sqrt(Raxis[1]^2 + Raxis[3]^2)))) * (((Raxis[3]<0)*-1) + ((Raxis[3]>0)*1))
  AZ_OQM <- (RadToDeg(acos( MEAN_QRSO[1]/sqrt(MEAN_QRSO[1]^2 + MEAN_QRSO[3]^2)))) * (((MEAN_QRSO[3]<0)*-1) + ((MEAN_QRSO[3]>0)*1))
  
  # % Azimuth of T: peak, area 
  AZ_OT <- (RadToDeg(acos(Taxis[1]/sqrt(Taxis[1]^2 + Taxis[3]^2)))) * (((Taxis[3]<0)*-1) + ((Taxis[3]>0)*1))
  AZ_OTM <- (RadToDeg(acos(MEAN_TO[1]/sqrt(MEAN_TO[1]^2 + MEAN_TO[3]^2)))) * (((MEAN_TO[3]<0)*-1) + ((MEAN_TO[3]>0)*1))
  
  # % AzCimuth of SVG: peak, area 
  AZ_SVG <- (RadToDeg(acos(SVG_axis[1]/sqrt(SVG_axis[1]^2 + SVG_axis[3]^2)))) * (((SVG_axis[3]<0)*-1) + ((SVG_axis[3]>0)*1))
  AZ_SVG_M <- (RadToDeg(acos(SVG_MO[1]/sqrt(SVG_MO[1]^2 + SVG_MO[3]^2)))) * (((SVG_MO[3]<0)*-1)+((SVG_MO[3]>0)*1))
  
  # % Elevation of QRS: peak, area 
  EL_OQ <- (RadToDeg(acos(Dot(Raxis,Ynew)/(sqrt(Raxis[1]^2+Raxis[2]^2+Raxis[3]^2)*sqrt(Ynew[1]^2+Ynew[2]^2+Ynew[3]^2)))))
  EL_OQM <- (RadToDeg(acos(Dot(MEAN_QRSO,Ynew)/(sqrt(MEAN_QRSO[1]^2 + MEAN_QRSO[2]^2 + MEAN_QRSO[3]^2)*sqrt(Ynew[1]^2+Ynew[2]^2+Ynew[3]^2)))))
  
  # % Elevation of T: peak, area 
  EL_OT <- (RadToDeg(acos(Dot(Taxis,Ynew)/(sqrt(Taxis[1]^2+Taxis[2]^2+Taxis[3]^2)*sqrt(Ynew[1]^2+Ynew[2]^2+Ynew[3]^2)))))
  EL_OTM <- (RadToDeg(acos(Dot(MEAN_TO,Ynew)/(sqrt(MEAN_TO[1]^2+MEAN_TO[2]^2+MEAN_TO[3]^2)*sqrt(Ynew[1]^2+Ynew[2]^2+Ynew[3]^2)))))
  
  # % Elevation of SVG: peak, area 
  EL_SVG <- (RadToDeg(acos(Dot(SVG_axis,Ynew)/(sqrt(SVG_axis[1]^2+SVG_axis[2]^2+SVG_axis[3]^2)*sqrt(Ynew[1]^2+Ynew[2]^2+Ynew[3]^2)))))
  EL_SVG_M <- (RadToDeg(acos(Dot(SVG_MO,Ynew)/(sqrt(SVG_MO[1]^2+SVG_MO[2]^2+SVG_MO[3]^2)*sqrt(Ynew[1]^2+Ynew[2]^2+Ynew[3]^2)))))
  
  
  # % ========================== Magnitudes Calculation ============================
  # % Magnitude of QRS: peak, area
  QRS_Mag <- sqrt(Raxis[1]^2+Raxis[2]^2+Raxis[3]^2)
  QRS_Mag_M <- sqrt(MEAN_QRSO[1]^2+MEAN_QRSO[2]^2+MEAN_QRSO[3]^2)
  
  # % Magnitude of T: peak, area
  T_Mag <- sqrt(Taxis[1]^2+Taxis[2]^2+Taxis[3]^2)
  T_Mag_M <- sqrt(MEAN_TO[1]^2+MEAN_TO[2]^2+MEAN_TO[3]^2)
  
  # % Magnitude of SVG: peak
  SVG_Mag <- sqrt(SVG_axis[1]^2+SVG_axis[2]^2+SVG_axis[3]^2)
  
  # % Magnitude of WVG: Wilson's Ventricular Gradient
  
  WVG <- sqrt((SumVGx^2) + (SumVGy^2) + (SumVGz^2))
  
  return(WVG)
  
}

# GEH QRS --------------------------------------------------------------------
find_QRS_intervals <- function(leads12, ann12) {
  # Find average QR, RS intervals across all leads 
  # Can specify specific leads rather than all 12
  # Input prep for GEH
  
  QR_interval <- array(0,dim(ann12)[[2]])
  RS_interval <- array(0,dim(ann12)[[2]])
  
  for (i in 1:dim(ann12)[[2]]) {
    Rwaves <- make_wave_table(ann12[, i], wave_value = 2)
    Rwaves$peak <- unlist(peak_isolation(leads12[, i], ann12[, i], wave_value = 2))
    
    if (Rwaves$wave_on[1] == 1) {
      start <- 2
    } else{
      start <- 1
    }
    
    QR_indv <- c()
    RS_indv <- c()
    
    for (j in start:length(Rwaves[, 1])) {
      QR_indv <- c(QR_indv, Rwaves$peak[j] - Rwaves$wave_on[j])
      RS_indv <- c(RS_indv, Rwaves$wave_off[j] - Rwaves$peak[j])
    }
    
    QR_interval[i] <- list(QR_indv)
    RS_interval[i] <- list(RS_indv)
  }
  
  QR_mean <- mean(unlist(QR_interval), na.rm = TRUE)
  RS_mean <- mean(unlist(RS_interval), na.rm = TRUE)
  
  return(data.frame(QR = QR_mean, RS = RS_mean))
  
}


# GEH Toff ----------------------------------------------------------------
find_RToff_interval <- function(leads12, ann12) {
  # Find average R - T-offset interval across all 12 leads. 
  # Can specify specific leads rather than all 12
  # GEH input prep
  
  Rtoff_interval <- array(0,dim(ann12)[[2]])
  # single lead
  for (i in 1:dim(ann12)[[2]]) {
    Rwaves <- make_wave_table(ann12[,i], wave_value = 2)
    Rwaves$peak <- unlist(peak_isolation(leads12[,i], ann12[,i], wave_value = 2))
    Twaves <- make_wave_table(ann12[,i], wave_value = 3)
    Twaves$peak <- unlist(peak_isolation(leads12[,i], ann12[,i], wave_value = 3))
    
    combined <- rbind(Rwaves, Twaves)
    combined <- combined[order(combined$wave_on),]
    
    # skip first t-wave if there's no preceding R-wave:
    
    t_ind <- which(combined$wave_type == "t")
    if (t_ind[1] == 1) {
      start <- 2
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
  
  return(mean(unlist(Rtoff_interval), na.rm = TRUE))
}

# GEH Tpeak --------------------------------------------------------------
find_RTpeak_interval <- function(leads12, ann12) {
  # Find R - T-peak
  # handles one sample
  Rt_interval <- array(0,dim(ann12)[[2]])
  for (i in 1:dim(ann12)[[2]]) {
  Rpeaks <- peak_isolation(leads12[,i], ann12[,i], wave_value = 2)
  Tpeaks <- peak_isolation(leads12[,i], ann12[,i], wave_value = 3)
  
  R_table <- data.frame(wave_type = "N", index = unlist(Rpeaks))
  T_table <- data.frame(wave_type = "t", index = unlist(Tpeaks))
  
  combined <- rbind(R_table,T_table)
  combined <- combined[order(combined$index),]
  
  t_ind <- which(combined$wave_type == "t")
  
  # skip first p-wave if there's no preceding t-wave:
  if (t_ind[1] == 1) {
    start <- 2
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
  
  return(mean(unlist(Rt_interval), na.rm = TRUE))
  
}

# GEH RPeak ----------------------------------------------------------------
find_geh_Rpeak <- function(XYZ_M, origin) {
  # Find R peak value across all 3D spatial leads (XYZ)
  # Input prep for GEH, added to the above GEH functions to find Q, S, T_peak, 
  # T_off indices
  
  Rpeaks <- array(0,3)
  for (i in 1:3){
    Rpeaks[i] <- which.max(abs(XYZ_M[,i] - XYZ_M[origin,i]))
  }
  mean <- mean(Rpeaks)
  return(mean)
}

# Median Beat Function ----------------------------------------------------
find_median_beat <- function(XYZ, Rpeaks, fs = 500) {
  # Inputs: 10 second XYZ, Rpeaks for a single lead
  # Output: one single average beat for all 3 spatial leads, centered on the 
  # Rpeak
  
  # Rx <- c(10,663,1343,2001,2644,3315,3970,4626) # Rpeaks on lead x
  
  # XYZ_O <- kors(input) # XYZ 10 sec ECG. 5000 x 3 ***
  XYZ_O <- XYZ
  
  if (class(Rpeaks) == "list") {
    Rpeaks <- unlist(Rpeaks)
  }
  
  window <- 0.1 # ms
  window_ind <- round(window*fs)
  
  Rx <- array(0,length(Rpeaks))
  for (i in 1:length(Rpeaks)) {
    
    start <- (Rpeaks[i] - window_ind)
    end <- (Rpeaks[i] + window_ind)
    
    start[start<1] <- 1
    end[end>length(XYZ_O[,1])] <- length(XYZ_O[,1])

    Rx[i] <- which.max(abs(XYZ_O[start:end,1])) + start - 1
  }
  
  R_no <- length(Rx)
  
  Total_samples <- dim(XYZ_O)[[1]] # verify dimension used for length
  Beat_length <- fs*1.2
  
  
  if (Rpeaks[1] < Beat_length/2) {
    c_int = 2
  } else {
    c_int = 1
  }
  
  if (abs(Rpeaks[length(Rpeaks - 1)] - Total_samples) < (Beat_length / 2) ) {
    c_fin = 2
  } else {
    c_fin = 1
  }
  
  I_dv=array(0,R_no); # I_dv: max slope change centered in each beat
  
  for (ii in c_int:(R_no-c_fin)) {
    
    I_dv_temp <- which.max(abs(diff(XYZ_O[ (Rx[ii] - 199) : (Rx[ii] + 200), 1])))
    # [value, index]
    # [~,I_dv_temp]=max(abs(diff(XYZ_O(     Rx(ii)-199:Rx(ii)+200,1    ))));
    
    I_dv[ii]=Rx[ii] + I_dv_temp - 199
  }
  
  
  # %% Separate the each beat in the X, Y, and Z leads. Each beat is centered on Maximum dv/dt detected.
  
  total_samples <- length(c_int:R_no-c_fin)
  
  Beats_x_T <- array(0, c(Beat_length, total_samples))# vector for centered beat ranges
  Beats_y_T <- array(0, c(Beat_length, total_samples))
  Beats_z_T <- array(0, c(Beat_length, total_samples))
  
  bct=1;
  for (ii in c_int:(R_no-c_fin)) {
    
    Beats_x_T[,bct] <- XYZ_O[(I_dv[ii] - (Beat_length/2 - 1)) : (I_dv[ii] + (Beat_length/2)), 1]
    Beats_y_T[,bct] <- XYZ_O[(I_dv[ii] - (Beat_length/2 - 1)) : (I_dv[ii] + (Beat_length/2)), 2]
    Beats_z_T[,bct] <- XYZ_O[(I_dv[ii] - (Beat_length/2 - 1)) : (I_dv[ii] + (Beat_length/2)), 3]
    
    bct=bct+1;
  }
  
  # %% Create the median beat of each X, Y, and Z leads 
  
  XYZ_M_T <- array(0, c(Beat_length, 3))
  for (i in 1 : Beat_length) {
    XYZ_M_T[i,1] <- median(Beats_x_T[i,])
    XYZ_M_T[i,2] <- median(Beats_y_T[i,])
    XYZ_M_T[i,3] <- median(Beats_z_T[i,]) 
  }
  
  VecMag_T <- array(0, Beat_length)
  # %% Calculate the Vector Magnitude (Euclidean norm)
  for (ii in 1:length(XYZ_M_T[,1]) ) {
    VecMag_T[ii] <- norm(XYZ_M_T[ii,], type = "2");
  }
  
  timeM=((1:length(VecMag_T))/fs)*1000;
  
  return(XYZ_M_T)
}

# Origin Function ---------------------------------------------------------
find_origin <- function(XYZ_M, Rpeaks, fs = 500) {
  # Finds the isoelectric point across all 3D spatial leads
  # Input: median XYZ, Rpeaks
  
  if (class(Rpeaks) == "list") {
    Rpeaks <- unlist(Rpeaks)
  }
  
  Rx <- Rpeaks
  
  # Rx <- which.max(abs(XYZ_M_T[,1] - median(XYZ_M_T[,1])))
  # Ry <- which.max(abs(XYZ_M_T[,2] - median(XYZ_M_T[,2])))
  # Rz <- which.max(abs(XYZ_M_T[,3] - median(XYZ_M_T[,3])))
  
  # Rx is the peak of the median beat** --> ~300

  
  # %% Calculate various windows analysis
  w <- floor( (mean(diff(Rx))*.8) - (372.69*(fs/1000)))
  
  # %% Is the RR interval < 600ms (0.6*fs)
  if (mean(diff(Rx)) < (600*(fs/1000)) ) { # (A.Label_beat~='S') ??? 
    cent <- (600-260)*(fs/1000)
  } else {
    cent <- (600-320)*(fs/1000)
  }
  
  
  w1_start <- cent-w
  w1_end <- cent+w
  
  w2_start <- cent-(160*(fs/1000))
  w2_end <- cent+(160*(fs/1000))
  
  if (w1_start < w2_start) {
    w1_start <- w2_start
  }
  
  if (w1_end > w2_end) {
    w1_end <- w2_end
  }
  
  # %% cut xyz median beat based on the window
  xyz_w1 = XYZ_M[w1_start:w1_end,]
  xyz_w2 = XYZ_M[w2_start:w2_end,]
  
  # %% calculate origin point with clustering algorithm ori_clustering
  output <- ori_clustering(xyz_w1)
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
  output <- abs_grad_fun(xyz_w1)
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
  
  Area_m1 <- AUC(x = x, y = VecMag_M1[w1_start:w1_end])
  Area_m2 <- AUC(x = x, y = VecMag_M2[w1_start:w1_end])
  
  
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
  diff_w1 <- data_processing(xyz_median,'gradsum',2)
  iter_num <- 0
  min_idx <- which.min(diff_w1)
  
  for (ii in 1:(length(diff_w1)-1)) {
    if (abs(diff_w1[length(diff_w1)-ii] - diff_w1[min_idx]) < 0.1) {
      iter_num <- length(diff_w1) - ii
    }
  }
  
  ori_pt_idx <- iter_num
  ori_val <- xyz_median[iter_num]
  
  output <- list(ori_pt_idx, ori_val)
  return(output)
}

ori_clustering <- function(xyz_median) {
  processed_data  <-  data_processing(xyz_median,'norm_mvvar',10);
  
  bindwidth <- 10
  bin_counts <- round((max(xyz_median) - min(xyz_median)) / bindwidth)
  breaks <- bin_counts + 1
  
  histogram <- hist(processed_data, breaks = breaks, plot = FALSE)
  bin_label <- cut(processed_data, histogram$breaks, labels = FALSE) # which bin each point is in
  N <- histogram$counts # number of points in each bin
  
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
    data_gradient_sum = data_processing(xyz_median,'gradsum',5)
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
  library(zoo)
  if (process == "norm_mvvar") {
    data_norm <- matrix(data = NA, dim(data)[[1]])
    for (ii in 1:length(data[, 1])) {
      data_norm[ii] <- norm(c(data[ii, 1], data[ii, 2], data[ii, 3]), type = "2")
      
    }
    processed_data <- rollapply(data_norm, width = process_window, FUN = var, fill = NA)
    
    end <- length(processed_data)
    for (i in 1: ((process_window - 1)/2)) {
      processed_data[i] <- var(data_norm[1: (i + (process_window - 1)/2) ])
      processed_data[end - i + 1] <- var(data_norm[(end - i + 1 - (process_window - 1)/2) : end ])
    }
    
    
    # gradient:
  } else if (process == "gradsum") {
    # **need to take averages while keeping them in their columns. Then take abs, gradient, and add up
    # pre_processed_data <- abs(gradient(rollmean(data[, 1], process_window, fill = NA))) + abs(gradient(rollmean(data[, 2], process_window, fill = NA))) + abs(gradient(rollmean(data[, 3], process_window, fill = NA)))
    
    pre_processed_data_x <- rollmean(data[, 1], process_window, fill = NA) 
    pre_processed_data_y <- rollmean(data[, 2], process_window, fill = NA) 
    pre_processed_data_z <- rollmean(data[, 3], process_window, fill = NA) 
    
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
korsMatrix <- array(c(0.38, -0.07, 0, 0, 0, 0, -0.13, 0.05, -0.01, 0.14, 0.06,0.54,  
                      -0.07, 0.93, 0, 0, 0, 0, 0.06, -0.02, -0.05, 0.06, -0.17, 0.13,  
                      0.11, -0.23, 0, 0, 0, 0, -0.43, -0.06,-0.14,-0.20,-0.11,0.31), 
                    c(3,12))

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