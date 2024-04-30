input_truth <- annotations12[test_rows, , 1]
input_predicted <- predictions_integer
wave_threshold <- 0.75

wave_accuracy <- array(0, c(dim(input_truth)[[1]], 3, 2))

for (sample_number in 1:dim(input_truth)[[1]]) {
  # sample_number <- 1
  for (wave in 1:3) {
    # wave <- 1 # 1 = p, 2 = qrs, 3 = t
    
    a <- which(input_truth[sample_number,] == wave)
    if (length(a) != 0) {
      b <- a[-1]
      len <- length(a)
      
      start_stop_index <-
        which(a[-len] + 1 != b) # Identify when the index increases by more than + 1
      start_stop_len <-
        length(start_stop_index) + 1 # add 1 length to account for 1st on, and last off times
      start_stop <- array(0, c(start_stop_len, 2))
      
      start_stop[1, 1] <- a[1]
      start_stop[start_stop_len, 2] <- tail(a, n = 1)
      
      for (i in 1:length(start_stop_index)) {
        start_stop[i, 2] <- a[start_stop_index[i]]
        start_stop[i + 1, 1] <- a[start_stop_index[i] + 1]
      }
      
      
      # test --------------------------------------------------------------------
      
      # n <- 4
      count <- 0
      for (n in 1:dim(start_stop)[[1]]) {
        sum <-
          sum(input_predicted[sample_number, start_stop[n, 1]:start_stop[n, 2]] == wave) / length(start_stop[n, 1]:start_stop[n, 2])
        if (sum > wave_threshold) {
          count <- count + 1
        }
      }
      wave_accuracy[sample_number, wave, ] <-
        c(count, dim(start_stop)[[1]])
    }
  }
}

p_accuracy <- sum(wave_accuracy[,1,1]) / sum(wave_accuracy[,1,2])
qrs_accuracy <- sum(wave_accuracy[,2,1]) / sum(wave_accuracy[,2,2])
t_accuracy <- sum(wave_accuracy[,3,1]) / sum(wave_accuracy[,3,2])

# for each sample:
# [ pwave_seen, pwave_total
#   qrs_seen, qrs_total
#   t_seen, t_total]