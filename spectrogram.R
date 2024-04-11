
# Filter ------------------------------------------------------------------


# **Decide raw vs. filtered data
# use full data without zero'd margins?

#**run bandpass first - remove baseline wander**

load("samples_full12.RData")

library(dplR)
f <- 500   # sample frequency
low <- 0.5 # high pass frequency threshold. Original: 0.5. 1 works well
high <- 40 # low  pass frequency threshold

dim <- dim(samples12)
samples_bp_full12 <- array(0,dim[1:3])

for (i in 1:dim[1]) {
  for (j in 1:dim[3]) {
    samples_bp_full12[i,,j] <- pass.filt(samples_full12[i,,j], c(f/high, f/low), "pass", method = c("Butterworth"), n = 4, Rp = 0.1)
  }
}


# Spectrogram -------------------------------------------------------------
library(signal)

input <- samples_bp_full12
nam <- "samples_filt_spec_i" # samples_spec_i

output <- array(0, c(170,5000,24))

for (i in 1:dim(output)[[1]]) {
a <- specgram(c(array(0, 64), input[i,,1], array(0, 64)), n = 128, Fs = 500, overlap = 127)

# a$f <- a$f[a$f < 45]
index = length(a[a$f < 45])

a$S <- a$S[1:index,]
a$f <- a$f[a$f < 45]

b <- (a$S - mean(a$S)) / sd(a$S)


real <- Re(b)
imaginary <- Im(b)

output[i,,] <- t(rbind(real,imaginary))

}

for (i in 1:dim(output)[[1]]) {
  # for (j in 1:dim(samples_full12)[[3]]) {
  output[i, 0:(which(annotations12[i, , 1] != 0)[[1]] - 1), ] <- 0
  output[i, (tail(which(annotations12[i, , 1] != 0), n = 1) + 1):dim(annotations12)[[2]], ] <- 0
  
}

assign(nam,output)

# a$S[i, 0:(which(annotations12[i, , 1] != 0)[[1]] - 1), ] <- 0
# a$S[i, (tail(which(annotations12[i, , 1] != 0), n = 1) + 1):dim(annotations12)[[2]], ] <- 0
# 
