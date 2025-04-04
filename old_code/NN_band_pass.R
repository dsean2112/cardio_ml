library(dplR)

f <- 500   # sample frequency
low <- 0.5 # high pass frequency threshold
high <- 40 # low  pass frequency threshold

dim <- dim(samples)
samples_bp <- array(0,dim[1:2])

for (i in 1:dim[1]) {
  samples_bp[i,] <- pass.filt(samples[i,,1], c(f/high, f/low), "pass", method = c("Butterworth"), n = 4, Rp = 0.1)
}

samples_bp <- array(c(samples_bp, samples[,,2]),c(170,5000,2))

# Padding signal prior to first and last annotation with zeros:
for (i in 1:dim(samples_bp)[[1]]){
  samples_bp[i,0:(which(samples_bp[i,,2] != 0 )[[1]]-1),1] <- 0 # can switch to other value
  samples_bp[i, (tail(which(samples_bp[i,,2] !=0), n = 1) + 1) :dim(samples_bp)[[2]],1] <- 0 # can switch to other value
}