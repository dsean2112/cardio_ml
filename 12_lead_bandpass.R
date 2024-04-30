# Filter ------------------------------------------------------------------
library(dplR)
library(ggplot2)
library(ggpubr)

f <- 500   # sample frequency
low <- 0.5 # high pass frequency threshold. Original: 0.5
high <- 40 # low  pass frequency threshold

dim <- dim(samples12)
samples_bp12 <- array(0,dim[1:3])

for (i in 1:dim[1]) {
  for (j in 1:dim[3]) {
    samples_bp12[i,,j] <- pass.filt(samples12[i,,j], c(f/high, f/low), "pass", method = c("Butterworth"), n = 4, Rp = 0.1)
  }
}

# samples_bp <- array(c(samples_bp, samples[,,2]),c(170,5000,2))

# Pad with zeros
for (i in 1:dim(samples_bp12)[[1]]) {
  for (j in 1:dim(samples_bp12)[[3]]) {
    samples_bp12[i, 0:(which(annotations12[i, , j] != 0)[[1]] - 1), j] <- 0
    samples_bp12[i, (tail(which(annotations12[i, , j] != 0), n = 1) + 1):dim(annotations12)[[2]], j] <- 0
  }
}

# Graph -------------------------------------------------------------------
sample_number <- 130
lead <- 1

original_frame <- data.frame(Time = 1:5000 / 500, Signal = samples12[sample_number, , lead])
filter_frame <- data.frame(Time = 1:5000 / 500, Signal = samples_bp12[sample_number, , lead])

colors <- annotations12[sample_number, , lead]
colors[colors == 1] <- 'p'
colors[colors == 2] <- 'N'
colors[colors == 3] <- 't'

original_plot <-
  ggplot(original_frame, aes(Time, Signal, color = colors)) + 
  geom_path(linewidth = 1, aes(group = 1)) + geom_point() + 
  scale_x_continuous(breaks = seq(0, 10, 1)) + 
  labs(title = "Original") + theme(legend.title=element_blank())

filter_plot <-
  ggplot(filter_frame, aes(Time, Signal, color = colors)) + 
  geom_path(linewidth = 1, aes(group = 1)) + geom_point() + 
  scale_x_continuous(breaks = seq(0, 10, 1)) + 
  labs(title = "Filter") + theme(legend.title=element_blank())


ggarrange(original_plot,filter_plot, nrow = 2,ncol = 1)



# other -------------------------------------------------------------------



x <- 1:5000 / 500

a <- 1:170
test_rows <- a[-train_rows]

leads <- c("i","ii","iii","avr","avl","avf","v1","v2","v3","v4","v5","v6")

predictions <- array(0, c(2,5000,12,4))

for (i in 1 : 12) {
  model <- load_model_tf(paste0("saved_model/lstm_", leads[i]))
  predictions <- model %>% predict(samples12[test_rows,,i])
  predictions <- predictions[1,,]
  
  
  predictions_integer <- max.col(predictions) - 1
  
  predictions_integer[predictions_integer == 1] <- 'p'
  predictions_integer[predictions_integer == 2] <- 'N'
  predictions_integer[predictions_integer == 3] <- 't'
  
  
  # frame <- data.frame(Time = x, Signal = samples12[test_rows[n],,i])
  # plot <-
  #   ggplot(frame, aes(Time, Signal, color = predictions_integer)) + 
  #   geom_path(linewidth = 0.25, aes(group = 1)) + geom_point() + 
  #   scale_x_continuous(breaks = seq(0, 10, 1)) + theme(legend.position="none") + theme(legend.title = element_blank())
  # nam <- paste("lead_", leads[i], sep = "")
  # assign(nam, plot)
}