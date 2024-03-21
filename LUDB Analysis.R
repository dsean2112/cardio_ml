library(keras)
library(caret)
model <- load_model_tf("saved_model/model_lstm1")



# Predict Samples ---------------------------------------------------------
if (exists("samples") == 0){
  load('LUDB_samples.RData')
  split <- 0.7
  train_size <- round(nrow(samples)*0.7)
  train_rows <- sample(nrow(samples), train_size)
  
  train_samples <- samples[train_rows,,]
  test_samples <- samples[-train_rows,,]
}

predictions <- model %>% predict(test_samples[,,1])

a <- dim(predictions)
predictions_integer <- array(0,a[1:2]) 

# use max.col to determine most probable value  for each time point
for (i in 1: a) {
  predictions_integer[i,] <- max.col(predictions[i,,])
}

#convert from dimension value 1,2,3,4 to 0,1,2,3
predictions_integer <- predictions_integer - 1

#binary true/false if prediciton matches annotations 
accuracy <- predictions_integer == test_samples[, ,2]


# Confusion Tables --------------------------------------------------------

confusion_matrix <- confusionMatrix(data = factor(predictions_integer), reference = factor(test_samples[,,2]))
confusion_table <- t(t(confusion_matrix$table) / colSums(confusion_matrix$table)) * 100
round(confusion_table, 1)


levels <- c(0,1,2,3) # number of value types must be the same when factoring
confusion_table_by_sample <- array(0,c(a[1],4,4))
for (i in 1:a[1]) {
  confusion_matrix_single  <- confusionMatrix(factor(predictions_integer[i,], levels = levels), factor(test_samples[i,,2], levels = levels))
  confusion_table_single <- t(t(confusion_matrix_single$table) / colSums(confusion_matrix_single$table)) * 100
  confusion_table_by_sample[i,1:4,1:4] <- confusion_table_single
}

confusion_by_sample_condensed <- array(0,c(a[1],4))
for (i in 1:4) {
  confusion_by_sample_condensed[,i] <- confusion_table_by_sample[,i,i]
}


# Plots -------------------------------------------------------------------
sample_number <- 24

library(ggpubr)

# Plot with Correct Labels:

correct_frame <- data.frame(Time = 1:5000 / 500, Signal = test_samples[sample_number, , 1])
#color code
colors <- test_samples[sample_number, , 2]
colors[colors == 1] <- 'p'
colors[colors == 2] <- 'N'
colors[colors == 3] <- 't'

correct_plot <-
  ggplot(correct_frame, aes(Time, Signal, color = colors)) + geom_path(linewidth =
                                                                      1, aes(group = 1)) + geom_point() + scale_x_continuous(breaks = seq(0, 10, 1))

# Plot with Predicted Labels:

title <- paste("Accuracy of Sample Number", sample_number)
color_frame <- data.frame(Time = 1:5000 / 500, Signal = test_samples[sample_number, , 1])
#color code
colors2 <- predictions_integer[sample_number, ]
colors2[colors2 == 1] <- 'p'
colors2[colors2 == 2] <- 'N'
colors2[colors2 == 3] <- 't'

# Plot for labeling predicted annotation
predicted_plot <-
  ggplot(color_frame, aes(Time, Signal, color = colors2)) + 
  geom_path(linewidth = 1, aes(group = 1)) + geom_point() + 
  scale_x_continuous(breaks = seq(0, 10, 1)) + 
  labs(title = title) + theme(legend.title=element_blank())

ggarrange(predicted_plot,correct_plot, nrow = 2,ncol = 1)
