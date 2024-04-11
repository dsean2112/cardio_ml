library(keras)
library(caret)

# Inputs:: 
range <- 1 # Lead number. 1:12, 1, etc. 1=i, 4=avr, 7=v1...
input <- samples_filt_spec_i #samples_bp12. sample inputs
model_name <- "spectro_i"



leads <- c("i","ii","iii","avr","avl","avf","v1","v2","v3","v4","v5","v6")

# for (lead in range) {
  model <- load_model_tf(paste0("saved_model/",model_name))
  
  
  a <- 1:170
  test_rows <- a[-train_rows]
  
  predictions <- model %>% predict(input[test_rows,,])
  
  a <- dim(predictions)
  predictions_integer <- array(0,a[1:2]) 
  
  # use max.col to determine most probable value  for each time point
  for (i in 1: a) {
    predictions_integer[i,] <- max.col(predictions[i,,])
  }
  
  #convert from dimension value 1,2,3,4 to 0,1,2,3
  predictions_integer <- predictions_integer - 1
  
  
  # Confusion Tables --------------------------------------------------------
  levels <- c(0,1,2,3) # number of value types must be the same when factoring
  
  confusion_matrix <- confusionMatrix(data = factor(predictions_integer, levels = levels), reference = factor(annotations12[test_rows,,lead], levels = levels))
  confusion_table <- t(t(confusion_matrix$table) / colSums(confusion_matrix$table)) * 100
  round(confusion_table, 1)
  assign(paste0("confusion_full_",model_name), confusion_matrix)
  assign(paste0("confusion_table_",model_name), confusion_table)
  
  
  # confusion_table_by_sample <- array(0,c(a[1],4,4))
  # for (i in 1:a[1]) {
  #   confusion_matrix_single  <- confusionMatrix(factor(predictions_integer[i,], levels = levels), factor(testing_data[i,,2], levels = levels))
  #   confusion_table_single <- t(t(confusion_matrix_single$table) / colSums(confusion_matrix_single$table)) * 100
  #   confusion_table_by_sample[i,1:4,1:4] <- confusion_table_single
  # }
  # 
  # confusion_by_sample_condensed <- array(0,c(a[1],4))
  # for (i in 1:4) {
  #   confusion_by_sample_condensed[,i] <- confusion_table_by_sample[,i,i]
  # }
# }