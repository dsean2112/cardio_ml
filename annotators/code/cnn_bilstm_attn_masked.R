# Prep --------------------------------------------------------------------
source('annotator_prep_functions.R')
out <- prep_ludb(lead = 1, annotator_style = 2)
#         1: 1 0 0 0 1 0 0 0 2 0 0 2 ...
#         2: 1 1 1 1 1 0 0 0 2 2 2 2 ...
#         3: 1 2 2 2 3 0 0 0 4 5 5 6 ...

training_signal <- out$training_signal
training_annotations <- out$training_annotations
testing_signal <- out$testing_signal
testing_annotations <- out$testing_annotations

# Build Model -------------------------------------------------------------------
library(keras)
remove(cnn_bilstm_attn_model)

dropout <- 0.5
bilstm_units <- 64
mask_value <- 0

num_classes <- length(unique(as.vector(training_annotations)))

# Input layer
input_shape <- c(5000, 1)

inputs <- layer_input(shape = input_shape)

# LSTM Branch (mask-aware)
masked_inputs <- inputs %>% 
  layer_masking(mask_value = mask_value)
lstm_branch <- masked_inputs %>% 
  bidirectional(layer_lstm(units = bilstm_units, return_sequences = TRUE)) %>%
  layer_dropout(rate = dropout)

# Convolution Branch (masking not propagated here)
conv_branch <- inputs %>% 
  layer_conv_1d(filters = 64, kernel_size = 5, activation = "relu", padding = "same") %>%
  layer_batch_normalization() %>%
  layer_max_pooling_1d(pool_size = 2, data_format = 'channels_first') %>%
  layer_dropout(rate = dropout)

# Depending on your design, you could merge these branches.
# For instance, if they’re spatially compatible after adjustments:
# Here we simply use the lstm_branch for the remaining network:
attn_block <- lstm_branch %>%
  layer_dense(units = 1, activation = "tanh") %>%
  layer_flatten() %>%
  layer_activation("softmax") %>%
  layer_repeat_vector(bilstm_units * 2) %>%
  layer_permute(c(2, 1))

mult_attn_block <- layer_multiply(list(lstm_branch, attn_block))

outputs <- mult_attn_block %>%
  time_distributed(layer_dense(units = num_classes, activation = "softmax"))

cnn_bilstm_attn_model <- keras_model(inputs = inputs, outputs = outputs)

cnn_bilstm_attn_model %>% compile(
  optimizer = 'rmsprop',
  loss = "sparse_categorical_crossentropy",
  # loss_weights = c(0.1, 1, 1, 1),  # Adjust weights for each class
  metrics = c("accuracy")
)

# Train model -------------------------------------------------------------
start <- Sys.time()
epochs <- 20
history <- cnn_bilstm_attn_model |> fit(training_signal, training_annotations, 
                                        epochs = epochs,
                                        # validation_data = list(val_inputs, val_targets),
                                        validation_data = list(testing_signal, testing_annotations),
                                        # loss_weights = c(0.1, 1, 1, 1), # experimenting with decreasing the weight of '0' annotations (ie no wave)
                                        verbose = 2)

# model_name <- paste0("../models/",output_name)
# save_model_tf(model, model_name)

end <- Sys.time()
time_spent <- end-start

# Test model --------------------------------------------------------------
# input_filtered <- filter_samples()

predictions <- cnn_bilstm_attn_model %>% predict(testing_signal)

predictions_integer <- array(0,c(nrow(predictions),ncol(predictions))) 

for (i in 1: nrow(predictions)) {
  predictions_integer[i,] <- max.col(predictions[i,,])
}

#convert from dimension value 1,2,3,4 to 0,1,2,3
predictions_integer <- predictions_integer - 1


# Analysis
confusion_analysis()

# plot --------------------------------------------------------------------

sample <- 60
plot_func(testing_signal[sample,],predictions_integer[sample,])