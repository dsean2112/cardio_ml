# Prep --------------------------------------------------------------------
source('code/annotator_prep_functions.R')
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
bilstm_layers <- 64 # original: 64

n_classes <- length(unique(as.vector(training_annotations)))

# Input layer
# 5000 steps and 1 channel at a time
input_shape <- c(5000, 1)
inputs <- layer_input(shape = input_shape)

# Convolutional Block
conv_block <-
  inputs |>
  layer_conv_1d(filters = 64, kernel_size = 5, activation = "relu", padding = "same") |>
  layer_batch_normalization() |>
  layer_max_pooling_1d(pool_size = 2, data_format = 'channels_first') |>
  layer_dropout(rate = dropout)  # Dropout after BiLSTM

# BiLSTM Block
bilstm_block <-
  conv_block |>
  bidirectional(layer_lstm(units = bilstm_layers, return_sequences = TRUE)) |>
  layer_dropout(rate = dropout)  # Dropout after BiLSTM

# Self Attention Block
attn_block <-
  bilstm_block |>
  layer_dense(units = 1, activation = "tanh") |>
  layer_flatten() |>
  layer_activation("softmax") |>
  layer_repeat_vector(bilstm_layers*2) |>
  layer_permute(c(2, 1))

mult_attn_block <- layer_multiply(list(bilstm_block, attn_block))

# Time Distributed Dense Layer
outputs <-
  mult_attn_block |>
  time_distributed(layer_dense(units = n_classes, activation = "sigmoid"))

# Create and compile the model
cnn_bilstm_attn_model <- keras_model(inputs = inputs, outputs = outputs)

cnn_bilstm_attn_model |> compile(
  # optimizer = optimizer_adam(),
  optimizer = 'rmsprop',
  loss = "sparse_categorical_crossentropy",
  # loss_weights = c(0.1, 1, 1, 1),  # Adjust weights for each class
  metrics = c("accuracy")
  
)


# Train model -------------------------------------------------------------
epochs <- 5
history <- cnn_bilstm_attn_model |> fit(training_signal, training_annotations, 
                                        epochs = epochs,
                                        # validation_data = list(val_inputs, val_targets),
                                        validation_data = list(testing_signal, testing_annotations),
                                        verbose = 2)


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

sample <- 10
plot_func(testing_signal[sample,],predictions_integer[sample,])