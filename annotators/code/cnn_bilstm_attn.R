# Prep --------------------------------------------------------------------
source('annotator_prep_functions.R')
annotator_style <- 2
lead <- 1
out <- prep_ludb(lead = lead, annotator_style = annotator_style)

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

dropout <- 0.3
bilstm_layers <- 64 # original: 64
mask_value <- 0

model_type <- 'CBA' # ie: cnn_bilstm_attn

date_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
model_name <- paste0(model_type,'_',date_time)

num_classes <- length(unique(as.vector(training_annotations)))

# Input layer
# 5000 steps and 1 channel at a time
input_shape <- c(5000, 1)
inputs <- layer_input(shape = input_shape) |>
  layer_masking(mask_value = mask_value) # Masking input values equal to 0
  # layer_embedding(input_dim = 5000, output_dim = 5000, mask_zero = TRUE)

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
  time_distributed(layer_dense(units = num_classes, activation = "sigmoid"))

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

# Print to command line
cat("model_type:", model_type, "\n",
    "model_name:", model_name, "\n",
    "bilstm_layers:", bilstm_layers, "\n",
    "epochs:", epochs, "\n",
    "dropout:", dropout, "\n",
    "annotator_style:", annotator_style, "\n")

start <- Sys.time()
history <- cnn_bilstm_attn_model |> fit(training_signal, training_annotations, 
                                        epochs = epochs,
                                        # validation_data = list(val_inputs, val_targets),
                                        validation_data = list(testing_signal, testing_annotations),
                                        # loss_weights = c(0.1, 1, 1, 1), # experimenting with decreasing the weight of '0' annotations (ie no wave)
                                        verbose = 2)

# output_name <- paste0("../models/",model_name)
# save_model_tf(model, output_name)

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
confusion <- confusion_analysis()


# Write to log ------------------------------------------------------------
# Add lock to avoid overwriting
library(filelock)
logFile <- '../models/model_log.RData'
lockFile <- paste0(logFile, ".lock")  # Create a lock file
# Acquire the lock before loading or saving
lock <- lock(lockFile)
load(logFile)

logFile <- '../models/model_log.RData'
load(logFile)


new_row <- data.frame(
  name = model_name,
  type = model_type,
  ann_style = annotator_style,
  lead = lead,
  bilstm_layers = bilstm_layers,
  dropout = dropout,
  filters = NA,
  epochs = epochs,
  time = round(time_spent, 2),
  training_samples = I(list(out$training_samples)),
  confusion = I(list(confusion))
)

model_log <- rbind(model_log, new_row)
save(model_log, file = logFile)

# Release the lock
unlock(lock)

# plot --------------------------------------------------------------------

# sample <- 60
# plot_func(testing_signal[sample,],predictions_integer[sample,])