# Prep --------------------------------------------------------------------
source('annotator_prep_functions.R')
annotator_style <- 2
lead <- 1
rounds <- 3

out <- prep_ludb(lead = lead, annotator_style = 2, rounds = rounds)
#         1: 1 0 0 0 1 0 0 0 2 0 0 2 ...
#         2: 1 1 1 1 1 0 0 0 2 2 2 2 ...
#         3: 1 2 2 2 3 0 0 0 4 5 5 6 ...

training_signal <- out$training_signal
training_annotations <- out$training_annotations
testing_signal <- out$testing_signal
testing_annotations <- out$testing_annotations



# Build Model -------------------------------------------------------------------
library(keras)
bilstm_layers <- 200 # original: 200
activation <- 'softmax' #'sigmoid': old version
mask_value <- 0
model_type <- 'bilstm'

date_time <- format(Sys.time(), "%Y%m%d_%H%M%S")

model_name <- paste0(model_type,'_',date_time)

num_classes <- length(unique(as.vector(training_annotations)))

inputs <- layer_input(shape = c(5000, 1), dtype = 'float32') |>
  layer_masking(mask_value = mask_value) # Masking input values equal to 0

outputs <-
  inputs |>
  layer_normalization() |> # normalize layers
  bidirectional(layer_lstm(units = bilstm_layers, return_sequences = 'True', activation = 'tanh')) |>
  layer_dense(units = num_classes, activation = activation, name = 'predictions')

model <- keras_model(inputs = inputs, outputs = outputs, name = 'mdl')

model |>
  compile(
    # optimizer = 'adam',
    optimizer = 'rmsprop',
    loss = 'sparse_categorical_crossentropy',
    metrics = 'accuracy'
  )

# Train model -------------------------------------------------------------
epochs <- 20

# Print to command line:
cat("model_type:", model_type, "\n",
    "model_name:", model_name, "\n",
    "bilstm_layers:", bilstm_layers, "\n",
    "epochs:", epochs, "\n",
    "activation:", activation, "\n",
    "annotator_style:", annotator_style, "\n")

# Train
start <- Sys.time()
history <- model |> fit(training_signal, training_annotations, 
                        epochs = epochs, 
                        validation_data = list(testing_signal, testing_annotations),
                        verbose = 2)

output_name <- paste0("../models/",model_name)
save_model_tf(model, output_name)

end <- Sys.time()
time_spent <- end-start

# bilstm: bilstm_layers, activation
# cnn_bilstm_attn: dropout, bilstm_layers
# UNET: dropout, filters
# **specify training/testing samples from out$

# Add to 

# Test model --------------------------------------------------------------
# input_filtered <- filter_samples()

predictions <- model %>% predict(testing_signal)

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


new_row <- data.frame(
  name = model_name,
  type = model_type,
  ann_style = annotator_style,
  lead = lead,
  bilstm_layers = bilstm_layers,
  dropout = NA,
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

# sample <- 12
# plot_func(testing_signal[sample,],predictions_integer[sample,])