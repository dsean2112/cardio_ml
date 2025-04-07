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

nputs <- layer_input(shape = c(5000, 24), dtype = 'float32')

outputs <-
  inputs |>
  layer_normalization() |>
  bidirectional(layer_lstm(units = 200, return_sequences = 'True', activation = 'tanh')) |>
  layer_dense(units = 4, activation = 'sigmoid', name = 'predictions')
# units = 200
model <- keras_model(inputs = inputs, outputs = outputs, name = 'mdl')

model |>
  compile(
    optimizer = 'adam',
    loss = 'sparse_categorical_crossentropy',
    metrics = 'accuracy'
  )

# Train model -------------------------------------------------------------
epochs <- 1
history <- model |> fit(training_signal, training_annotations, 
                        epochs = epochs, 
                        validation_data = list(testing_signal, testing_annotations),
                        verbose = 2)

# model_name <- paste0("models/",output_name)
# save_model_tf(model, model_name)