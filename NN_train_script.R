library(keras)

data <- samples_bp # Input Data: [i,5000,2]
name <- "saved_model/lstm_bp1"

split <- 0.7
train_size <- round(nrow(data)*0.7)
train_rows <- sample(nrow(data), train_size)

train_samples <- data[train_rows,,]
test_samples <- data[-train_rows,,]

inputs <- layer_input(shape = c(5000, 1), dtype = 'float32')

outputs <-
  inputs |>
  layer_normalization() |>
  layer_lstm(units = 200, return_sequences = 'True', activation = 'tanh') |>
  layer_dense(units = 4, activation = 'sigmoid', name = 'predictions')

model <- keras_model(inputs = inputs, outputs = outputs, name = 'mdl')

model |>
  compile(
    optimizer = 'adam',
    loss = 'sparse_categorical_crossentropy',
    metrics = 'accuracy'
  )

history <- model |> fit(train_samples[,,1], train_samples[,,2], epochs = 25, verbose = 2)

save_model_tf(model, name)
#load_model_tf("saved_model/")