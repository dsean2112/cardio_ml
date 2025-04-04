
output_name <- "spectro_ii"
lead_name <- "ii"
input_data <- samples_spec_ii[train_rows, , ]


lead_list <-
  c("i",
    "ii",
    "iii",
    "avr",
    "avl",
    "avf",
    "v1",
    "v2",
    "v3",
    "v4",
    "v5",
    "v6")
lead_index <- which(lead_list == lead_name)
labeled_data <- annotations12[train_rows, , lead_index]
  
  #**Handle split prior to function to train all 12 leads consistently 
  # Therefore, all input_data is used for training
  
  library(keras)
  
  inputs <- layer_input(shape = c(5000, 24), dtype = 'float32')
  
  outputs <-
    inputs |>
    layer_normalization() |>
    layer_lstm(units = 200, return_sequences = 'True', activation = 'tanh') |>
    layer_dense(units = 4, activation = 'sigmoid', name = 'predictions')
  # units = 200
  model <- keras_model(inputs = inputs, outputs = outputs, name = 'mdl')
  
  model |>
    compile(
      optimizer = 'adam',
      loss = 'sparse_categorical_crossentropy',
      metrics = 'accuracy'
    )
  
  history <- model |> fit(input_data, labeled_data, epochs = 25, verbose = 2)
  
  model_name <- paste0("saved_model/",output_name)
  save_model_tf(model, model_name)
  # return(history)
  
  # heatmap(foo,Rowv=NA,Colv=NA)