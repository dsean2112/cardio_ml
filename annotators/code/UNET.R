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

# Need to pad for UNET model
pad <- 8
half_pad <- array(0,c(nrow(testing_annotations),pad/2))
testing_annotations <- cbind(half_pad,testing_annotations,half_pad)
testing_signal <- cbind(half_pad,testing_signal,half_pad)

half_pad <- array(0,c(nrow(training_annotations),pad/2))
training_annotations <- cbind(half_pad,training_annotations,half_pad)
training_signal <- cbind(half_pad,training_signal,half_pad)

# Build model -------------------------------------------------------------
# INPUTS:
filters = 16 # typically 8, 16 or 32
dropout = 0.3 # typically between 0,3 - 0.5


library(keras)
rm(model)

# Find the number of classes in the annotation matrix (either 4 or 10 depending on annotator style)
num_classes <- length(unique(as.vector(training_annotations)))

input_shape = c(5008,1)
conv_kernal <- 9
deconv_kernal <- 8

inputs <- keras::layer_input(shape = input_shape)
normalized <- layer_batch_normalization(object = inputs)

conv1 = keras::layer_conv_1d(object = normalized, filters, conv_kernal, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
conv1 = keras::layer_conv_1d(object = conv1, filters, conv_kernal, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
pool1 = keras::layer_max_pooling_1d(conv1, pool_size = 2, strides = 2)
# pool1 = keras::layer_max_pooling_1d(conv1, pool_size = 2)

conv2 = keras::layer_conv_1d(object = pool1, filters*2, conv_kernal, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
conv2 = keras::layer_conv_1d(object = conv2, filters*2, conv_kernal, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
pool2 = keras::layer_max_pooling_1d(conv2, pool_size = 2, strides = 2)

conv3 = keras::layer_conv_1d(object = pool2, filters*4, conv_kernal, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
conv3 = keras::layer_conv_1d(object = conv3, filters*4, conv_kernal, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
pool3 = keras::layer_max_pooling_1d(conv3, pool_size = 2, strides = 2)

conv4 = keras::layer_conv_1d(object = pool3, filters*8, conv_kernal, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
conv4 = keras::layer_conv_1d(object = conv4, filters*8, conv_kernal, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
drop4 = keras::layer_dropout(conv4, rate = dropout)
pool4 = keras::layer_max_pooling_1d(drop4, pool_size = 2, strides = 2)

conv5 = keras::layer_conv_1d(object = pool4, filters*16, conv_kernal, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
conv5 = keras::layer_conv_1d(object = conv5, filters*16, conv_kernal, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
drop5 = keras::layer_dropout(conv5, rate = dropout)

upsample1 = keras::layer_upsampling_1d(object = drop5, size = 2)
up6 = keras::layer_conv_1d(object = upsample1, filters*8, 2, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
merge6 = keras::layer_concatenate(list(drop4,up6))
conv6 = keras::layer_conv_1d(object = merge6, filters*8, deconv_kernal, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
conv6 = keras::layer_conv_1d(object = conv6, filters*8, deconv_kernal, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')

upsample2 = keras::layer_upsampling_1d(object = conv6, size = 2)
up7 = keras::layer_conv_1d(object = upsample2, filters*4, 2, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
merge7 = keras::layer_concatenate(list(conv3,up7))
conv7 = keras::layer_conv_1d(object = merge7, filters*4, deconv_kernal, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
conv7 = keras::layer_conv_1d(object = conv7, filters*4, deconv_kernal, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')

upsample3 = keras::layer_upsampling_1d(object = conv7, size = 2)
up8 = keras::layer_conv_1d(object = upsample3, filters*2, 2, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
merge8 = keras::layer_concatenate(list(conv2,up8))
conv8 = keras::layer_conv_1d(object = merge8, filters*2, deconv_kernal, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
conv8 = keras::layer_conv_1d(object = conv8, filters*2, deconv_kernal, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')

upsample4 = keras::layer_upsampling_1d(object = conv8, size = 2)
up9 = keras::layer_conv_1d(object = upsample4, filters, 2, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
merge9 = keras::layer_concatenate(list(conv1,up9))
conv9 = keras::layer_conv_1d(object = merge9, filters, deconv_kernal, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
conv9 = keras::layer_conv_1d(object = conv9, filters, deconv_kernal, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')
conv9 = keras::layer_conv_1d(object = conv9, num_classes, 1, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')

outputs = layer_dense(object = conv9, units = num_classes, activation = 'sigmoid', name = 'predictions')
# outputs = keras::layer_conv_1d(object = conv9, 4, 1, activation = 'sigmoid')

model = keras::keras_model(inputs = inputs, outputs = outputs, name = 'UNet')

model |>
  compile(
    optimizer = 'rmsprop',
    loss = 'sparse_categorical_crossentropy',
    metrics = 'accuracy'
  )

rm(list = ls(pattern = 'conv[1-9]$|merge[1-9]|up[1-9]|upsample|outputs$|^inputs$|drop[1-9]|pool[1-9]|^normalized$'))


# Train -------------------------------------------------------------------
epochs <- 10 # ~5s per step
history <- model |> fit(
  training_signal,
  training_annotations,
  epochs = epochs,
  # validation_data = list(val_inputs, val_targets),
  validation_data = list(testing_signal, testing_annotations),
  # loss_weights = c(0.2, 1, 1, 1), # experimenting with decreasing the weight of '0' annotations (ie no wave)
  verbose = 2
  ) 


# Test --------------------------------------------------------------------
predictions <- model %>% predict(testing_signal)

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
