# File Sequence for Training: 
1. load LUDB_samples
2. Run NN_band_pass
3. Run NN_train_script

# For Analysis File/s: 
1. LUDB Analysis

Loading Models: 
library(keras)
model <- load_model_tf("saved_model/model_lstm2")
