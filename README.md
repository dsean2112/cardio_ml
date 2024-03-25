# File Sequence for Training: 
1. load LUDB_samples
2. Run NN_band_pass if desired
3. Run NN_train_script

# For Analysis File/s: 
1. LUDB Analysis

Loading Models: 
library(keras)
model <- load_model_tf("saved_model/model_lstm2")

# Models
model_lstm1: created using unfiltered data
lstm_bp1: created using a bandpass filter (cut-offs of 0.5 - 40 hz)
