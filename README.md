# File Sequence for Training: 
1. load LUDB_samples
2. Run NN_band_pass if desired
3. Run NN_train_script

# For Analysis File/s: 
1. LUDB Analysis  
   accuracy: simple binary yes/no if predictions match reference values  
   confusion_table: confusion table of all samples predicted  
   confusion_table_by_sample: creates matrix of confusion tables for each predicted sample  
   confusion_by_sample_condensed: for each sample, gives accuracy of predicting 0, P, QRS and T in respective columns  

# Models  
## Loading Models:  
library(keras)  
model <- load_model_tf("saved_model/model_lstm2")  

## List:  
model_lstm1: created using unfiltered data  
lstm_bp1: created using a bandpass filter (cut-offs of 0.5 - 40 hz)  
lstm_i, ii, iii etc: created using unfiltered data, for all 12 leads  

# Data:
samples12: 170 samples for all 12 leads. mV values.  
annotations12: 170 samples for all 12 leads. Annotations: 0 = 0, P = 1, QRS = 2, T = 3  
train_rows_12: Rows of samples12/annotations12 used to train the lstm_i, ii etc. Can use this to isolate testing samples  

samples_bp12: 170 samples for all 12 leads. mV values. Filtered through bandpass Butterworth filter, exclude <5 Hz, >40 Hz

LUDB_samples.RData: Variable name "samples" in R. Raw Samples of lead i only. 170 samples, 5000 time points. 
Third dimension:  
Value 1: mV values  
Value 2: annotation values: 0 = 0, P = 1, QRS = 2, T = 3  

# Misc
