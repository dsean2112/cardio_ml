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

# Misc
LUDB_samples.RData:  
Variable name "samples" in R. Raw Samples.  
First dimension: sample number (1-170).  
Second dimension: individual values over time (1-5000)  
Third dimension:  
   Value 1: mV values  
   Value 2: annotation values: 0 = 0, P = 1, QRS = 2, T = 3  
