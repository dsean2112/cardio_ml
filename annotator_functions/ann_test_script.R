# Sample script
library(shiva)
options(wfdb_path = 'wsl /usr/local/bin') 
model_path <- "saved_model/"

# test_wfdb <- load("test_wfdb.RData")
# test <- test_wfdb$signal$i

f <- 500   # sample frequency
high <- 40 # low  pass frequency threshold
low <- 0.5 # high pass frequency threshold. Original: 0.5. 1 works well
test_filt <- pass.filt(test, c(f/high, f/low), "pass", method = c("Butterworth"), n = 4, Rp = 0.1)

test_spectro <- build_spectrogram(input = test_filt)
test_ann <- predict_samples(signal = test_spectro, model_name = "spectro_bilstm_i", model_path = model_path)
plot_func(y = test, color = test_ann)

Rpeaks <- R_peak_isolation(test, test_ann) # change input var name to raw_signal
RR_table <- RR_table_make(Rpeaks)

test_spliced_spectro <- splice_signal(RR_table, test_spectro)
test_ann2 <- RR_waveform_prediction(test_spliced_spectro, model_name = "PT_bilstm_spectro_mask_i", model_path = model_path)
test_ann2_stitch <- RR_stitch_ann(test_ann2, RR_table)
plot_func(y = test_filt, color = test_ann2_stitch)
