
# Initialize --------------------------------------------------------------

library(shiva)
library(dplR)
options(wfdb_path = 'wsl /usr/local/bin') 
model_path <- "saved_model/"

# test_wfdb <- load("test_wfdb.RData")
# test <- test_wfdb$signal$i


# Inputs --------------------------------------------------------------

# xl_test <- read_wfdb(record = "00280_hr", record_dir = dir)
# test <- xl_test$signal$III

test <- samples12[1,,2]  # change as needed
model_10sec <- "spectro_bilstm_ii"
model_RR <- "PT_bilstm_spectro_mask_i"


# Calcs -------------------------------------------------------------------
sample <- 1


test_filt <- filter_samples(test)

test_spectro <- build_spectrogram(input = test_filt)
test_ann <- predict_samples(signal = test_spectro, model_name = model_10sec, model_path = model_path)
plot_func(y = test_filt[1,], color = test_ann[1,])

Rpeaks <- peak_isolation(test, test_ann, wave_value = 2) # change input var name to raw_signal
RR_table <- RR_table_make(Rpeaks)

test_spliced_spectro <- splice_signal(RR_table, test_spectro)
test_ann2 <- RR_waveform_prediction(test_spliced_spectro, model_name = model_RR, model_path = model_path)
test_ann2_stitch <- RR_stitch_ann(test_ann2, RR_table)

peaks <- peak_isolation(test_filt[1,], test_ann2_stitch[1,], wave_value = 1)
peaks
test_ann2_stitch[1,unlist(peaks)] <- "z"
# Plot --------------------------------------------------------------------

plot1 <- plot_func(y = test_filt[sample,], color = test_ann2_stitch[sample,])

isoelectric_line <- isoelec_find(signal = test_filt[sample,], annotations = test_ann2_stitch[sample,])

plot2 <- plot1 %>% add_trace(y = ~array(isoelectric_line, length(signal)), mode = 'line')
plot2
