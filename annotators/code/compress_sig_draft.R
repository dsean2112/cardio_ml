# Load required package

# Input: signal and annotations

# On input, if first indices are 0 (in signal), remove leading 0s from sig and ann
#     *Note length of leading 0s
# On input, if last  indices are 0 (in signal), remove lagging 0s from sig and ann
#     *Note length of lagging 0s

# Should apply same compression/dilation to sig and ann

# Ensure signal has 0s added on if compressed

# Add correct leading/lagging 0s on

# Output: list of sig and ann


# Should compress/dilate THEN translate across zeros vector

library(signal) # For resampling

# Example ECG voltage vector
ecg_voltage <- annotations # Replace with your ECG data

# Generate a random compression/dilation factor
compression_factor <- runif(1, min = 0.9, max = 1.1) # Random factor between 90% and 110%

# Define new length for the signal based on compression/dilation factor
new_length <- round(length(ecg_voltage) * compression_factor)

# Resample the ECG signal to the new length
compressed_ecg <- resample(ecg_voltage, p = new_length, q = length(ecg_voltage))

# Print the results
cat("Original ECG Length:", length(ecg_voltage), "\n")
cat("New ECG Length (Compressed/Dilated):", length(compressed_ecg), "\n")

# View the original and modified signals
# print(compressed_ecg)


# padding -----------------------------------------------------------------

# Original and compressed signal lengths
original_length <- length(ecg_voltage)
compressed_length <- length(compressed_ecg)

# Calculate the amount of padding required
padding <- original_length - compressed_length

# Add padding (zeros) evenly to both ends, or all to one end if preferred
if (padding > 0) {
  left_padding <- floor(padding / 2)
  right_padding <- ceiling(padding / 2)
  compressed_ecg <- c(rep(0, left_padding), compressed_ecg, rep(0, right_padding))
}


# testplot ----------------------------------------------------------------

truth_plot <- plot_func(annotations)
predicted_plot <- plot_func(ecg_voltage)

subplot(truth_plot,predicted_plot,nrows = 2)
