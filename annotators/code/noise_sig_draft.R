# Example ECG signal (replace with your actual data)
ecg_voltage <- testing_signal[sample,]

# Define the standard deviation of the Gaussian noise
sigma <- 50  # Adjust this value as needed

# Generate Gaussian random noise
noise <- rnorm(length(ecg_voltage), mean = 0, sd = sigma)

# Add the noise to the original ECG signal
noisy_ecg <- ecg_voltage + noise

# do filtering on singal snipped section, without padding?

# Plot --------------------------------------------------------------------

truth_plot <- plot_func(ecg_voltage,ann_wfdb2continuous2(testing_annotations[[sample]]))
predicted_plot <- plot_func(noisy_ecg,ann_wfdb2continuous2(testing_annotations[[sample]]))
filtered_plot <- plot_func(ecg_filter(noisy_ecg),ann_wfdb2continuous2(testing_annotations[[sample]]))
subplot(truth_plot,predicted_plot,filtered_plot,nrows = 3)
