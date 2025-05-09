import numpy as np
import matplotlib.pyplot as plt

# Load data from profiles_001_1.fit
file_path = 'profiles_001_1.fit'

# Use genfromtxt to load the data, skipping the first line
data = np.genfromtxt(file_path, skip_header=1, dtype=float)

# Extract columns
q_values = data[:, 0]  # First column: q (scattering angle)
iq_values = data[:, 1]  # Second column: I(q)
iq_errors = data[:, 2]  # Third column: Error
ensemble_fit = data[:, 3]  # Fourth column: Ensemble fit

# Filter out invalid values (zeros or negative values) before taking the logarithm
valid_indices = (iq_values > 0) & (ensemble_fit > 0)  # Avoid log of zero or negative numbers
q_values = q_values[valid_indices]
iq_values = iq_values[valid_indices]
iq_errors = iq_errors[valid_indices]
ensemble_fit = ensemble_fit[valid_indices]

# Calculate the range of log(I(q)) for the y-axis
log_iq_min = np.log(np.min(iq_values)) - 0.1  # Add a small buffer below the min log value
log_iq_max = np.log(np.max(iq_values)) + 0.1  # Add a small buffer above the max log value

# Define colors
exp_color = 'grey'  # Grey for experimental intensity
fit_color = 'blue'  # Blue for ensemble fit
errorbar_color = 'lightsteelblue'  # Very light blue for error bars

# Plotting log(I(q)) vs q
fig = plt.figure(figsize=(10, 6))
gs = fig.add_gridspec(2, 1, height_ratios=[5, 2], hspace=0.15)

# Subplot for log(I(q)) vs q
ax0 = fig.add_subplot(gs[0])
ax0.errorbar(
    q_values,
    np.log(iq_values),
    yerr=iq_errors / iq_values,  # Normalize errors for log scale
    fmt='o',
    markersize=3,
    label='Experimental',
    color=exp_color,
    ecolor=errorbar_color,
    capsize=2,
    alpha=0.7,
    zorder=1,
)
ax0.plot(
    q_values,
    np.log(ensemble_fit),
    color=fit_color,
    label='Ensemble Fit',
    linewidth=2,
    alpha=0.8,
    zorder=2,
)

# Set axis labels and limits
ax0.set_ylabel('log(I(q))', fontsize=14)
ax0.set_ylim(log_iq_min, log_iq_max)  # Adjust y-axis limits based on log data
ax0.legend()
ax0.grid(True, linestyle='-', alpha=0.5)

# Subplot for residuals
ax1 = fig.add_subplot(gs[1], sharex=ax0)
residuals = (iq_values - ensemble_fit) / iq_errors  # Error-weighted residuals
ax1.plot(
    q_values,
    residuals,
    'o',
    color=fit_color,  # Same color as the ensemble fit
    markersize=3,
#    label='Residuals',
)
ax1.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.7)  # Horizontal line at y=0

# Set axis labels
ax1.set_xlabel('q (\u00C5\u207B\u00B9)', fontsize=14)
ax1.set_ylabel('\u0394/\u03C3', fontsize=14)
ax1.grid(True, linestyle='-', alpha=0.5)

# Save the plot
plt.savefig('ensemble_saxs_fit.png', bbox_inches='tight', dpi=700)
plt.savefig('ensemble_saxs_fit.pdf')

