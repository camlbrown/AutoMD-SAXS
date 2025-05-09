import pandas as pd
import matplotlib.pyplot as plt

# Read the data from the file, skipping the first 5 lines
data = pd.read_csv('Rg_distr_001_1.txt', delim_whitespace=True, skiprows=5, header=None)

# Assign column names for clarity
data.columns = ['Rg', 'Pool_freq', 'Sel_freq']

# Extract the columns for plotting
Rg_values = data['Rg']  # Radius of gyration (Rg)
md_pool = data['Pool_freq']  # MD pool frequency
selected_pool = data['Sel_freq']  # SAXS-optimised pool frequency

# Plotting
plt.plot(Rg_values, md_pool, color='grey', label='MD pool')
plt.plot(Rg_values, selected_pool, color='blue', label='SAXS-optimised pool')
plt.xlabel('Rg (\u00C5)')
plt.ylabel('Frequency')
plt.legend()
plt.savefig('Rg_dist.png', dpi=700, bbox_inches='tight')
plt.savefig('Rg_dist.pdf')

