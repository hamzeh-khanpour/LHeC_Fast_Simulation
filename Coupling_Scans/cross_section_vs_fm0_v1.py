import numpy as np
import matplotlib.pyplot as plt

# Load data from the cross_sections.txt file
file_name = "cross_sections_FM0.txt"
data = []

# Read the data file
with open(file_name, "r") as file:
    for line in file.readlines()[1:]:  # Skip the header
        parts = line.strip().split()
        if len(parts) >= 2:
            fm0_value = float(parts[0])
            try:
                cross_section = float(parts[1])
            except ValueError:
                cross_section = None  # Handle missing/invalid values
            data.append((fm0_value, cross_section))

# Filter out rows with missing cross-section values
data = [(fm0, cs) for fm0, cs in data if cs is not None]

# Separate FM0 values and cross-sections
fm0_values, cross_sections = zip(*data)

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(fm0_values, cross_sections, 'o-', label='Cross-section')

# Add labels, grid, and title
plt.xscale('linear')
plt.yscale('log')
plt.xlabel('FM0 Value', fontsize=14)
plt.ylabel('Cross-section (pb)', fontsize=14)
plt.title('Cross-section vs FM0 Value', fontsize=16)
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.axhline(y=0.1, color='r', linestyle='--', label='Threshold Line (0.1 pb)')
plt.legend()

# Save the plot
plt.savefig('cross_section_vs_fm0.jpg', dpi=300)
plt.show()
