import numpy as np
import matplotlib.pyplot as plt

# Load data from the cross_sections_FM0.txt file
file_name_fm0 = "cross_sections_FM0.txt"
data_fm0 = []

# Read the FM0 data file
with open(file_name_fm0, "r") as file:
    for line in file.readlines()[1:]:  # Skip the header
        parts = line.strip().split()
        if len(parts) >= 2:
            fm0_value = float(parts[0])
            try:
                cross_section = float(parts[1])
            except ValueError:
                cross_section = None  # Handle missing/invalid values
            data_fm0.append((fm0_value, cross_section))

# Filter out rows with missing cross-section values for FM0
data_fm0 = [(fm0, cs) for fm0, cs in data_fm0 if cs is not None]

# Separate FM0 values and cross-sections
fm0_values, cross_sections_fm0 = zip(*data_fm0)

# Load data from the cross_sections_FM1.txt file
file_name_fm1 = "cross_sections_FM1.txt"
data_fm1 = []

# Read the FM1 data file
with open(file_name_fm1, "r") as file:
    for line in file.readlines()[1:]:  # Skip the header
        parts = line.strip().split()
        if len(parts) >= 2:
            fm1_value = float(parts[0])
            try:
                cross_section = float(parts[1])
            except ValueError:
                cross_section = None  # Handle missing/invalid values
            data_fm1.append((fm1_value, cross_section))

# Filter out rows with missing cross-section values for FM1
data_fm1 = [(fm1, cs) for fm1, cs in data_fm1 if cs is not None]

# Separate FM1 values and cross-sections
fm1_values, cross_sections_fm1 = zip(*data_fm1)

# Load data from the cross_sections_FM2.txt file
file_name_fm2 = "cross_sections_FM2.txt"
data_fm2 = []

# Read the FM2 data file
with open(file_name_fm2, "r") as file:
    for line in file.readlines()[1:]:  # Skip the header
        parts = line.strip().split()
        if len(parts) >= 2:
            fm2_value = float(parts[0])
            try:
                cross_section = float(parts[1])
            except ValueError:
                cross_section = None  # Handle missing/invalid values
            data_fm2.append((fm2_value, cross_section))

# Filter out rows with missing cross-section values for FM2
data_fm2 = [(fm2, cs) for fm2, cs in data_fm2 if cs is not None]

# Separate FM2 values and cross-sections
fm2_values, cross_sections_fm2 = zip(*data_fm2)

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(fm0_values, cross_sections_fm0, 'o-', label='$F_{M_0} / \Lambda^4$ (LHeC@1.2 TeV)', color='red', linewidth=2)
plt.plot(fm1_values, cross_sections_fm1, 's-', label='$F_{M_1} / \Lambda^4$ (LHeC@1.2 TeV)', color='orange', linewidth=2)
plt.plot(fm2_values, cross_sections_fm2, '^-', label='$F_{M_2} / \Lambda^4$ (LHeC@1.2 TeV)', color='green', linewidth=2)

# Add labels, grid, and title
plt.xscale('linear')
plt.yscale('log')

plt.xlabel('Coupling Values: $f_{M_i} / \Lambda^4$ (TeV$^{-4}$)', fontsize=14)
plt.ylabel('$\sigma$ (pb)', fontsize=14)

plt.title('The total cross-section of the process $e^- p \\to e^- W^+ W^- p$ as a function of $f_{M_i} / \Lambda^4$', fontsize=16)

plt.grid(True, which="both", linestyle="--", linewidth=0.5)
#plt.axhline(y=0.1, color='r', linestyle='--', label='Threshold Line (0.1 pb)')
plt.legend()

# Save the plot
plt.savefig('cross_section_vs_coupling.jpg', dpi=300)
plt.show()
