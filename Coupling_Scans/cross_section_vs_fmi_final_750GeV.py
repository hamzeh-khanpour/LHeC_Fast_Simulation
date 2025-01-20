import numpy as np
import matplotlib.pyplot as plt

# Matplotlib configuration for publication-quality plots
#import mplhep as hep

#hep.style.use("CMS")
#plt.style.use(hep.style.ROOT)

plt.rcParams["axes.linewidth"] = 1.8
plt.rcParams["xtick.major.width"] = 1.8
plt.rcParams["xtick.minor.width"] = 1.8
plt.rcParams["ytick.major.width"] = 1.8
plt.rcParams["ytick.minor.width"] = 1.8

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

plt.rcParams["xtick.labelsize"] = 15
plt.rcParams["ytick.labelsize"] = 15

plt.rcParams["legend.fontsize"] = 15

plt.rcParams['legend.title_fontsize'] = 'x-large'


# Define file names and labels for FM values
file_names = {
    "FM0": "cross_sections_FM0_750GeV.txt",
    "FM1": "cross_sections_FM1_750GeV.txt",
    "FM2": "cross_sections_FM2_750GeV.txt",
    "FM3": "cross_sections_FM3_750GeV.txt"
}
labels = {
    "FM0": "$F_{M_0} / \Lambda^4$ (LHeC@0.750 TeV)",
    "FM1": "$F_{M_1} / \Lambda^4$ (LHeC@0.750 TeV)",
    "FM2": "$F_{M_2} / \Lambda^4$ (LHeC@0.750 TeV)",
    "FM3": "$F_{M_3} / \Lambda^4$ (LHeC@0.750 TeV)"
}
colors = {
    "FM0": "red",
    "FM1": "orange",
    "FM2": "green",
    "FM3": "blue"
}

# Function to load data from a file
def load_data(file_name):
    data = []
    with open(file_name, "r") as file:
        for line in file.readlines()[1:]:  # Skip the header
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    fm_value = float(parts[0])
                    cross_section = float(parts[1])
                    data.append((fm_value, cross_section))
                except ValueError:
                    continue  # Skip rows with invalid data
    return zip(*data) if data else ([], [])

# Load data for each FM value
fm_data = {key: load_data(file_name) for key, file_name in file_names.items()}

# Plot the data
plt.figure(figsize=(10, 6))

for key, (fm_values, cross_sections) in fm_data.items():
    plt.plot(fm_values, cross_sections, marker='o', label=labels[key], color=colors[key], linewidth=2)

# Add labels, grid, and title
plt.xscale('linear')
plt.yscale('log')
#plt.ylim(1e-2, 1e3)  # Set y-axis limits

plt.xlabel('Coupling Values: $f_{M_i} / \Lambda^4$ (TeV$^{-4}$)', fontsize=16)
plt.ylabel('$\sigma$ (pb)  LHeC@0.750 TeV', fontsize=16)

plt.title('The total cross-section of the process $e^- p \\to e^- W^+ W^- p$ as a function of $f_{M_i} / \Lambda^4$', fontsize=16)

plt.grid(True, which="both", linestyle="--", linewidth=0.5)
#plt.axhline(y=0.1, color='r', linestyle='--', label='Threshold Line (0.1 pb)')
plt.legend()

# Save the plot
plt.savefig('cross_section_vs_coupling_750GeV.jpg', dpi=300)
plt.show()
