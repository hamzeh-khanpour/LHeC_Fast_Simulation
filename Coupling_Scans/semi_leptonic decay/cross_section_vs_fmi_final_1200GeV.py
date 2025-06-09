import numpy as np
import matplotlib.pyplot as plt

import matplotlib.ticker as ticker

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
#    "FM0": "cross_sections_FM0.txt",
#    "FM1": "cross_sections_FM1.txt",
    "FM2": "cross_sections_FM2.txt",
#    "FM3": "cross_sections_FM3.txt"
}
labels = {
#    "FM0": "$F_{M_0} / \Lambda^4$ (LHeC@1.2 TeV)",
#    "FM1": "$F_{M_1} / \Lambda^4$ (LHeC@1.2 TeV)",
    "FM2": "$F_{M_2} / \Lambda^4$ (LHeC@1.2 TeV)",
#    "FM3": "$F_{M_3} / \Lambda^4$ (LHeC@1.2 TeV)"
}
colors = {
#    "FM0": "red",
#    "FM1": "orange",
    "FM2": "green",
#    "FM3": "blue"
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
plt.figure(figsize=(10, 8))

for key, (fm_values, cross_sections) in fm_data.items():
    # Scale x-axis to 10^-12 units
    fm_values_scaled = [x * 1e12 for x in fm_values]
    plt.plot(fm_values_scaled, cross_sections, marker='o', label=labels[key], color=colors[key], linewidth=2)
# Axis scaling and labeling
plt.xscale('linear')
plt.yscale('log')

plt.xlabel(r'Coupling Values : $f_{M_i} / \Lambda^4$ [TeV$^{-4}$]', fontsize=16)
plt.ylabel(r'$\sigma$ (pb) LHeC@1.2 TeV', fontsize=16)
plt.title(r'Total cross-section : $p\,e^- \to p\,W^+W^-\,e^- \to p\,\ell^\pm \nu_{\ell} j j\,e^-$ as a function of $f_{M_i} / \Lambda^4$', fontsize=16)
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.legend()
# Save and show the plot
plt.savefig('cross_section_vs_FM2_coupling_1200GeV.pdf', dpi=600)
plt.show()
