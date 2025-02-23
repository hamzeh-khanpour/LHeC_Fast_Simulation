
# Hamzeh Khanpour -- 2025

import matplotlib.pyplot as plt
import numpy as np
import os  # ✅ Added for file existence checking

# Matplotlib configuration for publication-quality plots
import mplhep as hep


hep.style.use("CMS")
#plt.style.use(hep.style.ROOT)

'''plt.rcParams["axes.linewidth"] = 1.8
plt.rcParams["xtick.major.width"] = 1.8
plt.rcParams["xtick.minor.width"] = 1.8
plt.rcParams["ytick.major.width"] = 1.8
plt.rcParams["ytick.minor.width"] = 1.8

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

plt.rcParams["xtick.labelsize"] = 15
plt.rcParams["ytick.labelsize"] = 15

plt.rcParams["legend.fontsize"] = 15

plt.rcParams['legend.title_fontsize'] = 'x-large' '''



#======================================================================
#======================================================================


# ✅ Function to parse the LHE file and extract kinematic distributions
def parse_lhe_file(file_name):
    if not os.path.exists(file_name):
        print(f"❌❌❌❌ Error: File {file_name} not found!")
        return [], [], [], [], [], []

    pt_tau_plus = []
    eta_tau_plus = []
    pt_tau_minus = []
    eta_tau_minus = []
    rapidity_tau_pair = []
    invariant_mass_tau_pair = []

    with open(file_name, "r") as file:
        in_event = False
        tau_plus = None
        tau_minus = None

        for line in file:
            line = line.strip()
            if "<event>" in line:
                in_event = True
                tau_plus = None
                tau_minus = None
                continue
            if "</event>" in line:
                in_event = False
                if tau_plus and tau_minus:
                    px_pair = tau_plus["px"] + tau_minus["px"]
                    py_pair = tau_plus["py"] + tau_minus["py"]
                    pz_pair = tau_plus["pz"] + tau_minus["pz"]
                    energy_pair = tau_plus["energy"] + tau_minus["energy"]

                    if abs(energy_pair - abs(pz_pair)) > 1e-6:
                        rapidity = 0.5 * np.log((energy_pair + pz_pair) / (energy_pair - pz_pair))
                        rapidity_tau_pair.append(rapidity)

                    # ✅ More robust invariant mass calculation
                    mass_squared = energy_pair**2 - (px_pair**2 + py_pair**2 + pz_pair**2)
                    invariant_mass_tau_pair.append(np.sqrt(max(0, mass_squared)))  # Ensures no NaN values

                continue

            if in_event:
                if not line[0].isdigit() and not line[0] == '-':
                    continue

                parts = line.split()
                if len(parts) < 10:
                    continue

                try:
                    pdg_id = int(parts[0])
                    px = float(parts[6])
                    py = float(parts[7])
                    pz = float(parts[8])
                    energy = float(parts[9])
                    pt = np.sqrt(px**2 + py**2)

                    if abs(energy - abs(pz)) > 1e-6:
                        eta = 0.5 * np.log((energy + pz) / (energy - pz))

                        if pdg_id == 15:
                            pt_tau_plus.append(pt)
                            eta_tau_plus.append(eta)
                            tau_plus = {"px": px, "py": py, "pz": pz, "energy": energy}

                        if pdg_id == -15:
                            pt_tau_minus.append(pt)
                            eta_tau_minus.append(eta)
                            tau_minus = {"px": px, "py": py, "pz": pz, "energy": energy}

                except ValueError as e:
                    print(f"Error parsing line: {line}, error: {e}")
                    continue


    return pt_tau_plus, eta_tau_plus, pt_tau_minus, eta_tau_minus, rapidity_tau_pair, invariant_mass_tau_pair



#======================================================================
#======================================================================

# ✅ Function to load cross-section data from a txt file
def load_cross_section(file_name):
    if not os.path.exists(file_name):
        print(f"❌❌❌❌ Error: File {file_name} not found!")
        return [], []

    W_values = []
    Elastic_xsec = []

    with open(file_name, "r") as file:
        for line in file:
            if line.startswith("#"):  # Skip header lines
                continue

            parts = line.split()
            if len(parts) < 3:
                continue

            try:
                W_values.append(float(parts[0]))  # W [GeV]
                Elastic_xsec.append(float(parts[1]))  # Elastic cross-section [pb]
            except ValueError as e:
                print(f"Error parsing line: {line}, error: {e}")
                continue

    # ✅ Ensure data is sorted by W_values
    sorted_indices = np.argsort(W_values)
    W_values = np.array(W_values)[sorted_indices]
    Elastic_xsec = np.array(Elastic_xsec)[sorted_indices]

    return W_values, Elastic_xsec



#======================================================================


# ✅ Function to load rapidity cross-section data from a txt file
def load_rapidity_cross_section(file_name):
    if not os.path.exists(file_name):
        print(f"❌❌❌❌ Error: File {file_name} not found!")
        return [], []

    Yll_values = []
    Elastic_xsec = []

    with open(file_name, "r") as file:
        for line in file:
            if line.startswith("#"):  # Skip header lines
                continue

            parts = line.split()
            if len(parts) < 3:
                print(f"⚠️⚠️⚠️⚠️ Warning: Skipping malformed line: {line.strip()}")
                continue

            try:
                Yll_values.append(float(parts[0]))  # Yll_MPL (rapidity)
                Elastic_xsec.append(float(parts[1]))  # Elastic cross-section [pb]
            except ValueError as e:
                print(f"⚠️⚠️⚠️⚠️ Error parsing line: {line.strip()}, error: {e}")
                continue

    # ✅ Ensure data is sorted by Yll values for proper plotting
    sorted_indices = np.argsort(Yll_values)
    Yll_values = np.array(Yll_values)[sorted_indices]
    Elastic_xsec = np.array(Elastic_xsec)[sorted_indices]

    return Yll_values, Elastic_xsec



#======================================================================



# ✅ Function to plot weighted distributions with cross-section overlay
def plot_weighted_distribution_with_cross_section(data1, data2, cross_section_x, cross_section_y,
                                                  bins, range, color1, color2, color3, xlabel, ylabel, title, filename,
                                                  label1, label2, label3, integrated_cross_section_SM,
                                                  integrated_cross_section_a_tau, integrated_luminosity, log_scale=False):
    num_entries1 = len(data1)
    num_entries2 = len(data2)

    if num_entries1 == 0 or num_entries2 == 0:
        print(f"⚠️⚠️⚠️⚠️ Warning: No entries found for {label1} or {label2}. Skipping plot.")
        return

    bin_width = (range[1] - range[0]) / bins

    event_weight1 = (integrated_cross_section_SM * integrated_luminosity) / max(num_entries1, 1)
    event_weight2 = (integrated_cross_section_a_tau * integrated_luminosity) / max(num_entries2, 1)

    plt.hist(data1, bins=bins, range=range, color=color1, alpha=0.6, label=f"{label1}",
             weights=[event_weight1 / bin_width] * num_entries1, edgecolor="red", histtype="step", linewidth=2)

    plt.hist(data2, bins=bins, range=range, color=color2, alpha=0.6, label=f"{label2}",
             weights=[event_weight2 / bin_width] * num_entries2, edgecolor="blue", histtype="step", linestyle="dashed", linewidth=2)

    # ✅ Overlay elastic cross-section
    plt.plot(cross_section_x, cross_section_y, color=color3, linestyle="dashdot", linewidth=2, label=label3)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True, linestyle="--")
    plt.tight_layout()

    if log_scale:
        plt.xscale("log")
        plt.yscale("log")
        plt.ylim(1e-3, 20)
        plt.xlim(range[0], range[1])

    plt.savefig(filename, dpi=300)
    plt.show()



#======================================================================


# ✅ Function to plot weighted distributions with rapidity cross-section overlay
def plot_weighted_distribution_with_rapidity(data1, data2, cross_section_x, cross_section_y,
                                             bins, range, color1, color2, color3, xlabel, ylabel, title, filename,
                                             label1, label2, label3, integrated_cross_section_SM,
                                             integrated_cross_section_a_tau, integrated_luminosity, log_scale=False):
    num_entries1 = len(data1)
    num_entries2 = len(data2)

    if num_entries1 == 0 or num_entries2 == 0:
        print(f"⚠️⚠️⚠️⚠️ Warning: No entries found for {label1} or {label2}. Skipping plot.")
        return




    bin_width = (range[1] - range[0]) / bins

    event_weight1 = (integrated_cross_section_SM * integrated_luminosity) / max(num_entries1, 1)
    event_weight2 = (integrated_cross_section_a_tau * integrated_luminosity) / max(num_entries2, 1)



    plt.hist(data1, bins=bins, range=range, color=color1, alpha=0.6, label=f"{label1}",
             weights=[event_weight1 / bin_width] * num_entries1, edgecolor="red", histtype="step", linewidth=2)

    plt.hist(data2, bins=bins, range=range, color=color2, alpha=0.6, label=f"{label2}",
             weights=[event_weight2 / bin_width] * num_entries2, edgecolor="blue", histtype="step", linestyle="dashed", linewidth=2)

    # ✅ Overlay elastic cross-section
    plt.plot(cross_section_x, cross_section_y, color=color3, linestyle="dashdot", linewidth=2, label=label3)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True, linestyle="--")
    plt.tight_layout()

    if log_scale:
        plt.xscale("log")
        plt.yscale("log")
        plt.ylim(1e-3, 20)
        plt.xlim(range[0], range[1])

    plt.savefig(filename, dpi=300)
    plt.show()



#======================================================================



# ✅ Parse both files


file_name_SM = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_tautau_SM/merged_aa_tautau_SM.lhe" # 47.27


#file_name_a_tau = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO/merged_aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO_at_0_001_1TeV.lhe" # 50.975
#file_name_a_tau = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO/merged_aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO_at_0_001_2TeV.lhe" # 53.419


#file_name_a_tau = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO/merged_aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO_at_0_004_2TeV.lhe" # 53.419
file_name_a_tau = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO/merged_aa_tautau_SM_NP_2_SMEFTsim_top_alphaScheme_UFO_at_0_004_1TeV.lhe" # 53.419


cross_section_file = "cross_section_results.txt"
rapidity_cross_section_file = "Yll_elas_inel_data.txt"



pt_tau_plus_1, eta_tau_plus_1, pt_tau_minus_1, eta_tau_minus_1, rapidity_tau_pair_1, invariant_mass_tau_pair_1 = parse_lhe_file(file_name_SM)
pt_tau_plus_2, eta_tau_plus_2, pt_tau_minus_2, eta_tau_minus_2, rapidity_tau_pair_2, invariant_mass_tau_pair_2 = parse_lhe_file(file_name_a_tau)

# ✅ Load cross-section data
W_values, Elastic_xsec = load_cross_section(cross_section_file)

# ✅ Load rapidity cross-section data
Yll_values, Elastic_xsec_rapidity = load_rapidity_cross_section(rapidity_cross_section_file)


# ✅ Define integrated luminosity and separate cross-sections for each file

integrated_luminosity = 1.0  # fb^-1

integrated_cross_section_SM    = 47.27  # pb
integrated_cross_section_a_tau = 53.419  # pb





#======================================================================


# ✅ Plotting with correct notation for invariant mass
fig, ax = plt.subplots(figsize=(8.0, 9.0))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot comparisons with correct weighting + cross-section overlay (invariant mass)
plot_weighted_distribution_with_cross_section(invariant_mass_tau_pair_1, invariant_mass_tau_pair_2, W_values, Elastic_xsec,
                                              bins=490, range=(10, 500), color1="blue", color2="red", color3="teal",  # ✅ Changed from black for better contrast
                                              xlabel=r"$M_{\tau^+ \tau^-} \ \mathrm{[GeV]}$",
                                              ylabel=r"$d\sigma/dM_{\tau^+ \tau^-} \quad \mathrm{[pb/GeV]}$",
                                              title="LHeC @ 1.2 TeV", filename="Invariant_mass_tau_pair_SM_atau_EPA.jpg",
                                              label1=r"$\tau^+ \tau^-$ (SM)", label2=r"$\tau^+ \tau^- (a_{\tau} = 0.0042)$",
                                              label3=r"EPA",
                                              integrated_cross_section_SM=integrated_cross_section_SM,
                                              integrated_cross_section_a_tau=integrated_cross_section_a_tau,
                                              integrated_luminosity=integrated_luminosity, log_scale=True)






# ✅ Plotting with correct notation for invariant mass
fig, ax = plt.subplots(figsize=(8.0, 9.0))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot rapidity distribution with correct weighting + overlay
plot_weighted_distribution_with_rapidity(rapidity_tau_pair_1, rapidity_tau_pair_2, Yll_values, Elastic_xsec_rapidity,
                                         bins=20, range=(-10, 10), color1="purple", color2="green", color3="teal",  # ✅ Changed from black for better contrast
                                         xlabel=r"$Y_{\tau^+ \tau^-}$",
                                         ylabel=r"$d\sigma/dY_{\tau^+ \tau^-} \quad \mathrm{[pb]}$",  # ✅ Used subscript instead of superscript for consistency
                                         title="LHeC @ 1.2 TeV", filename="Rapidity_tau_pair_SM_atau_EPA.jpg",
                                         label1=r"$\tau^+ \tau^-$ (SM)", label2=r"$\tau^+ \tau^- (a_{\tau} = 0.0042)$",
                                         label3=r"EPA",
                                         integrated_cross_section_SM=integrated_cross_section_SM,
                                         integrated_cross_section_a_tau=integrated_cross_section_a_tau,
                                         integrated_luminosity=integrated_luminosity)





#======================================================================




# ✅ Function to plot the ratio of two distributions with statistical uncertainties
def plot_ratio(data1, data2, bins, range, xlabel, ylabel, title, filename):
    if len(data1) == 0 or len(data2) == 0:
        print("⚠️⚠️⚠️⚠️ Warning: One or both datasets are empty. Cannot compute ratio.")
        return

    # ✅ Compute histograms
    hist1, bin_edges = np.histogram(data1, bins=bins, range=range)
    hist2, _ = np.histogram(data2, bins=bins, range=range)

    # ✅ Compute ratio while avoiding division by zero
    ratio = np.divide(hist2, hist1, out=np.full_like(hist1, np.nan, dtype=float), where=hist1 != 0)

    # ✅ Compute bin centers
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # ✅ Compute Poisson statistical uncertainties
    err1 = np.sqrt(hist1)
    err2 = np.sqrt(hist2)
    ratio_err = ratio * np.sqrt((err1 / hist1) ** 2 + (err2 / hist2) ** 2)  # Error propagation

    # ✅ Replace NaN values (from division by zero) with np.nan instead of 1.0
    valid_idx = hist1 != 0
    ratio[~valid_idx] = np.nan

    # ✅ Interpolate missing values in ratio for smooth visualization
    if np.sum(valid_idx) > 1:
        ratio = np.interp(bin_centers, bin_centers[valid_idx], ratio[valid_idx])

    # ✅ Plot ratio with error bars
    fig, ax = plt.subplots(figsize=(14, 5))  # CMS-style wide ratio plot
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.25, top=0.95)  # More space at bottom

    ax.errorbar(bin_centers, ratio, yerr=ratio_err, fmt="o", color="blue", linewidth=3, markersize=8, markerfacecolor="magenta", markeredgecolor="red", label="Ratio $(a_{\\tau}/SM)$")

    # ✅ Add shaded uncertainty band centered on ratio
#    ax.fill_between(bin_centers, ratio - ratio_err, ratio + ratio_err, color="blue", alpha=0.3, label="Uncertainty Band")

    # ✅ Reference line at ratio = 1
    ax.axhline(y=1, color="green", linestyle="--", linewidth=2, label="y=1")

    # ✅ Format plot
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)


    # ✅ Adjust y-limits dynamically
#    ymin = max(0, np.nanpercentile(ratio, 5))  # Lower 5% percentile
#    ymax = np.nanpercentile(ratio, 95) * 1.5   # Upper 95% percentile

    ymin =   -0.1
    ymax =    2.0

    ax.set_ylim(ymin, ymax)

    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()

    # ✅ Save plot
    plt.savefig(filename, dpi=300)
    plt.show()


# ✅ Define parameters
bins = 10
range_limits = (10, 500)  # Adjust based on data
xlabel = r"$M_{\tau^+ \tau^-} \ \mathrm{[GeV]}$"
ylabel = "Ratio $(a_{\\tau}/SM)$"
title = "LHeC @ 1.2 TeV"
output_filename = "Ratio_NP_SM_Plot.jpg"

# ✅ Call function to plot ratio
plot_ratio(invariant_mass_tau_pair_1, invariant_mass_tau_pair_2, bins, range_limits, xlabel, ylabel, title, output_filename)




#======================================================================




# ✅ Function to plot the ratio of two distributions with statistical uncertainties
def plot_ratio(data1, data2, bins, range, xlabel, ylabel, title, filename):
    if len(data1) == 0 or len(data2) == 0:
        print("⚠️⚠️⚠️⚠️ Warning: One or both datasets are empty. Cannot compute ratio.")
        return

    # ✅ Compute histograms
    hist1, bin_edges = np.histogram(data1, bins=bins, range=range)
    hist2, _ = np.histogram(data2, bins=bins, range=range)

    # ✅ Compute ratio while avoiding division by zero
    ratio = np.divide(hist2, hist1, out=np.full_like(hist1, np.nan, dtype=float), where=hist1 != 0)

    # ✅ Compute complementary ratio (1 - ratio)
    complementary_ratio = 1 - ratio

    # ✅ Compute bin centers
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # ✅ Compute Poisson statistical uncertainties
    err1 = np.sqrt(hist1)
    err2 = np.sqrt(hist2)
    ratio_err = ratio * np.sqrt((err1 / hist1) ** 2 + (err2 / hist2) ** 2)  # Error propagation
    complementary_ratio_err = ratio_err  # Same error for complementary ratio

    # ✅ Replace NaN values (from division by zero) with np.nan instead of 1.0
    valid_idx = hist1 != 0
    ratio[~valid_idx] = np.nan
    complementary_ratio[~valid_idx] = np.nan

    # ✅ Plot ratio with error bars
    fig, ax = plt.subplots(figsize=(14, 5))  # CMS-style wide ratio plot
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.25, top=0.95)  # More space at bottom

#    ax.errorbar(bin_centers, ratio, yerr=ratio_err, fmt="o", color="blue", markersize=5, label="Ratio $(a_{\\tau}/SM)$")
    ax.errorbar(bin_centers, complementary_ratio, yerr=complementary_ratio_err, fmt="o", color="blue", linewidth=3, markersize=8, markerfacecolor="magenta", markeredgecolor="red", label="1 - $(a_{\\tau}/SM)$")

    # ✅ Add shaded uncertainty band centered on ratio
#    ax.fill_between(bin_centers, ratio - ratio_err, ratio + ratio_err, color="blue", alpha=0.3, label="Uncertainty Band")
#    ax.fill_between(bin_centers, complementary_ratio - complementary_ratio_err, complementary_ratio + complementary_ratio_err, color="green", alpha=0.3, label="Uncertainty Band (1 - Ratio)")

    # ✅ Reference line at ratio = 1
#    ax.axhline(y=1, color="red", linestyle="--", linewidth=2, label="Ref: Ratio = 1")
    ax.axhline(y=0, color="green", linestyle="--", linewidth=2, label="y=0")

    # ✅ Format plot
    ax.set_xlabel(xlabel)
#    ax.set_ylabel(ylabel + " & $1 - Ratio$")
    ax.set_ylabel("1 - $(a_{\\tau}/SM)$")
    ax.set_title(title)

    # ✅ Adjust y-limits dynamically
#    ymin = min(np.nanmin(ratio), np.nanmin(complementary_ratio)) * 0.8
#    ymax = max(np.nanmax(ratio), np.nanmax(complementary_ratio)) * 1.2

    ymin = -1.1
    ymax =  1.1

    ax.set_ylim(ymin, ymax)

    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()

    # ✅ Save plot
    plt.savefig(filename, dpi=300)
    plt.show()


# ✅ Define parameters
bins = 10
range_limits = (10, 500)  # Adjust based on data
xlabel = r"$M_{\tau^+ \tau^-} \ \mathrm{[GeV]}$"
ylabel = "Ratio $(a_{\\tau}/SM)$"
title = "LHeC @ 1.2 TeV"
output_filename = "Ratio_1_NP_SM_Plot.jpg"

# ✅ Call function to plot ratio
plot_ratio(invariant_mass_tau_pair_1, invariant_mass_tau_pair_2, bins, range_limits, xlabel, ylabel, title, output_filename)



#======================================================================





# ✅ Function to plot the ratio of two distributions with statistical uncertainties
def plot_ratio(data1, data2, bins, range, xlabel, ylabel, title, filename):
    if len(data1) == 0 or len(data2) == 0:
        print("⚠️ Warning: One or both datasets are empty. Cannot compute ratio.")
        return

    # ✅ Compute histograms
    hist1, bin_edges = np.histogram(data1, bins=bins, range=range)
    hist2, _ = np.histogram(data2, bins=bins, range=range)

    # ✅ Compute ratio while avoiding division by zero
    ratio = np.divide(hist2, hist1, out=np.full_like(hist1, np.nan, dtype=float), where=hist1 != 0)

    # ✅ Compute bin centers
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # ✅ Compute Poisson statistical uncertainties
    err1 = np.sqrt(hist1)
    err2 = np.sqrt(hist2)
    ratio_err = ratio * np.sqrt((err1 / hist1) ** 2 + (err2 / hist2) ** 2)  # Error propagation

    # ✅ Handle NaN values (from division by zero)
    valid_idx = hist1 != 0
    ratio[~valid_idx] = np.nan

    # ✅ Plot ratio with error bars
    fig, ax = plt.subplots(figsize=(14, 5))  # CMS-style wide ratio plot
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.25, top=0.95)  # More space at bottom

    # ✅ Shaded uncertainty band around ratio = 1
#    ax.fill_between(bin_centers, 1 - ratio_err, 1 + ratio_err, color="gray", alpha=0.3, label="Uncertainty Band")

    # ✅ Plot **continuous ratio line**
    ax.plot(bin_centers, ratio, drawstyle="steps-mid", linestyle="-", linewidth=3, label="Ratio $(a_{\\tau}/SM)$")

    # ✅ Overlay statistical error bars
#    ax.errorbar(bin_centers, ratio, yerr=ratio_err, fmt="o", color="black", markersize=5, label="Stat. Uncertainty")

    # ✅ Reference line at ratio = 1
    ax.axhline(y=1, color="green", linestyle="--", linewidth=2, label="y=1")

    # ✅ Format plot
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)


    # ✅ Adjust y-limits dynamically

    ymin =   -0.1
    ymax =    2.0

    ax.set_ylim(ymin, ymax)



    ax.legend(loc="upper right")
    ax.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()

    # ✅ Save plot
    plt.savefig(filename, dpi=300)
    plt.show()


# ✅ Define parameters
# ✅ Define parameters
bins = 10
range_limits = (10, 500)  # Adjust based on data
xlabel = r"$M_{\tau^+ \tau^-} \ \mathrm{[GeV]}$"
ylabel = "Ratio $(a_{\\tau}/SM)$"
title = "LHeC @ 1.2 TeV"
output_filename = "Ratio_Obs_Exp.jpg"

# ✅ Call function to plot ratio
plot_ratio(invariant_mass_tau_pair_1, invariant_mass_tau_pair_2, bins, range_limits, xlabel, ylabel, title, output_filename)




#======================================================================

# ✅ Function to plot the ratio of two distributions with statistical uncertainties and effective luminosity
def plot_ratio_with_luminosity(data1, data2, bins, range_limits, xlabel, ylabel, filename,
                               sigma_bsm, sigma_sm, n_events_bsm, n_events_sm):
    if len(data1) == 0 or len(data2) == 0:
        print("⚠️ Warning: One or both datasets are empty. Cannot compute ratio.")
        return

    # ✅ Compute histograms
    hist1, bin_edges = np.histogram(data1, bins=bins, range=range_limits)
    hist2, _ = np.histogram(data2, bins=bins, range=range_limits)

    # ✅ Compute ratio while avoiding division by zero
    ratio = np.divide(hist2, hist1, out=np.full_like(hist1, np.nan, dtype=float), where=hist1 != 0)

    # ✅ Compute bin centers
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # ✅ Compute Poisson statistical uncertainties
    err1 = np.sqrt(hist1)
    err2 = np.sqrt(hist2)
    ratio_err = ratio * np.sqrt((err1 / hist1) ** 2 + (err2 / hist2) ** 2)  # Error propagation

    # ✅ Compute effective luminosity for BSM and SM
    luminosity_bsm = n_events_bsm / sigma_bsm
    luminosity_sm = n_events_sm / sigma_sm

    # ✅ Compute the average effective luminosity
    avg_luminosity = (luminosity_bsm + luminosity_sm) / 2

    # ✅ Handle NaN values (from division by zero)
    valid_idx = hist1 != 0
    ratio[~valid_idx] = np.nan

    # ✅ Plot ratio with error bars
    fig, ax = plt.subplots(figsize=(14, 5))  # CMS-style wide ratio plot
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.25, top=0.95)  # More space at bottom

    # ✅ Plot **continuous ratio line**
    ax.plot(bin_centers, ratio, drawstyle="steps-mid", linestyle="-", linewidth=3, label="Ratio $(a_{\\tau}/SM)$")

    # ✅ Overlay statistical error bars
    ax.errorbar(bin_centers, ratio, yerr=ratio_err, fmt="o", color="red", markersize=8, label="Stat. Uncertainty")

    # ✅ Reference line at ratio = 1
    ax.axhline(y=1, color="green", linestyle="--", linewidth=2, label="y=1")

    # ✅ Format plot
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"LHeC @ 1.2 TeV ($a_t=0.0042, {{\cal L}} = 100.0$ fb$^{{-1}}$)")


    # ✅ Adjust y-limits dynamically
    ymin = -0.1
    ymax = 2.0
    ax.set_ylim(ymin, ymax)

    ax.legend(loc="upper right")
    ax.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()

    # ✅ Save plot
    plt.savefig(filename, dpi=300)
    plt.show()


# ✅ Given values for cross-sections and event counts
sigma_bsm = 53.419 * 1000  # Convert to fb
sigma_sm = 47.27 * 1000  # Convert to fb
n_events_bsm = 5000000  # Number of BSM events
n_events_sm = 5000000  # Number of SM events


# ✅ Define parameters
bins = 10
range_limits = (10, 500)  # Adjust based on data
xlabel = r"$M_{\tau^+ \tau^-} \ \mathrm{[GeV]}$"
ylabel = "Ratio $(a_{\\tau}/SM)$"
output_filename = "Ratio_Obs_Exp_with_Luminosity.jpg"


# ✅ Call function to plot ratio
plot_ratio_with_luminosity(invariant_mass_tau_pair_1, invariant_mass_tau_pair_2, bins, range_limits,
                           xlabel, ylabel, output_filename, sigma_bsm, sigma_sm,
                           n_events_bsm, n_events_sm)




#======================================================================


from scipy.optimize import curve_fit

# ✅ Linear function for fitting
def linear_fit(x, a, b):
    return a * x + b

# ✅ Function to plot the ratio of two distributions with statistical uncertainties and effective luminosity
def plot_ratio_with_luminosity(data1, data2, bins, range_limits, xlabel, ylabel, filename,
                               sigma_bsm, sigma_sm, n_events_bsm, n_events_sm):
    if len(data1) == 0 or len(data2) == 0:
        print("⚠️ Warning: One or both datasets are empty. Cannot compute ratio.")
        return

    # ✅ Compute histograms
    hist1, bin_edges = np.histogram(data1, bins=bins, range=range_limits)
    hist2, _ = np.histogram(data2, bins=bins, range=range_limits)

    # ✅ Compute ratio while avoiding division by zero
    ratio = np.divide(hist2, hist1, out=np.full_like(hist1, np.nan, dtype=float), where=hist1 != 0)

    # ✅ Compute bin centers
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # ✅ Compute Poisson statistical uncertainties
    err1 = np.sqrt(hist1)
    err2 = np.sqrt(hist2)
    ratio_err = ratio * np.sqrt((err1 / hist1) ** 2 + (err2 / hist2) ** 2)  # Error propagation

    # ✅ Compute effective luminosity for BSM and SM
    luminosity_bsm = n_events_bsm / sigma_bsm
    luminosity_sm = n_events_sm / sigma_sm

    # ✅ Compute the average effective luminosity
    avg_luminosity = (luminosity_bsm + luminosity_sm) / 2

    # ✅ Handle NaN values (from division by zero)
    valid_idx = hist1 != 0
    ratio[~valid_idx] = np.nan

    # ✅ Perform a linear fit to the ratio
    popt, pcov = curve_fit(linear_fit, bin_centers[valid_idx], ratio[valid_idx], sigma=ratio_err[valid_idx])
    slope, intercept = popt
    slope_err = np.sqrt(pcov[0, 0])
    intercept_err = np.sqrt(pcov[1, 1])

    # ✅ Plot ratio with error bars
    fig, ax = plt.subplots(figsize=(14, 5))  # CMS-style wide ratio plot
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.25, top=0.95)  # More space at bottom

    # ✅ Plot **continuous ratio line**
    ax.plot(bin_centers, ratio, drawstyle="steps-mid", linestyle="-", linewidth=3, label="Ratio $(a_{\\tau}/SM)$")

    # ✅ Overlay statistical error bars
    ax.errorbar(bin_centers, ratio, yerr=ratio_err, fmt="o", color="red", markersize=8, label="Stat. Uncertainty")

    # ✅ Plot linear fit line
    fit_x = np.linspace(range_limits[0], range_limits[1], 100)
    fit_y = linear_fit(fit_x, *popt)
    ax.plot(fit_x, fit_y, linestyle="--", color="magenta", linewidth=2, label=f"Fit: $y = ({slope:.4f} \pm {slope_err:.4f})x + ({intercept:.4f} \pm {intercept_err:.4f})$")

    # ✅ Reference line at ratio = 1
    ax.axhline(y=1, color="green", linestyle="--", linewidth=2, label="y=1")

    # ✅ Format plot
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"LHeC @ 1.2 TeV ($a_t=0.0042, {{\\cal L}} = {avg_luminosity:.2e}$ fb$^{{-1}}$)")

    # ✅ Adjust y-limits dynamically
    ymin = -0.1
    ymax = 2.0
    ax.set_ylim(ymin, ymax)

    ax.legend(loc="upper right")
    ax.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()

    # ✅ Save plot
    plt.savefig(filename, dpi=300)
    plt.show()

# ✅ Given values for cross-sections and event counts
sigma_bsm = 53.419 * 1000  # Convert to fb
sigma_sm = 47.27 * 1000  # Convert to fb
n_events_bsm = 5000000  # Number of BSM events
n_events_sm = 5000000  # Number of SM events

# ✅ Define parameters
bins = 10
range_limits = (10, 500)  # Adjust based on data
xlabel = r"$M_{\tau^+ \tau^-} \ \mathrm{[GeV]}$"
ylabel = "Ratio $(a_{\tau}/SM)$"
output_filename = "Ratio_Obs_Exp_with_Luminosity_Fit.jpg"

# ✅ Call function to plot ratio
plot_ratio_with_luminosity(invariant_mass_tau_pair_1, invariant_mass_tau_pair_2, bins, range_limits,
                           xlabel, ylabel, output_filename, sigma_bsm, sigma_sm,
                           n_events_bsm, n_events_sm)




#======================================================================


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ✅ Linear function for fitting
def linear_fit(x, a, b):
    return a * x + b

# ✅ Function to plot the ratio of two distributions with statistical uncertainties and effective luminosity
def plot_ratio_with_luminosity(data1, data2, bins, range_limits, xlabel, ylabel, filename,
                               sigma_bsm, sigma_sm, n_events_bsm, n_events_sm):
    if len(data1) == 0 or len(data2) == 0:
        print("⚠️ Warning: One or both datasets are empty. Cannot compute ratio.")
        return

    # ✅ Compute histograms
    hist1, bin_edges = np.histogram(data1, bins=bins, range=range_limits)
    hist2, _ = np.histogram(data2, bins=bins, range=range_limits)

    # ✅ Compute ratio while avoiding division by zero
    ratio = np.divide(hist2, hist1, out=np.full_like(hist1, np.nan, dtype=float), where=hist1 != 0)

    # ✅ Compute bin centers
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # ✅ Compute Poisson statistical uncertainties
    err1 = np.sqrt(hist1)
    err2 = np.sqrt(hist2)
    ratio_err = ratio * np.sqrt((err1 / hist1) ** 2 + (err2 / hist2) ** 2)  # Error propagation

    # ✅ Compute effective luminosity for BSM and SM
    luminosity_bsm = n_events_bsm / sigma_bsm
    luminosity_sm = n_events_sm / sigma_sm

    # ✅ Compute the average effective luminosity
    avg_luminosity = (luminosity_bsm + luminosity_sm) / 2

    # ✅ Handle NaN values (from division by zero)
    valid_idx = np.isfinite(ratio)  # Ensures only valid values are used

    # ✅ Perform a linear fit to the ratio
    popt, pcov = curve_fit(linear_fit, bin_centers[valid_idx], ratio[valid_idx], sigma=ratio_err[valid_idx])
    slope, intercept = popt
    slope_err = np.sqrt(pcov[0, 0])
    intercept_err = np.sqrt(pcov[1, 1])

    # ✅ Compute chi-squared statistic
    chi2 = np.sum(((ratio[valid_idx] - linear_fit(bin_centers[valid_idx], *popt)) / ratio_err[valid_idx])**2)
    dof = len(valid_idx) - 2  # Degrees of freedom (data points - fit parameters)
    reduced_chi2 = chi2 / dof

    # ✅ Plot ratio with error bars
    fig, ax = plt.subplots(figsize=(14, 5))  # CMS-style wide ratio plot
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.25, top=0.95)  # More space at bottom

    # ✅ Plot **continuous ratio line**
    ax.plot(bin_centers, ratio, drawstyle="steps-mid", linestyle="-", linewidth=3, label="Ratio $(a_{\\tau}/SM)$")

    # ✅ Overlay statistical error bars
    ax.errorbar(bin_centers, ratio, yerr=ratio_err, fmt="o", color="red", markersize=8, label="Stat. Uncertainty")

    # ✅ Plot linear fit line
    fit_x = np.linspace(range_limits[0], range_limits[1], 100)
    fit_y = linear_fit(fit_x, *popt)
    ax.plot(fit_x, fit_y, linestyle="--", color="magenta", linewidth=2, label=f"Fit: $y = ({slope:.4f} \pm {slope_err:.4f})x + ({intercept:.4f} \pm {intercept_err:.4f})$, $\chi^2/dof = {reduced_chi2:.2f}$")

    # ✅ Reference line at ratio = 1
    ax.axhline(y=1, color="green", linestyle="--", linewidth=2, label="y=1")

    # ✅ Format plot
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"LHeC @ 1.2 TeV ($a_t=0.0042, {{\\cal L}} = 100$ fb$^{{-1}}$)")

    # ✅ Adjust y-limits dynamically
    ymin = -0.1
    ymax = 2.0
    ax.set_ylim(ymin, ymax)

    ax.legend(loc="upper right")
    ax.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()

    # ✅ Save plot
    plt.savefig(filename, dpi=300)
    plt.show()

    # ✅ Return the fitted parameters and chi-squared value for further use
    return slope, intercept, reduced_chi2

# ✅ Given values for cross-sections and event counts
sigma_bsm = 53.419 * 1000  # Convert to fb
sigma_sm = 47.27 * 1000  # Convert to fb
n_events_bsm = 5000000  # Number of BSM events
n_events_sm = 5000000  # Number of SM events

# ✅ Define parameters
bins = 10
range_limits = (10, 500)  # Adjust based on data
xlabel = r"$M_{\tau^+ \tau^-} \ \mathrm{[GeV]}$"
ylabel = "Ratio $(a_{\tau}/SM)$"
output_filename = "Ratio_Obs_Exp_with_Luminosity_Fit_Final.jpg"

# ✅ Call function to plot ratio
slope, intercept, reduced_chi2 = plot_ratio_with_luminosity(invariant_mass_tau_pair_1, invariant_mass_tau_pair_2, bins, range_limits,
                           xlabel, ylabel, output_filename, sigma_bsm, sigma_sm,
                           n_events_bsm, n_events_sm)

print(f"Fit Results: Slope = {slope:.4f}, Intercept = {intercept:.4f}, Chi2/dof = {reduced_chi2:.2f}")





#======================================================================





