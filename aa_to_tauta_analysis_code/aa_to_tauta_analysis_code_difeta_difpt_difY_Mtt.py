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




# ✅ Function to parse the LHE file and extract kinematic distributions
def parse_lhe_file(file_name):
    if not os.path.exists(file_name):
        print(f"❌ Error: File {file_name} not found!")
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




# ✅ Function to plot weighted distributions (with correct separate cross-sections)
def plot_weighted_distribution(data1, data2, bins, range, color1, color2, xlabel, ylabel, title, filename, label1, label2,
                               integrated_cross_section_SM, integrated_cross_section_a_tau, integrated_luminosity, log_scale=False):
    num_entries1 = len(data1)
    num_entries2 = len(data2)

    if num_entries1 == 0 or num_entries2 == 0:
        print(f"⚠️ Warning: No entries found for {label1} or {label2}. Skipping plot.")
        return

    bin_width = (range[1] - range[0]) / bins
    event_weight1 = (integrated_cross_section_SM * integrated_luminosity) / num_entries1 if num_entries1 > 0 else 0
    event_weight2 = (integrated_cross_section_a_tau * integrated_luminosity) / num_entries2 if num_entries2 > 0 else 0


    plt.hist(data1, bins=bins, range=range, color=color1, alpha=0.6, label=f"{label1}",
             weights=[event_weight1] * num_entries1, edgecolor="red", histtype="step", linewidth=2)


    plt.hist(data2, bins=bins, range=range, color=color2, alpha=0.6, label=f"{label2}",
             weights=[event_weight2] * num_entries2, edgecolor="blue", histtype="step", linestyle="dashed", linewidth=2)


    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True, linestyle="--")
    plt.tight_layout()

    if log_scale:
        plt.xscale("log")
        plt.yscale("log")
        plt.ylim(1e-4, 1e2)
        plt.xlim(range[0], range[1])

    plt.savefig(filename, dpi=300)
    plt.show()



# ✅ Parse both files
file_name_SM = "/home/hamzeh-khanpour/aa_tautau_SM.lhe"
file_name_a_tau = "/home/hamzeh-khanpour/aa_tautau_a_tau.lhe"


pt_tau_plus_1, eta_tau_plus_1, pt_tau_minus_1, eta_tau_minus_1, rapidity_tau_pair_1, invariant_mass_tau_pair_1 = parse_lhe_file(file_name_SM)
pt_tau_plus_2, eta_tau_plus_2, pt_tau_minus_2, eta_tau_minus_2, rapidity_tau_pair_2, invariant_mass_tau_pair_2 = parse_lhe_file(file_name_a_tau)




# ✅ Define integrated luminosity and separate cross-sections for each file
integrated_luminosity = 1.0  # fb^-1
integrated_cross_section_SM = 47.27  # pb (for first file)
integrated_cross_section_a_tau = 50.72  # pb (for second file)




# Plotting
fig, ax = plt.subplots(figsize=(8.0, 9.0))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot comparisons with correct weighting
plot_weighted_distribution(rapidity_tau_pair_1, rapidity_tau_pair_2, bins=20, range=(-10, 10), color1="purple", color2="green",
                           xlabel=r"$Y_{\tau^+ \tau^-}$",
#                           ylabel=r"$d\sigma_{ep \to e(\gamma\gamma \to \tau^+ \tau^-)p^{(*)}} / dY^{\tau^+ \tau^-} \quad \mathrm{[pb]}$",
                           ylabel=r"$d\sigma/dY^{\tau^+ \tau^-} \quad \mathrm{[pb]}$",
                           title="LHeC@1.2 TeV", filename="rapidity_tau_pair_comparison.png",
                           label1=r"$\tau^+ \tau^-$ (SM)", label2=r"$\tau^+ \tau^- (a_{\tau} = 0.0042)$",
                           integrated_cross_section_SM=integrated_cross_section_SM,
                           integrated_cross_section_a_tau=integrated_cross_section_a_tau,
                           integrated_luminosity=integrated_luminosity)


# ✅ Plotting with correct notation for invariant mass
fig, ax = plt.subplots(figsize=(8.0, 9.0))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plot_weighted_distribution(invariant_mass_tau_pair_1, invariant_mass_tau_pair_2, bins=400, range=(10, 500), color1="blue", color2="red",
                           xlabel=r"$M_{\tau^+ \tau^-} \ \mathrm{[GeV]}$",
#                           ylabel=r"$d\sigma_{ep \to e(\gamma\gamma \to \tau^+ \tau^-)p^{(*)}} / dM_{\tau^+ \tau^-} \quad \mathrm{[pb/GeV]}$",
                           ylabel=r"$d\sigma/dM_{\tau^+ \tau^-} \quad \mathrm{[pb/GeV]}$",
                           title="LHeC@1.2 TeV", filename="invariant_mass_tau_pair_comparison.png",
                           label1=r"$\tau^+ \tau^-$ (SM)", label2=r"$\tau^+ \tau^- (a_{\tau} = 0.0042)$",
                           integrated_cross_section_SM=integrated_cross_section_SM,
                           integrated_cross_section_a_tau=integrated_cross_section_a_tau,
                           integrated_luminosity=integrated_luminosity, log_scale=True)








