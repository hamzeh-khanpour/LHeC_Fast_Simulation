import matplotlib.pyplot as plt
import numpy as np

# Function to parse the LHE file and extract Pt, eta of tau+, tau-, and rapidity of tau+ tau- pair
def parse_lhe_file(file_name):
    pt_tau_plus = []
    eta_tau_plus = []
    pt_tau_minus = []
    eta_tau_minus = []
    rapidity_tau_pair = []  # To store rapidity of tau+ tau- pair

    with open(file_name, "r") as file:
        in_event = False
        tau_plus = None
        tau_minus = None

        for line in file:
            line = line.strip()
            # Check if we are inside an <event> block
            if "<event>" in line:
                in_event = True
                tau_plus = None
                tau_minus = None
                continue
            if "</event>" in line:
                in_event = False

                # Calculate rapidity of the tau+ tau- pair if both are present
                if tau_plus and tau_minus:
                    px_pair = tau_plus["px"] + tau_minus["px"]
                    py_pair = tau_plus["py"] + tau_minus["py"]
                    pz_pair = tau_plus["pz"] + tau_minus["pz"]
                    energy_pair = tau_plus["energy"] + tau_minus["energy"]

                    if abs(energy_pair - abs(pz_pair)) > 1e-6:  # Avoid division by zero
                        rapidity = 0.5 * np.log((energy_pair + pz_pair) / (energy_pair - pz_pair))
                        rapidity_tau_pair.append(rapidity)

                continue

            if in_event:
                # Skip lines starting with non-numeric characters (metadata lines)
                if not line[0].isdigit() and not line[0] == '-':  # Include negative PDG IDs
                    continue

                # Split the line into components
                parts = line.split()
                # Check if the line has enough elements to extract particle data
                if len(parts) < 10:
                    continue

                try:
                    pdg_id = int(parts[0])
                    px = float(parts[6])
                    py = float(parts[7])
                    pz = float(parts[8])
                    energy = float(parts[9])

                    # Calculate transverse momentum Pt = sqrt(px^2 + py^2)
                    pt = (px**2 + py**2)**0.5

                    # Calculate pseudorapidity eta = 0.5 * ln((E + pz) / (E - pz)), avoid division by zero
                    if abs(energy - abs(pz)) > 1e-6:  # Add small epsilon to avoid precision issues
                        eta = 0.5 * np.log((energy + pz) / (energy - pz))

                        # Check if the particle is tau+ (PDG ID = 15)
                        if pdg_id == 15:
                            pt_tau_plus.append(pt)
                            eta_tau_plus.append(eta)
                            tau_plus = {"px": px, "py": py, "pz": pz, "energy": energy}

                        # Check if the particle is tau- (PDG ID = -15)
                        if pdg_id == -15:
                            pt_tau_minus.append(pt)
                            eta_tau_minus.append(eta)
                            tau_minus = {"px": px, "py": py, "pz": pz, "energy": energy}

                except ValueError as e:
                    print(f"Error parsing line: {line}, error: {e}")
                    continue

    return pt_tau_plus, eta_tau_plus, pt_tau_minus, eta_tau_minus, rapidity_tau_pair

# Function to streamline plotting
def plot_distribution(data, bins, range, color, xlabel, ylabel, title, filename, label, y_range=None):
    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=bins, range=range, color=color, alpha=0.7, label=label, edgecolor='black')
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.title(title, fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, linestyle="--", alpha=0.6)
    if y_range:
        plt.ylim(y_range)  # Set y-axis range if provided
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.show()


# Parse the file to get Pt, eta of tau+ and tau-, and rapidity of tau+ tau- pair
file_name = "unweighted_events.lhe"
pt_tau_plus, eta_tau_plus, pt_tau_minus, eta_tau_minus, rapidity_tau_pair = parse_lhe_file(file_name)


# Plot the Pt and eta distributions for tau+ and tau-
plot_distribution(pt_tau_plus, bins=50, range=(0, 100), color="blue",
                  xlabel=r"$P_T$ of $\tau^+$ (GeV)", ylabel="Number of Events",
                  title=r"$P_T$ Distribution of $\tau^+$", filename="pt_tau_plus_distribution.png",
                  label=r"$P_T$ of $\tau^+$", y_range=(0, 10000))

plot_distribution(eta_tau_plus, bins=50, range=(-10, 10), color="green",
                  xlabel=r"$\eta$ of $\tau^+$", ylabel="Number of Events",
                  title=r"$\eta$ Distribution of $\tau^+$", filename="eta_tau_plus_distribution.png",
                  label=r"$\eta$ of $\tau^+$", y_range=(0, 10000))

plot_distribution(pt_tau_minus, bins=50, range=(0, 100), color="red",
                  xlabel=r"$P_T$ of $\tau^-$ (GeV)", ylabel="Number of Events",
                  title=r"$P_T$ Distribution of $\tau^-$", filename="pt_tau_minus_distribution.png",
                  label=r"$P_T$ of $\tau^-$", y_range=(0, 10000))

plot_distribution(eta_tau_minus, bins=50, range=(-10, 10), color="orange",
                  xlabel=r"$\eta$ of $\tau^-$", ylabel="Number of Events",
                  title=r"$\eta$ Distribution of $\tau^-$", filename="eta_tau_minus_distribution.png",
                  label=r"$\eta$ of $\tau^-$", y_range=(0, 10000))


# Plot the rapidity distribution for tau+ tau- pair
plot_distribution(rapidity_tau_pair, bins=20, range=(-10, 10), color="purple",
                  xlabel=r"$Y_{\tau^+ \tau^-}$", ylabel="Number of Events",
                  title=r"Rapidity Distribution of $\tau^+ \tau^-$ Pair", filename="rapidity_tau_pair_distribution.png",
                  label=r"$Y_{\tau^+ \tau^-}$", y_range=(0, 10000))


