import matplotlib.pyplot as plt
import numpy as np

# Function to parse the LHE file and extract Pt and eta of tau+
def parse_lhe_file(file_name):
    pt_tau_plus = []
    eta_tau_plus = []

    with open(file_name, "r") as file:
        in_event = False
        for line in file:
            line = line.strip()
            # Check if we are inside an <event> block
            if "<event>" in line:
                in_event = True
                continue
            if "</event>" in line:
                in_event = False
                continue
            
            if in_event:
                # Skip lines starting with non-numeric characters (metadata lines)
                if not line[0].isdigit():
                    continue

                # Split the line into components
                parts = line.split()
                if len(parts) > 0:
                    # Check if the particle is tau+ (PDG ID = 15)
                    pdg_id = int(parts[0])
                    if pdg_id == 15:
                        px = float(parts[6])
                        py = float(parts[7])
                        pz = float(parts[8])
                        energy = float(parts[9])

                        # Calculate transverse momentum Pt = sqrt(px^2 + py^2)
                        pt = (px**2 + py**2)**0.5
                        pt_tau_plus.append(pt)

                        # Calculate pseudorapidity eta = 0.5 * ln((E + pz) / (E - pz))
                        if energy != abs(pz):  # Avoid division by zero
                            eta = 0.5 * np.log((energy + pz) / (energy - pz))
                            eta_tau_plus.append(eta)

    return pt_tau_plus, eta_tau_plus

# Parse the file to get Pt and eta of tau+
file_name = "unweighted_events.lhe"
pt_tau_plus, eta_tau_plus = parse_lhe_file(file_name)

# Plot the Pt distribution
plt.figure(figsize=(10, 6))
plt.hist(pt_tau_plus, bins=50, color="blue", alpha=0.7, label=r"$P_T$ of $\tau^+$")
plt.xlabel(r"$P_T$ of $\tau^+$ (GeV)", fontsize=14)
plt.ylabel("Number of Events", fontsize=14)
plt.title(r"$P_T$ Distribution of $\tau^+$", fontsize=16)
plt.legend(fontsize=12)
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# Save and show the Pt plot
plt.savefig("pt_tau_plus_distribution.png", dpi=300)
plt.show()




# Plot the eta distribution
plt.figure(figsize=(10, 6))
plt.hist(eta_tau_plus, bins=50, color="green", alpha=0.7, label=r"$\eta$ of $\tau^+$")
plt.xlabel(r"$\eta$ of $\tau^+$", fontsize=14)
plt.ylabel("Number of Events", fontsize=14)
plt.title(r"$\eta$ Distribution of $\tau^+$", fontsize=16)
plt.legend(fontsize=12)
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# Save and show the eta plot
plt.savefig("eta_tau_plus_distribution.png", dpi=300)
plt.show()



