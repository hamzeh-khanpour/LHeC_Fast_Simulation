


import ROOT
import matplotlib.pyplot as plt


# Matplotlib configuration for publication-quality plots
import mplhep as hep

hep.style.use("CMS")
#plt.style.use(hep.style.ROOT)

#plt.rcParams["axes.linewidth"] = 1.8
#plt.rcParams["xtick.major.width"] = 1.8
#plt.rcParams["xtick.minor.width"] = 1.8
#plt.rcParams["ytick.major.width"] = 1.8
#plt.rcParams["ytick.minor.width"] = 1.8

#plt.rcParams["xtick.direction"] = "in"
#plt.rcParams["ytick.direction"] = "in"

#plt.rcParams["xtick.labelsize"] = 15
#plt.rcParams["ytick.labelsize"] = 15

#plt.rcParams["legend.fontsize"] = 15

#plt.rcParams['legend.title_fontsize'] = 'x-large'  


# Function to parse the LHE file and extract lepton transverse momentum (pT)
def parse_lhe_file(file_name):
    pt_leptons = []

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
                if not line[0].isdigit() and not line[0] == "-":  # Include negative PDG IDs
                    continue

                # Split the line into components
                parts = line.split()
                # Check if the line has enough elements to extract particle data
                if len(parts) < 10:
                    continue

                pdg_id = int(parts[0])
                px = float(parts[6])
                py = float(parts[7])
                pz = float(parts[8])
                energy = float(parts[9])

                # Check if the particle is a lepton (e- or e+)
                if abs(pdg_id) in [11, 13, 15]:  # e, mu, tau
                    # Create TLorentzVector
                    lepton = ROOT.TLorentzVector()
                    lepton.SetPxPyPzE(px, py, pz, energy)

                    # Extract transverse momentum
                    pt_leptons.append(lepton.Pt())

    return pt_leptons


    fig, ax = plt.subplots(figsize = (12.0, 10.0))
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)


# Parse signal and background files
signal_file = "unweighted_events_signal.lhe"
background_file = "unweighted_events_background.lhe"

pt_leptons_signal = parse_lhe_file(signal_file)
pt_leptons_background = parse_lhe_file(background_file)

# Plot the transverse momentum distributions
# plt.figure(figsize=(8, 9))

plt.hist(pt_leptons_signal, bins=50, range=(0, 400), alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$]", color="blue", edgecolor="black")
plt.hist(pt_leptons_background, bins=50, range=(0, 400), alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="red", edgecolor="black")
plt.xlabel(r"$p_T^{\ell=\tau, \mu, e}$ [GeV]")
plt.ylabel("# of Events")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("pt_leptons_comparison.png", dpi=300)
plt.show()
