import ROOT
import matplotlib.pyplot as plt

# Load Delphes libraries
ROOT.gSystem.Load("libDelphes")


# Path to the Delphes ROOT file
file_path = "aa_ww_semi_leptonic_SM.root"

# Open the ROOT file
root_file = ROOT.TFile.Open(file_path)
if not root_file or root_file.IsZombie():
    raise FileNotFoundError(f"Cannot open ROOT file: {file_path}")

# Access the Delphes tree
tree = root_file.Get("Delphes")
if not tree:
    raise ValueError("Tree 'Delphes' not found in the ROOT file.")

# Create a histogram for lepton pT
h_pt = ROOT.TH1F("h_pt", "Transverse Momentum (p_{T}) of Leptons; p_{T} [GeV]; Entries", 50, 0, 200)

# Loop through the tree
for entry in tree:
    # Access electrons and fill their pT
    electrons = getattr(entry, "Electron")
    for i in range(electrons.GetEntries()):
        electron = electrons.At(i)
        h_pt.Fill(electron.PT)

    # Access muons and fill their pT
    muons = getattr(entry, "Muon")
    for i in range(muons.GetEntries()):
        muon = muons.At(i)
        h_pt.Fill(muon.PT)

# Convert the ROOT histogram to matplotlib for visualization
x_vals = [h_pt.GetBinCenter(i) for i in range(1, h_pt.GetNbinsX() + 1)]
y_vals = [h_pt.GetBinContent(i) for i in range(1, h_pt.GetNbinsX() + 1)]

# Plot the histogram using matplotlib
plt.bar(x_vals, y_vals, width=h_pt.GetBinWidth(1), color='blue', alpha=0.7, label="Leptons")
plt.xlabel("Transverse Momentum (p_{T}) [GeV]")
plt.ylabel("Entries")
plt.title("pT Distribution of Leptons")
plt.legend()
plt.grid()
plt.show()
