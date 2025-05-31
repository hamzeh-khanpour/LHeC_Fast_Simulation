import pyhf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#===========================
# Settings
#===========================
signal_cross_section_pb = 0.021425  # Cross section used to normalize signal weights [pb]
luminosity_fb = 100.0
bins = np.linspace(0, 1, 21)  # 20 bins over BDT score [0,1]

#===========================
# Load data
#===========================
signal_df = pd.read_csv("signal_after_cut.csv")
background_df = pd.read_csv("background_after_cut.csv")

#===========================
# Histogram BDT scores
#===========================
signal_hist, _ = np.histogram(signal_df["bdt_score"], bins=bins, weights=signal_df["weight"])
background_hist, _ = np.histogram(background_df["bdt_score"], bins=bins, weights=background_df["weight"])

# Avoid zero bins to stabilize likelihood
signal_hist += 1e-6
background_hist += 1e-6

#===========================
# Optional: plot input histograms
#===========================
plt.figure(figsize=(8, 6))
plt.hist(bins[:-1], bins=bins, weights=signal_hist, histtype='stepfilled', label='Signal', alpha=0.6, edgecolor='crimson')
plt.hist(bins[:-1], bins=bins, weights=background_hist, histtype='stepfilled', label='Background', alpha=0.6, edgecolor='navy')
plt.xlabel("BDT Score")
plt.ylabel("Weighted Events")
plt.legend()
plt.title("BDT Score Distribution (Post-Cut)")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("bdt_score_histograms_for_pyhf.pdf")
plt.close()

#===========================
# Build Workspace-compatible Model
#===========================
workspace_spec = {
    "channels": [{
        "name": "bdt_channel",
        "samples": [
            {
                "name": "signal",
                "data": signal_hist.tolist(),
                "modifiers": [{"name": "mu", "type": "normfactor", "data": None}]
            },
            {
                "name": "background",
                "data": background_hist.tolist(),
                "modifiers": [{
                    "name": "bkg_uncertainty",
                    "type": "shapesys",
                    "data": (0.1 * background_hist).tolist()
                }]
            }
        ]
    }],
    "parameters": [],
    "version": "1.0.0"
}

observations = [{"name": "bdt_channel", "data": background_hist.tolist()}]

workspace = pyhf.Workspace({
    "channels": workspace_spec["channels"],
    "observations": observations,
    "measurements": [{
        "name": "measurement",
        "config": {
            "poi": "mu",
            "parameters": []
        }
    }],
    "version": "1.0.0"
})

model = workspace.model(measurement_name="measurement")
data = pyhf.tensorlib.astensor(workspace.data(model))

#===========================
# Run hypothesis test
#===========================
CLs_obs, CLs_exp = pyhf.infer.hypotest(
    1.0,  # signal strength μ = 1
    data,
    model,
    test_stat="qtilde",
    return_expected=True
)

mu95 = pyhf.infer.intervals.upper_limit(model, data)

print("\n🔬 pyhf Inference Results:")
print(f"Observed CLs p-value: {CLs_obs:.3e}")
print(f"Expected 95% CL upper limit on signal strength μ: {mu95:.2f}")

#===========================
# Translate to cross section
#===========================
sigma95_pyhf = mu95 * signal_cross_section_pb
print(f"\n🔒 Corresponding cross section upper limit (95% CL): {sigma95_pyhf:.4f} pb")

#===========================
# Convert to FM2 Wilson coefficient limit
#===========================
sigma_SM = 0.0150  # [pb] from your SM-only simulation
A = (signal_cross_section_pb - sigma_SM) / (100.0**2)  # from FM2 = 100 TeV⁻⁴

if sigma95_pyhf <= sigma_SM:
    print("❌ No constraint on FM2 — σ_95 is below the SM cross section.")
else:
    fm2_limit = np.sqrt((sigma95_pyhf - sigma_SM) / A)
    print(f"✅ 95% CL upper limit on |fM2/Λ⁴|: {fm2_limit:.1f} TeV⁻⁴")
