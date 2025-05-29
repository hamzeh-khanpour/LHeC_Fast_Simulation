import ROOT
import numpy as np
import pandas as pd

# Open ROOT file
file = ROOT.TFile.Open("output_histograms.root")

# Mapping between logical observable names and actual ROOT histogram names
observable_map = {
    "lepton_pt": "hist_lepton_pt",
    "lepton_eta": "hist_lepton_eta",
    "leading_jet_pt": "hist_leading_jet_pt",
    "leading_jet_eta": "hist_leading_jet_eta",
    "subleading_jet_eta": "hist_subleading_jet_eta",
    "missing_et": "hist_missing_et",
    "delta_r": "hist_delta_r",
    "jet_centrality": "hist_jet_centrality",
    "delta_eta_jj": "hist_delta_eta_jj",
    "m_w_hadronic": "hist_m_w_hadronic",
    "m_w_leptonic": "hist_m_w_leptonic"
}

n_samples = 10000

def sample_from_hist(hist, n):
    values = []
    for i in range(hist.GetNbinsX()):
        count = int(hist.GetBinContent(i+1))
        xlow = hist.GetBinLowEdge(i+1)
        xhigh = hist.GetBinLowEdge(i+2)
        if count > 0:
            values.extend(np.random.uniform(xlow, xhigh, size=count))
    if len(values) == 0:
        return np.full(n, np.nan)
    return np.random.choice(values, size=n, replace=True)

# --- Signal ---
signal_data = {}
for obs, hist_root_name in observable_map.items():
    full_hist_name = f"{hist_root_name}_FM2_Lambda4"
    hist = file.Get(f"signal_FM2_Lambda4/{full_hist_name}")
    if hist:
        signal_data[obs] = sample_from_hist(hist, n_samples)
    else:
        print(f"❌ Missing: {full_hist_name} in signal")
        signal_data[obs] = np.full(n_samples, np.nan)
signal_df = pd.DataFrame(signal_data)
signal_df["label"] = 1

# --- Backgrounds ---
background_keys = [
    "aa_ww", "aa_ttbar", "aa_tautau", "aa_mumu",
    "inclusive_ttbar", "single_top", "w_production",
    "z_production", "wwj", "zzj", "wzj"
]

special_suffix = {"wwj", "zzj", "wzj"}
background_data = []

for bkg in background_keys:
    bkg_obs_data = {}
    suffix = "_production" if bkg in special_suffix else ""
    for obs, hist_root_name in observable_map.items():
        full_hist_name = f"{hist_root_name}_{bkg}{suffix}"
        hist = file.Get(f"background_{bkg}/{full_hist_name}")
        if hist:
            bkg_obs_data[obs] = sample_from_hist(hist, n_samples)
        else:
            print(f"⚠️ Missing: {full_hist_name} in {bkg}")
            bkg_obs_data[obs] = np.full(n_samples, np.nan)
    df_bkg = pd.DataFrame(bkg_obs_data)
    df_bkg["label"] = 0
    background_data.append(df_bkg)

background_df = pd.concat(background_data, ignore_index=True)
full_df = pd.concat([signal_df, background_df], ignore_index=True)

full_df.dropna(inplace=True)
full_df.to_csv("ml_input_from_histograms.csv", index=False)
print("✅ Saved: ml_input_from_histograms.csv with", full_df.shape[0], "events")
