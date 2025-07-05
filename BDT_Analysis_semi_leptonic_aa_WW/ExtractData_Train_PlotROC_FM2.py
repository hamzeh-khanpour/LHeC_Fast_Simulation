
## ================================================================================
##        Hamzeh Khanpour  --- June 2025
## ================================================================================


import ROOT
import numpy as np
import pandas as pd

#-------------------------------
# SETTINGS
#-------------------------------
n_samples = 100000
luminosity_fb = 100.0  # Target luminosity

# Cross sections (pb)
#signal_cross_section_fb = 0.015395800000000001 * 1000.0  #FM0
#signal_cross_section_fb = 0.015587299999999998 * 1000.0  #FM1
signal_cross_section_fb = 0.0173508 * 1000.0    # 0.022850099999999998 * 1000.0  #FM2
#signal_cross_section_fb = 0.017148100000000000 * 1000.0  #FM3

background_cross_sections_fb = {
    "aa_ww": 0.0149219 * 1000.0,
    "aa_ttbar": 4.824774e-05 * 1000.0,
    "aa_tautau": 1.806765e-01 * 1000.0,
    "aa_tautau_inel": 1.156165e-01 * 1000.0,
    "inclusive_ttbar": 0.00817326 * 1000.0,
    "single_top": 1.36211000   * 1000.0,
    "w_production": 1.965201542 * 1000.0,
    "z_production": 0.159347434 * 1000.0,
    "wwj": 0.02031491612402401   * 1000.0,
    "zzj": 8.106588466651764e-05 * 1000.0,
    "wzj": 0.0028587542003382592 * 1000.0
}



observable_map = {
    "lepton_pt": "hist_lepton_pt",
    "leading_jet_pt": "hist_leading_jet_pt",
    "lepton_eta": "hist_lepton_eta",
    "delta_r": "hist_delta_r",
    "missing_et": "hist_missing_et",
    "subleading_jet_eta": "hist_subleading_jet_eta",
    "leading_jet_eta": "hist_leading_jet_eta",
    "jet_centrality": "hist_jet_centrality",
    "delta_eta_jj": "hist_delta_eta_jj",
    "m_w_hadronic": "hist_m_w_hadronic",
    "m_w_leptonic": "hist_m_w_leptonic",

    # ✅ Newly added discriminating observables
    "pt_w_leptonic": "hist_pt_w_leptonic",
    "pt_w_hadronic": "hist_pt_w_hadronic",
    "delta_phi_lep_met": "hist_delta_phi_lep_met",
    "mt_w_leptonic": "hist_mt_w_leptonic",
    "ht_total": "hist_ht_total",
    "delta_phi_jj": "hist_delta_phi_jj",
    "delta_phi_wl_wh": "hist_delta_phi_wl_wh",
    "delta_eta_wl_wh": "hist_delta_eta_wl_wh",
    "m_jj": "hist_m_jj",
    "m_lvjj": "hist_m_lvjj"
}




special_suffix = {"wwj", "zzj", "wzj"}
background_keys = list(background_cross_sections_fb.keys())


#-------------------------------
# FUNCTIONS
#-------------------------------
def sample_from_hist(hist, n):
    values = []
    for i in range(hist.GetNbinsX()):
        count = int(hist.GetBinContent(i + 1))
        xlow = hist.GetBinLowEdge(i + 1)
        xhigh = hist.GetBinLowEdge(i + 2)
        if count > 0:
            values.extend(np.random.uniform(xlow, xhigh, size=count))
    if len(values) == 0:
        return np.full(n, np.nan)
    return np.random.choice(values, size=n, replace=True)


#-------------------------------
# PROCESS SIGNAL
#-------------------------------
file = ROOT.TFile.Open("output_histograms_FM2.root")
signal_data = {}
for obs, hist_base in observable_map.items():
    full_hist_name = f"{hist_base}_FM2_Lambda4"
    hist = file.Get(f"signal_FM2_Lambda4/{full_hist_name}")
    if hist:
        signal_data[obs] = sample_from_hist(hist, n_samples)
    else:
        print(f"❌ Missing: {full_hist_name} in signal")
        signal_data[obs] = np.full(n_samples, np.nan)

signal_df = pd.DataFrame(signal_data)
signal_df["label"] = 1
signal_df["weight"] = signal_cross_section_fb * luminosity_fb
signal_df["process"] = "FM2_Lambda4"


#-------------------------------
# PROCESS BACKGROUNDS
#-------------------------------
background_dataframes = []
for bkg in background_keys:
    bkg_data = {}
    suffix = "_production" if bkg in special_suffix else ""
    for obs, hist_base in observable_map.items():
        full_hist_name = f"{hist_base}_{bkg}{suffix}"
        hist = file.Get(f"background_{bkg}/{full_hist_name}")
        if hist:
            bkg_data[obs] = sample_from_hist(hist, n_samples)
        else:
            print(f"⚠️ Missing: {full_hist_name} in {bkg}")
            bkg_data[obs] = np.full(n_samples, np.nan)
    df_bkg = pd.DataFrame(bkg_data)
    df_bkg["label"] = 0
    df_bkg["weight"] = background_cross_sections_fb[bkg] * luminosity_fb
    df_bkg["process"] = bkg  # 🆕 Track process name
    background_dataframes.append(df_bkg)


#-------------------------------
# MERGE, CLEAN, SAVE
#-------------------------------
background_df = pd.concat(background_dataframes, ignore_index=True)
full_df = pd.concat([signal_df, background_df], ignore_index=True)
full_df.dropna(inplace=True)


# Save
full_df.to_csv("ml_input_from_histograms_FM2.csv", index=False)
print(f"✅ Saved: ml_input_from_histograms_FM2.csv with {full_df.shape[0]} events and weights for {luminosity_fb:.1f} fb^-1.")

