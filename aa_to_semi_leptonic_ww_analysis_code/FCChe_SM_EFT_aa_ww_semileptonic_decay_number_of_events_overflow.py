#!/usr/bin/env python3
"""
FCChe one-lepton + two-jet analysis with 'overflow -> last bin' folding.

- Inputs: two LHE files (SM, EFT)
- Observables: pT(lepton), pT(j1), pT(j2)
- Outputs: step plots (+ 'overflow' tag when enabled) and CSVs that match the plots
- Default yield mode: counts @ FCChe lumi = 1000 fb^-1
"""

import argparse
import gzip
import math
import re
from typing import Dict, List, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt
import csv

# -----------------------------
# LHE header parsing
# -----------------------------

RE_HEADER = re.compile(r"<header>(.*?)</header>", re.DOTALL | re.MULTILINE)
RE_MGRUNCARD = re.compile(r"<MGRunCard>\s*<!\[CDATA\[(.*?)\]\]>\s*</MGRunCard>", re.DOTALL | re.MULTILINE)
RE_MGPROCCARD = re.compile(r"<MG5ProcCard>\s*<!\[CDATA\[(.*?)\]\]>\s*</MG5ProcCard>", re.DOTALL | re.MULTILINE)
RE_SLHA = re.compile(r"<slha>\s*(.*?)\s*</slha>", re.DOTALL | re.MULTILINE)
RE_GENINFO = re.compile(r"<MGGenerationInfo>(.*?)</MGGenerationInfo>", re.DOTALL | re.MULTILINE)

def _open_maybe_gz(path: str):
    return gzip.open(path, "rt", encoding="utf-8", errors="ignore") if path.endswith(".gz") \
           else open(path, "rt", encoding="utf-8", errors="ignore")

def read_file_text(path: str) -> str:
    with _open_maybe_gz(path) as f:
        return f.read()

def _find_float(pattern: str, blob: str) -> Optional[float]:
    m = re.search(pattern, blob, re.MULTILINE)
    return float(m.group(1)) if m else None

def _find_int(pattern: str, blob: str) -> Optional[int]:
    m = re.search(pattern, blob, re.MULTILINE)
    return int(m.group(1)) if m else None

def _find_str(pattern: str, blob: str) -> Optional[str]:
    m = re.search(pattern, blob, re.MULTILINE)
    return m.group(1).strip() if m else None

def parse_header_blocks(text: str) -> Dict[str, str]:
    m_header = RE_HEADER.search(text)
    header = m_header.group(1) if m_header else text
    out = {}
    for name, regex in [("MGRunCard", RE_MGRUNCARD), ("MG5ProcCard", RE_MGPROCCARD),
                        ("SLHA", RE_SLHA), ("MGGenerationInfo", RE_GENINFO)]:
        m = regex.search(header)
        out[name] = m.group(1) if m else ""
    out["full_header"] = header
    return out

def parse_header_metadata(text: str) -> Dict[str, object]:
    blocks = parse_header_blocks(text)
    run = blocks.get("MGRunCard", "")
    gen = blocks.get("MGGenerationInfo", "")
    header_all = blocks.get("full_header", "")

    sigma_pb = _find_float(r"Integrated weight \(pb\)\s*:\s*([0-9Ee\+\-\.]+)", gen) or \
               _find_float(r"Integrated weight \(pb\)\s*:\s*([0-9Ee\+\-\.]+)", header_all)
    nevt = _find_int(r"Number of Events\s*:\s*([0-9]+)", gen) or \
           _find_int(r"Number of Events\s*:\s*([0-9]+)", header_all)

    lpp1 = _find_int(r"^\s*([\-0-9]+)\s*=\s*lpp1", run)
    lpp2 = _find_int(r"^\s*([\-0-9]+)\s*=\s*lpp2", run)
    ebeam1 = _find_float(r"^\s*([0-9Ee\+\-\.]+)\s*=\s*ebeam1", run)
    ebeam2 = _find_float(r"^\s*([0-9Ee\+\-\.]+)\s*=\s*ebeam2", run)
    pdlabel1 = _find_str(r"^\s*([A-Za-z0-9\-_]+)\s*=\s*pdlabel1", run)
    pdlabel2 = _find_str(r"^\s*([A-Za-z0-9\-_]+)\s*=\s*pdlabel2", run)

    return {
        "sigma_pb": sigma_pb, "nevents": nevt,
        "lpp1": lpp1, "lpp2": lpp2,
        "ebeam1_GeV": ebeam1, "ebeam2_GeV": ebeam2,
        "pdlabel1": pdlabel1, "pdlabel2": pdlabel2,
    }

# -----------------------------
# Event-level parsing
# -----------------------------

PDG_E = {11, -11}
PDG_MU = {13, -13}
PDG_CHARGED_LEPTONS = PDG_E | PDG_MU              # (e, mu) only as requested
PDG_NEUTRINOS = {12, 14, 16, -12, -14, -16}
PDG_JET_PARTONS = {1,2,3,4,5, -1,-2,-3,-4,-5, 21} # include b-quark jets + gluon

def parse_lhe_events_with_weights(path: str):
    """Yield tuples (w_ev, particles) where:
       - w_ev is XWGTUP if present, else None
       - particles is a list of dicts with keys: id, status, px, py, pz, E, M
    """
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="ignore") as f:
        in_event = False
        buf = []
        for line in f:
            if "<event" in line:
                in_event = True
                buf = []
                continue
            if in_event:
                if "</event" in line:
                    if not buf:
                        in_event = False
                        continue
                    header = buf[0].strip().split()
                    nup = int(header[0]) if header and header[0].strip("-").isdigit() else 0
                    w_ev = None
                    if len(header) >= 3:
                        try:
                            w_ev = float(header[2])     # XWGTUP
                        except Exception:
                            w_ev = None
                    rows = buf[1:1+nup]
                    particles = []
                    for row in rows:
                        cols = row.strip().split()
                        if len(cols) < 13:
                            continue
                        try:
                            pid = int(cols[0]); status = int(cols[1])
                            px, py, pz, E, M = map(float, cols[6:11])
                            particles.append({"id": pid, "status": status,
                                              "px": px, "py": py, "pz": pz, "E": E, "M": M})
                        except Exception:
                            continue
                    yield (w_ev, particles)
                    in_event = False
                else:
                    buf.append(line)

# Kinematics helpers
def pT(px: float, py: float) -> float:
    return math.hypot(px, py)

def eta(px: float, py: float, pz: float) -> float:
    p = math.sqrt(px*px + py*py + pz*pz)
    if p == abs(pz):
        return float("inf") if pz >= 0 else -float("inf")
    return 0.5 * math.log((p + pz) / (p - pz))

# Object pickers
def pick_hardest_lepton(event_particles) -> Optional[dict]:
    leptons = [p for p in event_particles if p["status"] == 1 and p["id"] in PDG_CHARGED_LEPTONS]
    return max(leptons, key=lambda x: pT(x["px"], x["py"])) if leptons else None

def get_sorted_jets(event_particles) -> List[dict]:
    jets = [p for p in event_particles if p["status"] == 1 and p["id"] in PDG_JET_PARTONS]
    jets.sort(key=lambda x: pT(x["px"], x["py"]), reverse=True)
    return jets

# -----------------------------
# Histogram helpers
# -----------------------------

def _hist_with_overflow_lastbin(data, weights, bins, rng):
    """Histogram where entries > max(range) are added to the last visible bin."""
    arr = np.asarray(data, dtype=float)
    wts = np.asarray(weights, dtype=float)
    lo, hi = float(rng[0]), float(rng[1])
    mask_in = arr <= hi
    hist, edges = np.histogram(arr[mask_in], bins=bins, range=(lo, hi), weights=wts[mask_in])
    if hist.size:
        hist[-1] += wts[~mask_in].sum()
    return hist, edges

def _step_with_last_plateau(edges, y):
    """Return (x,y) so the final bin is drawn as a plateau (not a spike)."""
    x = np.r_[edges[:-1], edges[-1]]
    yy = np.r_[y, y[-1] if len(y) else 0.0]
    return x, yy

def make_hist(data: List[float], weights: List[float], bins: int, rng: Tuple[float, float],
              overflow_lastbin: bool = False):
    """Return (centers, dsig, edges).  'dsig' is pb / unit (bin width).
       If overflow_lastbin=True, entries > max(range) are merged into the last bin *before* width division."""
    if overflow_lastbin:
        hist, edges = _hist_with_overflow_lastbin(data, weights, bins, rng)
    else:
        hist, edges = np.histogram(np.array(data), bins=bins, range=rng, weights=np.array(weights))
    widths = np.diff(edges)
    dsig = np.divide(hist, widths, out=np.zeros_like(hist, dtype=float), where=widths>0)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, dsig, edges

def _convert_for_mode(y_dsig, edges, mode: str, lumi_fb: float):
    """Convert dsig/dX to plotting array based on mode.
       - 'dsig': pb / unit
       - 'counts': expected events = (dsig * width [pb]) * (lumi_fb * 1e3 [pb^-1])
    """
    if mode == "counts":
        widths = np.diff(edges)
        sigma_bin_pb = y_dsig * widths
        return sigma_bin_pb * (lumi_fb * 1e3)  # 1 fb^-1 = 1e3 pb^-1
    return y_dsig

def save_plot_and_csv(x, y_sm, y_eft, edges, out_prefix: str, xlabel: str,
                      logy: bool, sm_label: str, eft_label: str,
                      sm_ls: str="--", eft_ls: str="-",
                      mode: str="counts", lumi_fb: float=1000.0,
                      annotate_overflow: bool=False):
    y_sm_plot  = _convert_for_mode(y_sm,  edges, mode, lumi_fb)
    y_eft_plot = _convert_for_mode(y_eft, edges, mode, lumi_fb)

    plt.figure()
    x_sm,  y_sm2  = _step_with_last_plateau(edges, y_sm_plot)
    x_eft, y_eft2 = _step_with_last_plateau(edges, y_eft_plot)
    plt.step(x_sm,  y_sm2,  where="post", label=sm_label,  linewidth=1.8, linestyle=sm_ls)
    plt.step(x_eft, y_eft2, where="post", label=eft_label, linewidth=1.8, linestyle=eft_ls)
    plt.xlabel(xlabel)
    plt.ylabel("Events" if mode == "counts" else r"d$\sigma$/dX  [pb / unit]")
    if logy:
        plt.yscale("log")
    plt.legend()
    plt.grid(True, which="both", alpha=0.3)
    plt.title(rf"EL–EL: $\gamma\gamma$ @ FCChe @ 1.2 TeV @ $L={lumi_fb:g}\,\mathrm{{fb}}^{{-1}}$")

    if annotate_overflow:
        ax = plt.gca()
        ymax = ax.get_ylim()[1]
        ax.text(edges[-1], 0.95*ymax, "overflow", rotation=90, va="top", ha="right", fontsize=18)

    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(f"{out_prefix}.{ext}")
    plt.close()

    with open(f"{out_prefix}.csv", "w", newline="") as cf:
        w = csv.writer(cf)
        if mode == "counts":
            w.writerow(["bin_lower_edge", "SM_events", "EFT_events"])
            for i in range(len(x)):
                w.writerow([f"{edges[i]:.6g}", f"{y_sm_plot[i]:.6g}", f"{y_eft_plot[i]:.6g}"])
        else:
            w.writerow(["bin_center", "SM_dsig", "EFT_dsig"])
            for i in range(len(x)):
                w.writerow([f"{x[i]:.6g}", f"{y_sm[i]:.6g}", f"{y_eft[i]:.6g}"])

# -------------------------------------------------
# Last-bin (overflow-inclusive) quick report helper
# -------------------------------------------------
def _report_last_bin(name, edges, y_sm, y_eft, mode, lumi_fb):
    """Print SM/EFT values and percent difference for the last bin (overflow folded),
       in the same units as the plot/CSV (counts or dsig)."""
    y_sm_plot  = _convert_for_mode(y_sm,  edges, mode, lumi_fb)
    y_eft_plot = _convert_for_mode(y_eft, edges, mode, lumi_fb)
    last = -1
    lo = edges[-2]
    denom = y_sm_plot[last] if y_sm_plot[last] != 0 else 1e-30
    delta_pct = (y_eft_plot[last] - y_sm_plot[last]) / denom * 100.0
    print(f"[last-bin] {name}  x>{lo:g} :  SM={y_sm_plot[last]:.3g} ,  EFT={y_eft_plot[last]:.3g} ,  Δ={delta_pct:.1f}%")

# -----------------------------
# One-sample analysis
# -----------------------------

def analyze_one_sample(path: str,
                       bins_lep=60, rng_lep=(0.0, 200.0),
                       bins_j =60, rng_j  =(0.0, 200.0),
                       use_per_event_weights=True):
    # header/meta
    text = read_file_text(path)
    meta = parse_header_metadata(text)
    sigma_pb = meta.get("sigma_pb") or 1.0
    nevt     = meta.get("nevents") or 1
    w_global = sigma_pb / nevt

    # observables
    pt_lep_vals, wts_lep = [], []
    pt_j1_vals,  wts_j1  = [], []
    pt_j2_vals,  wts_j2  = [], []

    n_total = 0
    for w_ev, ev in parse_lhe_events_with_weights(path):
        n_total += 1
        # weight choice
        w = (w_ev if (use_per_event_weights and (w_ev is not None)) else w_global)

        lep = pick_hardest_lepton(ev)
        if lep is not None:
            pt_lep_vals.append(pT(lep["px"], lep["py"]))
            wts_lep.append(w)

        jets = get_sorted_jets(ev)
        if len(jets) >= 1:
            pt_j1_vals.append(pT(jets[0]["px"], jets[0]["py"]))
            wts_j1.append(w)
        if len(jets) >= 2:
            pt_j2_vals.append(pT(jets[1]["px"], jets[1]["py"]))
            wts_j2.append(w)

    # histograms (overflow folding is applied by caller via flags)
    return {
        "meta": {**meta, "n_total": n_total, "sigma_pb": sigma_pb},
        "data": {
            "lep": (pt_lep_vals, wts_lep, bins_lep, rng_lep),
            "j1":  (pt_j1_vals,  wts_j1,  bins_j,   rng_j),
            "j2":  (pt_j2_vals,  wts_j2,  bins_j,   rng_j),
        }
    }

# -----------------------------
# CLI / main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(description="FCChe (γγ→WW→ℓ jj) spectra with overflow folding.")
    ap.add_argument("sm_lhe",  help="SM LHE path (.lhe or .lhe.gz)")
    ap.add_argument("eft_lhe", help="EFT LHE path (.lhe or .lhe.gz)")

    # binning
    ap.add_argument("--bins-lep", type=int, default=60, help="bins for pT(lepton)")
    ap.add_argument("--range-lep", type=float, nargs=2, default=(0.0, 200.0), help="range for pT(lepton) [GeV]")
    ap.add_argument("--bins-j",   type=int, default=60, help="bins for pT(j1) and pT(j2)")
    ap.add_argument("--range-j",  type=float, nargs=2, default=(0.0, 200.0), help="range for pT(j1/j2) [GeV]")

    # yields/units
    ap.add_argument("--yield-mode", choices=["dsig", "counts"], default="counts",
                    help="Plot dsig/dX (pb/unit) or expected counts (needs lumi).")
    ap.add_argument("--lumi-fb", type=float, default=1000.0,
                    help="Integrated luminosity in fb^-1 for --yield-mode counts (default: 1000).")
    ap.add_argument("--logy", action=argparse.BooleanOptionalAction, default=False, help="log-scale y")

    # weights
    ap.add_argument("--use-per-event-weights", action=argparse.BooleanOptionalAction, default=True,
                    help="Use XWGTUP per-event weights when available (fallback to σ/nevt).")

    # overflow controls (default ON as requested)
    ap.add_argument("--overflow-lastbin-lep", action=argparse.BooleanOptionalAction, default=True,
                    help="Fold pT(lepton) overflow into last bin and annotate (default: on).")
    ap.add_argument("--overflow-lastbin-j1",  action=argparse.BooleanOptionalAction, default=True,
                    help="Fold pT(j1) overflow into last bin and annotate (default: on).")
    ap.add_argument("--overflow-lastbin-j2",  action=argparse.BooleanOptionalAction, default=True,
                    help="Fold pT(j2) overflow into last bin and annotate (default: on).")
    ap.add_argument("--overflow-lastbin-all", action="store_true",
                    help="Shortcut: turn ON overflow folding for all three (overrides individual toggles).")

    args = ap.parse_args()

    if args.overflow_lastbin_all:
        args.overflow_lastbin_lep = args.overflow_lastbin_j1 = args.overflow_lastbin_j2 = True

    # analyze each sample
    sm = analyze_one_sample(
        args.sm_lhe,
        bins_lep=args.bins_lep, rng_lep=tuple(args.range_lep),
        bins_j=args.bins_j, rng_j=tuple(args.range_j),
        use_per_event_weights=args.use_per_event_weights
    )
    eft = analyze_one_sample(
        args.eft_lhe,
        bins_lep=args.bins_lep, rng_lep=tuple(args.range_lep),
        bins_j=args.bins_j, rng_j=tuple(args.range_j),
        use_per_event_weights=args.use_per_event_weights
    )

    # compact banner
    def label(meta, tag):
        s = f"{tag}"
        if meta.get('sigma_pb') is not None and meta.get('nevents') is not None:
            s += f" (σ={meta['sigma_pb']:.6g} pb, N={meta['nevents']})"
        return s

    for tag, R in (("SM", sm), ("EFT", eft)):
        m = R["meta"]
        print(f"== {tag} ==")
        print(f"sigma_pb: {m.get('sigma_pb')} | nevents: {m.get('nevents')} | "
              f"beams: (lpp1={m.get('lpp1')}, lpp2={m.get('lpp2')}), Ebeam(GeV): "
              f"({m.get('ebeam1_GeV')}, {m.get('ebeam2_GeV')}) | pdlabel: "
              f"({m.get('pdlabel1')}, {m.get('pdlabel2')}) | total events parsed: {m.get('n_total')}")

    # build histos (SM)
    (lep_vals_sm, lep_w_sm, nb_lep, rng_lep) = sm["data"]["lep"]
    (j1_vals_sm,  j1_w_sm,  nb_j,   rng_j)   = sm["data"]["j1"]
    (j2_vals_sm,  j2_w_sm,  _,      _)       = sm["data"]["j2"]
    c_lep, y_lep_sm, e_lep = make_hist(lep_vals_sm, lep_w_sm, nb_lep, rng_lep,
                                       overflow_lastbin=args.overflow_lastbin_lep)
    c_j1,  y_j1_sm,  e_j1  = make_hist(j1_vals_sm,  j1_w_sm,  nb_j,   rng_j,
                                       overflow_lastbin=args.overflow_lastbin_j1)
    c_j2,  y_j2_sm,  e_j2  = make_hist(j2_vals_sm,  j2_w_sm,  nb_j,   rng_j,
                                       overflow_lastbin=args.overflow_lastbin_j2)

    # build histos (EFT) with identical binning (re-histogram EFT into SM edges)
    (lep_vals_eft, lep_w_eft, _, _) = eft["data"]["lep"]
    (j1_vals_eft,  j1_w_eft,  _, _) = eft["data"]["j1"]
    (j2_vals_eft,  j2_w_eft,  _, _) = eft["data"]["j2"]

    # re-hist using same edges to guarantee identical binning
    def hist_with_edges(values, weights, edges, overflow_lastbin):
        lo, hi = edges[0], edges[-1]
        bins = len(edges) - 1
        return make_hist(values, weights, bins, (lo, hi), overflow_lastbin=overflow_lastbin)

    _, y_lep_eft, _ = hist_with_edges(lep_vals_eft, lep_w_eft, e_lep, args.overflow_lastbin_lep)
    _, y_j1_eft,  _ = hist_with_edges(j1_vals_eft,  j1_w_eft,  e_j1,  args.overflow_lastbin_j1)
    _, y_j2_eft,  _ = hist_with_edges(j2_vals_eft,  j2_w_eft,  e_j2,  args.overflow_lastbin_j2)

    # plot + CSV (counts by default @ FCChe lumi)
    mode = args.yield_mode
    Lfb  = args.lumi_fb

    save_plot_and_csv(c_lep, y_lep_sm, y_lep_eft, e_lep,
                      "pt_lepton_FCChe_SM_vs_EFT",
                      r"$p_T(\ell)$ [GeV]",
                      args.logy, label(sm["meta"], "SM"), label(eft["meta"], "EFT"),
                      sm_ls="--", eft_ls="-", mode=mode, lumi_fb=Lfb,
                      annotate_overflow=args.overflow_lastbin_lep)

    save_plot_and_csv(c_j1, y_j1_sm, y_j1_eft, e_j1,
                      "pt_j1_FCChe_SM_vs_EFT",
                      r"$p_T(j_1)$ [GeV]",
                      args.logy, label(sm["meta"], "SM"), label(eft["meta"], "EFT"),
                      sm_ls="--", eft_ls="-", mode=mode, lumi_fb=Lfb,
                      annotate_overflow=args.overflow_lastbin_j1)

    save_plot_and_csv(c_j2, y_j2_sm, y_j2_eft, e_j2,
                      "pt_j2_FCChe_SM_vs_EFT",
                      r"$p_T(j_2)$ [GeV]",
                      args.logy, label(sm["meta"], "SM"), label(eft["meta"], "EFT"),
                      sm_ls="--", eft_ls="-", mode=mode, lumi_fb=Lfb,
                      annotate_overflow=args.overflow_lastbin_j2)

    # after each plot: print last-bin (overflow-inclusive) summary
    _report_last_bin("pT(j1)", e_j1, y_j1_sm, y_j1_eft, mode, Lfb)
    _report_last_bin("pT(ℓ)",  e_lep, y_lep_sm, y_lep_eft, mode, Lfb)
    _report_last_bin("pT(j2)", e_j2, y_j2_sm, y_j2_eft, mode, Lfb)

    # quick σ summary
    sm_sig = sm["meta"]["sigma_pb"] or 0.0
    eft_sig = eft["meta"]["sigma_pb"] or 0.0
    if sm_sig > 0:
        rel = 100.0 * (eft_sig - sm_sig) / sm_sig
        print(f"[Summary] EFT total σ is {rel:.2f}% {'higher' if rel>=0 else 'lower'} than SM "
              f"({eft_sig:.6g} pb vs {sm_sig:.6g} pb).")

if __name__ == "__main__":
    main()
