#!/usr/bin/env python3
import argparse
import gzip
import io
import math
import os
import re
import sys
from typing import Dict, List, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt
import csv

# -----------------------------
# LHE parsing utilities
# -----------------------------

RE_HEADER = re.compile(r"<header>(.*?)</header>", re.DOTALL | re.MULTILINE)
RE_MGRUNCARD = re.compile(r"<MGRunCard>\s*<!\[CDATA\[(.*?)\]\]>\s*</MGRunCard>", re.DOTALL | re.MULTILINE)
RE_MGPROCCARD = re.compile(r"<MG5ProcCard>\s*<!\[CDATA\[(.*?)\]\]>\s*</MG5ProcCard>", re.DOTALL | re.MULTILINE)
RE_SLHA = re.compile(r"<slha>\s*(.*?)\s*</slha>", re.DOTALL | re.MULTILINE)
RE_GENINFO = re.compile(r"<MGGenerationInfo>(.*?)</MGGenerationInfo>", re.DOTALL | re.MULTILINE)

def _open_maybe_gz(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "rt", encoding="utf-8", errors="ignore")

def read_file_text(path: str) -> str:
    with _open_maybe_gz(path) as f:
        return f.read()

def parse_header_blocks(text: str) -> Dict[str, str]:
    m_header = RE_HEADER.search(text)
    header = m_header.group(1) if m_header else text
    out = {}
    for name, regex in [("MGRunCard", RE_MGRUNCARD), ("MG5ProcCard", RE_MGPROCCARD), ("SLHA", RE_SLHA), ("MGGenerationInfo", RE_GENINFO)]:
        m = regex.search(header)
        out[name] = m.group(1) if m else ""
    out["full_header"] = header
    return out

def _find_float(pattern: str, blob: str) -> Optional[float]:
    m = re.search(pattern, blob, re.MULTILINE)
    if not m:
        return None
    try:
        return float(m.group(1))
    except Exception:
        return None

def _find_int(pattern: str, blob: str) -> Optional[int]:
    m = re.search(pattern, blob, re.MULTILINE)
    if not m:
        return None
    try:
        return int(m.group(1))
    except Exception:
        return None

def _find_str(pattern: str, blob: str) -> Optional[str]:
    m = re.search(pattern, blob, re.MULTILINE)
    if not m:
        return None
    return m.group(1).strip()

def parse_header_metadata(text: str) -> Dict[str, object]:
    blocks = parse_header_blocks(text)
    run = blocks.get("MGRunCard", "")
    gen = blocks.get("MGGenerationInfo", "")
    proc = blocks.get("MG5ProcCard", "")
    slha = blocks.get("SLHA", "")
    header_all = blocks.get("full_header", "")

    sigma_pb = _find_float(r"Integrated weight \(pb\)\s*:\s*([0-9Ee\+\-\.]+)", gen) or \
               _find_float(r"Integrated weight \(pb\)\s*:\s*([0-9Ee\+\-\.]+)", header_all)
    nevt = _find_int(r"Number of Events\s*:\s*([0-9]+)", gen) or \
           _find_int(r"Number of Events\s*:\s*([0-9]+)", header_all)

    lpp1 = _find_int(r"^\s*([\-0-9]+)\s*=\s*lpp1", run)
    lpp2 = _find_int(r"^\s*([\-0-9]+)\s*=\s*lpp2", run)
    ebeam1 = _find_float(r"^\s*([0-9Ee\+\-\.]+)\s*=\s*ebeam1", run)
    ebeam2 = _find_float(r"^\s*([0-9Ee\+\-\.]+)\s*=\s*ebeam2", run)
    pdlabel = _find_str(r"^\s*([A-Za-z0-9\-_]+)\s*=\s*pdlabel", run)

    ptl = _find_float(r"^\s*([0-9Ee\+\-\.]+)\s*=\s*ptl", run)
    etal = _find_float(r"^\s*([0-9Ee\+\-\.]+)\s*=\s*etal", run)
    drll = _find_float(r"^\s*([0-9Ee\+\-\.]+)\s*=\s*drll", run)
    ptl1min = _find_float(r"^\s*([0-9Ee\+\-\.]+)\s*=\s*ptl1min", run)
    ptl2min = _find_float(r"^\s*([0-9Ee\+\-\.]+)\s*=\s*ptl2min", run)

    np_from_proc = _find_int(r"NP\s*=\s*([0-9]+)", proc)
    model_tag = None
    if re.search(r"\bSM_LM0123_UFO\b", proc) or re.search(r"\bSM_LM0123_UFO\b", header_all):
        model_tag = "SM_LM0123_UFO"

    eft_nonzero = []
    m_ano = re.search(r"Block\s+anoinputs(.*?)(?:\n\n|###################################|Block|\Z)", slha, re.DOTALL | re.IGNORECASE)
    if m_ano:
        for line in m_ano.group(1).splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("#")[0].split()
            if len(parts) >= 2:
                try:
                    val = float(parts[1])
                    if abs(val) > 0.0:
                        name = None
                        mname = re.search(r"#\s*([A-Za-z0-9_]+)", line)
                        if mname:
                            name = mname.group(1)
                        eft_nonzero.append((name or "c", val))
                except Exception:
                    pass

    return {
        "sigma_pb": sigma_pb, "nevents": nevt,
        "lpp1": lpp1, "lpp2": lpp2,
        "ebeam1_GeV": ebeam1, "ebeam2_GeV": ebeam2,
        "pdlabel": pdlabel,
        "ptl_GeV": ptl, "etal_max": etal, "drll_min": drll,
        "ptl1min_GeV": ptl1min, "ptl2min_GeV": ptl2min,
        "NP": np_from_proc, "model": model_tag,
        "eft_couplings": eft_nonzero,
    }

# -----------------------------
# Event-level parsing
# -----------------------------

PDG_E = {11, -11}
PDG_MU = {13, -13}
PDG_NU = {12, 14, 16, -12, -14, -16}
PDG_JET_PARTONS = {1,2,3,4, -1,-2,-3,-4, 21}
# NEW: leptons + (e,mu) neutrinos for M_WW reconstruction
PDG_LEPNU_EMU = {11, -11, 13, -13, 12, -12, 14, -14}

def parse_lhe_events(path: str):
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
                    nup = int(header[0]) if header and header[0].strip("-").isdigit() else None
                    rows = buf[1:1+(nup or 0)]
                    particles = []
                    for row in rows:
                        cols = row.strip().split()
                        if len(cols) < 13:  # needs 13 cols
                            continue
                        try:
                            pid = int(cols[0]); status = int(cols[1])
                            px, py, pz, E, M = map(float, cols[6:11])
                            particles.append({"id": pid, "status": status, "px": px, "py": py, "pz": pz, "E": E, "M": M})
                        except Exception:
                            continue
                    yield particles
                    in_event = False
                else:
                    buf.append(line)

def pT(px: float, py: float) -> float:
    return math.hypot(px, py)

def eta(px: float, py: float, pz: float) -> float:
    p = math.sqrt(px*px + py*py + pz*pz)
    if p == abs(pz):
        return float("inf") if pz >= 0 else -float("inf")
    return 0.5 * math.log((p + pz) / (p - pz))

def deltaR(a, b) -> float:
    phi1 = math.atan2(a["py"], a["px"])
    phi2 = math.atan2(b["py"], b["px"])
    dphi = (phi1 - phi2 + math.pi) % (2*math.pi) - math.pi
    eta1 = eta(a["px"], a["py"], a["pz"])
    eta2 = eta(b["px"], b["py"], b["pz"])
    return math.hypot(dphi, eta1 - eta2)

def pick_e_mu(event_particles):
    es = [p for p in event_particles if p["status"] == 1 and p["id"] in PDG_E]
    mus = [p for p in event_particles if p["status"] == 1 and p["id"] in PDG_MU]
    e = max(es, key=lambda x: pT(x["px"], x["py"])) if es else None
    mu = max(mus, key=lambda x: pT(x["px"], x["py"])) if mus else None
    return e, mu

def build_pmiss(event_particles):
    px = sum(p["px"] for p in event_particles if p["status"] == 1 and p["id"] in PDG_NU)
    py = sum(p["py"] for p in event_particles if p["status"] == 1 and p["id"] in PDG_NU)
    return math.hypot(px, py)

def has_cms_jet(event_particles, pt_thresh=30.0, eta_max=4.7) -> bool:
    for p in event_particles:
        if p["status"] != 1:
            continue
        if p["id"] in PDG_JET_PARTONS:
            if pT(p["px"], p["py"]) > pt_thresh and abs(eta(p["px"], p["py"], p["pz"])) < eta_max:
                return True
    return False

def cms_fiducial_pass(e, mu, pmiss, event_particles,
                      lead_pt=25.0, sublead_pt=15.0,
                      eta_max=2.4, dR_min=0.4,
                      jet_pt_veto=30.0, jet_eta_veto=4.7) -> bool:
    if e is None or mu is None:
        return False
    pts = sorted([pT(e["px"], e["py"]), pT(mu["px"], mu["py"])], reverse=True)
    if pts[0] <= lead_pt or pts[1] <= sublead_pt:
        return False
    if abs(eta(e["px"], e["py"], e["pz"])) >= eta_max:
        return False
    if abs(eta(mu["px"], mu["py"], mu["pz"])) >= eta_max:
        return False
    if deltaR(e, mu) <= dR_min:
        return False
    if pmiss <= 30.0:
        return False
    if has_cms_jet(event_particles, pt_thresh=jet_pt_veto, eta_max=jet_eta_veto):
        return False
    return True

# -----------------------------
# Histogram helpers
# -----------------------------

def make_hist(data: List[float], weights: List[float], bins: int, rng: Tuple[float, float]):
    hist, edges = np.histogram(np.array(data), bins=bins, range=rng, weights=np.array(weights))
    widths = np.diff(edges)
    dsig = np.divide(hist, widths, out=np.zeros_like(hist, dtype=float), where=widths>0)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, dsig, edges

# UPDATED: add linestyle control (SM dashed, EFT solid by default)
def save_plot_and_csv(x, y_sm, y_eft, edges, out_prefix: str, xlabel: str,
                      logy: bool, sm_label: str, eft_label: str,
                      sm_ls: str="--", eft_ls: str="-"):
    plt.figure()
    plt.step(edges[:-1], y_sm, where="post", label=sm_label, linewidth=1.8, linestyle=sm_ls)
    plt.step(edges[:-1], y_eft, where="post", label=eft_label, linewidth=1.8, linestyle=eft_ls)
    plt.xlabel(xlabel)
    plt.ylabel(r"d$\sigma$/dX  [pb / GeV]")
    if logy:
        plt.yscale("log")
    plt.legend()
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(f"{out_prefix}.{ext}")
    plt.close()

    with open(f"{out_prefix}.csv", "w", newline="") as cf:
        w = csv.writer(cf)
        w.writerow(["bin_center", "SM_dsig", "EFT_dsig"])
        for i in range(len(x)):
            w.writerow([f"{x[i]:.6g}", f"{y_sm[i]:.6g}", f"{y_eft[i]:.6g}"])


# -----------------------------
# One-sample analysis
# -----------------------------

def analyze_one_sample(path: str,
                       bins_e_mu=50, rng_e_mu=(0.0, 500.0),
                       bins_emu=50, rng_emu=(0.0, 600.0),
                       # NEW: M_WW binning controls
                       bins_mww=50, rng_mww=(160.0, 4000.0),
                       cms_fiducial=False,
                       fiducial_params=None):
    text = read_file_text(path)
    meta = parse_header_metadata(text)

    sigma_pb = meta.get("sigma_pb") or 1.0
    nevt = meta.get("nevents") or 1
    w = sigma_pb / nevt

    pt_e_vals, pt_mu_vals, pt_emu_vals = [], [], []
    wts_e, wts_mu, wts_emu = [], [], []
    # NEW: storage for M_WW
    mww_vals, wts_mww = [], []

    n_total = 0
    n_pass = 0

    for ev in parse_lhe_events(path):
        n_total += 1
        e, mu = pick_e_mu(ev)
        if e is None or mu is None:
            continue
        pmiss = build_pmiss(ev)
        passes = True
        if cms_fiducial:
            passes = cms_fiducial_pass(e, mu, pmiss, ev, **(fiducial_params or {}))
        if not passes:
            continue
        n_pass += 1

        # pT observables
        pe = pT(e["px"], e["py"])
        pm = pT(mu["px"], mu["py"])
        pemu = math.hypot(e["px"] + mu["px"], e["py"] + mu["py"])

        pt_e_vals.append(pe);     wts_e.append(w)
        pt_mu_vals.append(pm);    wts_mu.append(w)
        pt_emu_vals.append(pemu); wts_emu.append(w)

        # NEW: M_WW from status==1 leptons and (e,mu) neutrinos
        Ex = Ey = Ez = EE = 0.0
        for p in ev:
            if p["status"] == 1 and p["id"] in PDG_LEPNU_EMU:
                Ex += p["px"]; Ey += p["py"]; Ez += p["pz"]; EE += p["E"]
        m2 = EE*EE - (Ex*Ex + Ey*Ey + Ez*Ez)
        if m2 < 0.0:  # numerical safety
            m2 = 0.0
        mww = math.sqrt(m2)
        if mww > 0.0:
            mww_vals.append(mww)
            wts_mww.append(w)

    c_e,   dsig_e,   edges_e   = make_hist(pt_e_vals,  wts_e,   bins=bins_e_mu, rng=rng_e_mu)
    c_mu,  dsig_mu,  edges_mu  = make_hist(pt_mu_vals, wts_mu,  bins=bins_e_mu, rng=rng_e_mu)
    c_emu, dsig_emu, edges_emu = make_hist(pt_emu_vals,wts_emu, bins=bins_emu,  rng=rng_emu)
    # NEW: M_WW histogram
    c_mww, dsig_mww, edges_mww = make_hist(mww_vals,   wts_mww, bins=bins_mww, rng=rng_mww)

    return {
        "meta": meta,
        "n_total": n_total,
        "n_pass": n_pass,
        "sigma_pb": sigma_pb,
        "sigma_fid_pb": sigma_pb * (n_pass / nevt if nevt else 0.0),
        "hists": {
            "e":    (c_e,    dsig_e,    edges_e),
            "mu":   (c_mu,   dsig_mu,   edges_mu),
            "emu":  (c_emu,  dsig_emu,  edges_emu),
            "mww":  (c_mww,  dsig_mww,  edges_mww),  # NEW
        }
    }



# -----------------------------
# CLI / main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(description="Compare SM vs EFT pT spectra + M_WW from LHE files (aa->WW->eμ).")
    ap.add_argument("sm_lhe", help="SM LHE path (.lhe or .lhe.gz)")
    ap.add_argument("eft_lhe", help="EFT LHE path (.lhe or .lhe.gz)")
    ap.add_argument("--bins-e", type=int, default=50, help="bins for pT(e±) and pT(μ±)")
    ap.add_argument("--range-e", type=float, nargs=2, default=(0.0, 500.0), help="min max for pT(e±)/pT(μ±) GeV")
    ap.add_argument("--bins-emu", type=int, default=50, help="bins for pT(eμ)")
    ap.add_argument("--range-emu", type=float, nargs=2, default=(0.0, 600.0), help="min max for pT(eμ) GeV")
    # NEW: M_WW binning
    ap.add_argument("--bins-mww", type=int, default=50, help="bins for M_WW")
    ap.add_argument("--range-mww", type=float, nargs=2, default=(160.0, 4000.0), help="min max for M_WW [GeV]")
    ap.add_argument("--logy", action=argparse.BooleanOptionalAction, default=True, help="log-scale y")
    ap.add_argument("--cms-fiducial", action=argparse.BooleanOptionalAction, default=False, help="apply CMS-like selection")
    ap.add_argument("--lead-pt", type=float, default=25.0, help="leading lepton pT min (GeV)")
    ap.add_argument("--sublead-pt", type=float, default=15.0, help="subleading lepton pT min (GeV)")
    ap.add_argument("--eta-max", type=float, default=2.4, help="lepton |eta| max")
    ap.add_argument("--dR-min", type=float, default=0.4, help="ΔR(e,μ) min")
    ap.add_argument("--jet-pt-veto", type=float, default=30.0, help="jet veto pT (GeV)")
    ap.add_argument("--jet-eta-veto", type=float, default=4.7, help="jet veto |eta| max")
    args = ap.parse_args()

    fid = {
        "lead_pt": args.lead_pt,
        "sublead_pt": args.sublead_pt,
        "eta_max": args.eta_max,
        "dR_min": args.dR_min,
        "jet_pt_veto": args.jet_pt_veto,
        "jet_eta_veto": args.jet_eta_veto,
    }

    res_sm = analyze_one_sample(
        args.sm_lhe,
        bins_e_mu=args.bins_e, rng_e_mu=tuple(args.range_e),
        bins_emu=args.bins_emu, rng_emu=tuple(args.range_emu),
        bins_mww=args.bins_mww, rng_mww=tuple(args.range_mww),
        cms_fiducial=args.cms_fiducial, fiducial_params=fid,
    )
    res_eft = analyze_one_sample(
        args.eft_lhe,
        bins_e_mu=args.bins_e, rng_e_mu=tuple(args.range_e),
        bins_emu=args.bins_emu, rng_emu=tuple(args.range_emu),
        bins_mww=args.bins_mww, rng_mww=tuple(args.range_mww),
        cms_fiducial=args.cms_fiducial, fiducial_params=fid,
    )

    def label(meta, tag):
        s = f"{tag}"
        if meta.get('sigma_pb') is not None and meta.get('nevents') is not None:
            s += f" (σ={meta['sigma_pb']:.6g} pb, N={meta['nevents']})"
        return s

    # Print compact summaries
    for tag, R in (("SM", res_sm), ("EFT", res_eft)):
        m = R["meta"]
        print(f"== {tag} ==")
        print("sigma_pb:", m.get("sigma_pb"), "| nevents:", m.get("nevents"),
              "| beams:", (m.get("lpp1"), m.get("lpp2")), "Ebeam(GeV):", (m.get("ebeam1_GeV"), m.get("ebeam2_GeV")),
              "| pdlabel:", m.get("pdlabel"),
              "| cuts:", f"ptl≥{m.get('ptl_GeV')} GeV, |eta|≤{m.get('etal_max')}, drll≥{m.get('drll_min')} ;"
                         f" ptl1min={m.get('ptl1min_GeV')}, ptl2min={m.get('ptl2min_GeV')}",
              "| model:", m.get("model"), "NP:", m.get("NP"))
        if m.get("eft_couplings"):
            print("  EFT nonzero:", ", ".join([f"{n}={v:g}" for n,v in m["eft_couplings"]]))
        print(f"  Events: total={R['n_total']}, pass={R['n_pass']}, "
              f"sigma_tot={R['sigma_pb']} pb, sigma_fid={R['sigma_fid_pb']} pb")

    # Pull histos
    (x_e,  y_e_sm,  eedges)  = res_sm["hists"]["e"]
    (_,    y_e_eft, _)       = res_eft["hists"]["e"]
    (x_mu, y_mu_sm, muedges) = res_sm["hists"]["mu"]
    (_,    y_mu_eft, _)      = res_eft["hists"]["mu"]
    (x_em, y_em_sm, emedges) = res_sm["hists"]["emu"]
    (_,    y_em_eft, _)      = res_eft["hists"]["emu"]
    # NEW: M_WW
    (x_mww, y_mww_sm, mww_edges) = res_sm["hists"]["mww"]
    (_,     y_mww_eft, _ )       = res_eft["hists"]["mww"]

    # Save plots + CSVs (SM dashed, EFT solid)
    save_plot_and_csv(x_e,   y_e_sm,   y_e_eft,   eedges,   "pt_e_compare_SM_vs_EFT",
                      r"$p_T(e^\pm)$ [GeV]",  args.logy, label(res_sm["meta"], "SM"),  label(res_eft["meta"], "EFT"),
                      sm_ls="--", eft_ls="-")
    save_plot_and_csv(x_mu,  y_mu_sm,  y_mu_eft,  muedges,  "pt_mu_compare_SM_vs_EFT",
                      r"$p_T(\mu^\pm)$ [GeV]", args.logy, label(res_sm["meta"], "SM"),  label(res_eft["meta"], "EFT"),
                      sm_ls="--", eft_ls="-")
    save_plot_and_csv(x_em,  y_em_sm,  y_em_eft,  emedges,  "pt_emu_compare_SM_vs_EFT",
                      r"$p_T(e\mu)$ [GeV]",    args.logy, label(res_sm["meta"], "SM"),  label(res_eft["meta"], "EFT"),
                      sm_ls="--", eft_ls="-")
    # NEW: M_WW spectrum
    save_plot_and_csv(x_mww, y_mww_sm, y_mww_eft, mww_edges, "mww_compare_SM_vs_EFT",
                      r"$M_{WW}$ [GeV]",       args.logy, label(res_sm["meta"], "SM"),  label(res_eft["meta"], "EFT"),
                      sm_ls="--", eft_ls="-")

    # Quick textual ratio
    sm = res_sm["sigma_pb"] or 0.0
    eft = res_eft["sigma_pb"] or 0.0
    if sm > 0:
        rel = 100.0 * (eft-sm)/sm
        print(f"[Summary] EFT total σ is {rel:.2f}% {'higher' if rel>=0 else 'lower'} than SM ({eft:.6g} pb vs {sm:.6g} pb).")

if __name__ == "__main__":
    main()
