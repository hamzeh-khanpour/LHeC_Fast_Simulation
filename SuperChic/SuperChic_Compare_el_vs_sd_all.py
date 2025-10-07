#!/usr/bin/env python3
import argparse, re, math, sys, os
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple

Parsed = namedtuple("Parsed", "masses weights sigma_pb nev")

float_re = re.compile(r"[-+]?\d*\.?\d+(?:[Ee][+-]?\d+)?")

def read_init_sigma_pb(lines):
    """Return (sigma_pb, sigma_err_pb) from the <init> block, or (None,None)."""
    in_init = False
    for ln in lines:
        if "<init>" in ln:
            in_init = True
            continue
        if "</init>" in ln:
            break
        if in_init:
            # The single numeric line in <init> has: sigma, err, (often sigma again), 1
            nums = float_re.findall(ln)
            if len(nums) >= 2:
                sig = float(nums[0])
                err = float(nums[1])
                return sig, err
    return None, None

def parse_lhe_file(path):
    """Parse one SuperChic LHE file -> list of M_WW, per-event weights."""
    with open(path, "r") as f:
        lines = f.readlines()

    sigma_pb, _ = read_init_sigma_pb(lines)
    if sigma_pb is None:
        print(f"[WARN] Could not find sigma in <init> for {path}; assuming 0 pb")
        sigma_pb = 0.0

    masses = []
    # Count events robustly
    nev_total = sum(1 for ln in lines if "<event>" in ln)
    # Per-event weight in pb (convert to fb later)
    w_event_pb = sigma_pb / max(nev_total, 1)

    in_event = False
    wplus = None
    wminus = None
    for ln in lines:
        if "<event>" in ln:
            in_event = True
            wplus = None
            wminus = None
            continue
        if "</event>" in ln:
            if wplus is not None and wminus is not None:
                # invariant mass of the W pair
                E = wplus[3] + wminus[3]
                px = wplus[0] + wminus[0]
                py = wplus[1] + wminus[1]
                pz = wplus[2] + wminus[2]
                m2 = E*E - (px*px + py*py + pz*pz)
                m = math.sqrt(max(m2, 0.0))
                masses.append(m)
            in_event = False
            continue
        if not in_event:
            continue

        # Particle line: id, status, moth1, moth2, col1, col2, px,py,pz,E,m, ...
        parts = ln.strip().split()
        if len(parts) < 11:
            continue
        try:
            pid = int(parts[0])
        except ValueError:
            continue
        # Take W bosons regardless of status code (SuperChic uses 2)
        if pid == 24 or pid == -24:
            try:
                px = float(parts[6]); py = float(parts[7]); pz = float(parts[8]); E = float(parts[9])
            except ValueError:
                continue
            four = (px, py, pz, E)
            if pid == 24:
                wplus = four
            else:
                wminus = four

    # weights in fb (so hist is fb/GeV)
    weights_fb = np.full(len(masses), w_event_pb * 1e3)
    return Parsed(np.array(masses), weights_fb, sigma_pb, nev_total)

def auto_label(path):
    base = os.path.basename(path)
    if "_el_" in base or "elastic" in base.lower():
        return "elastic"
    if "_sd_" in base or "single" in base.lower():
        return "single-diss"
    return base

def load_group(paths):
    masses = []
    weights = []
    sigma_pb_sum = 0.0
    nev_sum = 0
    for p in paths:
        pr = parse_lhe_file(p)
        masses.append(pr.masses)
        weights.append(pr.weights)
        sigma_pb_sum += pr.sigma_pb
        nev_sum += pr.nev
    if masses:
        masses = np.concatenate(masses)
        weights = np.concatenate(weights)
    else:
        masses = np.array([])
        weights = np.array([])
    return Parsed(masses, weights, sigma_pb_sum, nev_sum)

def main():
    ap = argparse.ArgumentParser(description="Compare SuperChic γγ→W+W−: elastic vs single-diss")
    ap.add_argument("--el", nargs="+", required=True, help="elastic LHE files (eW+, μW+)")
    ap.add_argument("--sd", nargs="+", required=True, help="single-diss LHE files (eW+, μW+)")
    ap.add_argument("--bins", type=int, default=60)
    ap.add_argument("--mmin", type=float, default=160.0)
    ap.add_argument("--mmax", type=float, default=3000.0)
    ap.add_argument("--out", default="compare_el_vs_sd_MWW")
    args = ap.parse_args()

    el = load_group(args.el)
    sd = load_group(args.sd)

    # Bin edges
    bins = np.linspace(args.mmin, args.mmax, args.bins + 1)

    # Histograms (fb/GeV)
    H_el, _ = np.histogram(el.masses, bins=bins, weights=el.weights)
    H_sd, _ = np.histogram(sd.masses, bins=bins, weights=sd.weights)

    # Plot spectra
    centers = 0.5 * (bins[1:] + bins[:-1])
    width = np.diff(bins)

    plt.figure(figsize=(9,7))
    lbl_el = f"elastic (σ={el.sigma_pb*1e3:.3f} fb, N={el.nev})"
    lbl_sd = f"single-diss (σ={sd.sigma_pb*1e3:.3f} fb, N={sd.nev})"
    plt.step(centers, H_el/width, where="mid", label=lbl_el)
    plt.step(centers, H_sd/width, where="mid", label=lbl_sd)
    plt.xlabel(r"$M_{WW}$ [GeV]")
    plt.ylabel(r"$d\sigma/dM_{WW}$ [fb/GeV]")
    plt.title(r"SuperChic $\gamma\gamma\to W^+W^-$ @ 13 TeV: elastic vs single-diss (combined e$W^+$, $\mu W^+$)")
    if np.any(H_el>0) or np.any(H_sd>0):
        plt.yscale("log")
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.out + ".png", dpi=150)
    plt.savefig(args.out + ".pdf")

    # Ratio SD/EL
    ratio = np.divide(H_sd, H_el, out=np.zeros_like(H_sd), where=(H_el>0))
    plt.figure(figsize=(12,4))
    plt.step(centers, ratio, where="mid")
    plt.xlabel(r"$M_{WW}$ [GeV]")
    plt.ylabel("sd / el")
    plt.title("Ratio (single-diss / elastic)")
    plt.tight_layout()
    plt.savefig(args.out + "_ratio.png", dpi=150)
    plt.savefig(args.out + "_ratio.pdf")

    # CSVs
    with open(args.out + ".csv", "w") as g:
        g.write("bin_low,bin_high,center,dsig_dM_el_fb_per_GeV,dsig_dM_sd_fb_per_GeV,ratio_sd_over_el\n")
        for i in range(len(centers)):
            val_el = H_el[i]/width[i]
            val_sd = H_sd[i]/width[i]
            rat = ratio[i]
            g.write(f"{bins[i]},{bins[i+1]},{centers[i]},{val_el},{val_sd},{rat}\n")

    print("=== Group summary ===")
    print(f"ELASTIC:   files={args.el},  sigma = {el.sigma_pb*1e3:.3f} fb,  N = {el.nev}")
    print(f"SINGLE-D.: files={args.sd},  sigma = {sd.sigma_pb*1e3:.3f} fb,  N = {sd.nev}")
    print(f"Outputs: {args.out}.png/.pdf, {args.out}_ratio.png/.pdf, CSVs with same prefix")

if __name__ == "__main__":
    main()
