#!/usr/bin/env python3
import os, sys, gzip, math, re
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------
# DEFAULT INPUT FILES (edit here)
# -------------------------------
DEFAULT_FILES = [
    "/home/hamzeh-khanpour/SuperChic_1/BUILD/evrecs/evrecww_el_13TeV_elWplus.dat",      # SM edff
    "/home/hamzeh-khanpour/SuperChic_1/BUILD/evrecs/evrecww_sd_13TeV_elWplus.dat",     # EFT FM2 = 1 TeV^-4 edff
]

# Histogram settings (edit if needed)
MMIN   = 160.0
MMAX   = 5000.0
NBINS  = 60
LOG_Y  = True

# -------------------------------

def open_any(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'rt')

def parse_sigma_nevents(lhe_path):
    sigma_pb = None
    nevents = None
    with open_any(lhe_path) as f:
        for line in f:
            if '<MGGenerationInfo>' in line:
                break
        for line in f:
            if '</MGGenerationInfo>' in line:
                break
            if 'Integrated weight (pb)' in line:
                m = re.search(r'Integrated weight \(pb\)\s*:\s*([0-9Ee+\-\.]+)', line)
                if m: sigma_pb = float(m.group(1))
            if 'Number of Events' in line:
                m = re.search(r'Number of Events\s*:\s*([0-9]+)', line)
                if m: nevents = int(m.group(1))
    if sigma_pb is None or nevents is None:
        # Fallback scan (rarely needed)
        with open_any(lhe_path) as f:
            for line in f:
                if sigma_pb is None and 'Integrated weight (pb)' in line:
                    m = re.search(r'([0-9Ee+\-\.]+)', line)
                    if m: sigma_pb = float(m.group(1))
                if nevents is None and 'Number of Events' in line:
                    m = re.search(r'([0-9]+)', line)
                    if m: nevents = int(m.group(1))
                if sigma_pb is not None and nevents is not None:
                    break
    if sigma_pb is None or nevents is None:
        raise RuntimeError(f"Could not read sigma/nevents from {lhe_path}")
    return sigma_pb, nevents

def parse_MWW(lhe_path):
    MWW = []
    with open_any(lhe_path) as f:
        in_evt = False
        wplus = None
        wminus = None
        for line in f:
            s = line.strip()
            if s == '<event>':
                in_evt = True; wplus = None; wminus = None
                continue
            if s == '</event>':
                if wplus is not None and wminus is not None:
                    E  = wplus[0] + wminus[0]
                    px = wplus[1] + wminus[1]
                    py = wplus[2] + wminus[2]
                    pz = wplus[3] + wminus[3]
                    m2 = max(E*E - (px*px + py*py + pz*pz), 0.0)
                    MWW.append(math.sqrt(m2))
                in_evt = False
                continue
            if not in_evt or not s or (s[0] not in '-0123456789'):
                continue
            cols = s.split()
            if len(cols) < 10:
                continue
            pdg  = int(cols[0]); stat = int(cols[1])
            px   = float(cols[6]); py  = float(cols[7])
            pz   = float(cols[8]); E   = float(cols[9])
            if stat == 1 and pdg == 24:
                wplus  = [E, px, py, pz]
            elif stat == 1 and pdg == -24:
                wminus = [E, px, py, pz]
    return np.array(MWW, dtype=float)

def dsigma_hist(MWW, sigma_pb, nevents, mmin=MMIN, mmax=MMAX, nbins=NBINS):
    counts, edges = np.histogram(MWW, bins=nbins, range=(mmin, mmax))
    widths = np.diff(edges)
    factor = (sigma_pb / float(nevents))  # pb per event
    dsig = counts * factor / widths       # pb/GeV
    centers = 0.5*(edges[:-1] + edges[1:])
    return centers, dsig, edges

def auto_label(path):
    base = os.path.basename(path)
    for tag in ['iww','edff','chff','luxqed','lhapdf']:
        if tag in base.lower():
            return f"{tag.upper()} ({base})"
    return base

def main():
    # Use CLI files if provided, else defaults
    files = sys.argv[1:] if len(sys.argv) > 1 else DEFAULT_FILES

    # Basic existence check
    missing = [f for f in files if not os.path.isfile(f)]
    if missing:
        print("[ERROR] Missing files:")
        for m in missing:
            print("   ", m)
        print("\nEdit DEFAULT_FILES at the top of this script or pass paths on the command line.")
        sys.exit(1)

    plt.figure(figsize=(8,6))
    for fpath in files:
        sigma_pb, nevents = parse_sigma_nevents(fpath)
        MWW = parse_MWW(fpath)
        if len(MWW) == 0:
            print(f"[WARN] No W+W- found in {fpath}")
            continue
        x, y, edges = dsigma_hist(MWW, sigma_pb, nevents)
        label = auto_label(fpath) + f" (σ={sigma_pb:.4g} pb, N={nevents})"
        plt.step(x, y, where='mid', label=label)
        # Write CSV next to input
        csv = np.column_stack([x, y])
        out_csv = os.path.splitext(fpath)[0] + "_dSigmadM.csv"
        np.savetxt(out_csv, csv, delimiter=",",
                   header="M_WW[GeV], dSigma/dM[pb/GeV]", comments='')
        print(f"[OK] {fpath}: {len(MWW)} events with W+W-, wrote {out_csv}")



    ax = plt.gca()
    if LOG_Y:
        ax.set_yscale('log')
        ax.set_ylim(1e-8, 1e-3)


    plt.xlabel(r"$W_{\gamma\gamma}=M_{WW}\ \mathrm{[GeV]}$")
    plt.ylabel(r"$\mathrm{d}\sigma/\mathrm{d}M_{WW}\ \mathrm{[pb/GeV]}$")
    if LOG_Y: plt.yscale('log')
    plt.legend()
    plt.tight_layout()
    plt.savefig("gamma-UPC-aaww-SM-EFT-edff.png", dpi=200)
    plt.savefig("gamma-UPC-aaww-SM-EFT-edff.pdf")
    print("Wrote gamma-UPC-aaww-SM-EFT-edff.(png|pdf)")

if __name__ == "__main__":
    main()
