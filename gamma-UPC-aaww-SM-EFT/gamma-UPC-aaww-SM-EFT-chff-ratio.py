#!/usr/bin/env python3
import os, sys, gzip, math, re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# -------------------------------
# DEFAULT INPUT FILES (edit here)
# -------------------------------
DEFAULT_FILES = [
    "/home/hamzeh-khanpour/MG5_aMC_v3_6_4/gamma-UPC-aaww-SM-EFT/Events/run_01/gamma-UPC-aaww-SM-chff.lhe",
    "/home/hamzeh-khanpour/MG5_aMC_v3_6_4/gamma-UPC-aaww-SM-EFT/Events/run_02/gamma-UPC-aaww-EFT-chff.lhe",
]

# Histogram settings
MMIN   = 160.0
MMAX   = 5000.0
NBINS  = 60
LOG_Y  = True

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
        with open_any(lhe_path) as f:
            for line in f:
                if sigma_pb is None and 'Integrated weight (pb)' in line:
                    m = re.search(r'([0-9Ee+\-\.]+)', line);  sigma_pb = float(m.group(1)) if m else sigma_pb
                if nevents is None and 'Number of Events' in line:
                    m = re.search(r'([0-9]+)', line);        nevents = int(m.group(1)) if m else nevents
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

def hist_counts(MWW, mmin=MMIN, mmax=MMAX, nbins=NBINS):
    counts, edges = np.histogram(MWW, bins=nbins, range=(mmin, mmax))
    centers = 0.5*(edges[:-1] + edges[1:])
    widths  = np.diff(edges)
    return counts.astype(float), centers, widths, edges

def auto_role(path, fallback=None):
    base = os.path.basename(path).lower()
    if 'eft' in base and 'sm' not in base: return 'EFT'
    if 'sm'  in base and 'eft' not in base: return 'SM'
    return fallback

def dsigma_from_counts(counts, sigma_pb, nevents, widths):
    # pb/GeV per bin
    factor = sigma_pb / float(nevents)
    return (counts * factor) / widths

def main():
    files = sys.argv[1:] if len(sys.argv) > 1 else DEFAULT_FILES
    if len(files) != 2:
        print("[ERROR] Please provide exactly two LHE files (SM and EFT).")
        for f in files: print("  -", f)
        sys.exit(1)
    for f in files:
        if not os.path.isfile(f):
            print("[ERROR] Missing:", f); sys.exit(1)

    # Parse both
    meta = []
    for f in files:
        sigma_pb, nevents = parse_sigma_nevents(f)
        mww = parse_MWW(f)
        cts, ctrs, widths, edges = hist_counts(mww)
        meta.append({
            'path': f, 'sigma': sigma_pb, 'N': nevents,
            'MWW': mww, 'counts': cts, 'centers': ctrs, 'widths': widths, 'edges': edges
        })

    # Ensure identical binning
    if not np.allclose(meta[0]['edges'], meta[1]['edges']):
        raise RuntimeError("Input files produced different bin edges; use identical MMIN/MMAX/NBINS.")

    # Identify SM/EFT roles
    role0 = auto_role(meta[0]['path'], 'SM')
    role1 = auto_role(meta[1]['path'], 'EFT' if role0=='SM' else 'SM')
    if role0 == role1:
        # fallback based on order
        role0, role1 = 'SM', 'EFT'
    roles = [role0, role1]

    # Map to indices
    iSM  = roles.index('SM')
    iEFT = roles.index('EFT')

    # Top panel: spectra
    fig = plt.figure(figsize=(8,7.5))
    gs  = GridSpec(nrows=2, ncols=1, height_ratios=[3.0, 1.2], hspace=0.05)

    ax_top = fig.add_subplot(gs[0])
    for i in [iSM, iEFT]:
        y = dsigma_from_counts(meta[i]['counts'], meta[i]['sigma'], meta[i]['N'], meta[i]['widths'])
        lab = os.path.basename(meta[i]['path'])
        lab2 = f"{roles[i]} ({lab})  σ={meta[i]['sigma']:.5g} pb, N={meta[i]['N']}"
        ax_top.step(meta[i]['centers'], y, where='mid', label=lab2)
        # Save CSV per spectrum
        out_csv = os.path.splitext(meta[i]['path'])[0] + "_dSigmadM.csv"
        np.savetxt(out_csv, np.column_stack([meta[i]['centers'], y]),
                   delimiter=",", header="M_WW[GeV], dSigma/dM[pb/GeV]", comments='')
        print(f"[OK] Wrote {out_csv}")

    if LOG_Y:
        ax_top.set_yscale('log')
        ax_top.set_ylim(1e-8, 1e-3)
    ax_top.set_ylabel(r"$\mathrm{d}\sigma/\mathrm{d}M_{WW}\ \mathrm{[pb/GeV]}$")
    ax_top.legend(loc='best')
    ax_top.tick_params(labelbottom=False)  # share-x; hide x labels on top

    # Bottom panel: ratio EFT/SM with stat errors
    ax_bot = fig.add_subplot(gs[1], sharex=ax_top)
    c_sm   = meta[iSM]['counts']
    c_eft  = meta[iEFT]['counts']
    # Avoid divide-by-zero: mask bins with c_sm==0 or c_eft==0
    mask = (c_sm > 0) & (c_eft > 0)
    ratio = np.full_like(c_sm, np.nan, dtype=float)
    dr    = np.full_like(c_sm, np.nan, dtype=float)

    const = (meta[iEFT]['sigma']/meta[iEFT]['N']) / (meta[iSM]['sigma']/meta[iSM]['N'])
    ratio[mask] = (c_eft[mask] / c_sm[mask]) * const
    # Poisson (independent) error propagation
    dr[mask] = ratio[mask] * np.sqrt(1.0/c_eft[mask] + 1.0/c_sm[mask])

    ax_bot.errorbar(meta[iSM]['centers'][mask], ratio[mask], yerr=dr[mask],
                    fmt='o', ms=3, lw=1, capsize=2, label="EFT / SM (stat. only)")
    ax_bot.axhline(1.0, color='k', lw=1, alpha=0.6)
    ax_bot.set_xlabel(r"$W_{\gamma\gamma}=M_{WW}\ \mathrm{[GeV]}$")
    ax_bot.set_ylabel("ratio")
    ax_bot.set_ylim(0.0, 2.0)  # adjust as you like
    ax_bot.grid(alpha=0.25)
    ax_bot.legend(loc='best')

    # Save ratio CSV
    out_ratio = "EFT_over_SM_ratio.csv"
    np.savetxt(out_ratio,
               np.column_stack([meta[iSM]['centers'], ratio, dr]),
               delimiter=",",
               header="M_WW[GeV], ratio(EFT/SM), stat_err", comments='')
    print(f"[OK] Wrote {out_ratio}")

    fig.tight_layout()
    out_png = "gamma-UPC-aaww-SM-EFT-chff_with_ratio.png"
    out_pdf = "gamma-UPC-aaww-SM-EFT-chff_with_ratio.pdf"
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    print(f"[OK] Wrote {out_png} and {out_pdf}")

if __name__ == "__main__":
    main()
