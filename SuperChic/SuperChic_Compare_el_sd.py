#!/usr/bin/env python3
import argparse
import gzip
import math
import os
import sys
from typing import List, Tuple, Dict

import numpy as np
import matplotlib.pyplot as plt
import csv


def open_maybe_gzip(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def parse_init_block(lines_iter) -> Tuple[float, float]:
    """
    Parse the <init> block to get total cross section in pb (XSECUP) and error (XERRUP).
    SuperChic LHE writes:
      <init>
        2212 2212 ... [beams line]
        XSECUP  XERRUP  XMAXUP  LPRUP
      </init>
    Returns (xsec_pb, xerr_pb)
    """
    # Expect first numeric line (beams), then second numeric line with xsec.
    # We'll scan up to 10 lines to find two lines with >= 4 numbers.
    numeric_lines = []
    for _ in range(20):
        try:
            line = next(lines_iter)
        except StopIteration:
            break
        line = line.strip()
        if line.startswith('</init>'):
            break
        if not line or line.startswith('<'):
            continue
        parts = line.split()
        # Count how many parts parse as float
        ok = True
        floats = []
        for p in parts:
            try:
                floats.append(float(p.replace('D', 'E').replace('d', 'e')))
            except Exception:
                ok = False
                break
        if ok and len(floats) >= 4:
            numeric_lines.append(floats)
            if len(numeric_lines) >= 2:
                break
    if len(numeric_lines) < 2:
        raise ValueError("Could not parse <init> block for cross section.")
    xsec_pb = numeric_lines[1][0]
    xerr_pb = numeric_lines[1][1]
    return xsec_pb, xerr_pb


def invariant_mass(p4a, p4b) -> float:
    Ea, pxa, pya, pza = p4a
    Eb, pxb, pyb, pzb = p4b
    E = Ea + Eb
    px = pxa + pxb
    py = pya + pyb
    pz = pza + pzb
    m2 = E*E - (px*px + py*py + pz*pz)
    return math.sqrt(max(m2, 0.0))


def parse_event_block(lines_iter) -> Tuple[int, float, List[Dict]]:
    """
    Parse one <event> block.
    Returns (nup, xwgtup, particles)
    'particles' is a list of dicts with keys:
      id, status, moth1, moth2, col1, col2, px, py, pz, E, m
    """
    # Header line
    header = None
    for _ in range(10):
        line = next(lines_iter).strip()
        if line and not line.startswith('<'):
            header = line
            break
    if header is None:
        raise ValueError("Malformed <event> header.")
    h = header.split()
    if len(h) < 6:
        raise ValueError("Unexpected <event> header format: " + header)
    try:
        nup = int(float(h[0]))
        xwgtup = float(h[2].replace('D', 'E').replace('d', 'e'))
    except Exception as e:
        raise ValueError(f"Failed to parse NUP/XWGTUP from header: {header}") from e

    # Particle lines
    particles = []
    for _ in range(nup):
        line = next(lines_iter).strip()
        # Skip blank lines
        while not line:
            line = next(lines_iter).strip()
        parts = line.split()
        # Expected >= 13 columns, but SuperChic has at least 11 physics columns + extras.
        if len(parts) < 11:
            raise ValueError("Unexpected particle line: " + line)
        # Columns: IDUP ISTUP MOTH1 MOTH2 ICOL1 ICOL2 PX PY PZ E M ...
        def f(x): return float(x.replace('D', 'E').replace('d', 'e'))
        pid = int(parts[0])
        status = int(parts[1])
        moth1 = int(parts[2])
        moth2 = int(parts[3])
        col1 = int(parts[4])
        col2 = int(parts[5])
        px = f(parts[6])
        py = f(parts[7])
        pz = f(parts[8])
        E  = f(parts[9])
        m  = f(parts[10])
        particles.append({
            'id': pid, 'status': status, 'moth1': moth1, 'moth2': moth2,
            'col1': col1, 'col2': col2, 'px': px, 'py': py, 'pz': pz, 'E': E, 'm': m
        })
    # Seek to </event>
    for _ in range(1000):
        line = next(lines_iter)
        if line.strip().startswith('</event>'):
            break
    return nup, xwgtup, particles


def auto_label(path: str) -> str:
    b = os.path.basename(path)
    # Strip common "evrec" prefix and extension
    if b.lower().startswith('evrec'):
        b = b[5:]
    if b.lower().endswith('.dat'):
        b = b[:-4]
    if b.lower().endswith('.lhe'):
        b = b[:-4]
    return b


def parse_lhe(filepath: str, m_source: str = 'status2') -> Dict:
    """
    Parse a SuperChic LHE file.
    m_source: 'status2' (prefer W with ISTUP==2), or 'last' (take last two Ws).
    Returns dict with keys:
      'file', 'xsec_pb', 'xerr_pb', 'nevents', 'mww', 'weights_raw', 'sumw'
    """
    xsec_pb = None
    xerr_pb = None
    mww = []
    evt_weights = []
    sumw = 0.0
    nevents = 0

    with open_maybe_gzip(filepath) as f:
        lines = iter(f)
        for line in lines:
            s = line.strip()
            if s.startswith('<init>'):
                xsec_pb, xerr_pb = parse_init_block(lines)
            elif s.startswith('<event>'):
                nevents += 1
                nup, xwgtup, parts = parse_event_block(lines)
                sumw += xwgtup
                # pick W bosons
                ws = [p for p in parts if abs(p['id']) == 24]
                if not ws or len(ws) < 2:
                    # No explicit Ws? Skip MWW for this event
                    evt_weights.append(xwgtup)
                    continue
                if m_source == 'status2':
                    ws2 = [p for p in ws if p['status'] == 2]
                    if len(ws2) >= 2:
                        ws = ws2[:2]
                    else:
                        ws = ws[-2:]
                else:
                    ws = ws[-2:]
                p4a = (ws[0]['E'], ws[0]['px'], ws[0]['py'], ws[0]['pz'])
                p4b = (ws[1]['E'], ws[1]['px'], ws[1]['py'], ws[1]['pz'])
                m = invariant_mass(p4a, p4b)
                mww.append(m)
                evt_weights.append(xwgtup)

    if xsec_pb is None:
        raise RuntimeError(f"Did not find <init> block cross section in {filepath}")

    return {
        'file': filepath,
        'xsec_pb': xsec_pb,
        'xerr_pb': xerr_pb,
        'nevents': nevents,
        'mww': np.array(mww, dtype=float),
        'weights_raw': np.array(evt_weights, dtype=float),
        'sumw': float(sumw)
    }


def make_histogram(masses, weights_fb, bins):
    hist, edges = np.histogram(masses, bins=bins, weights=weights_fb)
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = (edges[1:] - edges[:-1])
    # differential xsec in fb/GeV
    with np.errstate(divide='ignore', invalid='ignore'):
        diff = np.where(widths > 0, hist / widths, 0.0)
    return centers, diff, edges, hist, widths


def main():
    parser = argparse.ArgumentParser(description="Compare SuperChic LHE files (M_WW spectrum & cross sections).")
    parser.add_argument('files', nargs='*', help="LHE files (.dat or .lhe, plain or .gz).")
    parser.add_argument('--bins', type=int, default=60, help="Number of M_WW bins (default: 60).")
    parser.add_argument('--mmin', type=float, default=160.0, help="Min M_WW (GeV). Default 160.")
    parser.add_argument('--mmax', type=float, default=3000.0, help="Max M_WW (GeV). Default 3000.")
    parser.add_argument('--out', default='compare_MWW', help="Output prefix for plots (default: compare_MWW).")
    parser.add_argument('--source', default='status2', choices=['status2','last'], help="How to pick Ws (default: status2).")
    args = parser.parse_args()

    # If no files passed, try common names in cwd
    files = args.files
    if not files:
        candidates = [
            'evrecww_el_13TeV_elWplus.dat',
            'evrecww_sd_13TeV_elWplus.dat',
            os.path.join('evrecs', 'evrecww_el_13TeV_elWplus.dat'),
            os.path.join('evrecs', 'evrecww_sd_13TeV_elWplus.dat'),
        ]
        files = [c for c in candidates if os.path.isfile(c)]
        if not files:
            print("No files given and none found in defaults; pass paths on the command line.", file=sys.stderr)
            sys.exit(1)

    parsed = []
    for fpath in files:
        try:
            rec = parse_lhe(fpath, m_source=args.source)
            parsed.append(rec)
        except Exception as e:
            print(f"[WARN] Failed to parse {fpath}: {e}", file=sys.stderr)

    if not parsed:
        print("No valid LHE files parsed.", file=sys.stderr)
        sys.exit(2)

    # Prepare bins
    bins = np.linspace(args.mmin, args.mmax, args.bins + 1)

    # Plot
    plt.figure(figsize=(8, 5.5))
    total_sigma_fb = 0.0

    for rec in parsed:
        xsec_pb = rec['xsec_pb']
        xerr_pb = rec['xerr_pb']
        nevents = rec['nevents']
        mww = rec['mww']
        raww = rec['weights_raw']
        sumw = rec['sumw'] if rec['sumw'] > 0 else len(raww)

        # Normalize weights to total xsec
        xsec_fb = xsec_pb * 1e3  # pb -> fb
        total_sigma_fb += xsec_fb
        # Per-event weights in fb
        scale = xsec_fb / sumw if sumw > 0 else 0.0
        weights_fb = raww * scale

        # Histogram
        centers, diff, edges, hist, widths = make_histogram(mww, weights_fb, bins)

        # Save CSV per file
        csv_name = os.path.splitext(os.path.basename(rec['file']))[0] + "_MWW.csv"
        with open(csv_name, 'w', newline='') as csvfile:
            w = csv.writer(csvfile)
            w.writerow(['bin_left', 'bin_right', 'bin_center', 'differential_xsec_fb_per_GeV', 'bin_integral_fb'])
            for i in range(len(centers)):
                w.writerow([edges[i], edges[i+1], centers[i], diff[i], hist[i]])

        label = f"{auto_label(rec['file'])} (σ={xsec_fb:.3f} fb, N={nevents})"
        plt.step(centers, diff, where='mid', label=label)

        print(f"File: {rec['file']}")
        print(f"  σ = {xsec_pb:.6e} pb = {xsec_fb:.6f} fb  (err {xerr_pb:.6e} pb)")
        print(f"  N events = {nevents}, sum(XWGTUP) = {sumw:.6e}")
        print(f"  CSV written: {csv_name}")

    plt.yscale('log')
    plt.xlabel(r"$M_{WW}$  [GeV]")
    plt.ylabel(r"$\mathrm{d}\sigma/\mathrm{d}M_{WW}$  [fb/GeV]")
    plt.title("SuperChic γγ → W⁺W⁻)")
    plt.legend()
    plt.tight_layout()
    png = args.out + ".png"
    pdf = args.out + ".pdf"
    plt.savefig(png, dpi=140)
    plt.savefig(pdf)
    print(f"\nSaved plots: {png}, {pdf}")
    print(f"Sum of σ over inputs = {total_sigma_fb:.6f} fb")


if __name__ == "__main__":
    main()
