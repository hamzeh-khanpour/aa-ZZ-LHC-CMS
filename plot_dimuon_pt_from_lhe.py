#!/usr/bin/env python3
"""
Read an LHE file produced by MadGraph/MadSpin, reconstruct the final-state
mu+ mu- system with ROOT.TLorentzVector, and plot:
  1) a unit-normalized pT(mu+mu-) distribution
  2) the differential distribution dSigma/dpT(mu+mu-)

Usage:
    python plot_dimuon_pt_from_lhe.py /path/to/decayed_events.lhe

Optional arguments:
    --nbins 40
    --xmin 0
    --xmax 2000
    --outdir plots_dimuon_pt

Notes
-----
* The script reads the total cross section after decay from the
  <MGGenerationInfo> block: "Integrated weight (pb)".
* For the differential distribution it uses:
      dSigma/dpT = Sigma_total * (1/N) * dN/dpT
  so the histogram integral reproduces Sigma_total (in pb).
* The unit-normalized shape is the probability density:
      (1/N) * dN/dpT
  and integrates to 1.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
from typing import Iterator, List, Optional, Tuple

try:
    import ROOT
except ImportError as exc:
    raise SystemExit(
        "ERROR: PyROOT is required for this script. Load your ROOT environment first, "
        "for example with 'source /path/to/thisroot.sh' or your local ROOT setup."
    ) from exc


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


Particle = Tuple[int, int, int, int, int, int, float, float, float, float, float, float, float]


def extract_total_xsec_pb(lhe_path: str) -> float:
    """Read the total post-decay cross section from the MGGenerationInfo block."""
    pattern = re.compile(r"Integrated weight \(pb\)\s*:\s*([+-]?[0-9]*\.?[0-9]+(?:[eE][+-]?\d+)?)")
    with open(lhe_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            match = pattern.search(line)
            if match:
                return float(match.group(1))
    raise RuntimeError(
        "Could not find 'Integrated weight (pb)' in the LHE header. "
        "Please check that the file contains the MGGenerationInfo block."
    )


def parse_event_blocks(lhe_path: str) -> Iterator[List[str]]:
    """Yield the raw lines inside each <event> ... </event> block."""
    in_event = False
    current: List[str] = []

    with open(lhe_path, "r", encoding="utf-8", errors="ignore") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if line == "<event>":
                in_event = True
                current = []
                continue
            if line == "</event>":
                if current:
                    yield current
                in_event = False
                current = []
                continue
            if in_event:
                if line and not line.startswith("#"):
                    current.append(line)


def parse_particle_line(line: str) -> Particle:
    parts = line.split()
    if len(parts) < 13:
        raise ValueError(f"Malformed particle line with {len(parts)} columns: {line}")
    return (
        int(parts[0]),   # IDUP
        int(parts[1]),   # ISTUP
        int(parts[2]),   # MOTHUP1
        int(parts[3]),   # MOTHUP2
        int(parts[4]),   # ICOLUP1
        int(parts[5]),   # ICOLUP2
        float(parts[6]), # PUP1 = px
        float(parts[7]), # PUP2 = py
        float(parts[8]), # PUP3 = pz
        float(parts[9]), # PUP4 = E
        float(parts[10]),# PUP5 = M
        float(parts[11]),# VTIMUP
        float(parts[12]) # SPINUP
    )


def tlv_from_particle(p: Particle) -> ROOT.TLorentzVector:
    vec = ROOT.TLorentzVector()
    vec.SetPxPyPzE(p[6], p[7], p[8], p[9])
    return vec


def extract_dimuon_pts(lhe_path: str) -> List[float]:
    """Loop over events, reconstruct mu+ mu- system, and return pT values."""
    dimuon_pts: List[float] = []

    for event_lines in parse_event_blocks(lhe_path):
        if not event_lines:
            continue

        header_cols = event_lines[0].split()
        if len(header_cols) < 1:
            continue

        try:
            nup = int(header_cols[0])
        except ValueError as exc:
            raise RuntimeError(f"Could not read NUP from event header: {event_lines[0]}") from exc

        particle_lines = event_lines[1:1 + nup]
        particles = [parse_particle_line(line) for line in particle_lines]

        mu_plus: Optional[ROOT.TLorentzVector] = None   # PDG ID = -13
        mu_minus: Optional[ROOT.TLorentzVector] = None  # PDG ID =  13

        for p in particles:
            pdg_id = p[0]
            status = p[1]
            if status != 1:
                continue
            if pdg_id == -13 and mu_plus is None:
                mu_plus = tlv_from_particle(p)
            elif pdg_id == 13 and mu_minus is None:
                mu_minus = tlv_from_particle(p)

        if mu_plus is None or mu_minus is None:
            # Skip events without a complete final-state dimuon pair.
            continue

        dimuon = mu_plus + mu_minus
        dimuon_pts.append(dimuon.Pt())

    if not dimuon_pts:
        raise RuntimeError(
            "No final-state mu+ mu- pairs were found in the LHE file. "
            "Please check that the file corresponds to the decayed sample."
        )

    return dimuon_pts


def choose_xmax(values: List[float], user_xmax: Optional[float]) -> float:
    if user_xmax is not None:
        return user_xmax
    vmax = max(values)
    if vmax <= 0:
        return 1.0
    return 1.05 * vmax


def make_histograms(
    pts: List[float],
    xsec_pb: float,
    nbins: int,
    xmin: float,
    xmax: float,
) -> Tuple[ROOT.TH1D, ROOT.TH1D]:
    """Create the normalized and differential histograms."""
    n_events = len(pts)
    if n_events == 0:
        raise ValueError("No pT values provided.")

    h_norm = ROOT.TH1D("h_norm", "", nbins, xmin, xmax)
    h_diff = ROOT.TH1D("h_diff", "", nbins, xmin, xmax)
    h_norm.Sumw2()
    h_diff.Sumw2()

    # Unit-normalized probability density: integral = 1
    weight_norm = 1.0 / n_events

    # Differential cross section density: integral = total cross section in pb
    weight_diff = xsec_pb / n_events

    for pt in pts:
        h_norm.Fill(pt, weight_norm)
        h_diff.Fill(pt, weight_diff)

    # Convert from per-bin weights to density by dividing by the bin width.
    h_norm.Scale(1.0, "width")
    h_diff.Scale(1.0, "width")

    return h_norm, h_diff


def style_hist(hist: ROOT.TH1D, x_title: str, y_title: str) -> None:
    hist.SetLineWidth(3)
    hist.SetTitle("")
    hist.GetXaxis().SetTitle(x_title)
    hist.GetYaxis().SetTitle(y_title)
    hist.GetXaxis().CenterTitle()
    hist.GetYaxis().CenterTitle()
    hist.GetXaxis().SetTitleSize(0.045)
    hist.GetYaxis().SetTitleSize(0.045)
    hist.GetXaxis().SetLabelSize(0.04)
    hist.GetYaxis().SetLabelSize(0.04)
    hist.GetYaxis().SetTitleOffset(1.3)


def draw_and_save(hist: ROOT.TH1D, out_base: str, extra_text: List[str]) -> None:
    canvas = ROOT.TCanvas(f"c_{hist.GetName()}", "", 900, 700)
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.05)
    canvas.SetBottomMargin(0.12)
    canvas.SetTopMargin(0.08)

    hist.Draw("hist e")

    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextSize(0.035)
    x0 = 0.18
    y0 = 0.88
    for i, text in enumerate(extra_text):
        latex.DrawLatex(x0, y0 - 0.045 * i, text)

    canvas.SaveAs(out_base + ".png")
    canvas.SaveAs(out_base + ".pdf")


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot pT(mu+mu-) from an LHE file.")
    parser.add_argument(
        "lhe_file",
        nargs="?",
        default="/home/hamzeh-khanpour/MG5_aMC_v3_6_6/AA_ZZ_madspin/Events/run_02_decayed_1/decayed_events.lhe",
        help="Path to the decayed LHE file.",
    )
    parser.add_argument("--nbins", type=int, default=40, help="Number of histogram bins.")
    parser.add_argument("--xmin", type=float, default=0.0, help="Minimum pT in GeV.")
    parser.add_argument("--xmax", type=float, default=None, help="Maximum pT in GeV.")
    parser.add_argument(
        "--outdir",
        type=str,
        default="plots_dimuon_pt",
        help="Directory where the output plots will be saved.",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.lhe_file):
        raise SystemExit(f"ERROR: File not found: {args.lhe_file}")

    os.makedirs(args.outdir, exist_ok=True)

    xsec_pb = extract_total_xsec_pb(args.lhe_file)
    pts = extract_dimuon_pts(args.lhe_file)
    xmax = choose_xmax(pts, args.xmax)

    h_norm, h_diff = make_histograms(
        pts=pts,
        xsec_pb=xsec_pb,
        nbins=args.nbins,
        xmin=args.xmin,
        xmax=xmax,
    )

    style_hist(
        h_norm,
        x_title="p_{T}(#mu^{+}#mu^{-}) [GeV]",
        y_title="1/N  dN/dp_{T}(#mu^{+}#mu^{-}) [GeV^{-1}]",
    )
    style_hist(
        h_diff,
        x_title="p_{T}(#mu^{+}#mu^{-}) [GeV]",
        y_title="d#sigma/dp_{T}(#mu^{+}#mu^{-}) [pb/GeV]",
    )

    info_lines = [
        f"Events with #mu^{{+}}#mu^{{-}} found: {len(pts)}",
        f"#sigma_{{after decay}} = {xsec_pb:.8e} pb",
    ]

    draw_and_save(h_norm, os.path.join(args.outdir, "pt_dimuon_normalized"), info_lines)
    draw_and_save(h_diff, os.path.join(args.outdir, "pt_dimuon_dsigma_dpT"), info_lines)

    # Also save the histograms in a ROOT file.
    root_out = ROOT.TFile(os.path.join(args.outdir, "pt_dimuon_histograms.root"), "RECREATE")
    h_norm.Write()
    h_diff.Write()
    root_out.Close()

    print("Done.")
    print(f"Input LHE file           : {args.lhe_file}")
    print(f"Events with dimuon pair  : {len(pts)}")
    print(f"Cross section after decay: {xsec_pb:.12e} pb")
    print(f"Output directory         : {os.path.abspath(args.outdir)}")
    print("Saved files:")
    print(f"  - {os.path.join(args.outdir, 'pt_dimuon_normalized.png')}")
    print(f"  - {os.path.join(args.outdir, 'pt_dimuon_normalized.pdf')}")
    print(f"  - {os.path.join(args.outdir, 'pt_dimuon_dsigma_dpT.png')}")
    print(f"  - {os.path.join(args.outdir, 'pt_dimuon_dsigma_dpT.pdf')}")
    print(f"  - {os.path.join(args.outdir, 'pt_dimuon_histograms.root')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
