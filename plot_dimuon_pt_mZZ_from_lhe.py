#!/usr/bin/env python3
"""
Read an LHE file produced by MadGraph/MadSpin, reconstruct:
  1) the final-state dimuon system mu+ mu-
  2) the diboson invariant mass m_ZZ from the two status-2 Z bosons
and plot for each observable:
  - a unit-normalized distribution
  - the differential distribution in pb/GeV

Usage:
    python plot_dimuon_pt_mZZ_from_lhe.py /path/to/decayed_events.lhe

Examples:
    python plot_dimuon_pt_mZZ_from_lhe.py
    python plot_dimuon_pt_mZZ_from_lhe.py /path/to/decayed_events.lhe \
        --pt-nbins 50 --mzz-nbins 50 --outdir plots_dimuon_pt_mZZ

Notes
-----
* The script reads the total cross section after decay from the
  <MGGenerationInfo> block: "Integrated weight (pb)".
* For the differential distributions it uses:
      dSigma/dX = Sigma_total * (1/N) * dN/dX
  so the histogram integral reproduces Sigma_total (in pb).
* The unit-normalized shapes are probability densities and integrate to 1.
"""

from __future__ import annotations

import argparse
import os
import re
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
            if in_event and line and not line.startswith("#"):
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


def extract_observables(lhe_path: str) -> Tuple[List[float], List[float]]:
    """Loop over events and return pT(mu+mu-) and mZZ values."""
    dimuon_pts: List[float] = []
    mzz_values: List[float] = []

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
        z_bosons: List[ROOT.TLorentzVector] = []

        for p in particles:
            pdg_id = p[0]
            status = p[1]

            if status == 1:
                if pdg_id == -13 and mu_plus is None:
                    mu_plus = tlv_from_particle(p)
                elif pdg_id == 13 and mu_minus is None:
                    mu_minus = tlv_from_particle(p)

            if status == 2 and pdg_id == 23:
                z_bosons.append(tlv_from_particle(p))

        if mu_plus is not None and mu_minus is not None:
            dimuon = mu_plus + mu_minus
            dimuon_pts.append(dimuon.Pt())

        if len(z_bosons) >= 2:
            zz_system = z_bosons[0] + z_bosons[1]
            mzz_values.append(zz_system.M())

    if not dimuon_pts:
        raise RuntimeError(
            "No final-state mu+ mu- pairs were found in the LHE file. "
            "Please check that the file corresponds to the decayed sample."
        )

    if not mzz_values:
        raise RuntimeError(
            "No status-2 Z boson pairs were found in the LHE file. "
            "Please check that the file corresponds to a decayed ZZ sample."
        )

    return dimuon_pts, mzz_values


def choose_xmax(values: List[float], user_xmax: Optional[float]) -> float:
    if user_xmax is not None:
        return user_xmax
    vmax = max(values)
    if vmax <= 0:
        return 1.0
    return 1.05 * vmax


def make_histograms(
    values: List[float],
    xsec_pb: float,
    nbins: int,
    xmin: float,
    xmax: float,
    stem: str,
) -> Tuple[ROOT.TH1D, ROOT.TH1D]:
    """Create normalized and differential histograms for one observable."""
    n_events = len(values)
    if n_events == 0:
        raise ValueError(f"No values provided for observable '{stem}'.")

    h_norm = ROOT.TH1D(f"h_{stem}_norm", "", nbins, xmin, xmax)
    h_diff = ROOT.TH1D(f"h_{stem}_diff", "", nbins, xmin, xmax)
    h_norm.Sumw2()
    h_diff.Sumw2()

    weight_norm = 1.0 / n_events
    weight_diff = xsec_pb / n_events

    for value in values:
        h_norm.Fill(value, weight_norm)
        h_diff.Fill(value, weight_diff)

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
    parser = argparse.ArgumentParser(description="Plot pT(mu+mu-) and mZZ from an LHE file.")
    parser.add_argument(
        "lhe_file",
        nargs="?",
        default="/home/hamzeh-khanpour/MG5_aMC_v3_6_6/AA_ZZ_madspin/Events/run_02_decayed_1/decayed_events.lhe",
        help="Path to the decayed LHE file.",
    )
    parser.add_argument("--pt-nbins", type=int, default=40, help="Number of bins for pT(mu+mu-).")
    parser.add_argument("--pt-xmin", type=float, default=0.0, help="Minimum pT(mu+mu-) in GeV.")
    parser.add_argument("--pt-xmax", type=float, default=None, help="Maximum pT(mu+mu-) in GeV.")
    parser.add_argument("--mzz-nbins", type=int, default=40, help="Number of bins for mZZ.")
    parser.add_argument("--mzz-xmin", type=float, default=0.0, help="Minimum mZZ in GeV.")
    parser.add_argument("--mzz-xmax", type=float, default=None, help="Maximum mZZ in GeV.")
    parser.add_argument(
        "--outdir",
        type=str,
        default="plots_dimuon_pt_mZZ",
        help="Directory where the output plots will be saved.",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.lhe_file):
        raise SystemExit(f"ERROR: File not found: {args.lhe_file}")

    os.makedirs(args.outdir, exist_ok=True)

    xsec_pb = extract_total_xsec_pb(args.lhe_file)
    dimuon_pts, mzz_values = extract_observables(args.lhe_file)

    pt_xmax = choose_xmax(dimuon_pts, args.pt_xmax)
    mzz_xmax = choose_xmax(mzz_values, args.mzz_xmax)

    h_pt_norm, h_pt_diff = make_histograms(
        values=dimuon_pts,
        xsec_pb=xsec_pb,
        nbins=args.pt_nbins,
        xmin=args.pt_xmin,
        xmax=pt_xmax,
        stem="pt_dimuon",
    )
    h_mzz_norm, h_mzz_diff = make_histograms(
        values=mzz_values,
        xsec_pb=xsec_pb,
        nbins=args.mzz_nbins,
        xmin=args.mzz_xmin,
        xmax=mzz_xmax,
        stem="mZZ",
    )

    style_hist(
        h_pt_norm,
        x_title="p_{T}(#mu^{+}#mu^{-}) [GeV]",
        y_title="1/N  dN/dp_{T}(#mu^{+}#mu^{-}) [GeV^{-1}]",
    )
    style_hist(
        h_pt_diff,
        x_title="p_{T}(#mu^{+}#mu^{-}) [GeV]",
        y_title="d#sigma/dp_{T}(#mu^{+}#mu^{-}) [pb/GeV]",
    )
    style_hist(
        h_mzz_norm,
        x_title="m_{ZZ} [GeV]",
        y_title="1/N  dN/dm_{ZZ} [GeV^{-1}]",
    )
    style_hist(
        h_mzz_diff,
        x_title="m_{ZZ} [GeV]",
        y_title="d#sigma/dm_{ZZ} [pb/GeV]",
    )

    info_lines = [
        f"Events with #mu^{{+}}#mu^{{-}} found: {len(dimuon_pts)}",
        f"Events with Z Z found: {len(mzz_values)}",
        f"#sigma_{{after decay}} = {xsec_pb:.8e} pb",
    ]

    draw_and_save(h_pt_norm, os.path.join(args.outdir, "pt_dimuon_normalized"), info_lines)
    draw_and_save(h_pt_diff, os.path.join(args.outdir, "pt_dimuon_dsigma_dpT"), info_lines)
    draw_and_save(h_mzz_norm, os.path.join(args.outdir, "mZZ_normalized"), info_lines)
    draw_and_save(h_mzz_diff, os.path.join(args.outdir, "mZZ_dsigma_dmZZ"), info_lines)

    root_out = ROOT.TFile(os.path.join(args.outdir, "dimuon_pt_mZZ_histograms.root"), "RECREATE")
    h_pt_norm.Write()
    h_pt_diff.Write()
    h_mzz_norm.Write()
    h_mzz_diff.Write()
    root_out.Close()

    print("Done.")
    print(f"Input LHE file             : {args.lhe_file}")
    print(f"Events with dimuon pair    : {len(dimuon_pts)}")
    print(f"Events with ZZ system      : {len(mzz_values)}")
    print(f"Cross section after decay  : {xsec_pb:.12e} pb")
    print(f"Output directory           : {os.path.abspath(args.outdir)}")
    print("Saved files:")
    print(f"  - {os.path.join(args.outdir, 'pt_dimuon_normalized.png')}")
    print(f"  - {os.path.join(args.outdir, 'pt_dimuon_normalized.pdf')}")
    print(f"  - {os.path.join(args.outdir, 'pt_dimuon_dsigma_dpT.png')}")
    print(f"  - {os.path.join(args.outdir, 'pt_dimuon_dsigma_dpT.pdf')}")
    print(f"  - {os.path.join(args.outdir, 'mZZ_normalized.png')}")
    print(f"  - {os.path.join(args.outdir, 'mZZ_normalized.pdf')}")
    print(f"  - {os.path.join(args.outdir, 'mZZ_dsigma_dmZZ.png')}")
    print(f"  - {os.path.join(args.outdir, 'mZZ_dsigma_dmZZ.pdf')}")
    print(f"  - {os.path.join(args.outdir, 'dimuon_pt_mZZ_histograms.root')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
