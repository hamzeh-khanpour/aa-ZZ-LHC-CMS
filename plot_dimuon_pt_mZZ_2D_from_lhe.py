#!/usr/bin/env python3
"""
Read an LHE file produced by MadGraph/MadSpin, reconstruct:
  1) the final-state dimuon system mu+ mu-
  2) the diboson invariant mass m_ZZ from the two status-2 Z bosons
and plot:
  - 1D unit-normalized distributions for pT(mu+mu-) and mZZ
  - 1D differential distributions dSigma/dpT and dSigma/dmZZ
  - a 2D correlation plot of pT(mu+mu-) versus mZZ
    with the kinematic boundary pT <= 0.5*sqrt(mZZ^2 - 4 mZ^2)

Usage:
    python plot_dimuon_pt_mZZ_2D_from_lhe.py /path/to/decayed_events.lhe

Examples:
    python plot_dimuon_pt_mZZ_2D_from_lhe.py
    python plot_dimuon_pt_mZZ_2D_from_lhe.py /path/to/decayed_events.lhe \
        --pt-nbins 60 --mzz-nbins 60 --outdir plots_dimuon_pt_mZZ_2D

Notes
-----
* The script reads the total cross section after decay from the
  <MGGenerationInfo> block: "Integrated weight (pb)".
* For the differential distributions it uses:
      dSigma/dX = Sigma_total * (1/N) * dN/dX
  so the histogram integral reproduces Sigma_total (in pb).
* The 2D differential histogram uses:
      d^2 Sigma/(dmZZ dpT) = Sigma_total * (1/N) * d^2N/(dmZZ dpT)
  and is expressed in pb/GeV^2.
"""

from __future__ import annotations

import argparse
import math
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

MZ_GEV = 91.1876


def extract_total_xsec_pb(lhe_path: str) -> float:
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
        float(parts[6]), # px
        float(parts[7]), # py
        float(parts[8]), # pz
        float(parts[9]), # E
        float(parts[10]),# M
        float(parts[11]),# VTIMUP
        float(parts[12]) # SPINUP
    )


def tlv_from_particle(p: Particle) -> ROOT.TLorentzVector:
    vec = ROOT.TLorentzVector()
    vec.SetPxPyPzE(p[6], p[7], p[8], p[9])
    return vec


def extract_observables(lhe_path: str) -> Tuple[List[float], List[float], List[Tuple[float, float]]]:
    dimuon_pts: List[float] = []
    mzz_values: List[float] = []
    correlated: List[Tuple[float, float]] = []

    for event_lines in parse_event_blocks(lhe_path):
        if not event_lines:
            continue

        header_cols = event_lines[0].split()
        if not header_cols:
            continue

        try:
            nup = int(header_cols[0])
        except ValueError as exc:
            raise RuntimeError(f"Could not read NUP from event header: {event_lines[0]}") from exc

        particle_lines = event_lines[1:1 + nup]
        particles = [parse_particle_line(line) for line in particle_lines]

        mu_plus: Optional[ROOT.TLorentzVector] = None
        mu_minus: Optional[ROOT.TLorentzVector] = None
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

        pt_val: Optional[float] = None
        mzz_val: Optional[float] = None

        if mu_plus is not None and mu_minus is not None:
            dimuon = mu_plus + mu_minus
            pt_val = dimuon.Pt()
            dimuon_pts.append(pt_val)

        if len(z_bosons) >= 2:
            zz_system = z_bosons[0] + z_bosons[1]
            mzz_val = zz_system.M()
            mzz_values.append(mzz_val)

        if pt_val is not None and mzz_val is not None:
            correlated.append((mzz_val, pt_val))

    if not dimuon_pts:
        raise RuntimeError("No final-state mu+ mu- pairs were found in the LHE file.")
    if not mzz_values:
        raise RuntimeError("No status-2 Z boson pairs were found in the LHE file.")
    if not correlated:
        raise RuntimeError("No events were found with both a dimuon pair and a ZZ system.")

    return dimuon_pts, mzz_values, correlated


def choose_xmax(values: List[float], user_xmax: Optional[float]) -> float:
    if user_xmax is not None:
        return user_xmax
    vmax = max(values)
    if vmax <= 0:
        return 1.0
    return 1.05 * vmax


def make_histograms_1d(
    values: List[float],
    xsec_pb: float,
    nbins: int,
    xmin: float,
    xmax: float,
    stem: str,
) -> Tuple[ROOT.TH1D, ROOT.TH1D]:
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


def make_histogram_2d(
    pairs: List[Tuple[float, float]],
    xsec_pb: float,
    xbins: int,
    xmin: float,
    xmax: float,
    ybins: int,
    ymin: float,
    ymax: float,
) -> Tuple[ROOT.TH2D, ROOT.TH2D]:
    n_events = len(pairs)
    if n_events == 0:
        raise ValueError("No correlated (mZZ, pT) pairs were provided.")

    h2_norm = ROOT.TH2D("h2_mZZ_vs_pt_norm", "", xbins, xmin, xmax, ybins, ymin, ymax)
    h2_diff = ROOT.TH2D("h2_mZZ_vs_pt_diff", "", xbins, xmin, xmax, ybins, ymin, ymax)
    h2_norm.Sumw2()
    h2_diff.Sumw2()

    weight_norm = 1.0 / n_events
    weight_diff = xsec_pb / n_events

    for mzz, pt in pairs:
        h2_norm.Fill(mzz, pt, weight_norm)
        h2_diff.Fill(mzz, pt, weight_diff)

    # Convert to density per GeV^2.
    for ix in range(1, h2_norm.GetNbinsX() + 1):
        dx = h2_norm.GetXaxis().GetBinWidth(ix)
        for iy in range(1, h2_norm.GetNbinsY() + 1):
            dy = h2_norm.GetYaxis().GetBinWidth(iy)
            area = dx * dy
            if area <= 0:
                continue
            h2_norm.SetBinContent(ix, iy, h2_norm.GetBinContent(ix, iy) / area)
            h2_norm.SetBinError(ix, iy, h2_norm.GetBinError(ix, iy) / area)
            h2_diff.SetBinContent(ix, iy, h2_diff.GetBinContent(ix, iy) / area)
            h2_diff.SetBinError(ix, iy, h2_diff.GetBinError(ix, iy) / area)

    return h2_norm, h2_diff


def style_hist_1d(hist: ROOT.TH1D, x_title: str, y_title: str) -> None:
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


def style_hist_2d(hist: ROOT.TH2D, x_title: str, y_title: str, z_title: str) -> None:
    hist.SetTitle("")
    hist.GetXaxis().SetTitle(x_title)
    hist.GetYaxis().SetTitle(y_title)
    hist.GetZaxis().SetTitle(z_title)
    hist.GetXaxis().CenterTitle()
    hist.GetYaxis().CenterTitle()
    hist.GetZaxis().CenterTitle()
    hist.GetXaxis().SetTitleSize(0.045)
    hist.GetYaxis().SetTitleSize(0.045)
    hist.GetZaxis().SetTitleSize(0.04)
    hist.GetXaxis().SetLabelSize(0.04)
    hist.GetYaxis().SetLabelSize(0.04)
    hist.GetZaxis().SetLabelSize(0.035)
    hist.GetYaxis().SetTitleOffset(1.25)
    hist.GetZaxis().SetTitleOffset(1.25)


def make_kinematic_boundary_graph(xmin: float, xmax: float, npoints: int = 800) -> ROOT.TGraph:
    x_start = max(xmin, 2.0 * MZ_GEV)
    x_end = max(x_start + 1.0, xmax)
    graph = ROOT.TGraph()
    ip = 0
    for i in range(npoints):
        x = x_start + (x_end - x_start) * i / max(1, npoints - 1)
        arg = x * x - 4.0 * MZ_GEV * MZ_GEV
        if arg < 0:
            continue
        y = 0.5 * math.sqrt(arg)
        graph.SetPoint(ip, x, y)
        ip += 1
    graph.SetLineColor(ROOT.kRed + 1)
    graph.SetLineWidth(3)
    graph.SetLineStyle(2)
    return graph


def draw_and_save_1d(hist: ROOT.TH1D, out_base: str, extra_text: List[str]) -> None:
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


def draw_and_save_2d(
    hist: ROOT.TH2D,
    out_base: str,
    extra_text: List[str],
    boundary: ROOT.TGraph,
    logz: bool = False,
) -> None:
    canvas = ROOT.TCanvas(f"c_{hist.GetName()}", "", 1000, 800)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.16)
    canvas.SetBottomMargin(0.12)
    canvas.SetTopMargin(0.08)
    if logz:
        canvas.SetLogz(True)

    hist.Draw("COLZ")
    boundary.Draw("L SAME")

    legend = ROOT.TLegend(0.54, 0.80, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(boundary, r"p_{T}^{max} = 0.5#sqrt{m_{ZZ}^{2}-4m_{Z}^{2}}", "l")
    legend.Draw()

    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextSize(0.032)
    x0 = 0.15
    y0 = 0.88
    for i, text in enumerate(extra_text):
        latex.DrawLatex(x0, y0 - 0.040 * i, text)

    canvas.SaveAs(out_base + ".png")
    canvas.SaveAs(out_base + ".pdf")


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot pT(mu+mu-), mZZ, and their 2D correlation from an LHE file.")
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
        default="plots_dimuon_pt_mZZ_2D",
        help="Directory where the output plots will be saved.",
    )
    parser.add_argument(
        "--linear-z",
        action="store_true",
        help="Use linear z scale for the 2D plots (default is log z).",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.lhe_file):
        raise SystemExit(f"ERROR: File not found: {args.lhe_file}")

    os.makedirs(args.outdir, exist_ok=True)

    xsec_pb = extract_total_xsec_pb(args.lhe_file)
    dimuon_pts, mzz_values, correlated = extract_observables(args.lhe_file)

    pt_xmax = choose_xmax(dimuon_pts, args.pt_xmax)
    mzz_xmax = choose_xmax(mzz_values, args.mzz_xmax)

    h_pt_norm, h_pt_diff = make_histograms_1d(
        values=dimuon_pts,
        xsec_pb=xsec_pb,
        nbins=args.pt_nbins,
        xmin=args.pt_xmin,
        xmax=pt_xmax,
        stem="pt_dimuon",
    )
    h_mzz_norm, h_mzz_diff = make_histograms_1d(
        values=mzz_values,
        xsec_pb=xsec_pb,
        nbins=args.mzz_nbins,
        xmin=args.mzz_xmin,
        xmax=mzz_xmax,
        stem="mZZ",
    )
    h2_norm, h2_diff = make_histogram_2d(
        pairs=correlated,
        xsec_pb=xsec_pb,
        xbins=args.mzz_nbins,
        xmin=args.mzz_xmin,
        xmax=mzz_xmax,
        ybins=args.pt_nbins,
        ymin=args.pt_xmin,
        ymax=pt_xmax,
    )

    style_hist_1d(
        h_pt_norm,
        x_title="p_{T}(#mu^{+}#mu^{-}) [GeV]",
        y_title="1/N  dN/dp_{T}(#mu^{+}#mu^{-}) [GeV^{-1}]",
    )
    style_hist_1d(
        h_pt_diff,
        x_title="p_{T}(#mu^{+}#mu^{-}) [GeV]",
        y_title="d#sigma/dp_{T}(#mu^{+}#mu^{-}) [pb/GeV]",
    )
    style_hist_1d(
        h_mzz_norm,
        x_title="m_{ZZ} [GeV]",
        y_title="1/N  dN/dm_{ZZ} [GeV^{-1}]",
    )
    style_hist_1d(
        h_mzz_diff,
        x_title="m_{ZZ} [GeV]",
        y_title="d#sigma/dm_{ZZ} [pb/GeV]",
    )
    style_hist_2d(
        h2_norm,
        x_title="m_{ZZ} [GeV]",
        y_title="p_{T}(#mu^{+}#mu^{-}) [GeV]",
        z_title="1/N  d^{2}N/(dm_{ZZ} dp_{T}) [GeV^{-2}]",
    )
    style_hist_2d(
        h2_diff,
        x_title="m_{ZZ} [GeV]",
        y_title="p_{T}(#mu^{+}#mu^{-}) [GeV]",
        z_title="d^{2}#sigma/(dm_{ZZ} dp_{T}) [pb/GeV^{2}]",
    )

    info_lines = [
        f"Events with #mu^{{+}}#mu^{{-}} found: {len(dimuon_pts)}",
        f"Events with Z Z found: {len(mzz_values)}",
        f"Correlated events used: {len(correlated)}",
        f"#sigma_{{after decay}} = {xsec_pb:.8e} pb",
    ]

    boundary = make_kinematic_boundary_graph(args.mzz_xmin, mzz_xmax)

    draw_and_save_1d(h_pt_norm, os.path.join(args.outdir, "pt_dimuon_normalized"), info_lines)
    draw_and_save_1d(h_pt_diff, os.path.join(args.outdir, "pt_dimuon_dsigma_dpT"), info_lines)
    draw_and_save_1d(h_mzz_norm, os.path.join(args.outdir, "mZZ_normalized"), info_lines)
    draw_and_save_1d(h_mzz_diff, os.path.join(args.outdir, "mZZ_dsigma_dmZZ"), info_lines)
    draw_and_save_2d(
        h2_norm,
        os.path.join(args.outdir, "mZZ_vs_pt_dimuon_2D_normalized"),
        info_lines,
        boundary,
        logz=not args.linear_z,
    )
    draw_and_save_2d(
        h2_diff,
        os.path.join(args.outdir, "mZZ_vs_pt_dimuon_2D_d2sigma"),
        info_lines,
        boundary,
        logz=not args.linear_z,
    )

    root_out = ROOT.TFile(os.path.join(args.outdir, "dimuon_pt_mZZ_2D_histograms.root"), "RECREATE")
    h_pt_norm.Write()
    h_pt_diff.Write()
    h_mzz_norm.Write()
    h_mzz_diff.Write()
    h2_norm.Write()
    h2_diff.Write()
    boundary.Write("kinematic_boundary")
    root_out.Close()

    print("Done.")
    print(f"Input LHE file               : {args.lhe_file}")
    print(f"Events with dimuon pair      : {len(dimuon_pts)}")
    print(f"Events with ZZ system        : {len(mzz_values)}")
    print(f"Correlated events            : {len(correlated)}")
    print(f"Cross section after decay    : {xsec_pb:.12e} pb")
    print(f"Output directory             : {os.path.abspath(args.outdir)}")
    print("Saved files:")
    for name in [
        "pt_dimuon_normalized.png",
        "pt_dimuon_normalized.pdf",
        "pt_dimuon_dsigma_dpT.png",
        "pt_dimuon_dsigma_dpT.pdf",
        "mZZ_normalized.png",
        "mZZ_normalized.pdf",
        "mZZ_dsigma_dmZZ.png",
        "mZZ_dsigma_dmZZ.pdf",
        "mZZ_vs_pt_dimuon_2D_normalized.png",
        "mZZ_vs_pt_dimuon_2D_normalized.pdf",
        "mZZ_vs_pt_dimuon_2D_d2sigma.png",
        "mZZ_vs_pt_dimuon_2D_d2sigma.pdf",
        "dimuon_pt_mZZ_2D_histograms.root",
    ]:
        print(f"  - {os.path.join(args.outdir, name)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
