#!/usr/bin/env python3
"""
Read an LHE file produced by MadGraph/MadSpin, reconstruct event-level observables,
and make a small EFT-validity diagnostic suite:

  1) 1D unit-normalized distributions for pT(mu+mu-), mZZ, and sqrt(shat)
  2) 1D differential distributions dSigma/dpT, dSigma/dmZZ, dSigma/dsqrt(shat)
  3) a 2D correlation plot of pT(mu+mu-) versus mZZ with the kinematic boundary
  4) a profile <mZZ> versus pT(mu+mu-)
  5) survival/fraction curves F(mZZ > Lambda_cut) and F(sqrt(shat) > Lambda_cut)
  6) optional truncated comparisons after imposing a cut on mZZ and/or sqrt(shat)
  7) standalone clipped-only reco-level shapes for mZZ and pT(mu+mu-) when truncation is active

Usage examples:
    python plot_dimuon_pt_mZZ_EFTchecks_from_lhe.py

    python plot_dimuon_pt_mZZ_EFTchecks_from_lhe.py \
        /home/hamzeh-khanpour/MG5_aMC_v3_6_6/AA_ZZ_madspin/Events/run_02_decayed_1/decayed_events.lhe \
        --pt-nbins 60 --mzz-nbins 60 --shat-nbins 60 \
        --pt-xmax 6000 --mzz-xmax 12000 --shat-xmax 12000 \
        --outdir plots_EFT_checks

    python plot_dimuon_pt_mZZ_EFTchecks_from_lhe.py \
        --mzz-cut 4000 --truncate-on mzz

    python plot_dimuon_pt_mZZ_EFTchecks_from_lhe.py \
        --shat-cut 5000 --truncate-on shat

    python plot_dimuon_pt_mZZ_EFTchecks_from_lhe.py \
        --mzz-cut 4000 --shat-cut 5000 --truncate-on both

Notes
-----
* The script reads the post-decay total cross section from the LHE header:
      Integrated weight (pb)
* It assumes unweighted events. Differential distributions are formed with
      dSigma/dX = Sigma_total * (1/N) * dN/dX
* Truncated differential distributions use the same per-event weight Sigma_total/N_all,
  applied only to passing events, so their integral is the estimated truncated cross section.
* For a pure 2->2 diboson final state at parton level, mZZ and sqrt(shat) should closely track,
  but both are kept explicitly for EFT checks.
"""

from __future__ import annotations

import argparse
import math
import os
import re
from dataclasses import dataclass
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


@dataclass
class EventObs:
    pt_dimuon: float
    mzz: float
    sqrt_shat: float


def extract_total_xsec_pb(lhe_path: str) -> float:
    pattern = re.compile(r"Integrated weight \(pb\)\s*:\s*([+-]?[0-9]*\.?[0-9]+(?:[eE][+-]?\d+)?)")
    with open(lhe_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            match = pattern.search(line)
            if match:
                return float(match.group(1))
    raise RuntimeError("Could not find 'Integrated weight (pb)' in the LHE header.")


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
    v = ROOT.TLorentzVector()
    v.SetPxPyPzE(p[6], p[7], p[8], p[9])
    return v


def extract_events(lhe_path: str) -> List[EventObs]:
    events: List[EventObs] = []

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
        incoming: List[ROOT.TLorentzVector] = []

        for p in particles:
            pdg_id = p[0]
            status = p[1]
            vec = tlv_from_particle(p)

            if status == 1:
                if pdg_id == -13 and mu_plus is None:
                    mu_plus = vec
                elif pdg_id == 13 and mu_minus is None:
                    mu_minus = vec

            if status == 2 and pdg_id == 23:
                z_bosons.append(vec)

            if status == -1:
                incoming.append(vec)

        if mu_plus is None or mu_minus is None:
            continue
        if len(z_bosons) < 2:
            continue
        if len(incoming) < 2:
            continue

        dimuon = mu_plus + mu_minus
        zz_system = z_bosons[0] + z_bosons[1]
        incoming_system = incoming[0] + incoming[1]

        events.append(
            EventObs(
                pt_dimuon=dimuon.Pt(),
                mzz=zz_system.M(),
                sqrt_shat=max(0.0, incoming_system.M()),
            )
        )

    if not events:
        raise RuntimeError("No events were found with dimuon, ZZ, and incoming-beam information.")

    return events


def choose_xmax(values: List[float], user_xmax: Optional[float]) -> float:
    if user_xmax is not None:
        return user_xmax
    vmax = max(values)
    if vmax <= 0:
        return 1.0
    return 1.05 * vmax


def make_hist_1d(values: List[float], nbins: int, xmin: float, xmax: float, name: str) -> ROOT.TH1D:
    h = ROOT.TH1D(name, "", nbins, xmin, xmax)
    h.Sumw2()
    return h


def fill_hist_shape(hist: ROOT.TH1D, values: List[float], weight: float) -> None:
    for value in values:
        hist.Fill(value, weight)


def density_scale_1d(hist: ROOT.TH1D) -> None:
    hist.Scale(1.0, "width")


def make_hist_2d(
    pairs: List[Tuple[float, float]],
    xbins: int,
    xmin: float,
    xmax: float,
    ybins: int,
    ymin: float,
    ymax: float,
    name: str,
) -> ROOT.TH2D:
    h = ROOT.TH2D(name, "", xbins, xmin, xmax, ybins, ymin, ymax)
    h.Sumw2()
    return h


def density_scale_2d(hist: ROOT.TH2D) -> None:
    for ix in range(1, hist.GetNbinsX() + 1):
        dx = hist.GetXaxis().GetBinWidth(ix)
        for iy in range(1, hist.GetNbinsY() + 1):
            dy = hist.GetYaxis().GetBinWidth(iy)
            area = dx * dy
            if area > 0:
                hist.SetBinContent(ix, iy, hist.GetBinContent(ix, iy) / area)
                hist.SetBinError(ix, iy, hist.GetBinError(ix, iy) / area)


def style_hist_1d(hist: ROOT.TH1D, x_title: str, y_title: str, color: int = ROOT.kBlue + 1) -> None:
    hist.SetLineColor(color)
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


def make_profile_mzz_vs_pt(events: List[EventObs], nbins_pt: int, pt_xmin: float, pt_xmax: float) -> ROOT.TProfile:
    prof = ROOT.TProfile("prof_mZZ_vs_pt", "", nbins_pt, pt_xmin, pt_xmax)
    for ev in events:
        prof.Fill(ev.pt_dimuon, ev.mzz)
    prof.SetLineColor(ROOT.kMagenta + 2)
    prof.SetMarkerColor(ROOT.kMagenta + 2)
    prof.SetMarkerStyle(20)
    prof.SetMarkerSize(0.9)
    prof.SetLineWidth(2)
    prof.GetXaxis().SetTitle("p_{T}(#mu^{+}#mu^{-}) [GeV]")
    prof.GetYaxis().SetTitle("<#it{m}_{ZZ}> [GeV]")
    prof.GetXaxis().CenterTitle()
    prof.GetYaxis().CenterTitle()
    prof.GetXaxis().SetTitleSize(0.045)
    prof.GetYaxis().SetTitleSize(0.045)
    prof.GetXaxis().SetLabelSize(0.04)
    prof.GetYaxis().SetLabelSize(0.04)
    prof.GetYaxis().SetTitleOffset(1.25)
    return prof


def make_survival_hist(values: List[float], nbins: int, xmin: float, xmax: float, name: str) -> ROOT.TH1D:
    if not values:
        raise ValueError("No values provided for survival histogram.")
    h = ROOT.TH1D(name, "", nbins, xmin, xmax)
    sorted_values = sorted(values)
    n = len(sorted_values)
    for ibin in range(1, nbins + 1):
        thr = h.GetBinLowEdge(ibin)
        # fraction with value > thr
        idx = 0
        lo, hi = 0, n
        while lo < hi:
            mid = (lo + hi) // 2
            if sorted_values[mid] <= thr:
                lo = mid + 1
            else:
                hi = mid
        idx = lo
        frac = (n - idx) / n
        h.SetBinContent(ibin, frac)
        h.SetBinError(ibin, math.sqrt(max(frac * (1.0 - frac) / n, 0.0)))
    h.SetLineWidth(3)
    h.SetLineColor(ROOT.kGreen + 2)
    h.SetTitle("")
    return h


def passes_truncation(ev: EventObs, truncate_on: str, mzz_cut: Optional[float], shat_cut: Optional[float]) -> bool:
    if truncate_on == "none":
        return True
    if truncate_on == "mzz":
        return (mzz_cut is None) or (ev.mzz <= mzz_cut)
    if truncate_on == "shat":
        return (shat_cut is None) or (ev.sqrt_shat <= shat_cut)
    if truncate_on == "both":
        ok_mzz = (mzz_cut is None) or (ev.mzz <= mzz_cut)
        ok_shat = (shat_cut is None) or (ev.sqrt_shat <= shat_cut)
        return ok_mzz and ok_shat
    raise ValueError(f"Unknown truncate_on option: {truncate_on}")


def draw_text_box(extra_text: List[str], x0: float = 0.16, y0: float = 0.88, size: float = 0.032) -> ROOT.TLatex:
    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextSize(size)
    for i, text in enumerate(extra_text):
        latex.DrawLatex(x0, y0 - 0.040 * i, text)
    return latex


def draw_and_save_1d(hist: ROOT.TH1D, out_base: str, extra_text: List[str]) -> None:
    canvas = ROOT.TCanvas(f"c_{hist.GetName()}", "", 900, 700)
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.05)
    canvas.SetBottomMargin(0.12)
    canvas.SetTopMargin(0.08)
    hist.Draw("hist e")
    draw_text_box(extra_text, x0=0.18, y0=0.88, size=0.035)
    canvas.SaveAs(out_base + ".png")
    canvas.SaveAs(out_base + ".pdf")


def draw_and_save_1d_overlay(
    hist_all: ROOT.TH1D,
    hist_cut: ROOT.TH1D,
    out_base: str,
    extra_text: List[str],
    label_all: str,
    label_cut: str,
) -> None:
    canvas = ROOT.TCanvas(f"c_{hist_all.GetName()}_{hist_cut.GetName()}", "", 900, 700)
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.05)
    canvas.SetBottomMargin(0.12)
    canvas.SetTopMargin(0.08)

    ymax = max(hist_all.GetMaximum(), hist_cut.GetMaximum()) * 1.25
    hist_all.SetMaximum(ymax)
    hist_all.Draw("hist e")
    hist_cut.Draw("hist e same")

    legend = ROOT.TLegend(0.58, 0.76, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hist_all, label_all, "l")
    legend.AddEntry(hist_cut, label_cut, "l")
    legend.Draw()

    draw_text_box(extra_text, x0=0.18, y0=0.72, size=0.032)
    canvas.SaveAs(out_base + ".png")
    canvas.SaveAs(out_base + ".pdf")


def draw_and_save_2d(hist: ROOT.TH2D, out_base: str, extra_text: List[str], boundary: ROOT.TGraph, logz: bool) -> None:
    canvas = ROOT.TCanvas(f"c_{hist.GetName()}", "", 1000, 800)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.16)
    canvas.SetBottomMargin(0.12)
    canvas.SetTopMargin(0.08)
    if logz:
        canvas.SetLogz(True)
    hist.Draw("COLZ")
    boundary.Draw("L SAME")
    legend = ROOT.TLegend(0.60, 0.80, 0.90, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(boundary, r"p_{T}^{max} = 0.5#sqrt{m_{ZZ}^{2}-4m_{Z}^{2}}", "l")
    legend.Draw()
    draw_text_box(extra_text, x0=0.15, y0=0.88, size=0.032)
    canvas.SaveAs(out_base + ".png")
    canvas.SaveAs(out_base + ".pdf")


def draw_and_save_profile(prof: ROOT.TProfile, out_base: str, extra_text: List[str]) -> None:
    canvas = ROOT.TCanvas(f"c_{prof.GetName()}", "", 900, 700)
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.05)
    canvas.SetBottomMargin(0.12)
    canvas.SetTopMargin(0.08)
    prof.Draw("E1")
    draw_text_box(extra_text, x0=0.18, y0=0.88, size=0.035)
    canvas.SaveAs(out_base + ".png")
    canvas.SaveAs(out_base + ".pdf")


def build_cut_label(truncate_on: str, mzz_cut: Optional[float], shat_cut: Optional[float]) -> str:
    if truncate_on == "none":
        return "No truncation"
    pieces: List[str] = []
    if truncate_on in ("mzz", "both") and mzz_cut is not None:
        pieces.append(f"truth-level clipping: m_{{ZZ}}^{{truth}} #leq #Lambda_{{clip}} = {mzz_cut:.0f} GeV")
    if truncate_on in ("shat", "both") and shat_cut is not None:
        pieces.append(f"truth-level clipping: #sqrt{{#hat{{s}}}} #leq #Lambda_{{clip}} = {shat_cut:.0f} GeV")
    return ", ".join(pieces) if pieces else "Truncated sample"


def main() -> int:
    parser = argparse.ArgumentParser(description="LHE EFT-validity diagnostic plots for pT(mu+mu-), mZZ, and sqrt(shat).")
    parser.add_argument(
        "lhe_file",
        nargs="?",
        default="/home/hamzeh-khanpour/MG5_aMC_v3_6_6/AA_ZZ_madspin/Events/run_02_decayed_1/decayed_events.lhe",
        help="Path to the decayed LHE file.",
    )
    parser.add_argument("--pt-nbins", type=int, default=40)
    parser.add_argument("--pt-xmin", type=float, default=0.0)
    parser.add_argument("--pt-xmax", type=float, default=None)
    parser.add_argument("--mzz-nbins", type=int, default=40)
    parser.add_argument("--mzz-xmin", type=float, default=0.0)
    parser.add_argument("--mzz-xmax", type=float, default=None)
    parser.add_argument("--shat-nbins", type=int, default=40)
    parser.add_argument("--shat-xmin", type=float, default=0.0)
    parser.add_argument("--shat-xmax", type=float, default=None)
    parser.add_argument("--survival-nbins", type=int, default=60)
    parser.add_argument("--mzz-cut", type=float, default=None, help="Optional truncation cut on mZZ in GeV.")
    parser.add_argument("--shat-cut", type=float, default=None, help="Optional truncation cut on sqrt(shat) in GeV.")
    parser.add_argument(
        "--truncate-on",
        choices=["none", "mzz", "shat", "both"],
        default="none",
        help="Which truncation criterion to apply when making comparison plots.",
    )
    parser.add_argument(
        "--lambda-cut",
        type=float,
        default=None,
        help="Optional threshold for a single-number report of fractions above Lambda_cut.",
    )
    parser.add_argument("--linear-z", action="store_true", help="Use linear z scale for the 2D plots.")
    parser.add_argument("--outdir", type=str, default="plots_dimuon_pt_mZZ_EFTchecks")
    args = parser.parse_args()

    if not os.path.isfile(args.lhe_file):
        raise SystemExit(f"ERROR: File not found: {args.lhe_file}")
    if args.truncate_on in ("mzz", "both") and args.mzz_cut is None:
        raise SystemExit("ERROR: --truncate-on mzz/both requires --mzz-cut.")
    if args.truncate_on in ("shat", "both") and args.shat_cut is None:
        raise SystemExit("ERROR: --truncate-on shat/both requires --shat-cut.")

    os.makedirs(args.outdir, exist_ok=True)

    xsec_pb = extract_total_xsec_pb(args.lhe_file)
    events = extract_events(args.lhe_file)
    n_all = len(events)

    pt_values = [ev.pt_dimuon for ev in events]
    mzz_values = [ev.mzz for ev in events]
    shat_values = [ev.sqrt_shat for ev in events]
    correlated = [(ev.mzz, ev.pt_dimuon) for ev in events]

    pt_xmax = choose_xmax(pt_values, args.pt_xmax)
    mzz_xmax = choose_xmax(mzz_values, args.mzz_xmax)
    shat_xmax = choose_xmax(shat_values, args.shat_xmax)

    # Full-sample 1D histograms
    h_pt_norm = make_hist_1d(pt_values, args.pt_nbins, args.pt_xmin, pt_xmax, "h_pt_norm")
    h_pt_diff = make_hist_1d(pt_values, args.pt_nbins, args.pt_xmin, pt_xmax, "h_pt_diff")
    h_mzz_norm = make_hist_1d(mzz_values, args.mzz_nbins, args.mzz_xmin, mzz_xmax, "h_mzz_norm")
    h_mzz_diff = make_hist_1d(mzz_values, args.mzz_nbins, args.mzz_xmin, mzz_xmax, "h_mzz_diff")
    h_shat_norm = make_hist_1d(shat_values, args.shat_nbins, args.shat_xmin, shat_xmax, "h_shat_norm")
    h_shat_diff = make_hist_1d(shat_values, args.shat_nbins, args.shat_xmin, shat_xmax, "h_shat_diff")

    fill_hist_shape(h_pt_norm, pt_values, 1.0 / n_all)
    fill_hist_shape(h_pt_diff, pt_values, xsec_pb / n_all)
    fill_hist_shape(h_mzz_norm, mzz_values, 1.0 / n_all)
    fill_hist_shape(h_mzz_diff, mzz_values, xsec_pb / n_all)
    fill_hist_shape(h_shat_norm, shat_values, 1.0 / n_all)
    fill_hist_shape(h_shat_diff, shat_values, xsec_pb / n_all)
    for h in [h_pt_norm, h_pt_diff, h_mzz_norm, h_mzz_diff, h_shat_norm, h_shat_diff]:
        density_scale_1d(h)

    # 2D histograms
    h2_norm = make_hist_2d(correlated, args.mzz_nbins, args.mzz_xmin, mzz_xmax, args.pt_nbins, args.pt_xmin, pt_xmax, "h2_norm")
    h2_diff = make_hist_2d(correlated, args.mzz_nbins, args.mzz_xmin, mzz_xmax, args.pt_nbins, args.pt_xmin, pt_xmax, "h2_diff")
    for mzz, pt in correlated:
        h2_norm.Fill(mzz, pt, 1.0 / n_all)
        h2_diff.Fill(mzz, pt, xsec_pb / n_all)
    density_scale_2d(h2_norm)
    density_scale_2d(h2_diff)

    # Profile and survival plots
    prof_mzz_vs_pt = make_profile_mzz_vs_pt(events, args.pt_nbins, args.pt_xmin, pt_xmax)
    h_surv_mzz = make_survival_hist(mzz_values, args.survival_nbins, args.mzz_xmin, mzz_xmax, "h_surv_mzz")
    h_surv_shat = make_survival_hist(shat_values, args.survival_nbins, args.shat_xmin, shat_xmax, "h_surv_shat")

    # Styles
    style_hist_1d(h_pt_norm, "p_{T}(#mu^{+}#mu^{-}) [GeV]", "1/N  dN/dp_{T}(#mu^{+}#mu^{-}) [GeV^{-1}]", ROOT.kBlue + 1)
    style_hist_1d(h_pt_diff, "p_{T}(#mu^{+}#mu^{-}) [GeV]", "d#sigma/dp_{T}(#mu^{+}#mu^{-}) [pb/GeV]", ROOT.kBlue + 1)
    style_hist_1d(h_mzz_norm, "m_{ZZ} [GeV]", "1/N  dN/dm_{ZZ} [GeV^{-1}]", ROOT.kBlue + 1)
    style_hist_1d(h_mzz_diff, "m_{ZZ} [GeV]", "d#sigma/dm_{ZZ} [pb/GeV]", ROOT.kBlue + 1)
    style_hist_1d(h_shat_norm, "#sqrt{#hat{s}} [GeV]", "1/N  dN/d#sqrt{#hat{s}} [GeV^{-1}]", ROOT.kBlue + 1)
    style_hist_1d(h_shat_diff, "#sqrt{#hat{s}} [GeV]", "d#sigma/d#sqrt{#hat{s}} [pb/GeV]", ROOT.kBlue + 1)
    style_hist_2d(h2_norm, "m_{ZZ} [GeV]", "p_{T}(#mu^{+}#mu^{-}) [GeV]", "1/N  d^{2}N/(dm_{ZZ} dp_{T}) [GeV^{-2}]")
    style_hist_2d(h2_diff, "m_{ZZ} [GeV]", "p_{T}(#mu^{+}#mu^{-}) [GeV]", "d^{2}#sigma/(dm_{ZZ} dp_{T}) [pb/GeV^{2}]")
    style_hist_1d(h_surv_mzz, "#Lambda_{cut} [GeV]", "Fraction with m_{ZZ} > #Lambda_{cut}", ROOT.kGreen + 2)
    style_hist_1d(h_surv_shat, "#Lambda_{cut} [GeV]", "Fraction with #sqrt{#hat{s}} > #Lambda_{cut}", ROOT.kOrange + 7)

    info_lines = [
        f"Correlated events used: {n_all}",
        f"#sigma_{{after decay}} = {xsec_pb:.8e} pb",
    ]

    boundary = make_kinematic_boundary_graph(args.mzz_xmin, mzz_xmax)

    # Save baseline plots
    draw_and_save_1d(h_pt_norm, os.path.join(args.outdir, "pt_dimuon_normalized"), info_lines)
    draw_and_save_1d(h_pt_diff, os.path.join(args.outdir, "pt_dimuon_dsigma_dpT"), info_lines)
    draw_and_save_1d(h_mzz_norm, os.path.join(args.outdir, "mZZ_normalized"), info_lines)
    draw_and_save_1d(h_mzz_diff, os.path.join(args.outdir, "mZZ_dsigma_dmZZ"), info_lines)
    draw_and_save_1d(h_shat_norm, os.path.join(args.outdir, "sqrtshat_normalized"), info_lines)
    draw_and_save_1d(h_shat_diff, os.path.join(args.outdir, "sqrtshat_dsigma_dsqrtshat"), info_lines)
    draw_and_save_2d(h2_norm, os.path.join(args.outdir, "mZZ_vs_pt_dimuon_2D_normalized"), info_lines, boundary, logz=not args.linear_z)
    draw_and_save_2d(h2_diff, os.path.join(args.outdir, "mZZ_vs_pt_dimuon_2D_d2sigma"), info_lines, boundary, logz=not args.linear_z)
    draw_and_save_profile(prof_mzz_vs_pt, os.path.join(args.outdir, "profile_mZZ_vs_pt_dimuon"), info_lines)
    draw_and_save_1d(h_surv_mzz, os.path.join(args.outdir, "fraction_mZZ_gt_Lambdacut"), info_lines)
    draw_and_save_1d(h_surv_shat, os.path.join(args.outdir, "fraction_sqrtshat_gt_Lambdacut"), info_lines)

    # Optional truncation comparisons
    summary_lines: List[str] = []
    h_pt_norm_cut = None
    h_pt_diff_cut = None
    h_mzz_norm_cut = None
    h_mzz_diff_cut = None
    cut_label = build_cut_label(args.truncate_on, args.mzz_cut, args.shat_cut)
    if args.truncate_on != "none":
        passed_events = [ev for ev in events if passes_truncation(ev, args.truncate_on, args.mzz_cut, args.shat_cut)]
        n_pass = len(passed_events)
        frac_pass = n_pass / n_all
        sigma_pass_pb = xsec_pb * frac_pass
        summary_lines.extend([
            f"Clipping mode: {args.truncate_on}",
            f"Truth-level clipping used: {cut_label}",
            f"Passing clipped events: {n_pass} / {n_all} = {frac_pass:.6f}",
            f"Estimated clipped cross section: {sigma_pass_pb:.8e} pb",
        ])

        if n_pass == 0:
            raise SystemExit("ERROR: No events survive the requested truncation cut.")

        pt_pass = [ev.pt_dimuon for ev in passed_events]
        mzz_pass = [ev.mzz for ev in passed_events]

        h_pt_norm_cut = make_hist_1d(pt_pass, args.pt_nbins, args.pt_xmin, pt_xmax, "h_pt_norm_cut")
        h_pt_diff_cut = make_hist_1d(pt_pass, args.pt_nbins, args.pt_xmin, pt_xmax, "h_pt_diff_cut")
        h_mzz_norm_cut = make_hist_1d(mzz_pass, args.mzz_nbins, args.mzz_xmin, mzz_xmax, "h_mzz_norm_cut")
        h_mzz_diff_cut = make_hist_1d(mzz_pass, args.mzz_nbins, args.mzz_xmin, mzz_xmax, "h_mzz_diff_cut")

        fill_hist_shape(h_pt_norm_cut, pt_pass, 1.0 / n_pass)
        fill_hist_shape(h_pt_diff_cut, pt_pass, xsec_pb / n_all)
        fill_hist_shape(h_mzz_norm_cut, mzz_pass, 1.0 / n_pass)
        fill_hist_shape(h_mzz_diff_cut, mzz_pass, xsec_pb / n_all)
        for h in [h_pt_norm_cut, h_pt_diff_cut, h_mzz_norm_cut, h_mzz_diff_cut]:
            density_scale_1d(h)

        style_hist_1d(h_pt_norm_cut, "p_{T}(#mu^{+}#mu^{-}) [GeV]", "1/N  dN/dp_{T}(#mu^{+}#mu^{-}) [GeV^{-1}]", ROOT.kRed + 1)
        style_hist_1d(h_pt_diff_cut, "p_{T}(#mu^{+}#mu^{-}) [GeV]", "d#sigma/dp_{T}(#mu^{+}#mu^{-}) [pb/GeV]", ROOT.kRed + 1)
        style_hist_1d(h_mzz_norm_cut, "m_{ZZ} [GeV]", "1/N  dN/dm_{ZZ} [GeV^{-1}]", ROOT.kRed + 1)
        style_hist_1d(h_mzz_diff_cut, "m_{ZZ} [GeV]", "d#sigma/dm_{ZZ} [pb/GeV]", ROOT.kRed + 1)

        trunc_info = info_lines + [
            f"Passing clipped fraction = {frac_pass:.6f}",
            f"#sigma_{{clip}} = {sigma_pass_pb:.8e} pb",
        ]
        draw_and_save_1d_overlay(
            h_pt_norm,
            h_pt_norm_cut,
            os.path.join(args.outdir, "pt_dimuon_normalized_truncated_overlay"),
            trunc_info,
            "Full sample",
            cut_label,
        )
        draw_and_save_1d_overlay(
            h_pt_diff,
            h_pt_diff_cut,
            os.path.join(args.outdir, "pt_dimuon_dsigma_dpT_truncated_overlay"),
            trunc_info,
            "Full sample",
            cut_label,
        )
        draw_and_save_1d_overlay(
            h_mzz_norm,
            h_mzz_norm_cut,
            os.path.join(args.outdir, "mZZ_normalized_truncated_overlay"),
            trunc_info,
            "Full sample",
            cut_label,
        )
        draw_and_save_1d_overlay(
            h_mzz_diff,
            h_mzz_diff_cut,
            os.path.join(args.outdir, "mZZ_dsigma_dmZZ_truncated_overlay"),
            trunc_info,
            "Full sample",
            cut_label,
        )

        # Standalone clipped-only reco-level shapes (CMS-style: clip first, then inspect reco distributions)
        clipped_info = [
            f"Correlated clipped events used: {n_pass}",
            f"{cut_label}",
            f"#sigma_{{clip}} = {sigma_pass_pb:.8e} pb",
        ]
        draw_and_save_1d(
            h_pt_norm_cut,
            os.path.join(args.outdir, "pt_dimuon_normalized_clipped_only"),
            clipped_info,
        )
        draw_and_save_1d(
            h_mzz_norm_cut,
            os.path.join(args.outdir, "mZZ_normalized_clipped_only"),
            clipped_info,
        )
    else:
        summary_lines.append("No truncation cut requested.")

    # Single-threshold fractions if requested
    if args.lambda_cut is not None:
        frac_mzz_gt = sum(1 for v in mzz_values if v > args.lambda_cut) / n_all
        frac_shat_gt = sum(1 for v in shat_values if v > args.lambda_cut) / n_all
        summary_lines.extend([
            f"Fraction with mZZ > {args.lambda_cut:.3f} GeV: {frac_mzz_gt:.8f}",
            f"Fraction with sqrt(shat) > {args.lambda_cut:.3f} GeV: {frac_shat_gt:.8f}",
        ])

    summary_path = os.path.join(args.outdir, "eft_truncation_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as fout:
        fout.write(f"Input LHE file: {args.lhe_file}\n")
        fout.write(f"Correlated events used: {n_all}\n")
        fout.write(f"Cross section after decay [pb]: {xsec_pb:.12e}\n")
        fout.write("\n")
        for line in summary_lines:
            fout.write(line + "\n")

    # ROOT output
    root_out = ROOT.TFile(os.path.join(args.outdir, "eft_diagnostic_histograms.root"), "RECREATE")
    root_objects = [
        h_pt_norm, h_pt_diff,
        h_mzz_norm, h_mzz_diff,
        h_shat_norm, h_shat_diff,
        h2_norm, h2_diff,
        prof_mzz_vs_pt,
        h_surv_mzz, h_surv_shat,
        boundary,
    ]
    for obj in [h_pt_norm_cut, h_pt_diff_cut, h_mzz_norm_cut, h_mzz_diff_cut]:
        if obj is not None:
            root_objects.append(obj)
    for obj in root_objects:
        obj.Write()
    root_out.Close()

    print("Done.")
    print(f"Input LHE file                : {args.lhe_file}")
    print(f"Correlated events             : {n_all}")
    print(f"Cross section after decay     : {xsec_pb:.12e} pb")
    print(f"Output directory              : {os.path.abspath(args.outdir)}")
    print(f"Summary file                  : {summary_path}")
    print("Saved baseline files include:")
    for name in [
        "pt_dimuon_normalized.pdf",
        "pt_dimuon_dsigma_dpT.pdf",
        "mZZ_normalized.pdf",
        "mZZ_dsigma_dmZZ.pdf",
        "sqrtshat_normalized.pdf",
        "sqrtshat_dsigma_dsqrtshat.pdf",
        "mZZ_vs_pt_dimuon_2D_normalized.pdf",
        "mZZ_vs_pt_dimuon_2D_d2sigma.pdf",
        "profile_mZZ_vs_pt_dimuon.pdf",
        "fraction_mZZ_gt_Lambdacut.pdf",
        "fraction_sqrtshat_gt_Lambdacut.pdf",
        "eft_diagnostic_histograms.root",
        "eft_truncation_summary.txt",
    ]:
        print(f"  - {os.path.join(args.outdir, name)}")
    if args.truncate_on != "none":
        print("Saved truth-level clipped files include:")
        for name in [
            "pt_dimuon_normalized_truncated_overlay.pdf",
            "pt_dimuon_dsigma_dpT_truncated_overlay.pdf",
            "mZZ_normalized_truncated_overlay.pdf",
            "mZZ_dsigma_dmZZ_truncated_overlay.pdf",
            "pt_dimuon_normalized_clipped_only.pdf",
            "mZZ_normalized_clipped_only.pdf",
        ]:
            print(f"  - {os.path.join(args.outdir, name)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
