# ============================================================
# WCharm Parton / Particle Ratio Plotter
# (With statistical error bars on ratios)
# ============================================================

import ROOT
import os
import glob
import re

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.AddDirectory(False)

selected_weight = 0

histo_titles = {
    "ptlepton": "Lepton Transverse Momentum",
    "etalepton": "Lepton Pseudorapidity",
    "met_et": "Missing Transverse Energy",
    "mtw": "W Boson Transverse Mass",
    "ptjet": "Jet Transverse Momentum",
    "philepton": "Lepton Azimuthal Angle",
    "phijet": "Jet Azimuthal Angle",
}

x_axis_labels = {
    "ptlepton": "p_{T}^{lepton} [GeV]",
    "etalepton": "#eta^{lepton}",
    "met_et": "E_{T}^{miss} [GeV]",
    "mtw": "m_{T}^{W} [GeV]",
    "ptjet": "p_{T}^{jet} [GeV]",
    "philepton": "#phi^{lepton}",
    "phijet": "#phi^{jet}",
}

valid_pairs = [
    ("WCH7minus",      "WCHPartonminus"),
    ("WCH7plus",       "WCHPartonplus"),
    ("WCPy8minus",     "WCPyPartonminus"),
    ("WCPy8plus",      "WCPyPartonplus"),
]

pair_colours = {
    ("WCH7minus", "WCHPartonminus"): ROOT.kBlue,
    ("WCH7plus", "WCHPartonplus"): ROOT.kGreen+2,
    ("WCPy8minus", "WCPyPartonminus"): ROOT.kRed,
    ("WCPy8plus", "WCPyPartonplus"): ROOT.kMagenta+1,
}

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUTS_DIR = os.path.join(BASE_DIR, "outputs_weights")
PLOTS_DIR = os.path.join(BASE_DIR, "plots_ratios")
os.makedirs(PLOTS_DIR, exist_ok=True)

pattern = os.path.join(OUTPUTS_DIR, f"output_w{selected_weight}_*.root")
files = glob.glob(pattern)

if not files:
    print("No files found for weight", selected_weight)
    raise SystemExit(1)

proc_files = {}
weight_name = None

for f in files:
    base = os.path.basename(f)
    match = re.match(rf"output_w{selected_weight}_(.+?)_([^_]+)\.root$", base)

    if match:
        wname = match.group(1)
        process = match.group(2)
        proc_files[process] = f
        if weight_name is None:
            weight_name = wname

if weight_name is None:
    weight_name = f"w{selected_weight}"

print(f"Using weight index {selected_weight} ({weight_name})")

def load_histo(path, histname):
    f = ROOT.TFile(path, "READ")
    if not f or f.IsZombie():
        return None
    h = f.Get(histname)
    if not h:
        f.Close()
        return None
    h_clone = h.Clone(histname + "_clone")
    ROOT.SetOwnership(h_clone, False)
    f.Close()
    return h_clone

# ============================================================
# Main Loop
# ============================================================

for histo in histo_titles:

    canvas = ROOT.TCanvas(f"c_{histo}", f"c_{histo}", 900, 700)
    ROOT.SetOwnership(canvas, False)
    canvas.SetRightMargin(0.35)

    legend = ROOT.TLegend(0.68, 0.65, 0.95, 0.88)
    ROOT.SetOwnership(legend, False)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    first_draw = True
    drawn_any = False
    keep_alive = []

    for particle, parton in valid_pairs:

        if particle not in proc_files or parton not in proc_files:
            continue

        h_particle = load_histo(proc_files[particle], histo)
        h_parton   = load_histo(proc_files[parton], histo)

        if not h_particle or not h_parton:
            continue

        if h_particle.Integral() == 0:
            continue

        # --------------------------------------------------------
        # Bin width normalisation (IMPORTANT for variable bins)
        # This ensures histograms represent density, not raw counts
        # Prevents bias when dividing histograms with different bin widths
        # --------------------------------------------------------
        h_particle.Scale(1.0, "width")
        h_parton.Scale(1.0, "width")
        # --------------------------------------------------------

        ratio = h_parton.Clone(f"ratio_{parton}_over_{particle}")
        ROOT.SetOwnership(ratio, False)

        # --------------------------------------------------------
        # Compute ratio: Parton / Particle
        # --------------------------------------------------------
        ratio.Divide(h_particle)

        ratio.SetLineWidth(2)
        ratio.SetLineColor(pair_colours[(particle, parton)])
        ratio.SetMarkerStyle(20)
        ratio.SetMarkerColor(pair_colours[(particle, parton)])
        ratio.SetMinimum(0.0)
        ratio.SetMaximum(1.2)

        # --------------------------------------------------------
        # Auto axis tightening (remove empty space)
        # --------------------------------------------------------
        ratio.GetXaxis().SetLimits(
            ratio.GetXaxis().GetXmin(),
            ratio.GetXaxis().GetXmax()
        )

        max_val = ratio.GetMaximum()
        if max_val > 0:
            ratio.SetMaximum(1.2 * max_val)
        # --------------------------------------------------------

        ratio.SetTitle(
            f"{histo_titles[histo]} "
            f"(Weight index = {selected_weight}, {weight_name});"
            f"{x_axis_labels[histo]};Parton / Particle"
        )

        if first_draw:
            ratio.Draw("E1")
            first_draw = False
        else:
            ratio.Draw("E1 SAME")

        legend.AddEntry(
            ratio,
            f"{parton.replace('WCH','Herwig ').replace('WCPy','Pythia ')}"
            f" / {particle.replace('WCH','Herwig ').replace('WCPy','Pythia ')}",
            "lep"
        )

        keep_alive.extend([h_particle, h_parton, ratio])
        drawn_any = True

    if not drawn_any:
        continue

    legend.Draw()

    out_base = os.path.join(PLOTS_DIR, f"ratio_{histo}_w{selected_weight}")
    canvas.SaveAs(out_base + ".png")
    canvas.SaveAs(out_base + ".root")

print("All ratio plots complete.")
