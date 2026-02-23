# ============================================================
# WCharm Histogram-to-Histogram Ratio Plotter
#
# Plots ONLY ratios between histogram values:
#   ratio(sample) = H_sample / H_reference
#
# Features:
#   - Same format as your ratio plotter
#   - Error bars shown (propagated by ROOT)
#   - Legend outside plot area
#   - Weight index AND weight name in title
# Each ratio plot will contain one line for every sample that is successfully loaded and divided by the chosen reference sample. Since there are eight possible samples in total, this means you would expect up to seven ratio curves (each sample divided by the reference), with the reference itself giving a flat line at 1 if it is included.
# So in general, the number of lines on each ratio histogram is equal to the number of available samples minus one (assuming one fixed reference).
# ============================================================

import ROOT
import os
import glob
import re

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.AddDirectory(False)

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------
selected_weight = 0

# Choose the denominator histogram (reference)
# Everything will be divided by this sample’s histogram.
reference_sample = "WCH7minus"

histo_titles = {
    "ptlepton": "Lepton Transverse Momentum",
    "etalepton": "Lepton Pseudorapidity",
    "met_et": "Missing Transverse Energy",
    "mtw": "W Boson Transverse Mass",
    "ptjet": "Jet Transverse Momentum",
}

x_axis_labels = {
    "ptlepton": "p_{T}^{lepton} [GeV]",
    "etalepton": "#eta^{lepton}",
    "met_et": "E_{T}^{miss} [GeV]",
    "mtw": "m_{T}^{W} [GeV]",
    "ptjet": "p_{T}^{jet} [GeV]",
}

# Which samples you want to consider (if present on disk)
samples = [
    "WCH7minus",
    "WCH7plus",
    "WCHPartonminus",
    "WCHPartonplus",
    "WCPy8minus",
    "WCPy8plus",
    "WCPyPartonminus",
    "WCPyPartonplus",
]

# Colours (kept consistent with your scheme)
sample_colours = {
    "WCH7minus": ROOT.kBlue,              # reference
    "WCH7plus": ROOT.kGreen+2,            # green (as requested earlier)
    "WCHPartonminus": ROOT.kBlue+3,
    "WCHPartonplus": ROOT.kGreen+3,
    "WCPy8minus": ROOT.kRed,
    "WCPy8plus": ROOT.kMagenta+1,
    "WCPyPartonminus": ROOT.kRed+2,
    "WCPyPartonplus": ROOT.kMagenta+3,
}

# ------------------------------------------------------------
# Directories
# ------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUTS_DIR = os.path.join(BASE_DIR, "outputs_weights")
PLOTS_DIR = os.path.join(BASE_DIR, "plots_uncerts")
os.makedirs(PLOTS_DIR, exist_ok=True)

pattern = os.path.join(OUTPUTS_DIR, f"output_w{selected_weight}_*.root")
files = glob.glob(pattern)

if not files:
    print("No files found for weight", selected_weight)
    raise SystemExit(1)

# ------------------------------------------------------------
# Map: process -> file, and extract weight name from filename
# ------------------------------------------------------------
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
print("Found processes:", sorted(proc_files.keys()))

# ------------------------------------------------------------
# Safe histogram loader
# ------------------------------------------------------------
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
# Main Loop: ratios vs reference_sample
# ============================================================
for histo in histo_titles:

    # Reference histogram must exist
    if reference_sample not in proc_files:
        print(f"Reference sample '{reference_sample}' not found for weight {selected_weight}.")
        continue

    h_ref = load_histo(proc_files[reference_sample], histo)
    if not h_ref or h_ref.Integral() == 0:
        print(f"Reference histogram missing/empty: {reference_sample} {histo}")
        continue

    canvas = ROOT.TCanvas(f"c_{histo}", f"c_{histo}", 900, 700)
    ROOT.SetOwnership(canvas, False)
    canvas.SetRightMargin(0.35)

    legend = ROOT.TLegend(0.68, 0.60, 0.95, 0.88)
    ROOT.SetOwnership(legend, False)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    first_draw = True
    keep_alive = [h_ref]

    for sample in samples:

        if sample not in proc_files:
            continue

        h = load_histo(proc_files[sample], histo)
        if not h or h.Integral() == 0:
            continue

        # ratio = sample / reference (propagates errors automatically)
        ratio = h.Clone(f"ratio_{sample}_over_{reference_sample}")
        ROOT.SetOwnership(ratio, False)
        ratio.Divide(h_ref)

        ratio.SetLineWidth(2)
        ratio.SetLineColor(sample_colours.get(sample, ROOT.kBlack))
        ratio.SetMarkerStyle(20)
        ratio.SetMarkerColor(sample_colours.get(sample, ROOT.kBlack))

        ratio.SetMinimum(0.0)
        ratio.SetMaximum(1.2)

        ratio.SetTitle(
            f"{histo_titles[histo]} "
            f"(Weight index = {selected_weight}, {weight_name});"
            f"{x_axis_labels[histo]};{sample} / {reference_sample}"
        )

        if first_draw:
            ratio.Draw("E1")
            first_draw = False
        else:
            ratio.Draw("E1 SAME")

        legend.AddEntry(ratio, f"{sample} / {reference_sample}", "lep")
        keep_alive.extend([h, ratio])

    if first_draw:
        # nothing drawn
        continue

    legend.Draw()

    out_base = os.path.join(PLOTS_DIR, f"ratio_{histo}_w{selected_weight}_ref_{reference_sample}")
    canvas.SaveAs(out_base + ".png")
    canvas.SaveAs(out_base + ".root")

print("All histogram-to-histogram ratio plots complete.")
