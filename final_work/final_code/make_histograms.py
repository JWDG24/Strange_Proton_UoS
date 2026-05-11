# ============================================================
# make_histograms.py
#
# Purpose:
#   Read WCharm ntuples and create one histogram ROOT file
#   for each valid (sample, weight index) combination.
#
# Behaviour:
#   - Missing weights are SKIPPED
#   - No fallback to weight 0
#   - No fake duplicate outputs
#   - One output file per real (sample, weight) pair
#
# Output histograms:
#   ptlepton, ptjet, etalepton, etajet,
#   philepton, phijet, met_et, mtw
#
# Notes:
#   - Uses ROOT-safe histogram ownership handling
#   - Uses explicit GetEntry(i) loops for better PyROOT behaviour
#   - Assumes weight vector length is consistent within a sample
# ============================================================

import ROOT
import sys
import os
from array import array

# ------------------------------------------------------------
# ROOT safety / batch mode
# ------------------------------------------------------------
# Prevent ROOT trying to pop up canvases
ROOT.gROOT.SetBatch(True)

# Prevent histograms being auto-registered in the current ROOT
# directory, which causes the repeated:
#   "Replacing existing TH1 ... (Potential memory leak)"
# warnings when making many histograms with the same names.
ROOT.TH1.AddDirectory(False)

# ------------------------------------------------------------
# User options
# ------------------------------------------------------------

# If True:
#   run all valid weights present in the sample
# If False:
#   run only the nominal/default weight 0
USE_ALL_WEIGHTS = True

# If True:
#   use variable eta bins
# If False:
#   use uniform eta bins
USE_VARIABLE_ETA_BINS = True

# ------------------------------------------------------------
# Locations
# ------------------------------------------------------------
DATA_DIR = "/mnt/c/Users/dugar/wcharm_analysis/data"
WEIGHT_DIR = "/mnt/c/Users/dugar/wcharm_analysis/weights"
OUTPUT_DIR = "/mnt/c/Users/dugar/wcharm_analysis/final_work/make_histograms.py_Outputs"

# Make sure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ------------------------------------------------------------
# Binning definitions
# ------------------------------------------------------------

# Lepton eta bins
eta_lepton_bins = sorted([0, 0.21, 0.42, 0.63, 0.84, 1.05, 1.37, 1.52, 1.74, 1.95, 2.18, 2.5])
eta_lepton_binArray = array('d', eta_lepton_bins)

# Missing transverse energy bins
met_bins = [
    0,
    10, 20, 25, 30, 35, 40, 45,
    50, 52, 55, 58, 60, 62, 65, 68, 70, 72, 75, 78, 80,
    85, 90, 100,
    120, 140, 160,
    200, 250, 300
]
met_binArray = array('d', met_bins)

# W transverse mass bins
mtw_bins = [
    0,
    10, 20, 30, 40,
    45, 50, 55, 58,
    60, 62, 65, 68,
    70, 72, 75, 78,
    80, 82, 85, 88, 90,
    95, 100, 110,
    120, 140, 160,
    200, 300
]
mtw_binArray = array('d', mtw_bins)

# Jet pT bins
ptjet_bins = [
    0,
    10, 20, 30, 40, 50, 60, 70, 80,
    90, 100, 110, 120,
    130, 140, 150,
    160, 170, 180,
    200, 220,
    240, 260,
    280, 300
]
ptjet_binArray = array('d', ptjet_bins)

# Lepton pT bins
ptlepton_bins = [
    0,
    6, 12, 18, 24, 30, 36, 42, 48, 54,
    60, 66, 72, 78, 84, 90, 96, 102, 108, 114,
    120,
    150, 180, 210, 240, 270, 300
]
ptlepton_binArray = array('d', ptlepton_bins)

# ------------------------------------------------------------
# Requested weights
# ------------------------------------------------------------
# The script will later keep only those that actually exist
# for a given sample.
ALL_REQUESTED_WEIGHTS = list(range(319))

# ------------------------------------------------------------
# Sample list
# ------------------------------------------------------------
SignalSamples_default = [
    "WCH7minus",
    "WCH7plus",
    "WCHPartonminus",
    "WCHPartonplus",
    "WCPy8minus",
    "WCPy8plus",
    "WCPyPartonminus",
    "WCPyPartonplus"
]

SignalSamples = sys.argv[1:] if len(sys.argv) > 1 else SignalSamples_default

# ============================================================
# Helper functions
# ============================================================

def load_weight_map_clean(path):
    """
    Read the cleaned weight mapping file.

    Expected format:
        0 Default
        1 MUR0.5_MUF0.5_PDF13300
        ...

    Returns:
        dict mapping weight index -> weight name
    """
    idx_to_name = {}

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            idx_str, name = line.split(None, 1)
            idx_to_name[int(idx_str)] = name.strip()

    return idx_to_name


def mapping_file_for_sample(sample):
    """
    Choose the correct weight mapping file based on the sample name.
    """
    if sample == "WCPyPartonplus":
        return os.path.join(WEIGHT_DIR, "WCPyPartonplus_weights_clean.txt")

    if sample == "WCPyPartonminus":
        return os.path.join(WEIGHT_DIR, "WCPyPartonminus_weights_clean.txt")

    if sample.startswith("WCPy8"):
        return os.path.join(WEIGHT_DIR, "weights.py.particle_clean.txt")

    if sample.startswith("WCH7"):
        return os.path.join(WEIGHT_DIR, "weights.hw.particle_clean.txt")

    if sample == "WCHPartonplus":
        return os.path.join(WEIGHT_DIR, "WCHPartonplus_weights_clean.txt")

    if sample == "WCHPartonminus":
        return os.path.join(WEIGHT_DIR, "WCHPartonminus_weights_clean.txt")

    raise RuntimeError(f"No mapping rule for sample '{sample}'")


def safe_name(text):
    """
    Make a string safe for use in filenames by replacing
    awkward characters with underscores.
    """
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in text)


def make_histograms():
    """
    Create a fresh set of histograms for one sample and one weight.

    Histograms are detached from ROOT directories via SetDirectory(0)
    to avoid ROOT ownership problems and repeated-name warnings.
    """
    hists = {}

    hists["ptlepton"] = ROOT.TH1F("ptlepton", "ptlepton", len(ptlepton_binArray) - 1, ptlepton_binArray)
    hists["ptjet"]    = ROOT.TH1F("ptjet", "ptjet", len(ptjet_binArray) - 1, ptjet_binArray)

    if USE_VARIABLE_ETA_BINS:
        hists["etalepton"] = ROOT.TH1F("etalepton", "etalepton", len(eta_lepton_binArray) - 1, eta_lepton_binArray)
        hists["etajet"]    = ROOT.TH1F("etajet", "etajet", len(eta_lepton_binArray) - 1, eta_lepton_binArray)
    else:
        hists["etalepton"] = ROOT.TH1F("etalepton", "etalepton", 12, 0, 2.5)
        hists["etajet"]    = ROOT.TH1F("etajet", "etajet", 12, 0, 2.5)

    hists["philepton"] = ROOT.TH1F("philepton", "philepton", 50, -1, 7)
    hists["phijet"]    = ROOT.TH1F("phijet", "phijet", 50, -1, 7)

    hists["met_et"] = ROOT.TH1F("met_et", "met_et", len(met_binArray) - 1, met_binArray)
    hists["mtw"]    = ROOT.TH1F("mtw", "mtw", len(mtw_binArray) - 1, mtw_binArray)

    for hist in hists.values():
        hist.Sumw2()
        hist.SetDirectory(0)

    return hists


def fill_histograms(tree, hists, weight_index, total_events):
    """
    Loop over all events and fill histograms using the chosen weight index.

    Important:
      We do NOT check len(weightvec) inside the event loop.
      The chosen weight_index has already been validated before this
      function is called, which is much faster in PyROOT.

    Returns:
      pos_eta_sum, neg_eta_sum
      purely as diagnostics.
    """
    pos_eta_sum = 0.0
    neg_eta_sum = 0.0

    nentries = tree.GetEntries()

    for ientry in range(nentries):
        tree.GetEntry(ientry)

        # Event weight:
        # 1 / total events for normalisation, multiplied by chosen weight variation
        weightxs = (1.0 / total_events) * tree.weightvec[weight_index]

        # Event selection cuts
        if tree.met_et <= 25:
            continue
        if tree.leptons_pt <= 20:
            continue
        if abs(tree.leptons_eta) >= 2.5:
            continue

        # Charge sign convention already used in your original script
        corr = -1.0 * tree.leptons_charge
        final_weight = corr * weightxs

        # Use absolute eta for the eta histograms
        eta_lep = abs(tree.leptons_eta)

        # Diagnostic positive/negative eta sums
        if tree.leptons_eta >= 0:
            pos_eta_sum += final_weight
        else:
            neg_eta_sum += final_weight

        # Fill lepton histograms
        hists["ptlepton"].Fill(tree.leptons_pt, final_weight)
        hists["philepton"].Fill(tree.leptons_phi, final_weight)
        hists["etalepton"].Fill(eta_lep, final_weight)

        # Fill MET histogram
        hists["met_et"].Fill(tree.met_et, final_weight)

        # Compute W transverse mass
        dphi = abs(tree.leptons_phi - tree.met_phi)
        mt2 = 2.0 * tree.leptons_pt * tree.met_et * (1.0 - ROOT.TMath.Cos(dphi))
        mtw_value = ROOT.TMath.Sqrt(mt2) if mt2 > 0 else 0.0
        hists["mtw"].Fill(mtw_value, final_weight)

        # Fill jet histograms
        njets = len(tree.jet_pt)
        for j in range(njets):
            hists["ptjet"].Fill(tree.jet_pt[j], final_weight)
            hists["etajet"].Fill(abs(tree.jet_eta[j]), final_weight)
            hists["phijet"].Fill(tree.jet_phi[j], final_weight)

    return pos_eta_sum, neg_eta_sum

# ============================================================
# Main loop over samples
# ============================================================

for SignalSample in SignalSamples:
    print(f"\nRunning sample: {SignalSample}")

    # --------------------------------------------------------
    # Load weight map for this sample
    # --------------------------------------------------------
    try:
        idx_to_name = load_weight_map_clean(mapping_file_for_sample(SignalSample))
    except Exception as e:
        print(f"  Failed to load weight map: {e}")
        continue

    # --------------------------------------------------------
    # Open ROOT input file
    # --------------------------------------------------------
    input_path = os.path.join(DATA_DIR, f"WCharm_{SignalSample}.root")
    tfile = ROOT.TFile(input_path, "READ")

    if not tfile or tfile.IsZombie():
        print(f"  Could not open file: {input_path}")
        continue

    # --------------------------------------------------------
    # Get tree
    # --------------------------------------------------------
    mytree = tfile.Get("WCharmTree")
    if not mytree:
        print("  Missing WCharmTree")
        tfile.Close()
        continue

    # --------------------------------------------------------
    # Count entries
    # --------------------------------------------------------
    total_events = mytree.GetEntries()
    if total_events <= 0:
        print("  Tree has no entries")
        tfile.Close()
        continue

    # --------------------------------------------------------
    # Read first entry to determine weight vector length
    # --------------------------------------------------------
    mytree.GetEntry(0)
    weight_vector_length = len(mytree.weightvec)

    # --------------------------------------------------------
    # Build list of real weights to run
    # --------------------------------------------------------
    if USE_ALL_WEIGHTS:
        weights_to_run = [
            w for w in ALL_REQUESTED_WEIGHTS
            if w in idx_to_name and w < weight_vector_length
        ]
    else:
        weights_to_run = [0] if (0 in idx_to_name and 0 < weight_vector_length) else []

    print(f"  Weight vector length: {weight_vector_length}")
    print(f"  Number of valid weights to run: {len(weights_to_run)}")

    # --------------------------------------------------------
    # Loop over valid weights for this sample
    # --------------------------------------------------------
    for weight_index in weights_to_run:
        weight_name = idx_to_name[weight_index]
        print(f"    Running weight {weight_index}: {weight_name}")

        # Create fresh histograms
        hists = make_histograms()

        # Fill them
        pos_eta_sum, neg_eta_sum = fill_histograms(
            tree=mytree,
            hists=hists,
            weight_index=weight_index,
            total_events=total_events
        )

        # ----------------------------------------------------
        # Prepare output file
        # ----------------------------------------------------
        out_name = f"output_w{weight_index}_{safe_name(weight_name)}_{SignalSample}.root"
        out_path = os.path.join(OUTPUT_DIR, out_name)

        outputfile = ROOT.TFile(out_path, "RECREATE")

        # Metadata stored in the ROOT file
        ROOT.TNamed("RequestedWeightIndex", str(weight_index)).Write()
        ROOT.TNamed("EffectiveWeightIndex", str(weight_index)).Write()
        ROOT.TNamed("EffectiveWeightName", weight_name).Write()
        ROOT.TNamed("WeightVecLength", str(weight_vector_length)).Write()
        ROOT.TNamed("SampleName", SignalSample).Write()

        ROOT.TParameter(float)("PosEtaWeightSum", pos_eta_sum).Write()
        ROOT.TParameter(float)("NegEtaWeightSum", neg_eta_sum).Write()

        # Write histograms
        for hist in hists.values():
            hist.Write()

        outputfile.Write()
        outputfile.Close()

    # Close input file for this sample
    tfile.Close()

print("\nAll valid histogram files created.")
