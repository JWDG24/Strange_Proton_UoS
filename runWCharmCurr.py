# ============================================================
# WCharm analysis script (ROOT/PyROOT)
#
# ROOT-safe version (no wildcard import)
#
# Behaviour:
#   • One output per (sample, requested weight index)
#   • If requested weight does NOT exist:
#         → use weight index 0
#         → label weight as "Default"
#   • No skipping
#   • No cross-generator mixing
# ============================================================

import ROOT
import sys
import os
from array import array

## Boolean Switches ##

# ------------------------------------------------------------
# --- Weight configuration switch ---
# If True  → run all weights in used_weights
# If False → run only default weight (index 0)
# ------------------------------------------------------------
USE_ALL_WEIGHTS = False

# ------------------------------------------------------------
# --- Eta binning switch ---
# True  → use variable (Kristin) bins
# False → use regular equidistant bins
# ------------------------------------------------------------
USE_VARIABLE_ETA_BINS = True


## Bins ##

# ------------------------------------------------------------
# Bin Lists
# (Phi works with uniform binning)
# ------------------------------------------------------------

# Lepton Pseudorapidity Bins (h_etalep)
eta_lepton_bins = sorted([0, 0.21, 0.42, 0.63, 0.84, 1.05, 1.37, 1.52, 1.74, 1.95, 2.18, 2.5])
eta_lepton_binArray = array('d', eta_lepton_bins)

# Missing Transverse Energy Bins (h_met)
met_bins = [
    0,
    10, 20, 25, 30, 35, 40, 45,
    50, 52, 55, 58, 60, 62, 65, 68, 70, 72, 75, 78, 80,
    85, 90, 100,
    120, 140, 160,
    200, 250, 300
]
met_binArray = array('d', met_bins)

# W Boson Transverse Mass Energy Bins (h_mtw)
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

# Jet Transverse Momentum Bins (h_ptjet)
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

# Letpon Transverse Momentum Bins (h_ptlepton)
ptlepton_bins = [
    0,
    6, 12, 18, 24, 30, 36, 42, 48, 54,
    60, 66, 72, 78, 84, 90, 96, 102, 108, 114,
    120,
    150, 180, 210, 240, 270, 300
]
ptlepton_binArray = array('d', ptlepton_bins)

# Jet and Lepton  Azimuthal Angles (h_phijet/philep)
    # These are fine for uniform binning as they are just angles so please see the binning code for the bins.

## Weighting Information ##
WEIGHT_DIR = "/mnt/c/Users/dugar/wcharm_analysis/weights"

# Requested weight indices
used_weights = [0, 292, 293, 314, 315]

## CODE ##

# ------------------------------------------------------------
# Load clean mapping file
# ------------------------------------------------------------
def load_weight_map_clean(path):

    idx_to_name = {}

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            idx_str, name = line.split(None, 1)
            idx_to_name[int(idx_str)] = name.strip()

    return idx_to_name


# ------------------------------------------------------------
# Choose mapping file based on sample name
# ------------------------------------------------------------
def mapping_file_for_sample(sample):

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
# MAIN LOOP OVER SAMPLES
# ============================================================

for SignalSample in SignalSamples:

    print("\nRunning on sample:", SignalSample)

    try:
        idx_to_name = load_weight_map_clean(mapping_file_for_sample(SignalSample))
    except Exception as e:
        print("Mapping load failed:", e)
        continue

    tfile = ROOT.TFile(
        f"/mnt/c/Users/dugar/wcharm_analysis/data/WCharm_{SignalSample}.root",
        "read"
    )

    if not tfile or tfile.IsZombie():
        print("File not found:", SignalSample)
        continue

    mytree = tfile.Get("WCharmTree")
    allevents = mytree.GetEntries()

    # ------------------------------------------------------------
    # --- Choose which weights to run based on boolean ---
    # ------------------------------------------------------------
    if USE_ALL_WEIGHTS:
        weights_to_run = used_weights
    else:
        weights_to_run = [0]

    # ============================================================
    # LOOP OVER REQUESTED WEIGHT INDICES
    # ============================================================

    for requested_widx in weights_to_run:

        mytree.GetEntry(0)
        wvec_len = len(mytree.weightvec)

        effective_widx = 0
        wname = "Default"

        if requested_widx in idx_to_name and requested_widx < wvec_len:
            effective_widx = requested_widx
            wname = idx_to_name[requested_widx]

        if wvec_len == 0:
            print(f"ERROR: weightvec empty for {SignalSample}")
            continue

        print(f"  Requested w{requested_widx} -> using index {effective_widx} ({wname})")

        outputfile = ROOT.TFile(
            f"/mnt/c/Users/dugar/wcharm_analysis/outputs_weights/"
            f"output_w{requested_widx}_{wname}_{SignalSample}.root",
            "recreate"
        )

        ROOT.TNamed("RequestedWeightIndex", str(requested_widx)).Write()
        ROOT.TNamed("EffectiveWeightIndex", str(effective_widx)).Write()
        ROOT.TNamed("EffectiveWeightName", wname).Write()
        ROOT.TNamed("WeightVecLength", str(wvec_len)).Write()

        ## --------------------------------------------------------
        ## Histograms
        ## --------------------------------------------------------

        h_ptlep = ROOT.TH1F("ptlepton","ptlepton",len(ptlepton_binArray)-1, ptlepton_binArray)
        h_ptjet = ROOT.TH1F("ptjet","ptjet",len(ptjet_binArray)-1, ptjet_binArray)

        if USE_VARIABLE_ETA_BINS:
            h_etalep = ROOT.TH1F("etalepton", "etalepton", len(eta_lepton_binArray)-1, eta_lepton_binArray)
            h_etajet = ROOT.TH1F("etajet", "etajet", len(eta_lepton_binArray)-1, eta_lepton_binArray)
        else:
            h_etalep = ROOT.TH1F("etalepton", "etalepton", 12, 0, 2.5)
            h_etajet = ROOT.TH1F("etajet", "etajet", 12, 0, 2.5)

        h_philep = ROOT.TH1F("philepton","philepton",50,-1,7)
        h_phijet = ROOT.TH1F("phijet","phijet",50,-1,7)

        h_met = ROOT.TH1F("met_et","met_et",len(met_binArray)-1, met_binArray)
        h_mtw = ROOT.TH1F("mtw","mtw",len(mtw_binArray)-1, mtw_binArray)

        for h in [h_ptlep, h_ptjet, h_etalep, h_etajet,
                  h_philep, h_phijet, h_met, h_mtw]:
            h.Sumw2()

        ## Double Value Checker ##
        pos_eta_sum = 0.0
        neg_eta_sum = 0.0

        # ============================================================
        # EVENT LOOP
        # ============================================================

        for event in mytree:

            weightxs = (1.0 / allevents) * event.weightvec[effective_widx]

            if event.met_et <= 25:
                continue
            if event.leptons_pt <= 20:
                continue
            if abs(event.leptons_eta) >= 2.5:
                continue

            corr = -1 * event.leptons_charge
            final_weight = corr * weightxs

            ## Consistent Eta Handling ##
            eta_lep = abs(event.leptons_eta)

            ## Double Value Checker ##
            if event.leptons_eta >= 0:
                pos_eta_sum += final_weight
            else:
                neg_eta_sum += final_weight

            h_ptlep.Fill(event.leptons_pt, final_weight)
            
            h_philep.Fill(event.leptons_phi, final_weight)

            ## Consistent Eta Handling ##
            h_etalep.Fill(eta_lep, final_weight)

            h_met.Fill(event.met_et, final_weight)

            dphi = abs(event.leptons_phi - event.met_phi)
            mt2 = 2.0 * event.leptons_pt * event.met_et * (1.0 - ROOT.TMath.Cos(dphi))
            mtw = ROOT.TMath.Sqrt(mt2) if mt2 > 0 else 0
            h_mtw.Fill(mtw, final_weight)

            for jet in range(len(event.jet_pt)):

                h_ptjet.Fill(event.jet_pt[jet], final_weight)

                ## Consistent Eta Handling ##
                eta_jet = abs(event.jet_eta[jet])
                h_etajet.Fill(eta_jet, final_weight)

                h_phijet.Fill(event.jet_phi[jet], final_weight)

        ## Double Value Checker ##
        print(f"  +eta weight sum: {pos_eta_sum}")
        print(f"  -eta weight sum: {neg_eta_sum}")
        if neg_eta_sum != 0:
            print(f"  Ratio (+/-): {pos_eta_sum / neg_eta_sum}")

        outputfile.Write()
        outputfile.Close()

    tfile.Close()

print("\nAll samples processed successfully.")
