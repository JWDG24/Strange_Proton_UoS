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

WEIGHT_DIR = "/mnt/c/Users/dugar/wcharm_analysis/weights"

# ------------------------------------------------------------
# --- Weight configuration switch ---
# If True  → run all weights in used_weights
# If False → run only default weight (index 0)
# ------------------------------------------------------------
USE_ALL_WEIGHTS = False

# Requested weight indices
used_weights = [0, 292, 293, 314, 315]


# ------------------------------------------------------------
# --- Kristin eta binning implementation ---
# Variable bin edges for |eta|
# Histogram will be filled with fabs(eta)
# ------------------------------------------------------------
bins = sorted([0, 0.21, 0.42, 0.63, 0.84, 1.05, 1.37, 1.52, 1.74, 1.95, 2.18, 2.5])
binArray = array('d', bins)


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

        # --------------------------------------------------------
        # Histograms
        # --------------------------------------------------------

        h_ptlep = ROOT.TH1F("ptlepton","ptlepton",30,0,300)
        h_ptjet = ROOT.TH1F("ptjet","ptjet",30,0,300)

        # --- Custom eta binning ---
        h_etalep = ROOT.TH1F("etalepton","etalepton",len(binArray)-1, binArray)
        h_etajet = ROOT.TH1F("etajet","etajet",len(binArray)-1, binArray)

        h_philep = ROOT.TH1F("philepton","philepton",100,-5,10)
        h_phijet = ROOT.TH1F("phijet","phijet",100,-5,10)
        h_met = ROOT.TH1F("met_et","met_et",30,0,300)
        h_mtw = ROOT.TH1F("mtw","mtw",30,0,300)

        for h in [h_ptlep, h_ptjet, h_etalep, h_etajet,
                  h_philep, h_phijet, h_met, h_mtw]:
            h.Sumw2()

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

            h_ptlep.Fill(event.leptons_pt, final_weight)

            # --- Fill using fabs(eta) ---
            h_etalep.Fill(ROOT.TMath.Abs(event.leptons_eta), final_weight)

            h_met.Fill(event.met_et, final_weight)

            dphi = abs(event.leptons_phi - event.met_phi)
            mt2 = 2.0 * event.leptons_pt * event.met_et * (1.0 - ROOT.TMath.Cos(dphi))
            mtw = ROOT.TMath.Sqrt(mt2) if mt2 > 0 else 0
            h_mtw.Fill(mtw, final_weight)

            for jet in range(len(event.jet_pt)):

                h_ptjet.Fill(event.jet_pt[jet], final_weight)

                # --- Fill jet eta with fabs ---
                h_etajet.Fill(ROOT.TMath.Abs(event.jet_eta[jet]), final_weight)

                h_phijet.Fill(event.jet_phi[jet], final_weight)

        outputfile.Write()
        outputfile.Close()

    tfile.Close()

print("\nAll samples processed successfully.")
