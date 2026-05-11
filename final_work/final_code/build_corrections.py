# ============================================================
# build_corrections.py
#
# Purpose:
#   Read the histogram ROOT files produced by make_histograms.py
#   and calculate:
#     - nominal correction factors
#     - uncertainty histograms by source
#     - total uncertainty
#     - optional generator-difference uncertainty
#
# Output:
#   1. ROOT files in final_corrections/
#   2. CSV tables in final_corrections_csv/
#
# Correction factor definition:
#     C(bin) = parton(bin) / particle(bin)
#
# Important:
#   - Only weights present in BOTH particle and parton samples
#     are used for a given generator pair.
#   - Missing weights are not faked.
# ============================================================

import ROOT
import os
import re
import math
import csv

# ------------------------------------------------------------
# ROOT settings
# ------------------------------------------------------------
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

# ------------------------------------------------------------
# Locations
# ------------------------------------------------------------
INPUT_DIR = "/mnt/c/Users/dugar/wcharm_analysis/final_work/make_histograms.py_Outputs"
OUTPUT_DIR = "/mnt/c/Users/dugar/wcharm_analysis/final_work/build_corrections.py_Outputs"
CSV_DIR = "/mnt/c/Users/dugar/wcharm_analysis/final_work/build_corrections.py_CSV"

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(CSV_DIR, exist_ok=True)

# ------------------------------------------------------------
# Observables expected in the histogram files
# ------------------------------------------------------------
OBSERVABLES = [
    "ptlepton",
    "ptjet",
    "etalepton",
    "etajet",
    "philepton",
    "phijet",
    "met_et",
    "mtw"
]

# ------------------------------------------------------------
# Particle / parton sample pairings
# ------------------------------------------------------------
SAMPLE_PAIRS = {
    "Pythia_plus":  ("WCPy8plus",      "WCPyPartonplus"),
    "Pythia_minus": ("WCPy8minus",     "WCPyPartonminus"),
    "Herwig_plus":  ("WCH7plus",       "WCHPartonplus"),
    "Herwig_minus": ("WCH7minus",      "WCHPartonminus"),
}

# ------------------------------------------------------------
# Weight groups
# ------------------------------------------------------------
SCALE_WEIGHTS = [1, 112, 215, 226, 237, 248, 259, 270]
MODEL_WEIGHTS = [292, 293, 314, 315]
SHOWER_WEIGHTS = list(range(294, 314)) + [316, 317]
EXCLUDED_WEIGHTS = [318]

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

def parse_output_filename(filename):
    pattern = (
        r"output_w(\d+)_(.+)_"
        r"(WCH7minus|WCH7plus|WCHPartonminus|WCHPartonplus|"
        r"WCPy8minus|WCPy8plus|WCPyPartonminus|WCPyPartonplus)\.root$"
    )

    match = re.match(pattern, filename)
    if not match:
        return None

    return {
        "weight_index": int(match.group(1)),
        "weight_name": match.group(2),
        "sample": match.group(3),
        "path": os.path.join(INPUT_DIR, filename),
    }


def scan_input_files():
    files_by_sample = {}

    for filename in os.listdir(INPUT_DIR):
        if not filename.endswith(".root"):
            continue

        info = parse_output_filename(filename)
        if info is None:
            continue

        sample = info["sample"]
        weight_index = info["weight_index"]

        files_by_sample.setdefault(sample, {})
        files_by_sample[sample][weight_index] = info

    return files_by_sample


def open_root_files_for_sample(sample_files):
    """
    Open all ROOT files for one sample once and keep them open.
    Returns:
      dict[weight_index] = TFile
    """
    opened = {}
    for w, info in sample_files.items():
        tf = ROOT.TFile.Open(info["path"], "READ")
        if not tf or tf.IsZombie():
            raise RuntimeError(f"Could not open ROOT file: {info['path']}")
        opened[w] = tf
    return opened


def close_root_files(opened_files):
    """
    Close all ROOT files in a dictionary.
    """
    for tf in opened_files.values():
        tf.Close()


def get_histogram_from_open_file(tf, hist_name):
    """
    Read a histogram from an already-open ROOT file and clone it into memory.
    """
    hist = tf.Get(hist_name)
    if not hist:
        raise RuntimeError(f"Histogram '{hist_name}' not found in file {tf.GetName()}")

    cloned = hist.Clone(f"{hist_name}_clone")
    cloned.SetDirectory(0)
    return cloned


def make_ratio_histogram(h_particle, h_parton, out_name):
    """
    Build the correction-factor histogram:
        parton / particle
    """
    ratio = h_parton.Clone(out_name)
    ratio.SetDirectory(0)
    ratio.Divide(h_parton, h_particle, 1.0, 1.0, "")
    return ratio


def empty_hist_like(reference_hist, name):
    out = reference_hist.Clone(name)
    out.Reset()
    out.SetDirectory(0)
    return out


def is_pdf_weight_name(weight_name):
    return weight_name.startswith("MUR1.0_MUF1.0_PDF")


def compute_envelope_uncertainty(nominal_hist, varied_hists, out_name):
    out = empty_hist_like(nominal_hist, out_name)

    for ibin in range(1, nominal_hist.GetNbinsX() + 1):
        central = nominal_hist.GetBinContent(ibin)
        max_dev = 0.0

        for hist in varied_hists:
            dev = abs(hist.GetBinContent(ibin) - central)
            if dev > max_dev:
                max_dev = dev

        out.SetBinContent(ibin, max_dev)
        out.SetBinError(ibin, 0.0)

    return out


def compute_rms_uncertainty(nominal_hist, varied_hists, out_name):
    out = empty_hist_like(nominal_hist, out_name)

    if len(varied_hists) == 0:
        return out

    for ibin in range(1, nominal_hist.GetNbinsX() + 1):
        central = nominal_hist.GetBinContent(ibin)
        sumsq = 0.0

        for hist in varied_hists:
            delta = hist.GetBinContent(ibin) - central
            sumsq += delta * delta

        rms = math.sqrt(sumsq / len(varied_hists))
        out.SetBinContent(ibin, rms)
        out.SetBinError(ibin, 0.0)

    return out


def compute_stat_uncertainty_from_ratio(ratio_hist, out_name):
    out = empty_hist_like(ratio_hist, out_name)

    for ibin in range(1, ratio_hist.GetNbinsX() + 1):
        out.SetBinContent(ibin, ratio_hist.GetBinError(ibin))
        out.SetBinError(ibin, 0.0)

    return out


def compute_total_uncertainty(nominal_hist, uncertainty_hists, out_name):
    out = empty_hist_like(nominal_hist, out_name)

    for ibin in range(1, nominal_hist.GetNbinsX() + 1):
        sumsq = 0.0
        for hist in uncertainty_hists:
            sigma = hist.GetBinContent(ibin)
            sumsq += sigma * sigma

        out.SetBinContent(ibin, math.sqrt(sumsq))
        out.SetBinError(ibin, 0.0)

    return out


def compute_generator_difference(hist_a, hist_b, out_name):
    out = empty_hist_like(hist_a, out_name)

    for ibin in range(1, hist_a.GetNbinsX() + 1):
        diff = abs(hist_a.GetBinContent(ibin) - hist_b.GetBinContent(ibin))
        out.SetBinContent(ibin, diff)
        out.SetBinError(ibin, 0.0)

    return out


def build_variation_histograms(weight_list, particle_open_files, parton_open_files, observable, label):
    varied = []

    for w in weight_list:
        h_particle = get_histogram_from_open_file(particle_open_files[w], observable)
        h_parton = get_histogram_from_open_file(parton_open_files[w], observable)

        h_corr = make_ratio_histogram(
            h_particle,
            h_parton,
            f"{observable}_corr_{label}_w{w}"
        )
        varied.append(h_corr)

    return varied


def write_hist_to_file(fout, hist):
    """
    Ensure the correct output file is the active ROOT directory before writing.
    """
    fout.cd()
    hist.Write()


def write_csv_summary(csv_path, pair_label, observable, h_corr, h_stat, h_scale, h_pdf, h_shower, h_model, h_total):
    with open(csv_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)

        writer.writerow([
            "pair_label",
            "observable",
            "bin",
            "bin_low_edge",
            "bin_up_edge",
            "correction_factor",
            "stat_unc",
            "scale_unc",
            "pdf_unc",
            "shower_unc",
            "model_unc",
            "total_unc"
        ])

        for ibin in range(1, h_corr.GetNbinsX() + 1):
            low_edge = h_corr.GetXaxis().GetBinLowEdge(ibin)
            up_edge = h_corr.GetXaxis().GetBinUpEdge(ibin)

            writer.writerow([
                pair_label,
                observable,
                ibin,
                low_edge,
                up_edge,
                h_corr.GetBinContent(ibin),
                h_stat.GetBinContent(ibin),
                h_scale.GetBinContent(ibin),
                h_pdf.GetBinContent(ibin),
                h_shower.GetBinContent(ibin),
                h_model.GetBinContent(ibin),
                h_total.GetBinContent(ibin)
            ])


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

files_by_sample = scan_input_files()
nominal_ratios = {}

for pair_label, (particle_sample, parton_sample) in SAMPLE_PAIRS.items():
    print(f"\nProcessing pair: {pair_label}")
    print(f"  Particle sample: {particle_sample}")
    print(f"  Parton sample:   {parton_sample}")

    if particle_sample not in files_by_sample:
        print("  Missing particle-level files")
        continue

    if parton_sample not in files_by_sample:
        print("  Missing parton-level files")
        continue

    particle_files = files_by_sample[particle_sample]
    parton_files = files_by_sample[parton_sample]

    common_weights = sorted(set(particle_files.keys()) & set(parton_files.keys()))
    common_weights = [w for w in common_weights if w not in EXCLUDED_WEIGHTS]

    print(f"  Common usable weights: {len(common_weights)}")

    if 0 not in common_weights:
        print("  No common nominal weight 0 found - skipping")
        continue

    pdf_weights = []
    for w in common_weights:
        weight_name = particle_files[w]["weight_name"]
        if (
            is_pdf_weight_name(weight_name)
            and w not in SCALE_WEIGHTS
            and w not in MODEL_WEIGHTS
            and w not in SHOWER_WEIGHTS
            and w != 281
        ):
            pdf_weights.append(w)

    scale_weights_common = [w for w in SCALE_WEIGHTS if w in common_weights]
    model_weights_common = [w for w in MODEL_WEIGHTS if w in common_weights]
    shower_weights_common = [w for w in SHOWER_WEIGHTS if w in common_weights]

    print(f"  Scale weights found:  {scale_weights_common}")
    print(f"  PDF weights found:    {len(pdf_weights)}")
    print(f"  Shower weights found: {shower_weights_common}")
    print(f"  Model weights found:  {model_weights_common}")

    particle_open_files = open_root_files_for_sample(particle_files)
    parton_open_files = open_root_files_for_sample(parton_files)

    out_root_path = os.path.join(OUTPUT_DIR, f"corrections_{pair_label}.root")
    fout = ROOT.TFile(out_root_path, "RECREATE")

    fout.cd()
    ROOT.TNamed("PairLabel", pair_label).Write()
    ROOT.TNamed("ParticleSample", particle_sample).Write()
    ROOT.TNamed("PartonSample", parton_sample).Write()

    for observable in OBSERVABLES:
        print(f"    Observable: {observable}")

        h_particle_nom = get_histogram_from_open_file(particle_open_files[0], observable)
        h_parton_nom = get_histogram_from_open_file(parton_open_files[0], observable)

        h_corr_nom = make_ratio_histogram(
            h_particle_nom,
            h_parton_nom,
            f"{observable}_corr_nominal"
        )
        write_hist_to_file(fout, h_corr_nom)

        nominal_ratios[(pair_label, observable)] = h_corr_nom.Clone(f"{pair_label}_{observable}_nominal")
        nominal_ratios[(pair_label, observable)].SetDirectory(0)

        h_unc_stat = compute_stat_uncertainty_from_ratio(
            h_corr_nom,
            f"{observable}_unc_stat"
        )
        write_hist_to_file(fout, h_unc_stat)

        varied_scale = build_variation_histograms(
            scale_weights_common,
            particle_open_files,
            parton_open_files,
            observable,
            "scale"
        )

        varied_pdf = build_variation_histograms(
            pdf_weights,
            particle_open_files,
            parton_open_files,
            observable,
            "pdf"
        )

        varied_shower = build_variation_histograms(
            shower_weights_common,
            particle_open_files,
            parton_open_files,
            observable,
            "shower"
        )

        varied_model = build_variation_histograms(
            model_weights_common,
            particle_open_files,
            parton_open_files,
            observable,
            "model"
        )

        h_unc_scale = compute_envelope_uncertainty(
            h_corr_nom,
            varied_scale,
            f"{observable}_unc_scale"
        )
        write_hist_to_file(fout, h_unc_scale)

        h_unc_pdf = compute_rms_uncertainty(
            h_corr_nom,
            varied_pdf,
            f"{observable}_unc_pdf"
        )
        write_hist_to_file(fout, h_unc_pdf)

        h_unc_shower = compute_envelope_uncertainty(
            h_corr_nom,
            varied_shower,
            f"{observable}_unc_shower"
        )
        write_hist_to_file(fout, h_unc_shower)

        h_unc_model = compute_envelope_uncertainty(
            h_corr_nom,
            varied_model,
            f"{observable}_unc_model"
        )
        write_hist_to_file(fout, h_unc_model)

        h_unc_total = compute_total_uncertainty(
            h_corr_nom,
            [h_unc_stat, h_unc_scale, h_unc_pdf, h_unc_shower, h_unc_model],
            f"{observable}_unc_total"
        )
        write_hist_to_file(fout, h_unc_total)

        csv_path = os.path.join(CSV_DIR, f"{pair_label}_{observable}.csv")
        write_csv_summary(
            csv_path=csv_path,
            pair_label=pair_label,
            observable=observable,
            h_corr=h_corr_nom,
            h_stat=h_unc_stat,
            h_scale=h_unc_scale,
            h_pdf=h_unc_pdf,
            h_shower=h_unc_shower,
            h_model=h_unc_model,
            h_total=h_unc_total
        )

    fout.Write()
    fout.Close()

    close_root_files(particle_open_files)
    close_root_files(parton_open_files)

# Optional generator comparison
generator_pairs = [
    ("plus", "Pythia_plus", "Herwig_plus"),
    ("minus", "Pythia_minus", "Herwig_minus"),
]

for charge_label, pythia_key, herwig_key in generator_pairs:
    out_root_path = os.path.join(OUTPUT_DIR, f"generator_difference_{charge_label}.root")
    fout = ROOT.TFile(out_root_path, "RECREATE")

    for observable in OBSERVABLES:
        key_py = (pythia_key, observable)
        key_hw = (herwig_key, observable)

        if key_py not in nominal_ratios or key_hw not in nominal_ratios:
            continue

        h_unc_generator = compute_generator_difference(
            nominal_ratios[key_py],
            nominal_ratios[key_hw],
            f"{observable}_unc_generator"
        )
        write_hist_to_file(fout, h_unc_generator)

    fout.Write()
    fout.Close()

print("\nFinished building correction factors, uncertainties, and CSV summaries.")
