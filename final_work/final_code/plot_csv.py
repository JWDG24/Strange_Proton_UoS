#!/usr/bin/env python3

# ============================================================
# plot_csv_E2E1.py
#
# Purpose:
#   Read all correction-factor CSV files from INPUT_DIR,
#   convert each one into a ROOT TH1D,
#   draw:
#     - top pad: correction factor with E2 band + E1 points
#     - bottom pad: all uncertainty curves
#   and save the plots to OUTPUT_DIR.
#
# Run with:
#   python3 plot_csv_E2E1.py
#
# Requirements:
#   - PyROOT installed
#   - CSV files contain columns:
#       bin_low_edge
#       bin_up_edge
#       correction_factor
#       stat_unc
#       scale_unc
#       pdf_unc
#       shower_unc
#       model_unc
#       total_unc
#
# Output:
#   For each CSV file:
#     - one PDF plot
#     - one PNG plot
#
# Notes:
#   - Processes all CSV files in INPUT_DIR
#   - Keeps output suffix: _E2E1
# ============================================================

import os
import csv
from array import array

import ROOT

# ------------------------------------------------------------
# ROOT settings
# ------------------------------------------------------------
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetErrorX(0.5)

# ------------------------------------------------------------
# Input / output directories
# ------------------------------------------------------------
INPUT_DIR = "/mnt/c/Users/dugar/wcharm_analysis/final_work/build_corrections.py_CSV"
OUTPUT_DIR = "/mnt/c/Users/dugar/wcharm_analysis/final_work/plot_csv_E2E1_Output"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ------------------------------------------------------------
# Plot settings
# ------------------------------------------------------------
SAVE_PDF = True
SAVE_PNG = True

CORR_COLUMN = "correction_factor"
TOTAL_COLUMN = "total_unc"
STAT_COLUMN = "stat_unc"
SCALE_COLUMN = "scale_unc"
PDF_COLUMN = "pdf_unc"
SHOWER_COLUMN = "shower_unc"
MODEL_COLUMN = "model_unc"

# Canvas / pad geometry
CANVAS_W = 1200
CANVAS_H = 850

LEFT_MARGIN = 0.12
RIGHT_MARGIN = 0.30   # room for outside legend
TOP_MARGIN = 0.08
BOTTOM_MARGIN = 0.12

TOP_PAD_YMIN = 0.34
TOP_PAD_YMAX = 1.00
BOT_PAD_YMIN = 0.08
BOT_PAD_YMAX = 0.34

# Font sizes
TOP_X_LABEL_SIZE = 0.0
TOP_X_TITLE_SIZE = 0.0

TOP_Y_TITLE_SIZE = 0.050
TOP_Y_LABEL_SIZE = 0.042
TOP_Y_TITLE_OFFSET = 1.15

BOT_X_TITLE_SIZE = 0.11
BOT_X_LABEL_SIZE = 0.095
BOT_X_TITLE_OFFSET = 1.05

BOT_Y_TITLE_SIZE = 0.10
BOT_Y_LABEL_SIZE = 0.085
BOT_Y_TITLE_OFFSET = 0.55

LEGEND_TEXT_SIZE = 0.030

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

def observable_to_x_title(observable_name: str) -> str:
    mapping = {
        "ptlepton": "Lepton p_{T} [GeV]",
        "ptjet": "Jet p_{T} [GeV]",
        "etalepton": "Lepton |#eta|",
        "etajet": "Jet |#eta|",
        "philepton": "Lepton #phi [rad]",
        "phijet": "Jet #phi [rad]",
        "met_et": "Missing E_{T} [GeV]",
        "mtw": "W transverse mass m_{T}^{W} [GeV]",
    }
    return mapping.get(observable_name, observable_name)


def observable_to_title(observable_name: str) -> str:
    mapping = {
        "ptlepton": "Lepton transverse momentum",
        "ptjet": "Jet transverse momentum",
        "etalepton": "Lepton pseudorapidity",
        "etajet": "Jet pseudorapidity",
        "philepton": "Lepton azimuthal angle",
        "phijet": "Jet azimuthal angle",
        "met_et": "Missing transverse energy",
        "mtw": "W transverse mass",
    }
    return mapping.get(observable_name, observable_name)


def parse_csv_filename(filename: str):
    stem = os.path.splitext(os.path.basename(filename))[0]
    parts = stem.rsplit("_", 1)
    if len(parts) != 2:
        return stem, stem
    return parts[0], parts[1]


def make_pretty_title(pair_label: str, observable: str) -> str:
    return f"{pair_label.replace('_', ' ')}: {observable_to_title(observable)}"


def load_csv_data(csv_path: str):
    low_edges = []
    up_edges = []

    corr_vals = []
    total_unc = []
    stat_unc = []
    scale_unc = []
    pdf_unc = []
    shower_unc = []
    model_unc = []

    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)

        required_columns = [
            "bin_low_edge",
            "bin_up_edge",
            CORR_COLUMN,
            TOTAL_COLUMN,
            STAT_COLUMN,
            SCALE_COLUMN,
            PDF_COLUMN,
            SHOWER_COLUMN,
            MODEL_COLUMN,
        ]
        for col in required_columns:
            if col not in reader.fieldnames:
                raise RuntimeError(f"Missing required column '{col}' in {csv_path}")

        for row in reader:
            low_edges.append(float(row["bin_low_edge"]))
            up_edges.append(float(row["bin_up_edge"]))
            corr_vals.append(float(row[CORR_COLUMN]))
            total_unc.append(float(row[TOTAL_COLUMN]))
            stat_unc.append(float(row[STAT_COLUMN]))
            scale_unc.append(float(row[SCALE_COLUMN]))
            pdf_unc.append(float(row[PDF_COLUMN]))
            shower_unc.append(float(row[SHOWER_COLUMN]))
            model_unc.append(float(row[MODEL_COLUMN]))

    if not low_edges:
        raise RuntimeError(f"No data rows found in {csv_path}")

    return {
        "low_edges": low_edges,
        "up_edges": up_edges,
        "corr_vals": corr_vals,
        "total_unc": total_unc,
        "stat_unc": stat_unc,
        "scale_unc": scale_unc,
        "pdf_unc": pdf_unc,
        "shower_unc": shower_unc,
        "model_unc": model_unc,
    }


def make_hist_from_csv(csv_path: str, data_dict):
    pair_label, observable = parse_csv_filename(csv_path)
    title = make_pretty_title(pair_label, observable)

    low_edges = data_dict["low_edges"]
    up_edges = data_dict["up_edges"]
    corr_vals = data_dict["corr_vals"]
    total_unc = data_dict["total_unc"]

    edges = [low_edges[0]] + up_edges
    edge_array = array("d", edges)

    hist_name = f"h_{pair_label}_{observable}"
    hist = ROOT.TH1D(hist_name, title, len(corr_vals), edge_array)
    hist.SetDirectory(0)

    for i, (y, e) in enumerate(zip(corr_vals, total_unc), start=1):
        hist.SetBinContent(i, y)
        hist.SetBinError(i, e)

    hist.SetTitle(title)
    return hist, pair_label, observable


def make_band_hist(source_hist):
    band = source_hist.Clone(source_hist.GetName() + "_band")
    band.SetDirectory(0)
    band.SetLineWidth(0)
    band.SetMarkerStyle(0)
    band.SetMarkerSize(0)
    band.SetFillStyle(1001)
    band.SetFillColorAlpha(ROOT.kBlue, 0.25)
    return band


def make_point_hist(source_hist):
    points = source_hist.Clone(source_hist.GetName() + "_points")
    points.SetDirectory(0)
    points.SetLineColor(ROOT.kBlue + 2)
    points.SetLineWidth(2)
    points.SetMarkerColor(ROOT.kBlue + 2)
    points.SetMarkerStyle(20)
    points.SetMarkerSize(0.9)
    return points


def make_uncertainty_graph(low_edges, up_edges, values, name, color, marker_style):
    n = len(values)
    x_vals = []
    y_vals = []
    for lo, hi, val in zip(low_edges, up_edges, values):
        x_vals.append(0.5 * (lo + hi))
        y_vals.append(val)

    graph = ROOT.TGraph(n, array("d", x_vals), array("d", y_vals))
    graph.SetName(name)
    graph.SetLineColor(color)
    graph.SetMarkerColor(color)
    graph.SetLineWidth(2)
    graph.SetMarkerStyle(marker_style)
    graph.SetMarkerSize(0.85)
    return graph


def max_of_list(values):
    return max(values) if values else 0.0


def style_top_axes(hist, observable):
    hist.GetXaxis().SetTitle(observable_to_x_title(observable))
    hist.GetYaxis().SetTitle("Correction factor")

    hist.GetXaxis().SetTitleSize(TOP_X_TITLE_SIZE)
    hist.GetXaxis().SetLabelSize(TOP_X_LABEL_SIZE)

    hist.GetYaxis().SetTitleSize(TOP_Y_TITLE_SIZE)
    hist.GetYaxis().SetLabelSize(TOP_Y_LABEL_SIZE)
    hist.GetYaxis().SetTitleOffset(TOP_Y_TITLE_OFFSET)


def style_bottom_axes(frame, observable):
    frame.GetXaxis().SetTitle(observable_to_x_title(observable))
    frame.GetYaxis().SetTitle("Uncertainty")

    frame.GetXaxis().SetTitleSize(BOT_X_TITLE_SIZE)
    frame.GetXaxis().SetLabelSize(BOT_X_LABEL_SIZE)
    frame.GetXaxis().SetTitleOffset(BOT_X_TITLE_OFFSET)

    frame.GetYaxis().SetTitleSize(BOT_Y_TITLE_SIZE)
    frame.GetYaxis().SetLabelSize(BOT_Y_LABEL_SIZE)
    frame.GetYaxis().SetTitleOffset(BOT_Y_TITLE_OFFSET)


def save_hist_plots(hist, csv_path: str, output_dir: str, data_dict):
    base = os.path.splitext(os.path.basename(csv_path))[0]
    _, observable = parse_csv_filename(csv_path)

    low_edges = data_dict["low_edges"]
    up_edges = data_dict["up_edges"]

    total_unc = data_dict["total_unc"]
    stat_unc = data_dict["stat_unc"]
    scale_unc = data_dict["scale_unc"]
    pdf_unc = data_dict["pdf_unc"]
    shower_unc = data_dict["shower_unc"]
    model_unc = data_dict["model_unc"]

    # Top pad y range for correction factor
    corr_high = max(hist.GetBinContent(i) + hist.GetBinError(i) for i in range(1, hist.GetNbinsX() + 1))
    corr_low = min(hist.GetBinContent(i) - hist.GetBinError(i) for i in range(1, hist.GetNbinsX() + 1))
    span = corr_high - corr_low
    if span <= 0:
        span = max(abs(corr_high), 1.0)

    top_y_min = corr_low - 0.10 * span
    top_y_max = corr_high + 0.10 * span

    # Bottom pad y range for uncertainties
    unc_max = max(
        max_of_list(total_unc),
        max_of_list(stat_unc),
        max_of_list(scale_unc),
        max_of_list(pdf_unc),
        max_of_list(shower_unc),
        max_of_list(model_unc),
    )
    bot_y_min = 0.0
    bot_y_max = 1.15 * unc_max if unc_max > 0 else 1.0

    canvas = ROOT.TCanvas(f"c_{base}", f"c_{base}", CANVAS_W, CANVAS_H)

    # -----------------------------
    # Top pad: correction factor
    # -----------------------------
    pad_top = ROOT.TPad(f"pad_top_{base}", "", 0.0, TOP_PAD_YMIN, 1.0, TOP_PAD_YMAX)
    pad_top.SetLeftMargin(LEFT_MARGIN)
    pad_top.SetRightMargin(RIGHT_MARGIN)
    pad_top.SetTopMargin(TOP_MARGIN)
    pad_top.SetBottomMargin(0.02)
    pad_top.Draw()
    pad_top.cd()

    band_hist = make_band_hist(hist)
    point_hist = make_point_hist(hist)

    style_top_axes(band_hist, observable)
    band_hist.SetMinimum(top_y_min)
    band_hist.SetMaximum(top_y_max)

    band_hist.Draw("E2")
    point_hist.Draw("E1 SAME")

    # Legend outside both plots, in right margin
    legend = ROOT.TLegend(0.73, 0.18, 0.98, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(LEGEND_TEXT_SIZE)

    # We'll add entries now and draw later on the top pad
    legend.AddEntry(point_hist, "Correction factor", "lep")

    canvas.cd()

    # -----------------------------
    # Bottom pad: uncertainties
    # -----------------------------
    pad_bot = ROOT.TPad(f"pad_bot_{base}", "", 0.0, BOT_PAD_YMIN, 1.0, BOT_PAD_YMAX)
    pad_bot.SetLeftMargin(LEFT_MARGIN)
    pad_bot.SetRightMargin(RIGHT_MARGIN)
    pad_bot.SetTopMargin(0.02)
    pad_bot.SetBottomMargin(BOTTOM_MARGIN + 0.12)
    pad_bot.Draw()
    pad_bot.cd()

    x_min = hist.GetXaxis().GetXmin()
    x_max = hist.GetXaxis().GetXmax()

    frame = ROOT.TH1F(f"frame_{base}", "", 1, x_min, x_max)
    frame.SetDirectory(0)
    frame.SetMinimum(bot_y_min)
    frame.SetMaximum(bot_y_max)
    style_bottom_axes(frame, observable)
    frame.Draw()

    g_total = make_uncertainty_graph(low_edges, up_edges, total_unc, "g_total", ROOT.kBlack, 24)
    g_stat = make_uncertainty_graph(low_edges, up_edges, stat_unc, "g_stat", ROOT.kRed + 1, 20)
    g_scale = make_uncertainty_graph(low_edges, up_edges, scale_unc, "g_scale", ROOT.kGreen + 2, 21)
    g_pdf = make_uncertainty_graph(low_edges, up_edges, pdf_unc, "g_pdf", ROOT.kMagenta + 1, 22)
    g_shower = make_uncertainty_graph(low_edges, up_edges, shower_unc, "g_shower", ROOT.kOrange + 7, 23)
    g_model = make_uncertainty_graph(low_edges, up_edges, model_unc, "g_model", ROOT.kCyan + 2, 33)

    g_total.Draw("LP SAME")
    g_stat.Draw("LP SAME")
    g_scale.Draw("LP SAME")
    g_pdf.Draw("LP SAME")
    g_shower.Draw("LP SAME")
    g_model.Draw("LP SAME")

    legend.AddEntry(g_total, "Total uncertainty", "lp")
    legend.AddEntry(g_stat, "Stat uncertainty", "lp")
    legend.AddEntry(g_scale, "Scale uncertainty", "lp")
    legend.AddEntry(g_pdf, "PDF uncertainty", "lp")
    legend.AddEntry(g_shower, "Shower uncertainty", "lp")
    legend.AddEntry(g_model, "Model uncertainty", "lp")

    # Draw legend on top pad so it sits outside both plot areas
    canvas.cd()
    pad_top.cd()
    legend.Draw()

    canvas.Update()

    saved_files = []

    if SAVE_PDF:
        pdf_path = os.path.join(output_dir, base + "_E2E1.pdf")
        canvas.SaveAs(pdf_path)
        saved_files.append(pdf_path)

    if SAVE_PNG:
        png_path = os.path.join(output_dir, base + "_E2E1.png")
        canvas.SaveAs(png_path)
        saved_files.append(png_path)

    return saved_files


def main():
    if not os.path.isdir(INPUT_DIR):
        raise RuntimeError(f"INPUT_DIR does not exist: {INPUT_DIR}")

    csv_files = sorted(
        os.path.join(INPUT_DIR, f)
        for f in os.listdir(INPUT_DIR)
        if f.endswith(".csv")
    )

    if not csv_files:
        raise RuntimeError(f"No CSV files found in INPUT_DIR: {INPUT_DIR}")

    print(f"Found {len(csv_files)} CSV files in:")
    print(f"  {INPUT_DIR}")
    print()

    saved_count = 0

    for csv_path in csv_files:
        try:
            data_dict = load_csv_data(csv_path)
            hist, pair_label, observable = make_hist_from_csv(csv_path, data_dict)

            saved_files = save_hist_plots(
                hist=hist,
                csv_path=csv_path,
                output_dir=OUTPUT_DIR,
                data_dict=data_dict
            )

            print(f"Plotted: {os.path.basename(csv_path)}")
            for path in saved_files:
                print(f"  -> {path}")

            saved_count += 1

        except Exception as e:
            print(f"Failed on {csv_path}")
            print(f"  Reason: {e}")

    print()
    print(f"Finished. Successfully plotted {saved_count} CSV files.")
    print("Output directory:")
    print(f"  {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
