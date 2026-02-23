import os
import ROOT

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")

samples = [
    "WCH7minus", "WCH7plus",
    "WCHPartonminus", "WCHPartonplus",
    "WCPy8minus", "WCPy8plus",
    "WCPyPartonminus", "WCPyPartonplus"
]

for s in samples:
    path = os.path.join(DATA_DIR, f"WCharm_{s}.root")
    f = ROOT.TFile(path, "READ")
    t = f.Get("WCharmTree")
    t.GetEntry(0)
    n = len(t.weightvec)
    print(f"{s:14s}  len(weightvec)={n}  has292={n>292}  has315={n>315}")
    f.Close()
