from ROOT import TMath, TStyle, TF1, TFile, TCanvas, gDirectory, TTree, TH1F, TH2F, TProfile2D , TProfile, TBrowser, gStyle, gMinuit, TText, TCut, TLorentzVector, kRed, kBlue, kGreen, TRandom3, TChain, TPad, TLegend, TColor, gStyle
from array import *
import sys
import os, glob, re

gStyle.SetOptStat(0)

doNorm=False
plotRatio=True

minimum=100000000
maximum=0

## Colours:

colors = [1,800,807,900,906,860,866,852,843,849,920,922,400,398,619,616,615,602,417,432,616, 800, 820,840,860,880,900,920,940]

#Histogram Lists:
histos=[]
histotitles=[]

histos.append("ptlepton")
histos.append("etalepton")

histotitles.append(" ; p_{T} (lepton); cross-section [nb] ")
histotitles.append(" ; #eta (lepton); cross-section [nb] ")

histos.append("met_et")
histotitles.append(" ; E_{T}^{miss} [GeV]; cross-section [nb] ")

histos.append("mtw")
histotitles.append(" ; m_{T}^{W} [GeV]; cross-section [nb] ")

histos.append("ptjet")
histotitles.append(" ; p_{T}^{jet} [GeV]; cross-section [nb] ")

name ="W+c"

## OLD (For no weightingg)
#settings=[]
#settings.append("Default")

# NEW (Automatically finds all the settings in ../outputs (e.g. Default, Var3cUp, Var3cDown))
settings = []
pat = re.compile(r"output_(.+?)_(.+)\.root$")
for f in glob.glob("../outputs/output_*.root"):
    m = pat.search(os.path.basename(f))
    if m:
        settings.append(m.group(1))
settings = sorted(list(set(settings)))
print("Found settings:", settings)


procs=[]
procs.append("WCH7minus")
procs.append("WCHPartonminus")
procs.append("WCH7plus")
procs.append("WCHPartonplus")
procs.append("WCPy8minus")
procs.append("WCPyPartonminus")
procs.append("WCPy8plus")
procs.append("WCPyPartonplus")

extra_string=""
extra_string_onplot=""

histograms=[]
legendname=[]

can=TCanvas()
h=TH1F()

histocounter=-1

for histo in histos:
    histocounter += 1
    histograms=[]
    legendname=[]

    for setting in settings:
        for proc in procs:
            ## OLD
            """
            outputfile = TFile("../outputs/output_"+setting+"_"+proc+".root", "open")

            h = TH1F()
            h = outputfile.Get(histo)
            h.SetDirectory(0)
            """

            ## NEW (prevetns crashing if the combo doesnt exist)
            fname = "../outputs/output_"+setting+"_"+proc+".root"
            if not os.path.exists(fname):
                continue

            outputfile = TFile(fname, "open")
            h = outputfile.Get(histo)
            if not h:
                outputfile.Close()
                continue

            h.SetDirectory(0)
            outputfile.Close()
            
            legendname.append(setting+"_"+proc)

            h.Draw()
            
            if (minimum>h.GetMinimum()):
                minimum=h.GetMinimum()

            if (maximum<h.GetMaximum()):
                maximum=h.GetMaximum()

            histograms.append(h)

    ##################################################################

    can = TCanvas("can_"+histo+extra_string,"can_"+histo+extra_string,800,600)
    can.Draw()
    can.cd(1)
    padNameUpper = "padUpper_"+histo
    padNameLower = "padLower_"+histo

    extraCan = 1.25
    padScaling = 1
    ratioPadScaling = 1

    padUpper = TPad()
    padLower = TPad()

    if plotRatio == True:
        ratioPadRatio = 1-(1-.16+.015) / extraCan
        padScaling = 1. / (1. - ratioPadRatio) / extraCan
        ratioPadScaling = (1. / ratioPadRatio) / extraCan

        padUpper = TPad(padNameUpper, padNameUpper, 0., 0., 1., 1.)
        padUpper.Draw()
    else:
        can.cd()
        padUpper = TPad(padNameUpper, padNameUpper, 0., 0., 1., 1.)
        padUpper.Draw()

    can.cd()
    can.SetWindowPosition(10,100)
    padUpper.cd()
    padUpper.SetBottomMargin(0.30)
    
    if histo != "etalepton":
        padUpper.SetLogx()

    can.Modified()
    can.Update()
    
    ##################################################################
    
    for colorcounter in range(0,len(histograms)):
        print(colorcounter)

        if colorcounter < len(colors):
            print("color")
            print(colors[colorcounter])
            histograms[colorcounter].SetLineColor(colors[colorcounter])
            histograms[colorcounter].SetMarkerColor(colors[colorcounter])

        histograms[colorcounter].SetLineWidth(2)
        if colorcounter>2:
            histograms[colorcounter].SetLineStyle(2)
        if colorcounter>5:
            histograms[colorcounter].SetLineStyle(5)
        if colorcounter>8:
            histograms[colorcounter].SetLineStyle(9)

        histograms[colorcounter].SetTitle(histotitles[histocounter])
        histograms[colorcounter].SetMarkerStyle(0)
        
        if (colorcounter==0):
            histograms[colorcounter].SetMinimum(minimum)
            histograms[colorcounter].SetMaximum(maximum*1.1)
            histograms[colorcounter].Draw("PE")

            if histo == "etalepton":
                histograms[colorcounter].GetXaxis().SetLimits(-3, 3)

            print("Drawing")
        elif (colorcounter>5):
            histograms[colorcounter].SetMarkerStyle(14+colorcounter)
            histograms[colorcounter].Draw("PE same")
        else:
            histograms[colorcounter].Draw("PE same")
        
        can.Modified()
        can.Update()
        can.Modified() 
        can.Update() 
        can.SaveAs("../plots/"+can.GetName()+".eps")
        
        if (colorcounter==0):
            x1 = 0.59
            y1 = 0.60
            x2 = 0.90
            y2 = 0.88
            y1 = y2 - (y2 - y1) * padScaling
            
            legendHeight = 1.5
            nLegendRows = 2
            legendHeight *= (y2 - y1) * (nLegendRows) / 5.
            y1 = y2 - legendHeight
            
            leg = TLegend(0.15, 0.02, 0.60, 0.22)
            leg.SetFillColor(0)
            leg.SetTextFont(42)
            leg.SetTextSize(0.02)
            leg.SetMargin(0.15)
            leg.SetBorderSize(1)
            leg.SetEntrySeparation(0.2)
            leg.SetMargin(0.12)
            
        leg.AddEntry(histograms[colorcounter],  legendname[colorcounter], "ple")
    
    if (extra_string_onplot!=""):
        leg.AddEntry(h,extra_string_onplot,"")
    
    leg.Draw()

    if plotRatio==True:

        canRatio = TCanvas("canRatio_"+histo+extra_string,"canRatio_"+histo+extra_string,800,600)
        canRatio.Draw()
        canRatio.cd()
        canRatio.SetBottomMargin(0.30)

        if histo != "etalepton":
            canRatio.SetLogx()
        
        ratios = []
        firstRatio = True

        legRatio = TLegend(0.15, 0.02, 0.60, 0.22)
        legRatio.SetFillColor(0)
        legRatio.SetBorderSize(1)
        legRatio.SetTextSize(0.018)
        legRatio.SetEntrySeparation(0.08)
        legRatio.SetMargin(0.12)

        for idx in range(1, len(histograms), 2):

            parton_idx   = idx
            particle_idx = idx - 1

            ratio = histograms[parton_idx].Clone("ratio_"+str(parton_idx))
            ratio.Divide(histograms[particle_idx])

            ratio.SetMinimum(0.00)
            ratio.SetMaximum(1.15)
            ratio.GetYaxis().SetTitle("parton / particle [nb]/[nb]")

            ratio.SetLineColor(histograms[parton_idx].GetLineColor())
            ratio.SetMarkerColor(histograms[parton_idx].GetMarkerColor())
            ratio.SetLineStyle(histograms[parton_idx].GetLineStyle())
            ratio.SetLineWidth(histograms[parton_idx].GetLineWidth())
            ratio.SetMarkerStyle(0)

            ratios.append(ratio)

            if firstRatio:
                ratio.Draw("PE")
                firstRatio = False
            else:
                ratio.Draw("PE same")

            label = legendname[parton_idx] + " / " + legendname[particle_idx]
            legRatio.AddEntry(ratio, label, "l")
            
            canRatio.Modified()
            canRatio.Update()

        legRatio.Draw()
        canRatio.Modified()
        canRatio.Update()
        canRatio.SaveAs("../plots/"+canRatio.GetName()+".png")
        canRatio.SaveAs("../plots/"+canRatio.GetName()+".root")

    print(can.GetName())
    can.cd()
    can.Modified()
    can.Update() 
    can.SaveAs("../plots/"+can.GetName()+".png")
    can.SaveAs("../plots/"+can.GetName()+".root")
