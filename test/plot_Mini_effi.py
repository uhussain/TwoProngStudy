'''
Usage:python plot.py RootFile.root label[optional]

Script to make some quick efficiency plots to test ntuple generation.


Author: L. Dodd, UW Madison

'''

from subprocess import Popen
from sys import argv, exit, stdout, stderr
import os
import ROOT

# So things don't look like crap.
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadTopMargin(0.1)
ROOT.gStyle.SetGridStyle(0)
ROOT.gStyle.SetTitleX(0.2)
ROOT.gStyle.SetTitleW(0.6)
######## File #########
if len(argv) < 2:
   print 'Usage:python plot.py RootFile.root label[optional]'
   exit()

infile = argv[1]
ntuple_file = ROOT.TFile(infile)

######## LABEL & SAVE WHERE #########

if len(argv)>2:
   saveWhere='~/private/output/tauAnalysis/7_6_5/'+argv[2]+'_'
else:
   saveWhere='~/private/output/tauAnalysis/7_6_5/'



#####################################
#Get Effi NTUPLE                 #
#####################################

byLooseCmbIso3 = ntuple_file.Get("byLooseCombinedIsolationDeltaBetaCorr3Hits/Ntuple")
byMedCmbIso3 = ntuple_file.Get("byMediumCombinedIsolationDeltaBetaCorr3Hits/Ntuple")
byTightCmbIso3 = ntuple_file.Get("byTightCombinedIsolationDeltaBetaCorr3Hits/Ntuple")

ntrlIsoPtSum = ntuple_file.Get("neutralIsoPtSum/Ntuple")
puCorrPtSum = ntuple_file.Get("puCorrPtSum/Ntuple")
MuLoose3 = ntuple_file.Get("againstMuonLoose3/Ntuple")
MuTight3 = ntuple_file.Get("againstMuonTight3/Ntuple")
EleVLooseMVA6 = ntuple_file.Get("againstElectronVLooseMVA6/Ntuple")
EleLooseMVA6 = ntuple_file.Get("againstElectronLooseMVA6/Ntuple")
EleMediumMVA6 = ntuple_file.Get("againstElectronMediumMVA6/Ntuple")
EleTightMVA6 = ntuple_file.Get("againstElectronTightMVA6/Ntuple")
EleVTightMVA6 = ntuple_file.Get("againstElectronVTightMVA6/Ntuple")
canvas = ROOT.TCanvas("asdf", "adsf", 1200, 800)

def make_plot(tree, variable, selection, binning, xaxis='', title=''):
    ''' Plot a variable using draw and return the histogram '''
    draw_string = "%s>>htemp(%s)" % (variable, ", ".join(str(x) for x in binning))
    tree.Draw(draw_string, selection, "goff")
    output_histo = ROOT.gDirectory.Get("htemp").Clone()
    output_histo.GetXaxis().SetTitle(xaxis)
    output_histo.SetTitle(title)
    return output_histo

def make_efficiency(denom, num):
    ''' Make an efficiency graph with the style '''
    eff = ROOT.TGraphAsymmErrors(num, denom)
    eff.SetMarkerStyle(20)
    eff.SetMarkerSize(1.5)
    eff.SetLineColor(ROOT.kBlack)
    return eff

def make_num(ntuple, variable,PtCut,binning):
    num = make_plot(
        ntuple, variable,
	"goodReco",
        binning
    )
    return num

def make_denom(ntuple, variable,PtCut,binning):
    denom = make_plot(
        ntuple, variable,
        "",
        binning
    )
    return denom

def produce_efficiency(ntuple, variable, PtCut,binning, filename,color):
    denom = make_denom(ntuple, variable,PtCut,binning)
    num = make_num(ntuple,variable,PtCut,binning)
    l1 = make_efficiency(denom,num)
    l1.SetMarkerColor(color)
    return l1

def compare_efficiencies(ntuple1,legend1,ntuple2, legend2, variable, PtCut, binning, filename,framemin,framemax,
                         title='', xaxis='',yaxis=''):
    frame = ROOT.TH1F("frame", "frame", *binning)
    l1 = produce_efficiency(ntuple1,variable, PtCut,binning, filename,ROOT.kMagenta-3)
    l2 = produce_efficiency(ntuple2,variable, PtCut,binning, filename,ROOT.kBlue-9)
    frame.SetMaximum(framemax)
    frame.SetMinimum(framemin)
    frame.SetTitle(title)
    frame.UseCurrentStyle()
    frame.GetXaxis().SetTitle(xaxis)
    frame.GetYaxis().SetTitle(yaxis)
    frame.GetYaxis().SetTitleOffset(1.2)
    frame.GetYaxis().CenterTitle()
    frame.GetXaxis().CenterTitle()
    frame.UseCurrentStyle()
    frame.Draw()
    l1.Draw('pe')
    l2.Draw('pesame')
    legend = ROOT.TLegend(0.6, 0.75, 0.9, 0.9, "", "brNDC")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1,legend1, "pe")
    legend.AddEntry(l2,legend2, "pe")
    legend.Draw()
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

def compare_3efficiencies(ntuple1,legend1,ntuple2, legend2,ntuple3, legend3, variable, PtCut, binning, filename,framemin,framemax,
                         title='', xaxis='',yaxis=''):
    frame = ROOT.TH1F("frame", "frame", *binning)
    l1 = produce_efficiency(ntuple1,variable, PtCut,binning, filename,ROOT.kMagenta-3)
    l2 = produce_efficiency(ntuple2,variable, PtCut,binning, filename,ROOT.kBlue-9)
    l3 = produce_efficiency(ntuple3,variable, PtCut,binning, filename,ROOT.kRed+3)
    frame.SetMinimum(framemin)
    frame.SetMaximum(framemax)
    frame.SetTitle(title)
    frame.UseCurrentStyle()
    frame.GetXaxis().SetTitle(xaxis)
    frame.GetYaxis().SetTitle(yaxis)
    frame.GetYaxis().SetTitleOffset(1.2)
    frame.GetYaxis().CenterTitle()
    frame.GetXaxis().CenterTitle()
    frame.Draw()
    l1.Draw('pe')
    l2.Draw('pesame')
    l3.Draw('pesame')
    
    legend = ROOT.TLegend(0.6, 0.7, 0.9, 0.9, "", "brNDC")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1,legend1, "pe")
    legend.AddEntry(l2,legend2, "pe")
    legend.AddEntry(l3,legend3, "pe")
    legend.Draw()
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

def compare_5efficiencies(ntuple1,legend1,ntuple2, legend2,ntuple3, legend3,ntuple4,legend4,ntuple5,legend5, variable, PtCut, binning, filename,framemin,framemax,
                         title='', xaxis='',yaxis=''):
    frame = ROOT.TH1F("frame", "frame", *binning)
    l1 = produce_efficiency(ntuple1,variable, PtCut,binning, filename,ROOT.kRed-3)
    l2 = produce_efficiency(ntuple2,variable, PtCut,binning, filename,ROOT.kRed+3)
    l3 = produce_efficiency(ntuple3,variable, PtCut,binning, filename,ROOT.kBlue+3)
    l4 = produce_efficiency(ntuple4,variable, PtCut,binning, filename,ROOT.kGreen-3)
    l5 = produce_efficiency(ntuple5,variable, PtCut,binning, filename,ROOT.kGreen+3)
    frame.SetMinimum(framemin)
    frame.SetMaximum(framemax)
    frame.SetTitle(title)
    frame.UseCurrentStyle()
    frame.GetXaxis().SetTitle(xaxis)
    frame.GetYaxis().SetTitle(yaxis)
    frame.GetYaxis().SetTitleOffset(1.2)
    frame.GetYaxis().CenterTitle()
    frame.GetXaxis().CenterTitle()
    frame.Draw()
    l1.Draw('pe')
    l2.Draw('pesame')
    l3.Draw('pesame')
    l4.Draw('pesame')
    l5.Draw('pesame')
    
    legend = ROOT.TLegend(0.6, 0.65, 0.9, 0.9, "", "brNDC")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1,legend1, "pe")
    legend.AddEntry(l2,legend2, "pe")
    legend.AddEntry(l3,legend3, "pe")
    legend.AddEntry(l4,legend4, "pe")
    legend.AddEntry(l5,legend5, "pe")
    legend.Draw()
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

################################################################################
# Efficiency for a 20 GeV cut on tau Pt 
################################################################################
## pT plots
compare_3efficiencies(byLooseCmbIso3, 'byLooseCombIsoDBCorr3Hits', byMedCmbIso3,'byMediumCombIsoDBCorr3Hits', byTightCmbIso3,'byTightCombIsoDBCorr3Hits','tauPt', 20, [20, 0, 120],#variable, ptcut, binning
                    'tau_iso_effi_pT', 0, 1,#filename
                    "Tau Efficiency",#title
                    "Tau p_{T} (GeV)",#xaxis
                    "Efficiency" #yaxis
             
)

compare_efficiencies(ntrlIsoPtSum,'neutralIsoPtSum',puCorrPtSum,'puCorrPtSum','tauPt',20,[20,0,120],
                    'tau_PtSum_effi_pT', 0, 1,
                    "Tau Efficiency",
                    "Tau p_{T} (GeV)",
                    "Efficiency"
)

compare_efficiencies(MuLoose3,'againstMuonLoose3',MuTight3,'againstMuonTight3','tauPt',20,[20,0,120],
                    'tau_Mu_effi_pT', 0, 1,
                    "Tau Efficiency",
                    "Tau p_{T} (GeV)",
                    "Efficiency"
)

compare_5efficiencies(EleVLooseMVA6,'againstElectronVLooseMVA6',EleLooseMVA6,'againstElectronLooseMVA6',EleMediumMVA6,'againstElectronMediumMVA6',EleTightMVA6,'againstElectronTightMVA6',EleVTightMVA6,'againstElectronVTightMVA6','tauPt',20,[20,0,120],
                    'tau_Ele_effi_pT', 0, 1,
                    "Tau Efficiency",
                    "Tau p_{T} (GeV)",
                    "Efficiency"
)

## eta plots
compare_3efficiencies(byLooseCmbIso3, 'byLooseCombIsoDBCorr3Hits', byMedCmbIso3,'byMediumCombIsoDBCorr3Hits', byTightCmbIso3,'byTightCombIsoDBCorr3Hits','tauEta', 20, [20,-2.4,2.4],#variable, ptcut, binning
                    'tau_iso_effi_eta', 0, 1,#filename
                    "Tau Efficiency",#title
                    "Tau #eta",#xaxis
                    "Efficiency" #yaxis             
)
compare_efficiencies(ntrlIsoPtSum,'neutralIsoPtSum',puCorrPtSum,'puCorrPtSum','tauEta',20,[20,-2.4,2.4],
                    'tau_PtSum_effi_eta',0, 1,
                    "Tau Efficiency",
                    "Tau #eta",
                    "Efficiency"
)
compare_efficiencies(MuLoose3,'againstMuonLoose3',MuTight3,'againstMuonTight3','tauEta',20,[20,-2.4,2.4],
                    'tau_Mu_effi_eta',0, 1,
                    "Tau Efficiency",
                    "Tau #eta",
                    "Efficiency"
)

compare_5efficiencies(EleVLooseMVA6,'againstElectronVLooseMVA6',EleLooseMVA6,'againstElectronLooseMVA6',EleMediumMVA6,'againstElectronMediumMVA6',EleTightMVA6,'againstElectronTightMVA6',EleVTightMVA6,'againstElectronVTightMVA6','tauEta',20,[20,-2.4,2.4],
                    'tau_Ele_effi_eta',0, 1,
                    "Tau Efficiency",
                    "Tau #eta",
                    "Efficiency"
)

## nvtx plots
compare_3efficiencies(byLooseCmbIso3, 'byLooseCombIsoDBCorr3Hits', byMedCmbIso3,'byMediumCombIsoDBCorr3Hits', byTightCmbIso3,'byTightCombIsoDBCorr3Hits','nvtx', 20, [20,0,35],#variable, ptcut, binning
                    'tau_iso_effi_nvtx',0, 1,#filename
                    "Tau Efficiency",#title
                    "N_{vtx}",#xaxis
                    "Efficiency" #yaxis             
)

compare_efficiencies(ntrlIsoPtSum,'neutralIsoPtSum',puCorrPtSum,'puCorrPtSum','nvtx',20,[20,0,35],
                    'tau_PtSum_effi_nvtx',0, 1,
                    "Tau Efficiency",
                    "N_{vtx}",
                    "Efficiency"
)

compare_efficiencies(MuLoose3,'againstMuonLoose3',MuTight3,'againstMuonTight3','nvtx',20,[20,0,35],
                    'tau_Mu_effi_nvtx',0, 1,
                    "Tau Efficiency",
                    "N_{vtx}",
                    "Efficiency"
)

compare_5efficiencies(EleVLooseMVA6,'againstElectronVLooseMVA6',EleLooseMVA6,'againstElectronLooseMVA6',EleMediumMVA6,'againstElectronMediumMVA6',EleTightMVA6,'againstElectronTightMVA6',EleVTightMVA6,'againstElectronVTightMVA6','nvtx',20,[20,0,35],
                    'tau_Ele_effi_nvtx',0, 1,
                    "Tau Efficiency",
                    "N_{vtx}",
                    "Efficiency"
)
