'''
Usage:python plot.py RootFile.root label[optional]

Script to make some quick efficiency plots to test ntuple generation.


Author: L. Dodd, UW Madison

'''

from subprocess import Popen
from sys import argv, exit, stdout, stderr

import ROOT

# So things don't look like crap.
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

######## File #########
if len(argv) < 2:
   print 'Usage:python plot.py RootFile.root label[optional]'
   exit()

infile = argv[1]
ntuple_file = ROOT.TFile(infile)

######## LABEL & SAVE WHERE #########

if len(argv)>2:
   saveWhere='~/private/output/tauAnalysis/'+argv[2]+'_'
else:
   saveWhere='~/private/output/tauAnalysis/'



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
EleVLooseMVA5 = ntuple_file.Get("againstElectronVLooseMVA5/Ntuple")
EleLooseMVA5 = ntuple_file.Get("againstElectronLooseMVA5/Ntuple")
EleMediumMVA5 = ntuple_file.Get("againstElectronMediumMVA5/Ntuple")

canvas = ROOT.TCanvas("asdf", "adsf", 800, 800)

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
        "goodReco==1",
        binning
    )
    return num

def make_denom(ntuple, variable,PtCut,binning):

    denom = make_plot(
        ntuple, variable,
        "genMatchedTau==1", #
        binning
    )
    return denom

def produce_efficiency(ntuple, variable, PtCut,binning, filename,color):
    denom = make_denom(ntuple, variable,PtCut,binning)
    num = make_num(ntuple,variable,PtCut,binning)
    l1 = make_efficiency(denom,num)
    l1.SetMarkerColor(color)
    return l1

def compare_efficiencies(ntuple1,legend1,ntuple2,legend2, variable, PtCut, binning, filename,
                         title='', xaxis='',yaxis=''):
    frame = ROOT.TH1F("frame", "frame", *binning)
    l1 = produce_efficiency(ntuple1,variable, PtCut,binning, filename,ROOT.kMagenta-3)
    l2 = produce_efficiency(ntuple2,variable, PtCut,binning, filename,ROOT.kBlue-9)
    frame.SetMaximum(1.2)
    frame.SetTitle(title)
    frame.GetXaxis().SetTitle(xaxis)
    frame.GetYaxis().SetTitle(yaxis)
    frame.Draw()
    l1.Draw('pe')
    l2.Draw('pesame')
    legend = ROOT.TLegend(0.5, 0.1, 0.89, 0.4, "", "brNDC")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1,legend1, "pe")
    legend.AddEntry(l2,legend2, "pe")
    legend.Draw()
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

def compare_3efficiencies(ntuple1,legend1,ntuple2, legend2,ntuple3, legend3, variable, PtCut, binning, filename,
                         title='', xaxis='',yaxis=''):
    frame = ROOT.TH1F("frame", "frame", *binning)
    l1 = produce_efficiency(ntuple1,variable, PtCut,binning, filename,ROOT.kMagenta-3)
    l2 = produce_efficiency(ntuple2,variable, PtCut,binning, filename,ROOT.kBlue-9)
    l3 = produce_efficiency(ntuple3,variable, PtCut,binning, filename,ROOT.kRed+3)
    frame.SetMaximum(1.2)
    frame.SetTitle(title)
    frame.GetXaxis().SetTitle(xaxis)
    frame.GetYaxis().SetTitle(yaxis)
    frame.Draw()
    l1.Draw('pe')
    l2.Draw('pesame')
    l3.Draw('pesame')
    legend = ROOT.TLegend(0.5, 0.7, 0.95, 0.9, "", "brNDC")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1,legend1, "pe")
    legend.AddEntry(l2,legend2, "pe")
    legend.AddEntry(l3,legend3, "pe")
    legend.Draw()
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

################################################################################
# Efficiency for a 20 GeV cut on tau Pt 
################################################################################

compare_3efficiencies(byLooseCmbIso3, 'byLooseCombIsoDBCorr3Hits', byMedCmbIso3,'byMediumCombIsoDBCorr3Hits', byTightCmbIso3,'byTightCombIsoDBCorr3Hits','tauPt', 20, [20, 0, 120],#variable, ptcut, binning
                    'tau_iso_effi_pT',#filename
                    "Tau Efficiency",#title
                    "pf Tau p_{T}",#xaxis
                    "efficiency" #yaxis             
)

compare_efficiencies(ntrlIsoPtSum,'neutralIsoPtSum',puCorrPtSum,'puCorrPtSum','tauPt',20,[20,0,120],
		    'tau_PtSum_effi_pT',
		    "Tau Efficiency",
		    "pf Tau p_{T} (GeV)"
		    "efficiency"
)

compare_efficiencies(MuLoose3,'againstMuonLoose3',MuTight3,'againstMuonTight3','tauPt',20,[20,0,120],
                    'tau_Mu_effi_pT',
                    "Tau Efficiency",
                    "pf Tau p_{T} (GeV)"
                    "efficiency"
)

compare_3efficiencies(EleVLooseMVA5,'againstElectronVLooseMVA5',EleLooseMVA5,'againstElectronLooseMVA5',EleMediumMVA5,'againstElectronMediumMVA5','tauPt',20,[20,0,120],
                    'tau_Ele_effi_pT',
                    "Tau Efficiency",
                    "pf Tau p_{T} (GeV)"
                    "efficiency"
)

## eta plots
compare_3efficiencies(byLooseCmbIso3, 'byLooseCombIsoDBCorr3Hits', byMedCmbIso3,'byMediumCombIsoDBCorr3Hits', byTightCmbIso3,'byTightCombIsoDBCorr3Hits','tauEta', 20, [20,-2.4,2.4],#variable, ptcut, binning
                    'tau_iso_effi_eta',#filename
                    "Tau Efficiency",#title
                    "Tau Eta",#xaxis
                    "efficiency" #yaxis             
)

compare_efficiencies(ntrlIsoPtSum,'neutralIsoPtSum',puCorrPtSum,'puCorrPtSum','tauEta',20,[20,-2.4,2.4],
                    'tau_PtSum_effi_eta',
                    "Tau Efficiency",
                    "Tau Eta",
                    "efficiency"
)

compare_efficiencies(MuLoose3,'againstMuonLoose3',MuTight3,'againstMuonTight3','tauEta',20,[20,-2.4,2.4],
                    'tau_Mu_effi_eta',
                    "Tau Efficiency",
                    "Tau Eta",
                    "efficiency"
)

compare_3efficiencies(EleVLooseMVA5,'againstElectronVLooseMVA5',EleLooseMVA5,'againstElectronLooseMVA5',EleMediumMVA5,'againstElectronMediumMVA5','tauEta',20,[20,-2.4,2.4],
                    'tau_Ele_effi_eta',
                    "Tau Efficiency",
                    "Tau Eta",
                    "efficiency"
)

## nvtx plots
compare_3efficiencies(byLooseCmbIso3, 'byLooseCombIsoDBCorr3Hits', byMedCmbIso3,'byMediumCombIsoDBCorr3Hits', byTightCmbIso3,'byTightCombIsoDBCorr3Hits','nvtx', 20, [20,0,35],#variable, ptcut, binning
                    'tau_iso_effi_nvtx',#filename
                    "Tau Efficiency",#title
                    "N_{vtx}",#xaxis
                    "efficiency" #yaxis             
)

compare_efficiencies(ntrlIsoPtSum,'neutralIsoPtSum',puCorrPtSum,'puCorrPtSum','nvtx',20,[20,0,35],
                    'tau_PtSum_effi_nvtx',
                    "Tau Efficiency",
                    "N_{vtx}",
                    "efficiency"
)

compare_efficiencies(MuLoose3,'againstMuonLoose3',MuTight3,'againstMuonTight3','nvtx',20,[20,0,35],
                    'tau_Mu_effi_nvtx',
                    "Tau Efficiency",
                    "N_{vtx}",
                    "efficiency"
)

compare_3efficiencies(EleVLooseMVA5,'againstElectronVLooseMVA5',EleLooseMVA5,'againstElectronLooseMVA5',EleMediumMVA5,'againstElectronMediumMVA5','nvtx',20,[20,0,35],
                    'tau_Ele_effi_nvtx',
                    "Tau Efficiency",
                    "N_{vtx}",
                    "efficiency"
)
