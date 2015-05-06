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
   saveWhere='~/www/Research/'+argv[2]+'_'
else:
   saveWhere='~/www/Research/'



#####################################
#Get Effi NTUPLE                 #
#####################################


LooseIso = ntuple_file.Get("byLooseIsolation/Ntuple")
VLooseIso = ntuple_file.Get("byVLooseIsolation/Ntuple")
byLooseCmbIso = ntuple_file.Get("byLooseCombinedIsolationDBSumPtCorr/Ntuple")
byVLooseCmbIso = ntuple_file.Get("byVLooseCombinedIsolationDBSumPtCorr/Ntuple")

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
        "pt[0] > %0.2f && dmf[0]>0 &&genMatchedTau[0]==1&&passDiscr[0]>0" % (PtCut),
        binning
    )
    return num
def make_denom(ntuple, variable,PtCut,binning):
    denom = make_plot(
        ntuple, variable,
        "dmf[0]>0&&genMatchedTau[0]==1&&pt[0]> %0.2f"%(PtCut), #
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
    legend = ROOT.TLegend(0.6, 0.2, 0.89, 0.5, "", "brNDC")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1,legend1, "pe")
    legend.AddEntry(l2,legend2, "pe")
    legend.Draw()
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

################################################################################
# Efficiency for a 20 GeV cut on tau Pt 
################################################################################

compare_efficiencies(LooseIso, "LooseIsolation",VLooseIso,"VLooseIso",'pt[0]', 0, [20, 0, 120],#variable, ptcut, binning
                    'tauPt_0',#filename
                    "Tau Efficiency",#title
                    "pf Tau p_{T} (GeV)",#xaxis
                    "(Pf Tau gen dmf tauID)/(dmf gen pftau)" #yaxis             
)

compare_efficiencies(byLooseCmbIso, 'byLooseCombIsoDBSumPtCorr', byVLooseCmbIso,'byVLooseCombIsoDBSumPtCorr','pt[0]', 0, [20, 0, 120],#variable, ptcut, binning
                    'tauPtcmb_0',#filename
                    "Tau Efficiency",#title
                    "pf Tau p_{T} (GeV)",#xaxis
                    "(Pf Tau gen dmf tauID)/(dmf gen pftau)" #yaxis             
)

compare_efficiencies(byLooseCmbIso, 'byLooseCombIsoDBSumPtCorr', LooseIso,'LooseIsolation','pt[0]', 0, [20, 0, 120],#variable, ptcut, binning
                    'tauPtloose_0',#filename
                    "Tau Efficiency",#title
                    "pf Tau p_{T} (GeV)",#xaxis
                    "(Pf Tau gen dmf tauID)/(dmf gen pftau)" #yaxis             
)

compare_efficiencies(byVLooseCmbIso, 'byVLooseCombIsoDBSumPtCorr', VLooseIso,'VLooseIsolation','pt[0]', 0, [20, 0, 120],#variable, ptcut, binning
                    'tauPtvloose_0',#filename
                    "Tau Efficiency",#title
                    "pf Tau p_{T} (GeV)",#xaxis
                    "(Pf Tau gen dmf tauID)/(dmf gen pftau)" #yaxis             
)
