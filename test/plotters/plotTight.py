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
if len(argv) < 1:
   print 'Usage:python plot.py RootFile.root label[optional]'
   exit()

infile = argv[1]
infile2= argv[2]
infile3= argv[3]
infile4= argv[4]
infile5= argv[5]
ntuple_file = ROOT.TFile(infile)
ntuple_file2 = ROOT.TFile(infile2)
ntuple_file3 = ROOT.TFile(infile3)
ntuple_file4 = ROOT.TFile(infile4)
ntuple_file5 = ROOT.TFile(infile5)

######## LABEL & SAVE WHERE #########

if len(argv)>6:
   saveWhere='~/www/Research/'+argv[6]+'_'
else:
   saveWhere='~/www/Research/'



#####################################
#Get Effi NTUPLE                 #
#####################################
byTightCmbIso3_03 = ntuple_file.Get("byTightCombinedIsolationDBSumPtCorr3Hits/Ntuple")
byTightCmbIso3_05 = ntuple_file2.Get("byTightCombinedIsolationDBSumPtCorr3Hits/Ntuple")
byTightCmbIso3_08 = ntuple_file3.Get("byTightCombinedIsolationDBSumPtCorr3Hits/Ntuple")
byTightCmbIso3_1 = ntuple_file4.Get("byTightCombinedIsolationDBSumPtCorr3Hits/Ntuple")
byTightCmbIso3_5 = ntuple_file5.Get("byTightCombinedIsolationDBSumPtCorr3Hits/Ntuple")

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
        "tauPt > %0.2f &&fabs(tauEta)<2.3&& dmf>0 &&genMatchedTau==1&&passDiscr>0" % (PtCut),
        binning
    )
    return num
def make_denom(ntuple, variable,PtCut,binning):
    denom = make_plot(
        ntuple, variable,
        "dmf>0&&fabs(tauEta)<2.3&&genMatchedTau==1&&tauPt> %0.2f"%(PtCut), #
        binning
    )
    return denom

def produce_efficiency(ntuple, variable, PtCut,binning, filename,color):
    denom = make_denom(ntuple, variable,PtCut,binning)
    num = make_num(ntuple,variable,PtCut,binning)
    l1 = make_efficiency(denom,num)
    l1.SetMarkerColor(color)
    return l1

def compare_efficiencies(ntuple1,legend1,ntuple2,legend2,ntuple3,legend3,ntuple4,legend4,ntuple5,legend5, variable, PtCut, binning, filename,
                         title='', xaxis='',yaxis=''):
    frame = ROOT.TH1F("frame", "frame", *binning)
    l1 = produce_efficiency(ntuple1,variable, PtCut,binning, filename,ROOT.kMagenta-3)
    l2 = produce_efficiency(ntuple2,variable, PtCut,binning, filename,ROOT.kBlue-9)
    l3 = produce_efficiency(ntuple3,variable, PtCut,binning, filename,ROOT.kOrange-9)
    l4 = produce_efficiency(ntuple4,variable, PtCut,binning, filename,ROOT.kBlue+2)
    l5 = produce_efficiency(ntuple5,variable, PtCut,binning, filename,ROOT.kGreen-3)
    frame.SetMaximum(1.2)
    frame.SetTitle(title)
    frame.GetXaxis().SetTitle(xaxis)
    frame.GetYaxis().SetTitle(yaxis)
    frame.Draw()
    l1.Draw('pe')
    l2.Draw('pesame')
    l3.Draw('pesame')
    l4.Draw('pesame')
    l5.Draw('pesame')
    legend = ROOT.TLegend(0.5, 0.1, 0.89, 0.4, "", "brNDC")
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
compare_efficiencies(
	            byTightCmbIso3_03, 'Tight dBIso 0.3', 
	            byTightCmbIso3_05, 'Tight dBIso 0.5', 
	            byTightCmbIso3_08, 'Tight dBIso 0.8', 
		    byTightCmbIso3_1,'Tight dBIso 1',
		    byTightCmbIso3_5,'Tight dBIso 5',
		    'tauPt', 20, [20, 0, 120],#variable, ptcut, binning
                    'tauPtTightCmbIso3',#filename
                    "byTightCombIsoDBSumPtCorr3Hits",#title
                    "pf Tau p_{T} (GeV)",#xaxis
                    "efficiency" #yaxis             
)
