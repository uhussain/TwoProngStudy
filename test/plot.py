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


ntuple = ntuple_file.Get("tauAnalyzer/Ntuple")

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
    eff.SetMarkerColor(ROOT.kMagenta+3)
    eff.SetMarkerSize(1.5)
    eff.SetLineColor(ROOT.kBlack)
    return eff

def compare_efficiencies(ntuple, variable, PtCut, binning, filename,
                         title='', xaxis='',yaxis=''):
    denom = make_plot(
        ntuple, variable,
        "", # No selection
        binning
    )

    num = make_plot(
        ntuple, variable,
        "pt > %0.2f && passDiscr>0" % (PtCut),
        binning
    )
    frame = ROOT.TH1F("frame", "frame", *binning)
    l1 = make_efficiency(denom, num)
    frame.SetMaximum(1.2)
    frame.SetTitle(title)
    frame.GetXaxis().SetTitle(xaxis)
    frame.GetYaxis().SetTitle(yaxis)
    frame.Draw()
    l1.Draw('pe')
    legend = ROOT.TLegend(0.7, 0.2, 0.89, 0.4, "", "brNDC")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1, "ByLooseIsolation", "pe")
    legend.Draw()
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

################################################################################
# Efficiency for a 20 GeV cut on tau Pt 
################################################################################

compare_efficiencies(ntuple, 'pt',  20, [20, 0, 120],#variable, ptcut, binning
                    'tauPt_20',#filename
                    "Tau Efficiency ByLooseIsolation",#title
                    "pf Tau p_{T} (GeV)",#xaxis
                    "Pf Tau pass tauID/(all Pf tau)" #yaxis             
)

