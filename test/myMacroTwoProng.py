import ROOT as rt
import CMS_lumi, tdrstyle
import array
from sys import argv,exit

rt.gROOT.SetBatch(True) #Running in batch mode prevents keeps plots from being drawn; much quicker
decayModes = ["h^{#pm}","h^{#pm}#pi_{0}","h^{#pm}h^{#mp}","h^{#pm}h^{#mp}#pi_{0}","h^{#pm}h^{#mp}h^{#pm}", "11"]
######## File #########
if len(argv) < 2:
   print 'Usage:python plot.py RootFile.root label[optional]'
   exit()

infile = argv[1]
ntuple_file = rt.TFile(infile,"READ")

byLooseCmbIso3 = ntuple_file.Get("byLooseCombinedIsolationDeltaBetaCorr3Hits/Ntuple")
byMedCmbIso3 = ntuple_file.Get("byMediumCombinedIsolationDeltaBetaCorr3Hits/Ntuple")
byTightCmbIso3 = ntuple_file.Get("byTightCombinedIsolationDeltaBetaCorr3Hits/Ntuple")
byNone = ntuple_file.Get("byNone/Ntuple")

######## LABEL & SAVE WHERE #########
if len(argv)>2:
   saveWhere='~/private/output/2prong/'+argv[2]+'_'
else:
   saveWhere='~/private/output/2prong/'

######## Style Choices ########
tdrstyle.ygrid = True
tdrstyle.logy = True
tdrstyle.logx = True
tdrstyle.setTDRStyle()
colors=[1,2,4,6,8,12]  #Add more colors for >5 sets of points in a single plot
styles=[20,21,22,23,33,34]  #Add more point styles if needed

#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation Preliminary"
#CMS_lumi.lumiText = "Simulation: Z* #rightarrow #tau#tau"
CMS_lumi.lumiText = "ggH #rightarrow #tau#tau"
iPos = 11

H_ref = 700; 
W_ref = 800; 
W = W_ref
H  = H_ref

# references for T, B, L, R
T = 0.08*H_ref
B = 0.12*H_ref 
L = 0.12*W_ref
R = 0.15*W_ref

def CompareEfficiencies(ntuple,filename,title,variable, maxdigit):
    canvas = rt.TCanvas("c2","c2",50,50,W,H)
    canvas.DrawFrame(1e-6,1e-6,1,1)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin( L/W )
    canvas.SetRightMargin( R/W )
    canvas.SetTopMargin( T/H )
    canvas.SetBottomMargin( B/H )
    canvas.SetTickx(0)
    canvas.SetTicky(1)
    h = rt.TH2F("h","hist", 100,0.1,1000,100,0.1,1000)
    h.SetStats(rt.kFALSE)
    ntuple.Draw("pT3:pT2>>h", "", "colz")
    h.SetMarkerStyle(1)
    xAxis = h.GetXaxis()
    #xAxis.SetNdivisions(6,5,0)
    xAxis.SetTitleOffset(1.1)
    xAxis.SetLimits(1e-1,1000)
    #xAxis.SetLimits(1e-4,1)
    xAxis.SetTitle("2^{nd} track: " + title + " isolation")

    yAxis = h.GetYaxis()
    #yAxis.SetNdivisions(6,5,0)
    yAxis.SetTitleOffset(1.1)
    yAxis.SetLimits(1e-1,1000)
    #yAxis.SetLimits(1e-4,1)
    yAxis.SetTitle("3^{rd} (failed) track")
    h.Draw("colz")
    rt.TGaxis.SetMaxDigits(maxdigit)
    #graph.Draw("pe")

    text = rt.TText(0.9,0.97,variable)
    text.SetTextAlign(22)
    text.SetTextColor(rt.kBlack)
    text.SetTextFont(43)
    text.SetTextSize(40)
    text.SetNDC()
    text.Draw()
    #draw the lumi text on the canvas
    CMS_lumi.CMS_lumi(canvas, 0, iPos)

    canvas.cd()
    canvas.Update()
    canvas.RedrawAxis()
    frame = canvas.GetFrame()
    frame.Draw()

    #save image of plot
    saveas = saveWhere + filename + '.png'
    canvas.SaveAs(saveas)

#Call function to create/save plots
CompareEfficiencies(byLooseCmbIso3, 'pTloose2', 'Loose', 'pT',6)
CompareEfficiencies(byMedCmbIso3, 'pTmed2', 'Medium', 'pT',6)
CompareEfficiencies(byTightCmbIso3, 'pTtight2', 'Tight', 'pT',3)
CompareEfficiencies(byNone,'pTnone2', 'No', 'pT',6)
