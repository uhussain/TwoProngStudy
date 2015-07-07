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
   saveWhere='~/myAnalysis/CMSSW_7_4_0/src/RecoTauTag/tauAnalysis/outputs/'+argv[2]+'_'
else:
   saveWhere='~/myAnalysis/CMSSW_7_4_0/src/RecoTauTag/tauAnalysis/outputs/'



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

################################################################################
# Efficiency for a 20 GeV cut on tau Pt 
################################################################################

compare_efficiencies(LooseIso, "LooseIsolation",VLooseIso,"VLooseIso",'tauPt', 20, [20, 0, 120],#variable, ptcut, binning
                    'tauPt',#filename
                    "Tau Efficiency",#title
                    "pf Tau p_{T} (GeV)",#xaxis
                    "efficiency" #yaxis             
)

compare_efficiencies(byLooseCmbIso, 'byLooseCombIsoDBSumPtCorr', byVLooseCmbIso,'byVLooseCombIsoDBSumPtCorr','tauPt', 20, [20, 0, 120],#variable, ptcut, binning
                    'tauPtcmb',#filename
                    "Tau Efficiency",#title
                    "pf Tau p_{T} (GeV)",#xaxis
                    "efficiency" #yaxis             
)

compare_efficiencies(byLooseCmbIso, 'byLooseCombIsoDBSumPtCorr', LooseIso,'LooseIsolation','tauPt', 20, [20, 0, 120],#variable, ptcut, binning
                    'tauPtloose',#filename
                    "Tau Efficiency",#title
                    "pf Tau p_{T} (GeV)",#xaxis
                    "efficiency" #yaxis             
)

compare_efficiencies(byVLooseCmbIso, 'byVLooseCombIsoDBSumPtCorr', VLooseIso,'VLooseIsolation','tauPt', 20, [20, 0, 120],#variable, ptcut, binning
                    'tauPtvloose',#filename
                    "Tau Efficiency",#title
                    "pf Tau p_{T} (GeV)",#xaxis
                    "efficiency" #yaxis             
)

#####################################################
  
#--------------------------------------------------------------------------------
# produce PAT-tuple
process.load("PhysicsTools/PatAlgos/patSequences_cff")
# configure pat::Jet production
# (enable L2L3Residual corrections in case running on Data)
jetCorrections = ( 'L1FastJet', 'L2Relative', 'L3Absolute')
from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection
switchJetCollection(
    process,
    jetSource = cms.InputTag('ak4PFJetsCHS'),
    jetCorrections = ( 'AK4PFchs', jetCorrections, "" ),
    outputModules = []
)


#process.patJets.addTagInfos = cms.bool(True)
#process.patJets.addBTagInfo = cms.bool(True)


#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff") #loading the configuration
# switch to HPS PFTaus (and disable all "cleaning" cuts)
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

#process.ak4PFJetsLegacyHPSPiZeros.stripPhiAssociationDistance = cms.double(0.9)

# switch on PAT trigger                                                                                                                      
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger( process )

process.makePatTrigger = cms.Sequence(process.patTrigger*process.patTriggerEvent)

#--------------------------------------------------------------------------------
# select "good" reconstructed vertices
#process.load("TauAnalysis/RecoTools/recoVertexSelection_cff")
process.selectedPrimaryVertexQuality = cms.EDFilter("VertexSelector",
    src = cms.InputTag('offlinePrimaryVertices'),
    cut = cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0"), # CV: cut >= 4 if using 'offlinePrimaryVertices',
                                                                         #         >= 7 if using 'offlinePrimaryVerticesWithBS' as input
    filter = cms.bool(False)                                          
)

process.selectedPrimaryVertexPosition = cms.EDFilter("VertexSelector",
    src = cms.InputTag('selectedPrimaryVertexQuality'),
    cut = cms.string("abs(z) < 24 & abs(position.Rho) < 2."),
    filter = cms.bool(False)                                           
)

process.selectedPrimaryVertexHighestPtTrackSum = cms.EDFilter("PATSingleVertexSelector",
    mode = cms.string('firstVertex'),
    vertices = cms.InputTag('selectedPrimaryVertexPosition'),
    filter = cms.bool(False)                                                    
)

process.selectPrimaryVertex = cms.Sequence(
    process.selectedPrimaryVertexQuality
   * process.selectedPrimaryVertexPosition
   * process.selectedPrimaryVertexHighestPtTrackSum
)


