'''
Usage: python plotEff.py RootFile1.root RootFile2.root RootFile3.root RootFile4.root RootFile5.root

Script to make Tau Id Efficiency Plots.
Requires Flat (non-vector) Trees

Author: L. Dodd, UW Madison

'''

from subprocess import Popen
from sys import argv, exit, stdout, stderr

import ROOT
from array import array
import plot_def


# So things don't look like crap.
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

######## File #########
if len(argv) < 6:
   print 'Usage:python plotEff.py RootFile1.root RootFile2.root RootFile3.root RootFile4.root RootFile5.root'
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
#This is defined in plot_def
#saveWhere='~/www/Research/'



#####################################
#Get Effi NTUPLE                 #
#####################################
byLooseCmbIso3_03 = ntuple_file.Get("byLooseCombinedIsolationDBSumPtCorr3Hits/Ntuple")
byLooseCmbIso3_05 = ntuple_file2.Get("byLooseCombinedIsolationDBSumPtCorr3Hits/Ntuple")
byLooseCmbIso3_08 = ntuple_file3.Get("byLooseCombinedIsolationDBSumPtCorr3Hits/Ntuple")
byLooseCmbIso3_1 = ntuple_file4.Get("byLooseCombinedIsolationDBSumPtCorr3Hits/Ntuple")
byLooseCmbIso3_5 = ntuple_file5.Get("byLooseCombinedIsolationDBSumPtCorr3Hits/Ntuple")
#####################################
byLooseCmbIso_03 = ntuple_file.Get("byLooseCombinedIsolationDBSumPtCorr/Ntuple")
byLooseCmbIso_05 = ntuple_file2.Get("byLooseCombinedIsolationDBSumPtCorr/Ntuple")
byLooseCmbIso_08 = ntuple_file3.Get("byLooseCombinedIsolationDBSumPtCorr/Ntuple")
byLooseCmbIso_1 = ntuple_file4.Get("byLooseCombinedIsolationDBSumPtCorr/Ntuple")
byLooseCmbIso_5 = ntuple_file5.Get("byLooseCombinedIsolationDBSumPtCorr/Ntuple")
#####################################
byMediumCmbIso3_03 = ntuple_file.Get("byMediumCombinedIsolationDBSumPtCorr3Hits/Ntuple")
byMediumCmbIso3_05 = ntuple_file2.Get("byMediumCombinedIsolationDBSumPtCorr3Hits/Ntuple")
byMediumCmbIso3_08 = ntuple_file3.Get("byMediumCombinedIsolationDBSumPtCorr3Hits/Ntuple")
byMediumCmbIso3_1 = ntuple_file4.Get("byMediumCombinedIsolationDBSumPtCorr3Hits/Ntuple")
byMediumCmbIso3_5 = ntuple_file5.Get("byMediumCombinedIsolationDBSumPtCorr3Hits/Ntuple")
#####################################
byTightCmbIso3_03 = ntuple_file.Get("byTightCombinedIsolationDBSumPtCorr3Hits/Ntuple")
byTightCmbIso3_05 = ntuple_file2.Get("byTightCombinedIsolationDBSumPtCorr3Hits/Ntuple")
byTightCmbIso3_08 = ntuple_file3.Get("byTightCombinedIsolationDBSumPtCorr3Hits/Ntuple")
byTightCmbIso3_1 = ntuple_file4.Get("byTightCombinedIsolationDBSumPtCorr3Hits/Ntuple")
byTightCmbIso3_5 = ntuple_file5.Get("byTightCombinedIsolationDBSumPtCorr3Hits/Ntuple")

################################################################################
# Efficiency for a 20 GeV cut on tau Pt 
################################################################################

bins=[0,20, 30, 40, 50, 60, 70, 80, 100, 120, 160,260.260]#change in defs
######################## LOOSE 3 Hits ##########################################
plot_def.compare_FR(
	            byLooseCmbIso3_03, 'Loose dBIso 0.2', 
	            byLooseCmbIso3_05, 'Loose dBIso 0.3', 
	            byLooseCmbIso3_08, 'Loose dBIso 0.35', 
		    byLooseCmbIso3_1,'Loose dBIso 0.4',
		    byLooseCmbIso3_5,'Loose dBIso 0.5',
		    'tauPt', 20, [36, 0, 360],#variable, ptcut, binning
                    'tauPtLooseCmbIso3FR',#filename
                    "byLooseCombIsoDBSumPtCorr3Hits",#title
                    "pf Tau p_{T} (GeV)",#xaxis
                    "Fake Rate" #yaxis             
)
######################## LOOSE #################################################
plot_def.compare_FR(
	            byLooseCmbIso_03, 'Loose dBIso 0.2', 
	            byLooseCmbIso_05, 'Loose dBIso 0.3', 
	            byLooseCmbIso_08, 'Loose dBIso 0.35', 
		    byLooseCmbIso_1,'Loose dBIso 0.4',
		    byLooseCmbIso_5,'Loose dBIso 0.5',
		    'tauPt', 20, [36, 0, 360],#variable, ptcut, binning
                    'tauPtLooseCmbIsoFR',#filename
                    "byLooseCombIsoDBSumPtCorr",#title
                    "pf Tau p_{T} (GeV)",#xaxis
                    "Fake Rate" #yaxis             
)

######################## MEDIUM #################################################
plot_def.compare_FR(
	            byMediumCmbIso3_03, 'Medium dBIso 0.2', 
	            byMediumCmbIso3_05, 'Medium dBIso 0.3', 
	            byMediumCmbIso3_08, 'Medium dBIso 0.35', 
		    byMediumCmbIso3_1,'Medium dBIso 0.4',
		    byMediumCmbIso3_5,'Medium dBIso 0.5',
		    'tauPt', 20, [36, 0, 360],#variable, ptcut, binning
                    'tauPtMediumCmbIso3FR',#filename
                    "byMediumCombIsoDBSumPtCorr3Hits",#title
                    "pf Tau p_{T} (GeV)",#xaxis
                    "Fake Rate" #yaxis             
)

######################## TIGHT #################################################
plot_def.compare_FR(
	            byTightCmbIso3_03, 'Tight dBIso 0.2', 
	            byTightCmbIso3_05, 'Tight dBIso 0.3', 
	            byTightCmbIso3_08, 'Tight dBIso 0.35', 
		    byTightCmbIso3_1,'Tight dBIso 0.4',
		    byTightCmbIso3_5,'Tight dBIso 0.5',
		    'tauPt', 20, [36, 0, 360],#variable, ptcut, binning
                    'tauPtTightCmbIso3FR',#filename
                    "byTightCombIsoDBSumPtCorr3Hits",#title
                    "pf Tau p_{T} (GeV)",#xaxis
                    "Fake Rate" #yaxis             
)
