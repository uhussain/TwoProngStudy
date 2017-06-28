#Calculate if ChHadronPtFrac and PhotonPtFraction add up.
from ROOT import TFile,TTree,TDirectory,gROOT
from array import array
from ROOT import TH1F, TCanvas
f = TFile("Zprime_AOD_dR01_j1PFCons.root")
t = f.Get("Reco/RecoTree")

nWeirdEvents = 0.0
nSigEvents = 0.0
TrkEvents=0.0
TotalEvents=0.0
OneTrkEvent=0.0
NoTrkEvent=0.0
HadrEvents=0.0
OneHadrEvent=0.0
Hadr12Ptfrac=Trk12Ptfrac=Trk12Pt=Pho12Ptfrac=Pho1Ptfrac=Neutral12Ptfrac=Neutral1Ptfrac=TrkPtfrac=j1Pt=0.0
OnePosHadr=OneNegHadr=OnePosTrk=OneNegTrk=0.0
TwoPhotonEvt = OnePhotonEvt =NoPhotonEvt=0.0
TwoNeutralEvt = OneNeutralEvt =NoNeutralEvt=0.0
ExactlyOneTrk=0.0
metcut = 0.0
#Root file to store all histograms
f_out = TFile("PostAnalysis.root", "recreate")

#general histograms
OnePosTrkPt= TH1F("OnePosTrkPt","Positive Track 1",100,0,1000)
OneNegTrkPt= TH1F("OneNegTrkPt","Negative Track 1",100,0,1000)
OnePosHadrPt= TH1F("OnePosHadrPt","Positive Hadr 1",100,0,1000)
OneNegHadrPt= TH1F("OneNegHadrPt","Negative Hadr 1",100,0,1000)
dRPosTrk1NegTrk1=TH1F("dRPosTrk1NegTrk1","DeltaR between positive and negative Track",50,0,0.1)
dRPosHadr1NegHadr1=TH1F("dRPosHadr1NegHadr1","DeltaR between positive and negative Hadron",50,0,0.1)
#histograms for events having atleast TwoTracks
j1Pt_TwoTracks = TH1F("j1Pt_TwoTracks","Leading Jet P_{T}",100,0,2000)
j1PF1ConsPt_TwoTracks = TH1F("j1PF1ConsPt_TwoTracks","P_{T} of Leading PFCons of the leading Jet",100,0,1000)
j1PF1ConsEt_TwoTracks = TH1F("j1PF1ConsEt_TwoTracks","E_{t} of Leading PFCons of the leading Jet",100,0,1000)
j1PF1ConsPID_TwoTracks = TH1F("j1PF1ConsPID_TwoTracks","PID of Leading PFCons of the leading Jet",100,-300,300)
#second PFcons
j1PF2ConsPt_TwoTracks = TH1F("j1PF2ConsPt_TwoTracks","P_{T} of Second PFCons of the leading Jet",100,0,2000)
j1PF2ConsEt_TwoTracks = TH1F("j1PF2ConsEt_TwoTracks","E_{t} of Second PFCons of the leading Jet",100,0,2000)
j1PF2ConsPID_TwoTracks = TH1F("j1PF2ConsPID_TwoTracks","PID of Second PFCons of the leading Jet",100,-300,300)

Trk12Ptfrac_plot = TH1F("Trk12Ptfrac_plot","Trk12Pt/j1Pt in events with at least 2 tracks",100,0,1.5)
Trk12Ptfrac_plot.GetXaxis().SetTitle("P_{T} fraction")
Trk12Ptfrac_plot.GetYaxis().SetTitle("Events")

#histograms for Exactly OneTrack events
PosTrk1Pt_OneTracks = TH1F("PosTrk1Pt_OneTracks","Leading Positive Track P_{T}",100,0,1000)
NegTrk1Pt_OneTracks = TH1F("NegTrk1Pt_OneTracks","Leading Negative Track P_{T}",100,0,1000)
PosHadr1Pt_OneTracks = TH1F("PosHadr1Pt_OneTracks","Leading Positive Hadron P_{T}",100,0,1000)
NegHadr1Pt_OneTracks = TH1F("NegHadr1Pt_OneTracks","Leading Negative Hadron P_{T}",100,0,1000)
j1Pt_OneTracks = TH1F("j1Pt_OneTracks","Leading Jet P_{T}",100,0,3000)
j1PF1ConsPt_OneTracks = TH1F("j1PF1ConsPt_OneTracks","P_{T} of Leading PFCons of the leading Jet",100,0,1500)
j1PF1ConsEt_OneTracks = TH1F("j1PF1ConsEt_OneTracks","E_{t} of Leading PFCons of the leading Jet",100,0,1500)
j1PF1ConsPID_OneTracks = TH1F("j1PF1ConsPID_OneTracks","PID of Leading PFCons of the leading Jet",100,-300,300)
#second PFcons
j1PF2ConsPt_OneTracks = TH1F("j1PF2ConsPt_OneTracks","P_{T} of Second PFCons of the leading Jet",100,0,2000)
j1PF2ConsEt_OneTracks = TH1F("j1PF2ConsEt_OneTracks","E_{t} of Second PFCons of the leading Jet",100,0,2000)
j1PF2ConsPID_OneTracks = TH1F("j1PF2ConsPID_OneTracks","PID of Second PFCons of the leading Jet",100,-300,300)

#histograms for NoTrack events
PosTrk1Pt_NoTracks = TH1F("PosTrk1Pt_NoTracks","Leading Positive Track P_{T}",100,0,1000)
NegTrk1Pt_NoTracks = TH1F("NegTrk1Pt_NoTracks","Leading Negative Track P_{T}",100,0,1000)
PosHadr1Pt_NoTracks = TH1F("PosHadr1Pt_NoTracks","Leading Positive Hadron P_{T}",100,0,1000)
NegHadr1Pt_NoTracks = TH1F("NegHadr1Pt_NoTracks","Leading Negative Hadron P_{T}",100,0,1000)
j1Pt_NoTracks = TH1F("j1Pt_NoTracks","Leading Jet P_{T}",100,0,3000)
j1PF1ConsPt_NoTracks = TH1F("j1PF1ConsPt_NoTracks","P_{T} of Leading PFCons of the leading Jet",100,0,2000)
j1PF1ConsEt_NoTracks = TH1F("j1PF1ConsEt_NoTracks","E_{t} of Leading PFCons of the leading Jet",100,0,2000)
j1PF1ConsPID_NoTracks = TH1F("j1PF1ConsPID_NoTracks","PID of Leading PFCons of the leading Jet",100,-300,300)
#second PFcons
j1PF2ConsPt_NoTracks = TH1F("j1PF2ConsPt_NoTracks","P_{T} of Second PFCons of the leading Jet",100,0,2000)
j1PF2ConsEt_NoTracks = TH1F("j1PF2ConsEt_NoTracks","E_{t} of Second PFCons of the leading Jet",100,0,2000)
j1PF2ConsPID_NoTracks = TH1F("j1PF2ConsPID_NoTracks","PID of Second PFCons of the leading Jet",100,-300,300)

#c1 = TCanvas()
for event in t:
    metcut = (abs(event.pfMET-event.caloMET))/event.pfMET
    if(event.pfMET>300 and event.j1Pt>100 and event.j1Eta<2.4 and event.j1NHF<0.8 and event.j1CHF>0.1 and metcut<0.5):
        TotalEvents+=1
        TotalPt= event.posTrk1Pt+event.negTrk1Pt+event.Pho1_pt+event.Pho2_pt
        Ptfrac = (TotalPt/event.j1Pt)
        #how many events have j1Pt largely occupied by Tracks and Photons
        if(Ptfrac>0.90):
            nSigEvents+=1
        else:
            nWeirdEvents+=1
        #how many events have atleast 2 chargedHadrons
        if(event.nChHadrPosj1>0 and event.nChHadrNegj1>0):
            HadrEvents+=1
            Hadr12Ptfrac = Hadr12Ptfrac + ((event.ChHadrPos1_pt+event.ChHadrNeg1_pt)/event.j1Pt)
            dRPosHadr1NegHadr1.Fill(event.dRPosHadr1NegHadr1)
        if(event.nChHadrPosj1>0):
            OnePosHadr+=1
            OnePosHadrPt.Fill(event.ChHadrPos1_pt)
        if(event.nChHadrNegj1>0):
            OneNegHadr+=1
            OneNegHadrPt.Fill(event.ChHadrNeg1_pt)
        #how many events have atleast 1 charged Hadron
        if(event.nChHadrPosj1>0 or event.nChHadrNegj1>0):
            OneHadrEvent+=1
        #how many events have 2 tracks
        if(event.nPosTrkj1>0 and event.nNegTrkj1>0):
            TrkEvents+=1
            dRPosTrk1NegTrk1.Fill(event.dRPosTrk1negTrk1)
            Trk12Ptfrac= Trk12Ptfrac + ((event.posTrk1Pt+event.negTrk1Pt)/event.j1Pt)
            Trk12Ptfrac_plot.Fill((event.posTrk1Pt+event.negTrk1Pt)/event.j1Pt)
            j1Pt_TwoTracks.Fill(event.j1Pt)
            j1PF1ConsPt_TwoTracks.Fill(event.j1PFConsPt[0])
            j1PF1ConsEt_TwoTracks.Fill(event.j1PFConsEt[0])
            j1PF1ConsPID_TwoTracks.Fill(event.j1PFConsPID[0])
            #Second PFConstituent of the leading Jet
            j1PF2ConsPt_TwoTracks.Fill(event.j1PFConsPt[1])
            j1PF2ConsEt_TwoTracks.Fill(event.j1PFConsEt[1])
            j1PF2ConsPID_TwoTracks.Fill(event.j1PFConsPID[1])
        #No.of events with 1 positive Trk
        if(event.nPosTrkj1>0):
            OnePosTrk+=1
            OnePosTrkPt.Fill(event.posTrk1Pt)
        if(event.nNegTrkj1>0):
            OneNegTrk+=1
            OneNegTrkPt.Fill(event.negTrk1Pt)
        #how many events have atleast 1 track
        if(event.nTotTrkj1>0):
            OneTrkEvent+=1
        #what happens in events with no tracks    
        if(event.nTotTrkj1==0):
            NoTrkEvent+=1
            #print event.posTrk1Pt
            Trk12Pt = Trk12Pt + (event.posTrk1Pt + event.negTrk1Pt)
            j1Pt = j1Pt + event.j1Pt
            PosTrk1Pt_NoTracks.Fill(event.posTrk1Pt)
            NegTrk1Pt_NoTracks.Fill(event.negTrk1Pt)
            PosHadr1Pt_NoTracks.Fill(event.ChHadrPos1_pt)
            NegHadr1Pt_NoTracks.Fill(event.ChHadrNeg1_pt)
            j1Pt_NoTracks.Fill(event.j1Pt)
            j1PF1ConsPt_NoTracks.Fill(event.j1PFConsPt[0])
            j1PF1ConsEt_NoTracks.Fill(event.j1PFConsEt[0])
            j1PF1ConsPID_NoTracks.Fill(event.j1PFConsPID[0])
            #Second PFCons in NoTrack events
            j1PF2ConsPt_NoTracks.Fill(event.j1PFConsPt[1])
            j1PF2ConsEt_NoTracks.Fill(event.j1PFConsEt[1])
            j1PF2ConsPID_NoTracks.Fill(event.j1PFConsPID[1])
            if(event.nPhoj1>1):
                TwoPhotonEvt+=1
                Pho12Ptfrac = Pho12Ptfrac + ((event.Pho1_pt+event.Pho2_pt)/event.j1Pt)
            elif(event.nPhoj1>0):
                OnePhotonEvt+=1
                Pho1Ptfrac = Pho1Ptfrac + event.Pho1_ptfrac
            else:
                NoPhotonEvt+=1 
            if(event.nNeutralj1>1):
                TwoNeutralEvt+=1
                Neutral12Ptfrac = Neutral12Ptfrac + ((event.Neutral1_pt+event.Neutral2_pt)/event.j1Pt)
            elif(event.nNeutralj1>0):
                OneNeutralEvt+=1
                Neutral1Ptfrac = Neutral1Ptfrac + (event.Neutral1_pt/event.j1Pt)
            else:
                NoNeutralEvt+=1
        if(event.nTotTrkj1==1):
            ExactlyOneTrk+=1
            PosTrk1Pt_OneTracks.Fill(event.posTrk1Pt)
            NegTrk1Pt_OneTracks.Fill(event.negTrk1Pt)
            PosHadr1Pt_OneTracks.Fill(event.ChHadrPos1_pt)
            NegHadr1Pt_OneTracks.Fill(event.ChHadrNeg1_pt)
            j1Pt_OneTracks.Fill(event.j1Pt)
            #First PFCons
            j1PF1ConsPt_OneTracks.Fill(event.j1PFConsPt[0])
            j1PF1ConsEt_OneTracks.Fill(event.j1PFConsEt[0])
            j1PF1ConsPID_OneTracks.Fill(event.j1PFConsPID[0])
            #second PFCons
            j1PF2ConsPt_OneTracks.Fill(event.j1PFConsPt[1])
            j1PF2ConsEt_OneTracks.Fill(event.j1PFConsEt[1])
            j1PF2ConsPID_OneTracks.Fill(event.j1PFConsPID[1])
            #Pho1_pt.Fill(event.Pho1_pt)
            if(event.nPosTrkj1>0):
                TrkPtfrac= TrkPtfrac + (event.posTrk1Pt/event.j1Pt)
            if(event.nNegTrkj1>0):
                TrkPtfrac= TrkPtfrac + (event.negTrk1Pt/event.j1Pt)
f_out.Write()
f_out.Close()
#Trk12Ptfrac_plot.Draw()
#c1.SaveAs('Trk12Ptfrac_plot.pdf','pdf')
print "All tracks and Photons within dR<0.1 of the leading Jet"    
print "Total Events: ", TotalEvents
#print "Frac.of events with good Ptfrac in leading charged Trks and/or photons: ", (nSigEvents/TotalEvents)
#print ""
print "No of events with at least 1 positive and 1 negative hadron: ",HadrEvents
print "Frac. of events with 1 positive hadron:",(OnePosHadr/TotalEvents)
print "Frac. of events with 1 negative hadron:",(OneNegHadr/TotalEvents)
print "Frac.of events with at least 1 good positive and 1 good negative Hadron: ",(HadrEvents/TotalEvents)
print "In these events, average Hadr12Pt/j1Pt: ",(Hadr12Ptfrac/HadrEvents)

print "Frac. of events with atleast 1 good hadron: ",(OneHadrEvent/TotalEvents)
print ""
print "No.of events with atleast 1 positive and 1 negative Track: ",TrkEvents
print "Frac. of events with 1 positive Track:",(OnePosTrk/TotalEvents)
print "Frac. of events with 1 negative Track:",(OneNegTrk/TotalEvents)
print "Frac.of events with at least 1 good positive and 1 good negative Track: ",(TrkEvents/TotalEvents)
print "In these events, average Trk12Pt/j1Pt: ",(Trk12Ptfrac/TrkEvents)
print ""
print "Fraction of events with exactly one Track:",(ExactlyOneTrk/TotalEvents)
print "In these events, TrkPt/j1Pt: ",(TrkPtfrac/ExactlyOneTrk)

print "Frac. of events with atleast 1 good track: ",(OneTrkEvent/TotalEvents)
print "Frac. of events with no tracks: ",(1-(OneTrkEvent/TotalEvents))
print "In events with no tracks, Trk12Pt: ",(Trk12Pt/NoTrkEvent)
print "In events with no tracks, average j1Pt: ",(j1Pt/NoTrkEvent)
print "No.of events with no tracks: ",NoTrkEvent
print ""
print "Of these events with no tracks, events with 2 Photons: ",TwoPhotonEvt
print "In events with no tracks, Pho12Ptfrac: ",(Pho12Ptfrac/TwoPhotonEvt)

print "Of these events with no tracks, events with 1 Photon: ",OnePhotonEvt
print "In events with no tracks, Pho1Ptfrac: ",(Pho1Ptfrac/OnePhotonEvt)
print "Of these no track events, events with 0 Photons: ",NoPhotonEvt

print "Of these events with no tracks, events with 2 Neutrals: ",TwoNeutralEvt
print "In events with no tracks, Neutral12Ptfrac: ",(Neutral12Ptfrac/TwoNeutralEvt)

print "Of these events with no tracks, events with 1 Neutral: ",OneNeutralEvt
print "In events with no tracks, Neutral1Ptfrac: ",(Neutral1Ptfrac/OneNeutralEvt)
print "Of these no track events, events with 0 Neutrals: ",NoNeutralEvt
