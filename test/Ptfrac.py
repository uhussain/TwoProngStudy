#Calculate if ChHadronPtFrac and PhotonPtFraction add up.
from ROOT import TFile
from array import array
from ROOT import TH1F, TCanvas
f = TFile("Zprime_AOD_dR02_Test.root")
t = f.Get("Reco/RecoTree")

nWeirdEvents = 0.0
nSigEvents = 0.0
TrkEvents=0.0
TotalEvents=0.0
OneTrkEvent=0.0
HadrEvents=0.0
OneHadrEvent=0.0
TwoPhotonEvt = OnePhotonEvt =NoPhotonEvt=0.0
TwoNeutralEvt = OneNeutralEvt =NoNeutralEvt=0.0
metcut = 0.0
dRPosHadr1NegHadr1 = TH1F("dRPosHadr1NegHadr1","DeltaR between PF Positive and Negatively Charged Hadron",200,0,0.1)
dRPosHadr1NegHadr1.GetXaxis().SetTitle("#Delta R")
dRPosHadr1NegHadr1.GetYaxis().SetTitle("Events")
c1 = TCanvas()
for event in t:
    metcut = (abs(event.pfMET-event.caloMET))/event.pfMET
    if(event.pfMET>300 and event.j1Pt>100 and event.j1Eta<2.4 and event.j1NHF<0.8 and event.j1CHF>0.1 and metcut<0.5):
        TotalEvents+=1
        TotalPt= event.posTrk1Pt+event.negTrk1Pt+event.Pho1_pt+event.Pho2_pt
        Ptfrac = TotalPt/event.j1Pt
        #how many events have j1Pt largely occupied by Tracks and Photons
        if(Ptfrac>0.90):
            nSigEvents+=1
        else:
            nWeirdEvents+=1
        #how many events have atleast 2 chargedHadrons
        if(event.nChHadrPosj1>0 and event.nChHadrNegj1>0):
            HadrEvents+=1
            Hadr12Ptfrac = ((event.ChHadrPos1_pt+event.ChHadrNeg1_pt)/event.j1Pt) 
        #how many events have atleast 1 charged Hadron
        if(event.nChHadrPosj1>0 or event.nChHadrNegj1>0):
            OneHadrEvent+=1
        #how many events have 2 tracks
        if(event.nPosTrkj1>0 and event.nNegTrkj1>0):
            TrkEvents+=1
            Trk12Ptfrac = ((event.posTrk1Pt+event.negTrk1Pt)/event.j1Pt) 
        #how many events have atleast 1 track
        if(event.nTotTrkj1>0):
            OneTrkEvent+=1
        #what happens in events with no tracks    
        if(event.nTotTrkj1==0):
            #print event.posTrk1Pt
            Trk12Pt = event.posTrk1Pt + event.negTrk1Pt
            j1Pt = event.j1Pt
            if(event.nPhoj1>1):
                TwoPhotonEvt+=1
                Pho12Ptfrac = ((event.Pho1_pt+event.Pho2_pt)/event.j1Pt)
            elif(event.nPhoj1>0):
                OnePhotonEvt+=1
                Pho1Ptfrac = event.Pho1_ptfrac
            else:
                NoPhotonEvt+=1 
            if(event.nNeutralj1>1):
                TwoNeutralEvt+=1
                Neutral12Ptfrac = ((event.Neutral1_pt+event.Neutral2_pt)/event.j1Pt)
            elif(event.nNeutralj1>0):
                OneNeutralEvt+=1
                Neutral1Ptfrac = event.Neutral1_ptfrac
            else:
                NoNeutralEvt+=1
        if(event.ChHadrPos1_pt>0 and event.ChHadrNeg1_pt>0):
            dRPosHadr1NegHadr1.Fill(event.dRPosHadr1NegHadr1)

dRPosHadr1NegHadr1.Draw()
c1.SaveAs('dRPosHadr1NegHadr1.pdf','pdf')
print "All tracks and Photons within dR<0.2 of the leading Jet"    
print "Total Events: ", TotalEvents
print "Frac.of events with good Ptfrac in leading charged Trks and/or photons: ", (nSigEvents/TotalEvents)
print ""
print "No of events with at least 1 positive and 1 negative hadron: ",HadrEvents
print "Frac.of events with at least 1 good positive and 1 good negative Hadron: ",(HadrEvents/TotalEvents)
print "In these events, Hadr12Pt/j1Pt: ",Hadr12Ptfrac

print "Frac. of events with atleast 1 good hadron: ",(OneHadrEvent/TotalEvents)
print ""
print "No.of events with atleast 1 positive and 1 negative Track: ",TrkEvents
print "Frac.of events with at least 1 good positive and 1 good negative Track: ",(TrkEvents/TotalEvents)
print "In these events, Trk12Pt/j1Pt: ",Trk12Ptfrac

print "Frac. of events with atleast 1 good track: ",(OneTrkEvent/TotalEvents)
print "Frac. of events with no tracks: ",(1-(OneTrkEvent/TotalEvents))
print "In events with no tracks, Trk12Ptfrac: ",Trk12Pt/j1Pt
print "No.of events with no tracks: ",(1-(OneTrkEvent/TotalEvents))*TotalEvents
print ""
print "Of these events with no tracks, events with 2 Photons: ",TwoPhotonEvt
print "In events with no tracks, Pho12Ptfrac: ",Pho12Ptfrac

print "Of these events with no tracks, events with 1 Photon: ",OnePhotonEvt
print "In events with no tracks, Pho1Ptfrac: ",Pho1Ptfrac
print "Of these no track events, events with 0 Photons: ",NoPhotonEvt

print "Of these events with no tracks, events with 2 Neutrals: ",TwoNeutralEvt
print "In events with no tracks, Pho12Ptfrac: ",Pho12Ptfrac

print "Of these events with no tracks, events with 1 Neutral: ",OneNeutralEvt
print "In events with no tracks, Pho1Ptfrac: ",Pho1Ptfrac
print "Of these no track events, events with 0 Neutrals: ",NoNeutralEvt
