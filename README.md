Study to determine whether a track based analysis is possible for Z' PencilJet similar to taus.
Example installation:
```
cmsrel CMSSW_8_0_26_patch1
cd CMSSW_8_0_26_patch1/src
cmsenv
git cms-init
git cms-addpkg RecoTauTag/RecoTau
git clone https://github.com/uhussain/TwoProngStudy.git RecoTauTag/tauAnalysis
scram b -j 8
```
Generator level and Reco study with MiniAOD signal samples
```
cd CMSSW_8_0_26_patch1/src/RecoTauTag/RecoTau/TwoProngStudy/test
cmsRun runMINIAODtwoprong.py
```
For closer look at Tracks through reco::Track collection in AOD signal samples
```
cd CMSSW_8_0_26_patch1/src/RecoTauTag/RecoTau/TwoProngStudy/test
cmsRun runAODtwoprong.py
```
