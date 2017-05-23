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
This is an attempt to study TwoProng nature of the Z' PencilJet using the tau infrastructure. We are interested in TwoProng0Pi0 and TwoProng1Pi0 based on our signal sample.

