```
cd test 
cmsRun runAODSIM_fakeRate.py
cmsRun runAODSIM_effi.py
python plot.py [effi].root label[optional]
python plotFR.py [fakeRate].root label[optional]
```

note: in plot[FR].py change the save directory to your favortie place. e.g. not ~/www/Research 

note: I hope I pattuplized/reran the hps taus correctly- ressembles AN2014-008
Working on MiniAod package!


for isolation cone changes run 
```
cd test
cmsRun runAODSIM_fakeRate_iso.py isoDBCone=[double]
cmsRun runAODSIM_effi.py isoDBCone=[double]
python plot.py [effi].root label[optional]
python plotFR.py [fakeRate].root label[optional]
```


At Wisconsin for farmout- the files submitted are in WJets.txt and HTauTau.txt
```
cd test
source farmout.sh [JOB_NAME] #for standard isolation cones
source farmout_iso.sh [JOB_NAME] #Change Iso cone size in the file before submitting
python plot.py [effi].root label[optional]
python plotFR.py [fakeRate].root label[optional]
```


For MiniAOD- the efficiency script in production  in plugins/MiniAODtester.cc can be run
Gen Level matching needs to be implemented here!!!
```
cd test
cmsRun runMiniAOD.py
```
