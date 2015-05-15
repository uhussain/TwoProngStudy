#!/bin/bash

EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME"
fi

farmoutAnalysisJobs $1-WJets \
  --infer-cmssw-path \
  --input-file-list=WJets.txt \
  --input-dir=root://cmsxrootd.fnal.gov/ \
  --assume-input-files-exist \
  ./runAODSIM_fakeRate_iso.py isoDBCone=0.3 \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'

farmoutAnalysisJobs $1-ggHtautau \
  --infer-cmssw-path \
  --input-file-list=HTauTau.txt \
  --input-dir=root://cmsxrootd.fnal.gov/ \
  --assume-input-files-exist \
  ./runAODSIM_effi_iso.py isoDBCone=0.3 \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'
