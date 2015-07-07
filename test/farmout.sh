#!/bin/bash

EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME"
fi

farmoutAnalysisJobs $1-WJets \
  --job-generates-output-name \
  --input-files-per-job=50 \
  --infer-cmssw-path \
  --input-file-list=WJets.txt \
  --input-dir=root://cmsxrootd.fnal.gov/ \
  --assume-input-files-exist \
  ./runMINIAOD_FR.py  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'

farmoutAnalysisJobs $1-ggHtautau \
  --job-generates-output-name \
  --input-files-per-job=50 \
  --infer-cmssw-path \
  --input-file-list=HTauTau.txt \
  --input-dir=root://cmsxrootd.fnal.gov/ \
  --assume-input-files-exist \
  ./runMINIAOD_effi.py  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'

farmoutAnalysisJobs $1-QCD \
  --job-generates-output-name \
  --input-files-per-job=50 \
  --infer-cmssw-path \
  --input-file-list=QCD_multijet.txt \
  --input-dir=root://cmsxrootd.fnal.gov/ \
  --assume-input-files-exist \
  ./runMINIAOD_FR_QCD.py  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'
