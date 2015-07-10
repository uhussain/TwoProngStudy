#!/bin/bash

EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME"
fi
farmoutAnalysisJobs $1-QCD \
  --job-generates-output-name \
  --infer-cmssw-path \
  --input-file-list=QCD_multijet.txt \
  --input-dir=root://cmsxrootd.fnal.gov/ \
  --assume-input-files-exist \
  ./runMINIAOD_FR_QCD.py  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'
