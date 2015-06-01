#!/bin/bash

EXPECTED_ARGS=2
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME ISO_CONE"
fi

echo "Submitting Jobs with isolation cone $2 and job name $1-iso$2"

farmoutAnalysisJobs $1-iso$2-FR \
  --infer-cmssw-path \
  --input-file-list=WJets.txt \
  --input-dir=root://cmsxrootd.fnal.gov/ \
  --assume-input-files-exist \
  ./testFR.py isoConeSize=$2 isoDBFactor=$3\
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'

farmoutAnalysisJobs $1-iso$2-EFF \
  --infer-cmssw-path \
  --input-file-list=HTauTau.txt \
  --input-dir=root://cmsxrootd.fnal.gov/ \
  --assume-input-files-exist \
  ./testEFF.py isoConeSize=$2 isoDBFactor=$3\
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'
