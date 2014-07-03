#!/bin/bash

if [ "$CMSSW_VERSION" != "CMSSW_5_3_16_patch1" ];then
    echo "[ERROR] Invalid CMSSW release: only CMSSW_5_3_16_patch1 is admitted" >> /dev/stderr
    exit 1
fi

git clone -b v1-2-3 https://github.com/amarini/QuarkGluonTagger.git

git clone https://github.com/gdimperi/JECGammaJet.git JetMETCorrections/GammaJetFilter/
git clone https://github.com/gdimperi/GBRLikelihoodEGToolsValueMap.git Calibration/EleNewEnergiesProducer

git clone -b legacyCompatibility https://github.com/bendavid/GBRLikelihood.git HiggsAnalysis/GBRLikelihood
git clone -b hggpaperV6 https://github.com/bendavid/GBRLikelihoodEGTools.git HiggsAnalysis/GBRLikelihoodEGTools

scram b -j10
