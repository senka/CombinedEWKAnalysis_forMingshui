#!/bin/bash

limType=$1

python buildWVworkspace.py muon ${limType}
python buildWVworkspace.py electron ${limType}

combineCards.py *boosted_${limType}.txt > wv_semil_combined_${limType}.txt
text2workspace.py -m 126 wv_semil_combined_${limType}.txt -o ATGC_WV_elmu_workspace_${limType}.root -P CombinedEWKAnalysis.CommonTools.HagiwaraAndZeppenfeldTwoDimensionalModel:${limType}Model