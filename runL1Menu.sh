#!/bin/bash


make; 
#./bin/doPrescalingL1.exe /eos/cms/store/user/mguilbau/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_2018/MBHydjet_HI_L1_HLT_20181019/181019_150623/0000/openHLT_Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_2018.root inputs/L1prescalesV4_0_0.txt inputs/L1Menu_CollisionsHeavyIons2018_v4.xml >& outL1MenuratesAll.log
./bin/doPrescalingL1.exe /eos/cms/store/user/mguilbau/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_2018/MBHydjet_HI_L1_HLT_20181019/181019_150623/0000/openHLT_Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_2018.root inputs/L1prescalesV4_0_0_noHF2Not.txt inputs/L1Menu_CollisionsHeavyIons2018_v4.xml >& outL1MenuratesNoMBHF2Not.log
