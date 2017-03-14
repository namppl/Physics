#ifndef MuonMomentumCorrection_H
#define MuonMomentumCorrection_H

#include "Physics/NtupleMaker/interface/NtupleClasses.h"
#include "Physics/Rochester/interface/RoccoR.h"
// #include "Physics/Correction/interface/rochcor2016.h"

void MuonMomentumCorrection( std::vector<NtupleMuon>& muons, RoccoR& rc, const bool& isData, const bool& doSystematics );

#endif
