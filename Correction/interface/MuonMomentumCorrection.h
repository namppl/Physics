#ifndef MuonMomentumCorrection_H
#define MuonMomentumCorrection_H

#include "Physics/NtupleMaker/interface/NtupleClasses.h"
#include "Physics/Correction/interface/rochcor2016.h"

void MuonMomentumCorrection( std::vector<NtupleMuon>& muons, rochcor2016* rmcor, const bool isData );

#endif
