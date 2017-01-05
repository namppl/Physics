#ifndef EgammaEnergyCorrection_H
#define EgammaEnergyCorrection_H

#include "Physics/NtupleMaker/interface/NtupleClasses.h"
#include "EgammaAnalysis/ElectronTools/interface/EnergyScaleCorrection_class.h"

void EgammaEnergyCorrection( NtupleEvent* event, EnergyScaleCorrection_class& electron_correction, EnergyScaleCorrection_class& photon_correction, const bool& isData );

#endif
