#ifndef EgammaEnergyCorrection_H
#define EgammaEnergyCorrection_H

#include "Physics/NtupleMaker/interface/NtupleClasses.h"
#include "Physics/Correction/interface/EnergyScaleCorrection_class.h"
#include "Physics/Correction/interface/EgammaEnergyCorrection.h"

void EgammaEnergyCorrection( NtupleEvent* event, EnergyScaleCorrection_class& electron_correction, EnergyScaleCorrection_class& photon_correction, const bool& isData );

#endif
