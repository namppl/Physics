#include "Physics/Correction/interface/EgammaEnergyCorrection.h"

void EgammaEnergyCorrection( NtupleEvent* event, EnergyScaleCorrection_class& electron_correction, EnergyScaleCorrection_class& photon_correction, const bool& isData ) {

	float sigma = 0;
	float factor = 0;
	TRandom3 myRandom_(13);

	if( isData ) {
		electron_correction.doScale = true;
		photon_correction.doScale = true;
		electron_correction.doSmearings = false;
		photon_correction.doSmearings = false;

		for( auto& el : event->electrons ) {
			if( el.pt > 8 && fabs(el.eta) < 2.4  && (fabs(el.eta) < 1.4442 || fabs(el.eta) > 1.566) ) {
				//cout << "< Electron correction >" << endl;
				//cout << el.pt << " : ";
				factor = electron_correction.ScaleCorrection(event->run, fabs(el.eta)<1.479, el.r9, el.eta, el.pt);	
				el.pt = el.pt * factor; 
				//cout << el.pt << endl;
			}
		}
		for( auto& pho : event->photons ) {
			if( pho.pt > 30 && fabs(pho.eta) < 2.5  && (fabs(pho.eta) < 1.4442 || fabs(pho.eta) > 1.566) ) {
				//cout << "< Photon correction >" << endl;
				//cout << pho.pt << " : ";
				factor = photon_correction.ScaleCorrection(event->run, fabs(pho.eta)<1.479, pho.r9, pho.eta, pho.pt);	
				pho.pt = pho.pt * factor; 
				//cout << pho.pt << endl;
			}
		}
	}
	else {
		electron_correction.doScale = false;
		photon_correction.doScale = false;
		electron_correction.doSmearings = true;
		photon_correction.doSmearings = true;

		for( auto& el : event->electrons ) {
			if( el.pt > 8 && fabs(el.eta) < 2.4  && (fabs(el.eta) < 1.4442 || fabs(el.eta) > 1.566) ) {
				//cout << "< Electron correction >" << endl;
				//cout << el.pt << " : ";
				sigma = electron_correction.getSmearingSigma(event->run, fabs(el.eta)<1.479, el.r9, el.eta, el.pt, 0., 0.);	
				factor = myRandom_.Gaus( 1., sigma );
				el.pt = el.pt * factor; 
				//cout << el.pt << endl;
			}
		}
		for( auto& pho : event->photons ) {
			if( pho.pt > 30 && fabs(pho.eta) < 2.5  && (fabs(pho.eta) < 1.4442 || fabs(pho.eta) > 1.566) ) {
				//cout << "< Photon correction >" << endl;
				//cout << pho.pt << " : ";
				sigma = photon_correction.getSmearingSigma(event->run, fabs(pho.eta)<1.479, pho.r9, pho.eta, pho.pt, 0., 0.);	
				factor = myRandom_.Gaus( 1., sigma );
				pho.pt = pho.pt * factor; 
				//cout << pho.pt << endl;
			}
		}
	}
}
