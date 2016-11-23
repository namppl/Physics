#include <vector>
#include <iostream>
#include "Physics/Correction/interface/MuonMomentumCorrection.h"

void MuonMomentumCorrection( std::vector<NtupleMuon>& muons, rochcor2016* rmcor, const bool isData ) {

	TLorentzVector mom;
	float qter1 = 1.0;

	for( auto& mu : muons ) {

		// mu.miniIso *= mu.pt;
		float pt = mu.pt, eta = mu.eta, phi = mu.phi, mass = muon_mass;
		mom.SetPtEtaPhiM(pt, eta, phi, mass);
		qter1 = 1.0;

		if(isData) rmcor->momcor_data(mom, mu.charge, 0, qter1);
        else if( (mu.isGlobalMuon||mu.isTrackerMuon) && mu.pt > 9 && mu.nTrackerLayers>0 && mu.nTrackerLayers<20 ) {
			//cout << mu.isGlobalMuon << " : " << mu.isTrackerMuon << endl;
			//cout << mu.pt << " : " << mu.eta << " : " << mu.nTrackerLayers << endl;
			rmcor->momcor_mc(mom, mu.charge, mu.nTrackerLayers, qter1);
		}

		mu.pt = mom.Pt();
		//mu.miniIso /= mu.pt;

	}
}
