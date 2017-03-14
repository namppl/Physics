#include <vector>
#include <iostream>
#include "Physics/Correction/interface/MuonMomentumCorrection.h"
#include "TH1D.h"

void MuonMomentumCorrection( std::vector<NtupleMuon>& muons, RoccoR& rc, const bool& isData, const bool& doSystematics ) {

	// TLorentzVector mom;

	// float qter1 = 1.0;
	double scaleFactor = 1.0;
	double u1 = gRandom->Rndm(), u2 = gRandom->Rndm();

	for( auto& mu : muons ) {

		// mu.miniIso *= mu.pt;
		double charge = mu.charge, pt = mu.pt, eta = mu.eta, phi = mu.phi;
		int nl = mu.nTrackerLayers; //, mass = muon_mass;
		// std::cout << charge << " : " << pt << " : " << eta << " : " << phi << " : " << nl << std::endl;
		// mom.SetPtEtaPhiM(pt, eta, phi, mass);
		// qter1 = 1.0;

		if(isData) scaleFactor = rc.kScaleDT(charge, pt, eta, phi, 0, 0);
        else if( (mu.isGlobalMuon||mu.isTrackerMuon) && mu.pt > 5 && nl>0 && nl<20 ) {
		// else {
        	if(!doSystematics) scaleFactor = rc.kScaleAndSmearMC(charge, pt, eta, phi, nl, u1, u2, 0, 0);
        	else {
        		TH1D* deviations = new TH1D("deviations", "", 1000, 0., 2.);
        		for(auto i=0; i<100; ++i) deviations->Fill( rc.kScaleAndSmearMC(charge, pt, eta, phi, nl, u1, u2, 1, i) );
        		scaleFactor = 1 + deviations->GetRMS() / deviations->GetMean();
        		delete deviations;
        		// double CorDm, FitDm, deviation, maxDeviation = 1.;
        		// for(auto i=0; i<5; ++i) {
        			// CorDm = rc.kScaleAndSmearMC(charge, pt, eta, phi, nl, u1, u2, 4, i);
        			// FitDm = rc.kScaleAndSmearMC(charge, pt, eta, phi, nl, u1, u2, 5, i);
    				// deviation = (fabs(CorDm-1) > fabs(FitDm-1)) ? CorDm : FitDm;
        			// if( fabs(deviation-1) > fabs(maxDeviation-1) ) maxDeviation = deviation;
        		// }
        		// scaleFactor = maxDeviation;
        	}
		}

		mu.pt *= scaleFactor;
		//mu.miniIso /= mu.pt;

	}
}
