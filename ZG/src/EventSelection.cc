#include "Physics/ZG/interface/EventSelection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <assert.h>

EventSelector::EventSelector( const bool& _verbose): verbose{_verbose}, elPtCut{65.}, elExtraPtCut{10.}, elVetoPtCut{10.}, phoPtCut_ee{65.}, muPtCut{26.}, muExtraPtCut{10.}, muVetoPtCut{10.}, phoPtCut_mm{40.}, lowerMassCut{50}, upperMassCut{130}, phoPtOverMZG{40./150.}, multiPho{false}, norm{-1.0}
{}

EventSelector::EventSelector( const bool& _verbose, const double& _elPtCut, const double& _elExtraPtCut, const double& _elVetoPtCut, const double& _phoPtCut_ee, const double& _muPtCut, const double& _muExtraPtCut, const double& _muVetoPtCut, const double& _phoPtCut_mm, const double& _lowerMassCut, const double& _upperMassCut, const double& _phoPtOverMZG, const bool& _multiPho, const double& _norm ):
verbose{_verbose}, elPtCut{_elPtCut}, elExtraPtCut{_elExtraPtCut}, elVetoPtCut{_elVetoPtCut}, phoPtCut_ee{_phoPtCut_ee}, muPtCut{_muPtCut}, muExtraPtCut{_muExtraPtCut}, muVetoPtCut{_muVetoPtCut}, phoPtCut_mm{_phoPtCut_mm}, lowerMassCut{_lowerMassCut}, upperMassCut{_upperMassCut}, phoPtOverMZG{_phoPtOverMZG}, multiPho{_multiPho}, norm{_norm}
{
	std::cout << std::endl;
	std::cout << ":: selection rule ::" << std::endl;
	std::cout << "lower mass threshold: " << lowerMassCut << std::endl; 
	std::cout << "upper mass threshold: " << upperMassCut << std::endl; 
	std::cout << "photon pt over ZG mass: " << phoPtOverMZG << std::endl;
	std::cout << "multiple photon: "<< multiPho << std::endl;
	ee = new ElectronSelector(elPtCut, 2.5, elExtraPtCut, elVetoPtCut);
	mm = new MuonSelector(muPtCut, 2.4, muExtraPtCut, muVetoPtCut);
	gamma_ee = new PhotonSelector(phoPtCut_ee, 2.5);
	gamma_mm = new PhotonSelector(phoPtCut_mm, 2.5);
	scaleFactors = new ScaleFactor();
	std::cout << std::endl;
}

EventSelector::~EventSelector() {
	delete mm;
	delete ee;
	delete gamma_ee;
	delete gamma_mm;
}

void EventSelector::setEvent(NtupleEvent* _event) {
	event = _event;
}

void EventSelector::selectICHEP16() {
	
	run = event->run;
	lumi = event->lumi;
	eventN = event->event;
	eventNum = event->event;

	// std::cout << "event: " << event->event << std::endl;
	// std::cout << "eventN: " << eventN << std::endl;
	// std::cout << "eventNum: " << eventNum << std::endl;

	weight = event->weight;	
	nPU = event->nPU;
	nVert = event->nVertices;

	muon = false;
	electron = false;

	ee->select(event->electrons);
	mm->selectICHEP16(event->muons);

	if( ee->passSelection() && mm->nMuons()==0 ) {	

		if(multiPho) gamma_ee->select(event->photons, ee->at(0), ee->at(1));
		else gamma_ee->select(event->photons);

		if(gamma_ee->passSelection()) selectDielGamma();

	}
	else if( mm->passSelection() && ee->nElectrons()==0 ) {

		if(multiPho) gamma_mm->select(event->photons, mm->at(0), mm->at(1));
		else gamma_mm->select(event->photons);

		if(gamma_mm->passSelection()) selectDimuGamma();

	}

}

void EventSelector::fillEventVariables() {
	run = event->run;
	lumi = event->lumi;
	eventN = event->event;
	eventNum = event->event;
	// std::cout << "event: " << event->event << std::endl;
    // std::cout << "eventN: " << eventN << std::endl;
    // std::cout << "eventNum: " << eventNum << std::endl;
	nPU = event->nPU;
	aMCNLO = event->weight;
	weight = event->weight;
	if( norm > 0. ) weight *= PUReweight(nPU) * norm;	// norm = 0. for data
	nVert = event->nVertices;

	Flag_badMuons = isTriggered("Flag_badMuons");
	Flag_duplicatedMuons = isTriggered("Flag_duplicatedMuons");
	Flag_noBadMuons = isTriggered("Flag_noBadMuons");

	Flag_dupECALClusters = isTriggered("Flag_dupECALClusters");
	Flag_noHitsNotReplaced = isTriggered("Flag_noHitsNotReplaced");

}

void EventSelector::selectMoriond17(const int& option) {

	muon = false;
	electron = false;
	
	switch(option) {
		case 0:
			ee->select(event->electrons);
			break;
		case 1:
			ee->selectEE(event->electrons);
			break;
		case 5:
			ee->select(event->electrons);
			break;
		case 6:
			ee->selectEE(event->electrons);
			break;
		case 7:
			ee->selectEE(event->electrons);
			break;
		case 8:
			ee->select(event->electrons);
			break;
	}
	
	if(option < 8) mm->selectMoriond17(event->muons);
	else mm->selectICHEP16(event->muons); 

	if( ee->passSelection() && mm->nMuons()==0 ) {
		switch(option) {
			case 0:
				gamma_ee->select(event->photons);
				break;
			case 1:
				gamma_ee->selectEE(event->photons);
				break;
			case 5:
				gamma_ee->selectMVA(event->photons);
				break;
			case 6:
				gamma_ee->selectMVAEE(event->photons);
				break;
			case 7:
				gamma_ee->selectMVAEE(event->photons, ee->at(0), ee->at(1));
				break;
			case 8:
				gamma_ee->selectMVA(event->photons, ee->at(0), ee->at(1));
				break;
		}
		if( gamma_ee->passSelection() ) {
			selectDielGamma();
		}
	}
	else if( mm->passSelection() && ee->nElectrons()==0 ) {
		switch(option) {
			case 0:
				gamma_mm->select(event->photons);
				break;
			case 1:
				gamma_mm->select(event->photons);
				break;
			case 5:
				gamma_mm->selectMVA(event->photons);
				break;
			case 6:
				gamma_mm->selectMVA(event->photons);
				break;
			case 7:
				gamma_mm->selectMVA(event->photons, mm->at(0), mm->at(1));
				break;
			case 8:
				gamma_mm->selectMVA(event->photons, mm->at(0), mm->at(1));
				break;
		}
		if( gamma_mm->passSelection() ) {
			selectDimuGamma();
		}
	}

}

void EventSelector::selectDielGamma() {

	assert( ee->passSelection() && mm->nMuons()==0 && gamma_ee->passSelection() );	

	const NtupleElectron& el1{ee->at(0)};
	const NtupleElectron& el2{ee->at(1)};
	const NtuplePhoton& pho{gamma_ee->at(0)};

	double z_mass = ( momentum(el1,electron_mass) + momentum(el2,electron_mass) ).M();
	double boss_mass = ( momentum(el1,electron_mass) + momentum(el2,electron_mass) + momentum(pho,0.) ).M();
	double dR1 = deltaR(el1.eta,el1.phi,pho.eta,pho.phi); 
	double dR2 = deltaR(el2.eta,el2.phi,pho.eta,pho.phi);

	if( /*el1.charge*el2.charge < 0 &&*/ dR1 > 0.4 && dR2 > 0.4 && z_mass > lowerMassCut && z_mass < upperMassCut && pho.pt/boss_mass > phoPtOverMZG ) {
		electron = true;
	}
		
}

void EventSelector::selectDimuGamma() {

	assert( mm->passSelection() && ee->nElectrons()==0 && gamma_mm->passSelection() ); 

	const NtupleMuon& mu1{mm->at(0)};
	const NtupleMuon& mu2{mm->at(1)};
	const NtuplePhoton& pho{gamma_mm->at(0)};

	double z_mass = ( momentum(mu1,muon_mass) + momentum(mu2,muon_mass) ).M();
	double boss_mass = ( momentum(mu1,muon_mass) + momentum(mu2,muon_mass) + momentum(pho,0.) ).M();
	double dR1 = deltaR(mu1.eta,mu1.phi,pho.eta,pho.phi); 
	double dR2 = deltaR(mu2.eta,mu2.phi,pho.eta,pho.phi);

	if( mu1.charge*mu2.charge < 0 && dR1 > 0.4 && dR2 > 0.4 && z_mass > lowerMassCut && z_mass < upperMassCut && pho.pt/boss_mass > phoPtOverMZG ) {
		muon = true;
	}

}

bool EventSelector::isTriggered( const TString& triggerName ) {
	for( auto& trigger : event->triggers ) {
	    if( trigger.name.find(triggerName) == 0 ) {
	    	if(verbose) std::cout << "Fired: " << trigger.name << std::endl;
	    	return trigger.isFired;
	    }
	}
    return false;
}

void EventSelector::setGenMass() {
	genWeight = event->weight;
	TLorentzVector mom_system;
	TLorentzVector mom_Z;
	TLorentzVector mom_lept0;
	TLorentzVector mom_lept1;
	TLorentzVector mom_pho;
	int nGenEl = 0;
	int nGenMu = 0;
	int nGenPho = 0;

	genLept.clear();
	genPho = nullptr;

	for( auto& par : event->genparticles ) {
		//if( par.fromHardProcessFinalState ) {
		if( par.isHardProcess ) {
			if( fabs(par.id) == 11 ) {
				genLept.push_back(&par);
				++nGenEl;
			}
			else if( fabs(par.id) == 13 ) {
				genLept.push_back(&par);
				++nGenMu;
			}
			else if( fabs(par.id) == 22 ) {
				genPho = &par;
				++nGenPho;
			}
		} 
		if( par.fromHardProcessFinalState ) {
			if( fabs(par.id) == 11 ) {
				if(par.id>0) mom_lept0.SetPtEtaPhiM(par.pt, par.eta, par.phi, electron_mass);
				else mom_lept1.SetPtEtaPhiM(par.pt, par.eta, par.phi, electron_mass);
			}
			else if( fabs(par.id) == 13 ) {
				if(par.id>0) mom_lept0.SetPtEtaPhiM(par.pt, par.eta, par.phi, muon_mass);
				else mom_lept1.SetPtEtaPhiM(par.pt, par.eta, par.phi, muon_mass);
			}
			else if( fabs(par.id) == 22 ) {
				mom_pho.SetPtEtaPhiM(par.pt, par.eta, par.phi, 0.);
			}
		} 
	}
	assert( (nGenEl==2 && nGenMu==0 && nGenPho==1) || (nGenEl==0 && nGenMu==2 && nGenPho==1) || (nGenEl==0 && nGenMu==0 && nGenPho==1) );
	
	// std::cout << "evaluate kinematics" << std::endl;
	mom_system = mom_lept0 + mom_lept1 + mom_pho;
	mom_Z = mom_lept0 + mom_lept1;
	accepted = false;
	double lept0_pt = mom_lept0.Pt();
	double lept0_eta = mom_lept0.Eta();
	double lept0_phi = mom_lept0.Phi();
	double lept1_pt = mom_lept1.Pt();
	double lept1_eta = mom_lept1.Eta();
	double lept1_phi = mom_lept1.Phi();
	double pho_pt = mom_pho.Pt();
	// double pho_pt = genPho->pt;
	double pho_eta = mom_pho.Eta();
	double pho_phi = mom_pho.Phi();

	if( mom_Z.M() > lowerMassCut && mom_Z.M() < upperMassCut && fabs(lept0_eta) < 2.4 && fabs(lept1_eta) < 2.4 && fabs(pho_eta) < 2.5 && (fabs(pho_eta)<1.4442||fabs(pho_eta)>1.566) ) accepted = true;

	if(nGenEl==2) {
		channel = "ee";
		genMass = mom_system.M();
		if(accepted) {
			double dR1 = deltaR(lept0_eta,lept0_phi,pho_eta,pho_phi);
			double dR2 = deltaR(lept1_eta,lept1_phi,pho_eta,pho_phi);
			accepted = ( ((lept0_pt>elPtCut && lept1_pt>elExtraPtCut)||(lept0_pt>elExtraPtCut && lept1_pt>elPtCut)) && (fabs(lept0_eta)<1.4442||fabs(lept0_eta)>1.566) && (fabs(lept1_eta)<1.4442||fabs(lept1_eta)>1.566) && pho_pt > phoPtCut_ee && dR1 > 0.4 && dR2 > 0.4 && pho_pt/genMass > phoPtOverMZG  );
		}
	}
	else if(nGenMu==2) {
		channel = "mm";
		genMass = mom_system.M();
		if(accepted) {
			double dR1 = deltaR(lept0_eta,lept0_phi,pho_eta,pho_phi);
			double dR2 = deltaR(lept1_eta,lept1_phi,pho_eta,pho_phi);
			accepted = ( ((lept0_pt>muPtCut && lept1_pt>muExtraPtCut)||(lept0_pt>muExtraPtCut && lept1_pt>muPtCut)) && pho_pt > phoPtCut_mm && dR1 > 0.4 && dR2 > 0.4 && pho_pt/genMass > phoPtOverMZG  );
		}
	}
	else {
		channel = "tt";
		genMass = -999;
	}
}

bool EventSelector::matchedAcceptance() {

	TLorentzVector mom_ZG;
	TLorentzVector mom_Z;
	TLorentzVector mom0;
	TLorentzVector mom1;
	TLorentzVector mom_pho;

	double dR0;
	double dR1;
	double min0 = 0.1;
	double min1 = 0.1;

	for(auto& pho:event->photons) {
		dR0 = deltaR(genPho->eta,genPho->phi,pho.eta,pho.phi);
		if( dR0 < 0.1 ) {
			mom_pho.SetPtEtaPhiM(pho.pt, pho.eta, pho.phi, 0.);
			break;
		}
	}

	if(channel=="ee") {
		for(auto& el:event->electrons) {
			dR0 = deltaR(genLept[0]->eta,genLept[0]->phi,el.eta,el.phi);
			dR1 = deltaR(genLept[1]->eta,genLept[1]->phi,el.eta,el.phi);
			if( dR0 < min0 ) {
				min0 = dR0;
				mom0.SetPtEtaPhiM(el.pt, el.eta, el.phi, electron_mass);
			} 
			else if( dR1 < min1 ) {
				min1 = dR1;
				mom1.SetPtEtaPhiM(el.pt, el.eta, el.phi, electron_mass);
			} 
		}
	}
	else if(channel=="mm") {
		for(auto& mu:event->muons) {
			dR0 = deltaR(genLept[0]->eta,genLept[0]->phi,mu.eta,mu.phi);
			dR1 = deltaR(genLept[1]->eta,genLept[1]->phi,mu.eta,mu.phi);
			if( dR0 < min0 ) {
				min0 = dR0;
				mom0.SetPtEtaPhiM(mu.pt, mu.eta, mu.phi, muon_mass);
			} 
			else if( dR1 < min1 ) {
				min1 = dR1;
				mom1.SetPtEtaPhiM(mu.pt, mu.eta, mu.phi, muon_mass);
			} 
		}
	}
	else return false;

	mom_ZG = mom0 + mom1 + mom_pho;
	mom_Z = mom0 + mom1;
	double lept0_pt = mom0.Pt();
	double lept0_eta = mom0.Eta();
	double lept0_phi = mom0.Phi();
	double lept1_pt = mom1.Pt();
	double lept1_eta = mom1.Eta();
	double lept1_phi = mom1.Phi();
	double pho_pt = mom_pho.Pt();
	double pho_eta = mom_pho.Eta();
	double pho_phi = mom_pho.Phi();
	dR0 = deltaR(lept0_eta,lept0_phi,pho_eta,pho_phi);
	dR1 = deltaR(lept1_eta,lept1_phi,pho_eta,pho_phi);

	if( mom_Z.M() > lowerMassCut && mom_Z.M() < upperMassCut && fabs(lept0_eta) < 2.4 && fabs(lept1_eta) < 2.4 && fabs(pho_eta) < 2.5 && (fabs(pho_eta)<1.4442||fabs(pho_eta)>1.566) ) {
		if( channel=="ee" )  return ( ((lept0_pt>elPtCut && lept1_pt>elExtraPtCut)||(lept0_pt>elExtraPtCut && lept1_pt>elPtCut)) && (fabs(lept0_eta)<1.4442||fabs(lept0_eta)>1.566) && (fabs(lept1_eta)<1.4442||fabs(lept1_eta)>1.566) && pho_pt > phoPtCut_ee && dR0 > 0.4 && dR1 > 0.4 && pho_pt/genMass > phoPtOverMZG  );
		else if( channel == "mm" ) return ( ((lept0_pt>muPtCut && lept1_pt>muExtraPtCut)||(lept0_pt>muExtraPtCut && lept1_pt>muPtCut)) && pho_pt > phoPtCut_mm && dR0 > 0.4 && dR1 > 0.4 && pho_pt/genMass > phoPtOverMZG  );
	}

	return false;

}

double PUReweight(const int& nPU) {
    int index = nPU-1;
	double ReReco_Moriond17_65[] = {0.385248, 1.23178, 1.2836, 1.14691, 1.30067, 1.3103, 0.978232, 0.886271, 1.24131, 1.54662, 1.7043, 1.70436, 1.63711, 1.62373, 1.60381, 1.52131, 1.42025, 1.3345, 1.25322, 1.1859, 1.13768, 1.09294, 1.05956, 1.03098, 1.00485, 0.984915, 0.968839, 0.949522, 0.929222, 0.903614, 0.853359, 0.803061, 0.733009, 0.657447, 0.574058, 0.486491, 0.395985, 0.310427, 0.232546, 0.167488, 0.112897, 0.0725097, 0.0448537, 0.0264411, 0.0152595, 0.0084306, 0.00440492, 0.00231261, 0.00118615, 0.000622782, 0.000345365, 0.0002281, 0.0001919, 0.000200888, 0.000262581, 0.000372703, 0.000517395, 0.000743546, 0.00108736, 0.00150363, 0.00238769, 0.00307876, 0.00354461, 0.00369242, 0.00392163, 0.00360387, 0.00314731, 0.00267078, 0.00233614, 0.00198402, 0.00166782, 0.0013905, 0.00115199, 0.000955533, 0.000777455};

    return (index < 75) ? ReReco_Moriond17_65[index] : 1. ;
}

void EventSelector::applyScaleFactors() {
	assert( electron || muon );
	if( electron ) {
		double scaleFactor = scaleFactors->getPhotonScaleFactor(gamma_ee->at(0));
		scaleFactor *= scaleFactors->getElectronScaleFactor(ee->at(0)) * scaleFactors->getElectronScaleFactor(ee->at(1));
		assert(scaleFactor > 0.);
		// std::cout << "photon(ee): " << scaleFactor << std::endl;
		weight *= scaleFactor;
	} 
	else if( muon ) {
		double scaleFactor = scaleFactors->getPhotonScaleFactor(gamma_mm->at(0));
		assert(scaleFactor > 0.);
		// if(scaleFactor < 0.9) std::cout << "photon(mm): " << scaleFactor << std::endl;
		weight *= scaleFactor;
		scaleFactor = scaleFactors->getMuonScaleFactors(mm);
		// if(scaleFactor < 0.9) std::cout << "muon: " << scaleFactor << std::endl;
		weight *= scaleFactor;
	} 
}

ScaleFactor::ScaleFactor() {
	// ELECTRON SFS
	ELECTRON_SF_ID[0][0] = 0.921569;
	ELECTRON_SF_ID[0][1] = 0.947436;
	ELECTRON_SF_ID[0][2] = 0.988277;
	ELECTRON_SF_ID[0][3] = 0.997753;
	ELECTRON_SF_ID[0][4] = 1.00333;
	ELECTRON_SF_ID[0][5] = 1.00879;
	ELECTRON_SF_ID[0][6] = 1.03122;
	ELECTRON_SF_ID[0][7] = 1.01438;
	ELECTRON_SF_ID[1][0] = 1.24183;
	ELECTRON_SF_ID[1][1] = 0.978548;
	ELECTRON_SF_ID[1][2] = 1;
	ELECTRON_SF_ID[1][3] = 0.992941;
	ELECTRON_SF_ID[1][4] = 0.978923;
	ELECTRON_SF_ID[1][5] = 1.05313;
	ELECTRON_SF_ID[1][6] = 1.05942;
	ELECTRON_SF_ID[1][7] = 1.06619;
	ELECTRON_SF_ID[2][0] = 0.987342;
	ELECTRON_SF_ID[2][1] = 0.967933;
	ELECTRON_SF_ID[2][2] = 0.978022;
	ELECTRON_SF_ID[2][3] = 0.97983;
	ELECTRON_SF_ID[2][4] = 0.98;
	ELECTRON_SF_ID[2][5] = 0.989496;
	ELECTRON_SF_ID[2][6] = 0.98218;
	ELECTRON_SF_ID[2][7] = 0.975;
	ELECTRON_SF_ID[3][0] = 1.00417;
	ELECTRON_SF_ID[3][1] = 0.970484;
	ELECTRON_SF_ID[3][2] = 0.976999;
	ELECTRON_SF_ID[3][3] = 0.981992;
	ELECTRON_SF_ID[3][4] = 0.982143;
	ELECTRON_SF_ID[3][5] = 0.990546;
	ELECTRON_SF_ID[3][6] = 0.993717;
	ELECTRON_SF_ID[3][7] = 1.00947;
	ELECTRON_SF_ID[4][0] = 1.14989;
	ELECTRON_SF_ID[4][1] = 0.902479;
	ELECTRON_SF_ID[4][2] = 0.984148;
	ELECTRON_SF_ID[4][3] = 0.98014;
	ELECTRON_SF_ID[4][4] = 0.973099;
	ELECTRON_SF_ID[4][5] = 0.962834;
	ELECTRON_SF_ID[4][6] = 0.9718;
	ELECTRON_SF_ID[4][7] = 1.01222;
	ELECTRON_SF_ID[5][0] = 0.944268;
	ELECTRON_SF_ID[5][1] = 0.926209;
	ELECTRON_SF_ID[5][2] = 0.974389;
	ELECTRON_SF_ID[5][3] = 0.98771;
	ELECTRON_SF_ID[5][4] = 0.992274;
	ELECTRON_SF_ID[5][5] = 1.00444;
	ELECTRON_SF_ID[5][6] = 1.01774;
	ELECTRON_SF_ID[5][7] = 1.00887;

	ELECTRON_SF_RECO[0] = 1.3176;
	ELECTRON_SF_RECO[1] = 1.11378;
	ELECTRON_SF_RECO[2] = 1.02463;
	ELECTRON_SF_RECO[3] = 1.01364;
	ELECTRON_SF_RECO[4] = 1.00728;
	ELECTRON_SF_RECO[5] = 0.994819;
	ELECTRON_SF_RECO[6] = 0.994786;
	ELECTRON_SF_RECO[7] = 0.991632;
	ELECTRON_SF_RECO[8] = 0.963129;
	ELECTRON_SF_RECO[9] = 0.989701;
	ELECTRON_SF_RECO[10] = 0.985656;
	ELECTRON_SF_RECO[11] = 0.981595;
	ELECTRON_SF_RECO[12] = 0.984678;
	ELECTRON_SF_RECO[13] = 0.981614;
	ELECTRON_SF_RECO[14] = 0.980433;
	ELECTRON_SF_RECO[15] = 0.984552;
	ELECTRON_SF_RECO[16] = 0.988764;
	ELECTRON_SF_RECO[17] = 0.987743;
	ELECTRON_SF_RECO[18] = 0.987743;
	ELECTRON_SF_RECO[19] = 0.987743;
	ELECTRON_SF_RECO[20] = 0.98768;
	ELECTRON_SF_RECO[21] = 0.967598;
	ELECTRON_SF_RECO[22] = 0.989627;
	ELECTRON_SF_RECO[23] = 0.992761;
	ELECTRON_SF_RECO[24] = 0.991761;
	ELECTRON_SF_RECO[25] = 0.99794;
	ELECTRON_SF_RECO[26] = 1.00104;
	ELECTRON_SF_RECO[27] = 0.989507;
	ELECTRON_SF_RECO[28] = 0.970519;
	ELECTRON_SF_RECO[29] = 0.906667;

	// ELECTRON SF UNCERTAINTIES
	ELECTRON_SF_UNC_ID[0][0] = 0.0121155;
	ELECTRON_SF_UNC_ID[0][1] = 0.0152835;
	ELECTRON_SF_UNC_ID[0][2] = 0.00750628;
	ELECTRON_SF_UNC_ID[0][3] = 0.00350003;
	ELECTRON_SF_UNC_ID[0][4] = 0.00416162;
	ELECTRON_SF_UNC_ID[0][5] = 0.00836258;
	ELECTRON_SF_UNC_ID[0][6] = 0.0283599;
	ELECTRON_SF_UNC_ID[0][7] = 0.0289233;
	ELECTRON_SF_UNC_ID[1][0] = 0.0266196;
	ELECTRON_SF_UNC_ID[1][1] = 0.0335694;
	ELECTRON_SF_UNC_ID[1][2] = 0.00856674;
	ELECTRON_SF_UNC_ID[1][3] = 0.0101916;
	ELECTRON_SF_UNC_ID[1][4] = 0.0341891;
	ELECTRON_SF_UNC_ID[1][5] = 0.0465936;
	ELECTRON_SF_UNC_ID[1][6] = 0.0688819;
	ELECTRON_SF_UNC_ID[1][7] = 0.100967;
	ELECTRON_SF_UNC_ID[2][0] = 0.00833912;
	ELECTRON_SF_UNC_ID[2][1] = 0.00167066;
	ELECTRON_SF_UNC_ID[2][2] = 0.00279245;
	ELECTRON_SF_UNC_ID[2][3] = 0.00211752;
	ELECTRON_SF_UNC_ID[2][4] = 0.00267665;
	ELECTRON_SF_UNC_ID[2][5] = 0.00277478;
	ELECTRON_SF_UNC_ID[2][6] = 0.0116938;
	ELECTRON_SF_UNC_ID[2][7] = 0.00788104;
	ELECTRON_SF_UNC_ID[3][0] = 0.00833912;
	ELECTRON_SF_UNC_ID[3][1] = 0.00167066;
	ELECTRON_SF_UNC_ID[3][2] = 0.00279245;
	ELECTRON_SF_UNC_ID[3][3] = 0.00211752;
	ELECTRON_SF_UNC_ID[3][4] = 0.00267665;
	ELECTRON_SF_UNC_ID[3][5] = 0.00277478;
	ELECTRON_SF_UNC_ID[3][6] = 0.0116938;
	ELECTRON_SF_UNC_ID[3][7] = 0.00809816;
	ELECTRON_SF_UNC_ID[4][0] = 0.0269229;
	ELECTRON_SF_UNC_ID[4][1] = 0.0332945;
	ELECTRON_SF_UNC_ID[4][2] = 0.00856674;
	ELECTRON_SF_UNC_ID[4][3] = 0.0101916;
	ELECTRON_SF_UNC_ID[4][4] = 0.0341891;
	ELECTRON_SF_UNC_ID[4][5] = 0.046389;
	ELECTRON_SF_UNC_ID[4][6] = 0.068783;
	ELECTRON_SF_UNC_ID[4][7] = 0.0989201;
	ELECTRON_SF_UNC_ID[5][0] = 0.0121155;
	ELECTRON_SF_UNC_ID[5][1] = 0.0152417;
	ELECTRON_SF_UNC_ID[5][2] = 0.00750628;
	ELECTRON_SF_UNC_ID[5][3] = 0.00350003;
	ELECTRON_SF_UNC_ID[5][4] = 0.00416162;
	ELECTRON_SF_UNC_ID[5][5] = 0.00851604;
	ELECTRON_SF_UNC_ID[5][6] = 0.028257;
	ELECTRON_SF_UNC_ID[5][7] = 0.0298121;

	ELECTRON_SF_UNC_RECO[0] = 0.0182387;
	ELECTRON_SF_UNC_RECO[1] = 0.0110667;
	ELECTRON_SF_UNC_RECO[2] = 0.00815835;
	ELECTRON_SF_UNC_RECO[3] = 0.00713343;
	ELECTRON_SF_UNC_RECO[4] = 0.00420281;
	ELECTRON_SF_UNC_RECO[5] = 0.00649267;
	ELECTRON_SF_UNC_RECO[6] = 0.00516608;
	ELECTRON_SF_UNC_RECO[7] = 0.00551198;
	ELECTRON_SF_UNC_RECO[8] = 0.0260302;
	ELECTRON_SF_UNC_RECO[9] = 0.00359897;
	ELECTRON_SF_UNC_RECO[10] = 0.00506368;
	ELECTRON_SF_UNC_RECO[11] = 0.00331157;
	ELECTRON_SF_UNC_RECO[12] = 0.0061287;
	ELECTRON_SF_UNC_RECO[13] = 0.00635848;
	ELECTRON_SF_UNC_RECO[14] = 0.00530156;
	ELECTRON_SF_UNC_RECO[15] = 0.00530156;
	ELECTRON_SF_UNC_RECO[16] = 0.00635848;
	ELECTRON_SF_UNC_RECO[17] = 0.0061287;
	ELECTRON_SF_UNC_RECO[18] = 0.00331157;
	ELECTRON_SF_UNC_RECO[19] = 0.00506368;
	ELECTRON_SF_UNC_RECO[20] = 0.00359897;
	ELECTRON_SF_UNC_RECO[21] = 0.0260302;
	ELECTRON_SF_UNC_RECO[22] = 0.00551198;
	ELECTRON_SF_UNC_RECO[23] = 0.00516608;
	ELECTRON_SF_UNC_RECO[24] = 0.00649267;
	ELECTRON_SF_UNC_RECO[25] = 0.00420281;
	ELECTRON_SF_UNC_RECO[26] = 0.00713343;
	ELECTRON_SF_UNC_RECO[27] = 0.00815835;
	ELECTRON_SF_UNC_RECO[28] = 0.0110667;
	ELECTRON_SF_UNC_RECO[29] = 0.0182387;

	// PHOTONS SFS
	PHOTON_SF[0][0] = 0.884007;
	PHOTON_SF[0][1] = 0.896628;
	PHOTON_SF[0][2] = 0.90556;
	PHOTON_SF[0][3] = 0.961499;
	PHOTON_SF[0][4] = 0.961499;
	PHOTON_SF[1][0] = 0.884401;
	PHOTON_SF[1][1] = 0.900318;
	PHOTON_SF[1][2] = 0.892832;
	PHOTON_SF[1][3] = 0.907278;
	PHOTON_SF[1][4] = 0.907278;
	PHOTON_SF[2][0] = 0.9938;
	PHOTON_SF[2][1] = 0.9938;
	PHOTON_SF[2][2] = 0.9938;
	PHOTON_SF[2][3] = 0.9938;
	PHOTON_SF[2][4] = 0.9938;
	PHOTON_SF[3][0] = 0.946261;
	PHOTON_SF[3][1] = 0.948876;
	PHOTON_SF[3][2] = 0.943398;
	PHOTON_SF[3][3] = 0.989422;
	PHOTON_SF[3][4] = 0.989422;
	PHOTON_SF[4][0] = 0.945988;
	PHOTON_SF[4][1] = 0.952929;
	PHOTON_SF[4][2] = 0.946219;
	PHOTON_SF[4][3] = 0.98613;
	PHOTON_SF[4][4] = 0.98613;
	PHOTON_SF[5][0] = 0.947127;
	PHOTON_SF[5][1] = 0.959271;
	PHOTON_SF[5][2] = 0.952527;
	PHOTON_SF[5][3] = 0.990491;
	PHOTON_SF[5][4] = 0.990491;
	PHOTON_SF[6][0] = 0.956532;
	PHOTON_SF[6][1] = 0.954528;
	PHOTON_SF[6][2] = 0.950164;
	PHOTON_SF[6][3] = 1.00478;
	PHOTON_SF[6][4] = 1.00478;
	PHOTON_SF[7][0] = 0.9875;
	PHOTON_SF[7][1] = 0.9875;
	PHOTON_SF[7][2] = 0.9875;
	PHOTON_SF[7][3] = 0.9875;
	PHOTON_SF[7][4] = 0.9875;
	PHOTON_SF[8][0] = 0.908546;
	PHOTON_SF[8][1] = 0.921225;
	PHOTON_SF[8][2] = 0.910669;
	PHOTON_SF[8][3] = 0.946951;
	PHOTON_SF[8][4] = 0.946951;
	PHOTON_SF[9][0] = 0.897727;
	PHOTON_SF[9][1] = 0.921013;
	PHOTON_SF[9][2] = 0.927619;
	PHOTON_SF[9][3] = 0.99187;
	PHOTON_SF[9][4] = 0.99187;

	// PHOTON SF UNCERTAINTIES
	PHOTON_SF_UNC[0][0] = 0.0134801;
	PHOTON_SF_UNC[0][1] = 0.015658;
	PHOTON_SF_UNC[0][2] = 0.016519;
	PHOTON_SF_UNC[0][3] = 0.0168844;
	PHOTON_SF_UNC[0][4] = 0.0168846;
	PHOTON_SF_UNC[1][0] = 0.0153668;
	PHOTON_SF_UNC[1][1] = 0.0137329;
	PHOTON_SF_UNC[1][2] = 0.0159758;
	PHOTON_SF_UNC[1][3] = 0.0173759;
	PHOTON_SF_UNC[1][4] = 0.0173761;
	PHOTON_SF_UNC[2][0] = 0.993871;
	PHOTON_SF_UNC[2][1] = 0.993871;
	PHOTON_SF_UNC[2][2] = 0.993871;
	PHOTON_SF_UNC[2][3] = 0.993871;
	PHOTON_SF_UNC[2][4] = 0.993871;
	PHOTON_SF_UNC[3][0] = 0.0154168;
	PHOTON_SF_UNC[3][1] = 0.01479;
	PHOTON_SF_UNC[3][2] = 0.018913;
	PHOTON_SF_UNC[3][3] = 0.0145483;
	PHOTON_SF_UNC[3][4] = 0.0145485;
	PHOTON_SF_UNC[4][0] = 0.0201885;
	PHOTON_SF_UNC[4][1] = 0.0154222;
	PHOTON_SF_UNC[4][2] = 0.0150579;
	PHOTON_SF_UNC[4][3] = 0.0131747;
	PHOTON_SF_UNC[4][4] = 0.0131751;
	PHOTON_SF_UNC[5][0] = 0.0201962;
	PHOTON_SF_UNC[5][1] = 0.0154784;
	PHOTON_SF_UNC[5][2] = 0.0151148;
	PHOTON_SF_UNC[5][3] = 0.0132216;
	PHOTON_SF_UNC[5][4] = 0.0132219;
	PHOTON_SF_UNC[6][0] = 0.0155074;
	PHOTON_SF_UNC[6][1] = 0.014842;
	PHOTON_SF_UNC[6][2] = 0.0189615;
	PHOTON_SF_UNC[6][3] = 0.0147876;
	PHOTON_SF_UNC[6][4] = 0.0147879;
	PHOTON_SF_UNC[7][0] = 0.98751;
	PHOTON_SF_UNC[7][1] = 0.98751;
	PHOTON_SF_UNC[7][2] = 0.98751;
	PHOTON_SF_UNC[7][3] = 0.98751;
	PHOTON_SF_UNC[7][4] = 0.98751;
	PHOTON_SF_UNC[8][0] = 0.0117818;
	PHOTON_SF_UNC[8][1] = 0.00939707;
	PHOTON_SF_UNC[8][2] = 0.0124745;
	PHOTON_SF_UNC[8][3] = 0.01412;
	PHOTON_SF_UNC[8][4] = 0.0141203;
	PHOTON_SF_UNC[9][0] = 0.00920787;
	PHOTON_SF_UNC[9][1] = 0.0120458;
	PHOTON_SF_UNC[9][2] = 0.0130544;
	PHOTON_SF_UNC[9][3] = 0.0130826;
	PHOTON_SF_UNC[9][4] = 0.0130829;

	// MUON EFFICIENCIES
	MUON_ID_LOOSE_DATA[0][0] = 1;
	MUON_ID_LOOSE_DATA[0][1] = 0.997894;
	MUON_ID_LOOSE_DATA[0][2] = 0.996786;
	MUON_ID_LOOSE_DATA[0][3] = 0.998304;
	MUON_ID_LOOSE_DATA[0][4] = 0.996966;
	MUON_ID_LOOSE_DATA[0][5] = 0.996025;
	MUON_ID_LOOSE_DATA[0][6] = 1;
	MUON_ID_HIGHPT_DATA[0][0] = 0.99199;
	MUON_ID_HIGHPT_DATA[0][1] = 0.963946;
	MUON_ID_HIGHPT_DATA[0][2] = 0.958723;
	MUON_ID_HIGHPT_DATA[0][3] = 0.960721;
	MUON_ID_HIGHPT_DATA[0][4] = 0.959969;
	MUON_ID_HIGHPT_DATA[0][5] = 0.960132;
	MUON_ID_HIGHPT_DATA[0][6] = 0.969657;
	MUON_ISO_LOOSE_DATA[0][0] = 0.879075;
	MUON_ISO_LOOSE_DATA[0][1] = 0.917952;
	MUON_ISO_LOOSE_DATA[0][2] = 0.946867;
	MUON_ISO_LOOSE_DATA[0][3] = 0.973483;
	MUON_ISO_LOOSE_DATA[0][4] = 0.988399;
	MUON_ISO_LOOSE_DATA[0][5] = 0.991894;
	MUON_ISO_LOOSE_DATA[0][6] = 0.993687;
	MUON_ISO_HIGHPT_DATA[0][0] = 0.87812;
	MUON_ISO_HIGHPT_DATA[0][1] = 0.917361;
	MUON_ISO_HIGHPT_DATA[0][2] = 0.94659;
	MUON_ISO_HIGHPT_DATA[0][3] = 0.973363;
	MUON_ISO_HIGHPT_DATA[0][4] = 0.988378;
	MUON_ISO_HIGHPT_DATA[0][5] = 0.9919;
	MUON_ISO_HIGHPT_DATA[0][6] = 0.993711;
	MUON_TRIGGER_LOOSE_DATA[0][0] = 0.000241997;
	MUON_TRIGGER_LOOSE_DATA[0][1] = 0.0164698;
	MUON_TRIGGER_LOOSE_DATA[0][2] = 0.856545;
	MUON_TRIGGER_LOOSE_DATA[0][3] = 0.91209;
	MUON_TRIGGER_LOOSE_DATA[0][4] = 0.917525;
	MUON_TRIGGER_LOOSE_DATA[0][5] = 0.91802;
	MUON_TRIGGER_LOOSE_DATA[0][6] = 0.914173;
	MUON_TRIGGER_LOOSE_DATA[0][7] = 0.90638;
	MUON_TRIGGER_LOOSE_DATA[0][8] = 0.888707;
	MUON_TRIGGER_HIGHPT_DATA[0][0] = 0.000241275;
	MUON_TRIGGER_HIGHPT_DATA[0][1] = 0.016697;
	MUON_TRIGGER_HIGHPT_DATA[0][2] = 0.875062;
	MUON_TRIGGER_HIGHPT_DATA[0][3] = 0.931227;
	MUON_TRIGGER_HIGHPT_DATA[0][4] = 0.936236;
	MUON_TRIGGER_HIGHPT_DATA[0][5] = 0.936195;
	MUON_TRIGGER_HIGHPT_DATA[0][6] = 0.93235;
	MUON_TRIGGER_HIGHPT_DATA[0][7] = 0.926049;
	MUON_TRIGGER_HIGHPT_DATA[0][8] = 0.911697;
	MUON_ID_LOOSE_MC[0][0] = 0.999507;
	MUON_ID_LOOSE_MC[0][1] = 0.99852;
	MUON_ID_LOOSE_MC[0][2] = 0.99839;
	MUON_ID_LOOSE_MC[0][3] = 0.998733;
	MUON_ID_LOOSE_MC[0][4] = 0.998767;
	MUON_ID_LOOSE_MC[0][5] = 0.998828;
	MUON_ID_LOOSE_MC[0][6] = 0.999069;
	MUON_ID_HIGHPT_MC[0][0] = 0.969871;
	MUON_ID_HIGHPT_MC[0][1] = 0.972887;
	MUON_ID_HIGHPT_MC[0][2] = 0.973377;
	MUON_ID_HIGHPT_MC[0][3] = 0.973653;
	MUON_ID_HIGHPT_MC[0][4] = 0.973631;
	MUON_ID_HIGHPT_MC[0][5] = 0.974453;
	MUON_ID_HIGHPT_MC[0][6] = 0.974497;
	MUON_ISO_LOOSE_MC[0][0] = 0.89471;
	MUON_ISO_LOOSE_MC[0][1] = 0.926319;
	MUON_ISO_LOOSE_MC[0][2] = 0.952042;
	MUON_ISO_LOOSE_MC[0][3] = 0.976816;
	MUON_ISO_LOOSE_MC[0][4] = 0.990825;
	MUON_ISO_LOOSE_MC[0][5] = 0.993989;
	MUON_ISO_LOOSE_MC[0][6] = 0.995229;
	MUON_ISO_HIGHPT_MC[0][0] = 0.893834;
	MUON_ISO_HIGHPT_MC[0][1] = 0.92589;
	MUON_ISO_HIGHPT_MC[0][2] = 0.951948;
	MUON_ISO_HIGHPT_MC[0][3] = 0.976794;
	MUON_ISO_HIGHPT_MC[0][4] = 0.990713;
	MUON_ISO_HIGHPT_MC[0][5] = 0.993995;
	MUON_ISO_HIGHPT_MC[0][6] = 0.995186;
	MUON_TRIGGER_LOOSE_MC[0][0] = 8.77849e-05;
	MUON_TRIGGER_LOOSE_MC[0][1] = 0.00829161;
	MUON_TRIGGER_LOOSE_MC[0][2] = 0.917587;
	MUON_TRIGGER_LOOSE_MC[0][3] = 0.939189;
	MUON_TRIGGER_LOOSE_MC[0][4] = 0.939964;
	MUON_TRIGGER_LOOSE_MC[0][5] = 0.940834;
	MUON_TRIGGER_LOOSE_MC[0][6] = 0.938727;
	MUON_TRIGGER_LOOSE_MC[0][7] = 0.933991;
	MUON_TRIGGER_LOOSE_MC[0][8] = 0.911105;
	MUON_TRIGGER_HIGHPT_MC[0][0] = 8.76155e-05;
	MUON_TRIGGER_HIGHPT_MC[0][1] = 0.00838101;
	MUON_TRIGGER_HIGHPT_MC[0][2] = 0.932104;
	MUON_TRIGGER_HIGHPT_MC[0][3] = 0.953411;
	MUON_TRIGGER_HIGHPT_MC[0][4] = 0.953739;
	MUON_TRIGGER_HIGHPT_MC[0][5] = 0.954368;
	MUON_TRIGGER_HIGHPT_MC[0][6] = 0.952725;
	MUON_TRIGGER_HIGHPT_MC[0][7] = 0.95148;
	MUON_TRIGGER_HIGHPT_MC[0][8] = 0.928014;
	MUON_ID_LOOSE_DATA[1][0] = 0.998829;
	MUON_ID_LOOSE_DATA[1][1] = 0.995886;
	MUON_ID_LOOSE_DATA[1][2] = 0.994204;
	MUON_ID_LOOSE_DATA[1][3] = 0.995501;
	MUON_ID_LOOSE_DATA[1][4] = 0.995282;
	MUON_ID_LOOSE_DATA[1][5] = 0.994512;
	MUON_ID_LOOSE_DATA[1][6] = 0.99618;
	MUON_ID_HIGHPT_DATA[1][0] = 0.964333;
	MUON_ID_HIGHPT_DATA[1][1] = 0.952049;
	MUON_ID_HIGHPT_DATA[1][2] = 0.954987;
	MUON_ID_HIGHPT_DATA[1][3] = 0.959168;
	MUON_ID_HIGHPT_DATA[1][4] = 0.960021;
	MUON_ID_HIGHPT_DATA[1][5] = 0.96149;
	MUON_ID_HIGHPT_DATA[1][6] = 0.962634;
	MUON_ISO_LOOSE_DATA[1][0] = 0.890724;
	MUON_ISO_LOOSE_DATA[1][1] = 0.930753;
	MUON_ISO_LOOSE_DATA[1][2] = 0.952599;
	MUON_ISO_LOOSE_DATA[1][3] = 0.976171;
	MUON_ISO_LOOSE_DATA[1][4] = 0.990396;
	MUON_ISO_LOOSE_DATA[1][5] = 0.993229;
	MUON_ISO_LOOSE_DATA[1][6] = 0.994868;
	MUON_ISO_HIGHPT_DATA[1][0] = 0.889531;
	MUON_ISO_HIGHPT_DATA[1][1] = 0.929915;
	MUON_ISO_HIGHPT_DATA[1][2] = 0.952341;
	MUON_ISO_HIGHPT_DATA[1][3] = 0.976095;
	MUON_ISO_HIGHPT_DATA[1][4] = 0.990344;
	MUON_ISO_HIGHPT_DATA[1][5] = 0.993236;
	MUON_ISO_HIGHPT_DATA[1][6] = 0.99501;
	MUON_TRIGGER_LOOSE_DATA[1][0] = 0.000833003;
	MUON_TRIGGER_LOOSE_DATA[1][1] = 0.0192174;
	MUON_TRIGGER_LOOSE_DATA[1][2] = 0.84156;
	MUON_TRIGGER_LOOSE_DATA[1][3] = 0.917703;
	MUON_TRIGGER_LOOSE_DATA[1][4] = 0.922287;
	MUON_TRIGGER_LOOSE_DATA[1][5] = 0.923559;
	MUON_TRIGGER_LOOSE_DATA[1][6] = 0.918514;
	MUON_TRIGGER_LOOSE_DATA[1][7] = 0.903209;
	MUON_TRIGGER_LOOSE_DATA[1][8] = 0.885367;
	MUON_TRIGGER_HIGHPT_DATA[1][0] = 0.000842043;
	MUON_TRIGGER_HIGHPT_DATA[1][1] = 0.0193012;
	MUON_TRIGGER_HIGHPT_DATA[1][2] = 0.853633;
	MUON_TRIGGER_HIGHPT_DATA[1][3] = 0.928499;
	MUON_TRIGGER_HIGHPT_DATA[1][4] = 0.933267;
	MUON_TRIGGER_HIGHPT_DATA[1][5] = 0.934552;
	MUON_TRIGGER_HIGHPT_DATA[1][6] = 0.929864;
	MUON_TRIGGER_HIGHPT_DATA[1][7] = 0.91739;
	MUON_TRIGGER_HIGHPT_DATA[1][8] = 0.899688;
	MUON_ID_LOOSE_MC[1][0] = 0.997363;
	MUON_ID_LOOSE_MC[1][1] = 0.99791;
	MUON_ID_LOOSE_MC[1][2] = 0.997684;
	MUON_ID_LOOSE_MC[1][3] = 0.998098;
	MUON_ID_LOOSE_MC[1][4] = 0.998487;
	MUON_ID_LOOSE_MC[1][5] = 0.998757;
	MUON_ID_LOOSE_MC[1][6] = 0.999684;
	MUON_ID_HIGHPT_MC[1][0] = 0.977565;
	MUON_ID_HIGHPT_MC[1][1] = 0.982473;
	MUON_ID_HIGHPT_MC[1][2] = 0.9818;
	MUON_ID_HIGHPT_MC[1][3] = 0.98312;
	MUON_ID_HIGHPT_MC[1][4] = 0.982535;
	MUON_ID_HIGHPT_MC[1][5] = 0.98349;
	MUON_ID_HIGHPT_MC[1][6] = 0.982766;
	MUON_ISO_LOOSE_MC[1][0] = 0.892923;
	MUON_ISO_LOOSE_MC[1][1] = 0.935223;
	MUON_ISO_LOOSE_MC[1][2] = 0.954136;
	MUON_ISO_LOOSE_MC[1][3] = 0.976164;
	MUON_ISO_LOOSE_MC[1][4] = 0.991295;
	MUON_ISO_LOOSE_MC[1][5] = 0.994055;
	MUON_ISO_LOOSE_MC[1][6] = 0.995817;
	MUON_ISO_HIGHPT_MC[1][0] = 0.892516;
	MUON_ISO_HIGHPT_MC[1][1] = 0.935019;
	MUON_ISO_HIGHPT_MC[1][2] = 0.953824;
	MUON_ISO_HIGHPT_MC[1][3] = 0.976132;
	MUON_ISO_HIGHPT_MC[1][4] = 0.991314;
	MUON_ISO_HIGHPT_MC[1][5] = 0.994047;
	MUON_ISO_HIGHPT_MC[1][6] = 0.99578;
	MUON_TRIGGER_LOOSE_MC[1][0] = 0.000444135;
	MUON_TRIGGER_LOOSE_MC[1][1] = 0.0115531;
	MUON_TRIGGER_LOOSE_MC[1][2] = 0.924932;
	MUON_TRIGGER_LOOSE_MC[1][3] = 0.968141;
	MUON_TRIGGER_LOOSE_MC[1][4] = 0.968743;
	MUON_TRIGGER_LOOSE_MC[1][5] = 0.970193;
	MUON_TRIGGER_LOOSE_MC[1][6] = 0.971123;
	MUON_TRIGGER_LOOSE_MC[1][7] = 0.969557;
	MUON_TRIGGER_LOOSE_MC[1][8] = 0.939054;
	MUON_TRIGGER_HIGHPT_MC[1][0] = 0.000435311;
	MUON_TRIGGER_HIGHPT_MC[1][1] = 0.0115903;
	MUON_TRIGGER_HIGHPT_MC[1][2] = 0.930451;
	MUON_TRIGGER_HIGHPT_MC[1][3] = 0.97381;
	MUON_TRIGGER_HIGHPT_MC[1][4] = 0.974611;
	MUON_TRIGGER_HIGHPT_MC[1][5] = 0.976113;
	MUON_TRIGGER_HIGHPT_MC[1][6] = 0.975809;
	MUON_TRIGGER_HIGHPT_MC[1][7] = 0.978771;
	MUON_TRIGGER_HIGHPT_MC[1][8] = 0.948755;
	MUON_ID_LOOSE_DATA[2][0] = 1;
	MUON_ID_LOOSE_DATA[2][1] = 0.997041;
	MUON_ID_LOOSE_DATA[2][2] = 0.997307;
	MUON_ID_LOOSE_DATA[2][3] = 0.998397;
	MUON_ID_LOOSE_DATA[2][4] = 0.99865;
	MUON_ID_LOOSE_DATA[2][5] = 0.995921;
	MUON_ID_LOOSE_DATA[2][6] = 0.998732;
	MUON_ID_HIGHPT_DATA[2][0] = 0.999958;
	MUON_ID_HIGHPT_DATA[2][1] = 0.980702;
	MUON_ID_HIGHPT_DATA[2][2] = 0.984571;
	MUON_ID_HIGHPT_DATA[2][3] = 0.984735;
	MUON_ID_HIGHPT_DATA[2][4] = 0.985623;
	MUON_ID_HIGHPT_DATA[2][5] = 0.983197;
	MUON_ID_HIGHPT_DATA[2][6] = 0.987808;
	MUON_ISO_LOOSE_DATA[2][0] = 0.905713;
	MUON_ISO_LOOSE_DATA[2][1] = 0.943181;
	MUON_ISO_LOOSE_DATA[2][2] = 0.963979;
	MUON_ISO_LOOSE_DATA[2][3] = 0.980984;
	MUON_ISO_LOOSE_DATA[2][4] = 0.992715;
	MUON_ISO_LOOSE_DATA[2][5] = 0.995036;
	MUON_ISO_LOOSE_DATA[2][6] = 0.995881;
	MUON_ISO_HIGHPT_DATA[2][0] = 0.9047;
	MUON_ISO_HIGHPT_DATA[2][1] = 0.942782;
	MUON_ISO_HIGHPT_DATA[2][2] = 0.96371;
	MUON_ISO_HIGHPT_DATA[2][3] = 0.980898;
	MUON_ISO_HIGHPT_DATA[2][4] = 0.992533;
	MUON_ISO_HIGHPT_DATA[2][5] = 0.995026;
	MUON_ISO_HIGHPT_DATA[2][6] = 0.995919;
	MUON_TRIGGER_LOOSE_DATA[2][0] = 0.0010555;
	MUON_TRIGGER_LOOSE_DATA[2][1] = 0.0281372;
	MUON_TRIGGER_LOOSE_DATA[2][2] = 0.788449;
	MUON_TRIGGER_LOOSE_DATA[2][3] = 0.874352;
	MUON_TRIGGER_LOOSE_DATA[2][4] = 0.881563;
	MUON_TRIGGER_LOOSE_DATA[2][5] = 0.885286;
	MUON_TRIGGER_LOOSE_DATA[2][6] = 0.884138;
	MUON_TRIGGER_LOOSE_DATA[2][7] = 0.882944;
	MUON_TRIGGER_LOOSE_DATA[2][8] = 0.873511;
	MUON_TRIGGER_HIGHPT_DATA[2][0] = 0.00104391;
	MUON_TRIGGER_HIGHPT_DATA[2][1] = 0.0280684;
	MUON_TRIGGER_HIGHPT_DATA[2][2] = 0.794699;
	MUON_TRIGGER_HIGHPT_DATA[2][3] = 0.880748;
	MUON_TRIGGER_HIGHPT_DATA[2][4] = 0.887591;
	MUON_TRIGGER_HIGHPT_DATA[2][5] = 0.890361;
	MUON_TRIGGER_HIGHPT_DATA[2][6] = 0.889464;
	MUON_TRIGGER_HIGHPT_DATA[2][7] = 0.888038;
	MUON_TRIGGER_HIGHPT_DATA[2][8] = 0.881197;
	MUON_ID_LOOSE_MC[2][0] = 1;
	MUON_ID_LOOSE_MC[2][1] = 0.998649;
	MUON_ID_LOOSE_MC[2][2] = 0.998903;
	MUON_ID_LOOSE_MC[2][3] = 0.999125;
	MUON_ID_LOOSE_MC[2][4] = 0.999336;
	MUON_ID_LOOSE_MC[2][5] = 0.999383;
	MUON_ID_LOOSE_MC[2][6] = 0.999179;
	MUON_ID_HIGHPT_MC[2][0] = 0.995347;
	MUON_ID_HIGHPT_MC[2][1] = 0.995144;
	MUON_ID_HIGHPT_MC[2][2] = 0.995095;
	MUON_ID_HIGHPT_MC[2][3] = 0.995066;
	MUON_ID_HIGHPT_MC[2][4] = 0.995325;
	MUON_ID_HIGHPT_MC[2][5] = 0.995007;
	MUON_ID_HIGHPT_MC[2][6] = 0.995643;
	MUON_ISO_LOOSE_MC[2][0] = 0.911568;
	MUON_ISO_LOOSE_MC[2][1] = 0.945881;
	MUON_ISO_LOOSE_MC[2][2] = 0.962769;
	MUON_ISO_LOOSE_MC[2][3] = 0.980191;
	MUON_ISO_LOOSE_MC[2][4] = 0.992551;
	MUON_ISO_LOOSE_MC[2][5] = 0.995341;
	MUON_ISO_LOOSE_MC[2][6] = 0.996667;
	MUON_ISO_HIGHPT_MC[2][0] = 0.911241;
	MUON_ISO_HIGHPT_MC[2][1] = 0.945816;
	MUON_ISO_HIGHPT_MC[2][2] = 0.962724;
	MUON_ISO_HIGHPT_MC[2][3] = 0.980177;
	MUON_ISO_HIGHPT_MC[2][4] = 0.992666;
	MUON_ISO_HIGHPT_MC[2][5] = 0.995482;
	MUON_ISO_HIGHPT_MC[2][6] = 0.996669;
	MUON_TRIGGER_LOOSE_MC[2][0] = 0.000275576;
	MUON_TRIGGER_LOOSE_MC[2][1] = 0.0160394;
	MUON_TRIGGER_LOOSE_MC[2][2] = 0.835014;
	MUON_TRIGGER_LOOSE_MC[2][3] = 0.892034;
	MUON_TRIGGER_LOOSE_MC[2][4] = 0.892813;
	MUON_TRIGGER_LOOSE_MC[2][5] = 0.894989;
	MUON_TRIGGER_LOOSE_MC[2][6] = 0.897931;
	MUON_TRIGGER_LOOSE_MC[2][7] = 0.890065;
	MUON_TRIGGER_LOOSE_MC[2][8] = 0.912288;
	MUON_TRIGGER_HIGHPT_MC[2][0] = 0.000274416;
	MUON_TRIGGER_HIGHPT_MC[2][1] = 0.0160097;
	MUON_TRIGGER_HIGHPT_MC[2][2] = 0.836504;
	MUON_TRIGGER_HIGHPT_MC[2][3] = 0.892962;
	MUON_TRIGGER_HIGHPT_MC[2][4] = 0.893937;
	MUON_TRIGGER_HIGHPT_MC[2][5] = 0.895715;
	MUON_TRIGGER_HIGHPT_MC[2][6] = 0.898596;
	MUON_TRIGGER_HIGHPT_MC[2][7] = 0.892444;
	MUON_TRIGGER_HIGHPT_MC[2][8] = 0.912583;
	MUON_ID_LOOSE_DATA[3][0] = 0.999673;
	MUON_ID_LOOSE_DATA[3][1] = 0.992308;
	MUON_ID_LOOSE_DATA[3][2] = 0.991898;
	MUON_ID_LOOSE_DATA[3][3] = 0.991747;
	MUON_ID_LOOSE_DATA[3][4] = 0.992971;
	MUON_ID_LOOSE_DATA[3][5] = 0.992314;
	MUON_ID_LOOSE_DATA[3][6] = 0.993223;
	MUON_ID_HIGHPT_DATA[3][0] = 0.973049;
	MUON_ID_HIGHPT_DATA[3][1] = 0.964256;
	MUON_ID_HIGHPT_DATA[3][2] = 0.961626;
	MUON_ID_HIGHPT_DATA[3][3] = 0.95978;
	MUON_ID_HIGHPT_DATA[3][4] = 0.959157;
	MUON_ID_HIGHPT_DATA[3][5] = 0.957948;
	MUON_ID_HIGHPT_DATA[3][6] = 0.958645;
	MUON_ISO_LOOSE_DATA[3][0] = 0.929878;
	MUON_ISO_LOOSE_DATA[3][1] = 0.961148;
	MUON_ISO_LOOSE_DATA[3][2] = 0.976787;
	MUON_ISO_LOOSE_DATA[3][3] = 0.988362;
	MUON_ISO_LOOSE_DATA[3][4] = 0.995224;
	MUON_ISO_LOOSE_DATA[3][5] = 0.996369;
	MUON_ISO_LOOSE_DATA[3][6] = 0.997034;
	MUON_ISO_HIGHPT_DATA[3][0] = 0.928773;
	MUON_ISO_HIGHPT_DATA[3][1] = 0.960429;
	MUON_ISO_HIGHPT_DATA[3][2] = 0.976338;
	MUON_ISO_HIGHPT_DATA[3][3] = 0.988262;
	MUON_ISO_HIGHPT_DATA[3][4] = 0.995232;
	MUON_ISO_HIGHPT_DATA[3][5] = 0.996419;
	MUON_ISO_HIGHPT_DATA[3][6] = 0.997041;
	MUON_TRIGGER_LOOSE_DATA[3][0] = 0.00673737;
	MUON_TRIGGER_LOOSE_DATA[3][1] = 0.0689954;
	MUON_TRIGGER_LOOSE_DATA[3][2] = 0.63076;
	MUON_TRIGGER_LOOSE_DATA[3][3] = 0.757008;
	MUON_TRIGGER_LOOSE_DATA[3][4] = 0.78507;
	MUON_TRIGGER_LOOSE_DATA[3][5] = 0.801005;
	MUON_TRIGGER_LOOSE_DATA[3][6] = 0.804341;
	MUON_TRIGGER_LOOSE_DATA[3][7] = 0.797834;
	MUON_TRIGGER_LOOSE_DATA[3][8] = 0.785411;
	MUON_TRIGGER_HIGHPT_DATA[3][0] = 0.0068072;
	MUON_TRIGGER_HIGHPT_DATA[3][1] = 0.0698809;
	MUON_TRIGGER_HIGHPT_DATA[3][2] = 0.645432;
	MUON_TRIGGER_HIGHPT_DATA[3][3] = 0.77488;
	MUON_TRIGGER_HIGHPT_DATA[3][4] = 0.802788;
	MUON_TRIGGER_HIGHPT_DATA[3][5] = 0.816415;
	MUON_TRIGGER_HIGHPT_DATA[3][6] = 0.821003;
	MUON_TRIGGER_HIGHPT_DATA[3][7] = 0.814951;
	MUON_TRIGGER_HIGHPT_DATA[3][8] = 0.774534;
	MUON_ID_LOOSE_MC[3][0] = 1;
	MUON_ID_LOOSE_MC[3][1] = 0.996517;
	MUON_ID_LOOSE_MC[3][2] = 0.995512;
	MUON_ID_LOOSE_MC[3][3] = 0.997213;
	MUON_ID_LOOSE_MC[3][4] = 0.998159;
	MUON_ID_LOOSE_MC[3][5] = 0.998186;
	MUON_ID_LOOSE_MC[3][6] = 0.999899;
	MUON_ID_HIGHPT_MC[3][0] = 0.987032;
	MUON_ID_HIGHPT_MC[3][1] = 0.9888;
	MUON_ID_HIGHPT_MC[3][2] = 0.984009;
	MUON_ID_HIGHPT_MC[3][3] = 0.986752;
	MUON_ID_HIGHPT_MC[3][4] = 0.986902;
	MUON_ID_HIGHPT_MC[3][5] = 0.986838;
	MUON_ID_HIGHPT_MC[3][6] = 0.990606;
	MUON_ISO_LOOSE_MC[3][0] = 0.935102;
	MUON_ISO_LOOSE_MC[3][1] = 0.960554;
	MUON_ISO_LOOSE_MC[3][2] = 0.974836;
	MUON_ISO_LOOSE_MC[3][3] = 0.987728;
	MUON_ISO_LOOSE_MC[3][4] = 0.995222;
	MUON_ISO_LOOSE_MC[3][5] = 0.996178;
	MUON_ISO_LOOSE_MC[3][6] = 0.996701;
	MUON_ISO_HIGHPT_MC[3][0] = 0.9345;
	MUON_ISO_HIGHPT_MC[3][1] = 0.960529;
	MUON_ISO_HIGHPT_MC[3][2] = 0.974598;
	MUON_ISO_HIGHPT_MC[3][3] = 0.987671;
	MUON_ISO_HIGHPT_MC[3][4] = 0.99518;
	MUON_ISO_HIGHPT_MC[3][5] = 0.996131;
	MUON_ISO_HIGHPT_MC[3][6] = 0.996694;
	MUON_TRIGGER_LOOSE_MC[3][0] = 0.00224537;
	MUON_TRIGGER_LOOSE_MC[3][1] = 0.0440228;
	MUON_TRIGGER_LOOSE_MC[3][2] = 0.733557;
	MUON_TRIGGER_LOOSE_MC[3][3] = 0.838855;
	MUON_TRIGGER_LOOSE_MC[3][4] = 0.853953;
	MUON_TRIGGER_LOOSE_MC[3][5] = 0.860836;
	MUON_TRIGGER_LOOSE_MC[3][6] = 0.857747;
	MUON_TRIGGER_LOOSE_MC[3][7] = 0.840479;
	MUON_TRIGGER_LOOSE_MC[3][8] = 0.934433;
	MUON_TRIGGER_HIGHPT_MC[3][0] = 0.00212545;
	MUON_TRIGGER_HIGHPT_MC[3][1] = 0.0437604;
	MUON_TRIGGER_HIGHPT_MC[3][2] = 0.740227;
	MUON_TRIGGER_HIGHPT_MC[3][3] = 0.846157;
	MUON_TRIGGER_HIGHPT_MC[3][4] = 0.860103;
	MUON_TRIGGER_HIGHPT_MC[3][5] = 0.866098;
	MUON_TRIGGER_HIGHPT_MC[3][6] = 0.863988;
	MUON_TRIGGER_HIGHPT_MC[3][7] = 0.84842;
	MUON_TRIGGER_HIGHPT_MC[3][8] = 0.930904;

	// MUON EFFICIENCY UNCERTAINTIES
	MUON_ID_UNC_LOOSE_DATA[0][0] = 6.94613e-08;
	MUON_ID_UNC_LOOSE_DATA[0][1] = 0.00151717;
	MUON_ID_UNC_LOOSE_DATA[0][2] = 0.00201308;
	MUON_ID_UNC_LOOSE_DATA[0][3] = 0.000313195;
	MUON_ID_UNC_LOOSE_DATA[0][4] = 0.000612682;
	MUON_ID_UNC_LOOSE_DATA[0][5] = 0.00226443;
	MUON_ID_UNC_LOOSE_DATA[0][6] = 1.36565e-09;
	MUON_ID_UNC_HIGHPT_DATA[0][0] = 0.0058325;
	MUON_ID_UNC_HIGHPT_DATA[0][1] = 0.00579926;
	MUON_ID_UNC_HIGHPT_DATA[0][2] = 0.00148282;
	MUON_ID_UNC_HIGHPT_DATA[0][3] = 0.000848865;
	MUON_ID_UNC_HIGHPT_DATA[0][4] = 0.000924437;
	MUON_ID_UNC_HIGHPT_DATA[0][5] = 0.000580344;
	MUON_ID_UNC_HIGHPT_DATA[0][6] = 0.0025401;
	MUON_ISO_UNC_LOOSE_DATA[0][0] = 0.00278236;
	MUON_ISO_UNC_LOOSE_DATA[0][1] = 0.00360292;
	MUON_ISO_UNC_LOOSE_DATA[0][2] = 0.00189865;
	MUON_ISO_UNC_LOOSE_DATA[0][3] = 0.000460735;
	MUON_ISO_UNC_LOOSE_DATA[0][4] = 0.000192883;
	MUON_ISO_UNC_LOOSE_DATA[0][5] = 0.000178093;
	MUON_ISO_UNC_LOOSE_DATA[0][6] = 0.00021287;
	MUON_ISO_UNC_HIGHPT_DATA[0][0] = 0.00322658;
	MUON_ISO_UNC_HIGHPT_DATA[0][1] = 0.0036953;
	MUON_ISO_UNC_HIGHPT_DATA[0][2] = 0.00192867;
	MUON_ISO_UNC_HIGHPT_DATA[0][3] = 0.000466109;
	MUON_ISO_UNC_HIGHPT_DATA[0][4] = 0.000185791;
	MUON_ISO_UNC_HIGHPT_DATA[0][5] = 0.000102084;
	MUON_ISO_UNC_HIGHPT_DATA[0][6] = 0.000214082;
	// MUON_TRIGGER_UNC_LOOSE_DATA[0][0] = 1.96732e-05;
	// MUON_TRIGGER_UNC_LOOSE_DATA[0][1] = 0.000364793;
	// MUON_TRIGGER_UNC_LOOSE_DATA[0][2] = 0.000616222;
	// MUON_TRIGGER_UNC_LOOSE_DATA[0][3] = 0.000715149;
	// MUON_TRIGGER_UNC_LOOSE_DATA[0][4] = 0.000824014;
	// MUON_TRIGGER_UNC_LOOSE_DATA[0][5] = 0.000488275;
	// MUON_TRIGGER_UNC_LOOSE_DATA[0][6] = 0.000860823;
	// MUON_TRIGGER_UNC_LOOSE_DATA[0][7] = 0.00171415;
	// MUON_TRIGGER_UNC_LOOSE_DATA[0][8] = 0.00366344;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[0][0] = 7.70321e-06;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[0][1] = 0.00479116;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[0][2] = 0.000493172;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[0][3] = 0.000372831;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[0][4] = 0.000394522;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[0][5] = 0.000338051;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[0][6] = 0.000730853;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[0][7] = 0.00148436;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[0][8] = 0.0033648;
	MUON_ID_UNC_LOOSE_MC[0][0] = 0.000187578;
	MUON_ID_UNC_LOOSE_MC[0][1] = 0.000830821;
	MUON_ID_UNC_LOOSE_MC[0][2] = 0.000150059;
	MUON_ID_UNC_LOOSE_MC[0][3] = 5.16007e-05;
	MUON_ID_UNC_LOOSE_MC[0][4] = 2.64581e-05;
	MUON_ID_UNC_LOOSE_MC[0][5] = 0.000828951;
	MUON_ID_UNC_LOOSE_MC[0][6] = 0.000517748;
	MUON_ID_UNC_HIGHPT_MC[0][0] = 0.00127528;
	MUON_ID_UNC_HIGHPT_MC[0][1] = 0.000826834;
	MUON_ID_UNC_HIGHPT_MC[0][2] = 0.000682959;
	MUON_ID_UNC_HIGHPT_MC[0][3] = 0.000499467;
	MUON_ID_UNC_HIGHPT_MC[0][4] = 0.000263423;
	MUON_ID_UNC_HIGHPT_MC[0][5] = 0.00111686;
	MUON_ID_UNC_HIGHPT_MC[0][6] = 0.000669284;
	MUON_ISO_UNC_LOOSE_MC[0][0] = 0.00242987;
	MUON_ISO_UNC_LOOSE_MC[0][1] = 0.00332098;
	MUON_ISO_UNC_LOOSE_MC[0][2] = 0.00199579;
	MUON_ISO_UNC_LOOSE_MC[0][3] = 0.000195305;
	MUON_ISO_UNC_LOOSE_MC[0][4] = 0.000140876;
	MUON_ISO_UNC_LOOSE_MC[0][5] = 0.000134217;
	MUON_ISO_UNC_LOOSE_MC[0][6] = 0.000180631;
	MUON_ISO_UNC_HIGHPT_MC[0][0] = 0.00247862;
	MUON_ISO_UNC_HIGHPT_MC[0][1] = 0.00333309;
	MUON_ISO_UNC_HIGHPT_MC[0][2] = 0.00199584;
	MUON_ISO_UNC_HIGHPT_MC[0][3] = 0.000194271;
	MUON_ISO_UNC_HIGHPT_MC[0][4] = 0.000273527;
	MUON_ISO_UNC_HIGHPT_MC[0][5] = 0.00013478;
	MUON_ISO_UNC_HIGHPT_MC[0][6] = 0.000182703;
	MUON_TRIGGER_UNC_LOOSE_MC[0][0] = 1.30234e-05;
	MUON_TRIGGER_UNC_LOOSE_MC[0][1] = 0.000195184;
	MUON_TRIGGER_UNC_LOOSE_MC[0][2] = 0.000630042;
	MUON_TRIGGER_UNC_LOOSE_MC[0][3] = 0.000706813;
	MUON_TRIGGER_UNC_LOOSE_MC[0][4] = 0.000614214;
	MUON_TRIGGER_UNC_LOOSE_MC[0][5] = 0.000545957;
	MUON_TRIGGER_UNC_LOOSE_MC[0][6] = 0.00114231;
	MUON_TRIGGER_UNC_LOOSE_MC[0][7] = 0.00257207;
	MUON_TRIGGER_UNC_LOOSE_MC[0][8] = 0.00535755;
	MUON_TRIGGER_UNC_HIGHPT_MC[0][0] = 4.57307e-05;
	MUON_TRIGGER_UNC_HIGHPT_MC[0][1] = 0.000213643;
	MUON_TRIGGER_UNC_HIGHPT_MC[0][2] = 0.000567509;
	MUON_TRIGGER_UNC_HIGHPT_MC[0][3] = 0.000588149;
	MUON_TRIGGER_UNC_HIGHPT_MC[0][4] = 0.000507257;
	MUON_TRIGGER_UNC_HIGHPT_MC[0][5] = 0.000471968;
	MUON_TRIGGER_UNC_HIGHPT_MC[0][6] = 0.00106998;
	MUON_TRIGGER_UNC_HIGHPT_MC[0][7] = 0.00197248;
	MUON_TRIGGER_UNC_HIGHPT_MC[0][8] = 0.0054568;
	MUON_ID_UNC_LOOSE_DATA[1][0] = 0.00117129;
	MUON_ID_UNC_LOOSE_DATA[1][1] = 0.00327549;
	MUON_ID_UNC_LOOSE_DATA[1][2] = 0.0013365;
	MUON_ID_UNC_LOOSE_DATA[1][3] = 0.000411608;
	MUON_ID_UNC_LOOSE_DATA[1][4] = 0.000238296;
	MUON_ID_UNC_LOOSE_DATA[1][5] = 0.0007724;
	MUON_ID_UNC_LOOSE_DATA[1][6] = 0.00342998;
	MUON_ID_UNC_HIGHPT_DATA[1][0] = 0.0100371;
	MUON_ID_UNC_HIGHPT_DATA[1][1] = 0.00577705;
	MUON_ID_UNC_HIGHPT_DATA[1][2] = 0.00157254;
	MUON_ID_UNC_HIGHPT_DATA[1][3] = 0.000788025;
	MUON_ID_UNC_HIGHPT_DATA[1][4] = 0.000572921;
	MUON_ID_UNC_HIGHPT_DATA[1][5] = 0.000929799;
	MUON_ID_UNC_HIGHPT_DATA[1][6] = 0.00325309;
	MUON_ISO_UNC_LOOSE_DATA[1][0] = 0.00252035;
	MUON_ISO_UNC_LOOSE_DATA[1][1] = 0.00359154;
	MUON_ISO_UNC_LOOSE_DATA[1][2] = 0.00208887;
	MUON_ISO_UNC_LOOSE_DATA[1][3] = 0.000430911;
	MUON_ISO_UNC_LOOSE_DATA[1][4] = 0.000257928;
	MUON_ISO_UNC_LOOSE_DATA[1][5] = 0.000204708;
	MUON_ISO_UNC_LOOSE_DATA[1][6] = 0.000377402;
	MUON_ISO_UNC_HIGHPT_DATA[1][0] = 0.0025507;
	MUON_ISO_UNC_HIGHPT_DATA[1][1] = 0.00366345;
	MUON_ISO_UNC_HIGHPT_DATA[1][2] = 0.00210372;
	MUON_ISO_UNC_HIGHPT_DATA[1][3] = 0.000429083;
	MUON_ISO_UNC_HIGHPT_DATA[1][4] = 0.000180012;
	MUON_ISO_UNC_HIGHPT_DATA[1][5] = 0.000212362;
	MUON_ISO_UNC_HIGHPT_DATA[1][6] = 0.00283321;
	// MUON_TRIGGER_UNC_LOOSE_DATA[1][0] = 2.99985e-05;
	// MUON_TRIGGER_UNC_LOOSE_DATA[1][1] = 0.000317854;
	// MUON_TRIGGER_UNC_LOOSE_DATA[1][2] = 0.000907689;
	// MUON_TRIGGER_UNC_LOOSE_DATA[1][3] = 0.000995557;
	// MUON_TRIGGER_UNC_LOOSE_DATA[1][4] = 0.00088836;
	// MUON_TRIGGER_UNC_LOOSE_DATA[1][5] = 0.000770469;
	// MUON_TRIGGER_UNC_LOOSE_DATA[1][6] = 0.00149059;
	// MUON_TRIGGER_UNC_LOOSE_DATA[1][7] = 0.00369728;
	// MUON_TRIGGER_UNC_LOOSE_DATA[1][8] = 0.00956495;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[1][0] = 3.04514e-05;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[1][1] = 0.000324651;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[1][2] = 0.00100692;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[1][3] = 0.00106928;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[1][4] = 0.000699363;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[1][5] = 0.000629975;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[1][6] = 0.00138993;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[1][7] = 0.00346677;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[1][8] = 0.00885234;
	MUON_ID_UNC_LOOSE_MC[1][0] = 0.00138712;
	MUON_ID_UNC_LOOSE_MC[1][1] = 0.000501587;
	MUON_ID_UNC_LOOSE_MC[1][2] = 0.000277977;
	MUON_ID_UNC_LOOSE_MC[1][3] = 0.000115235;
	MUON_ID_UNC_LOOSE_MC[1][4] = 0.000743016;
	MUON_ID_UNC_LOOSE_MC[1][5] = 0.00015394;
	MUON_ID_UNC_LOOSE_MC[1][6] = 0.000241467;
	MUON_ID_UNC_HIGHPT_MC[1][0] = 0.0017869;
	MUON_ID_UNC_HIGHPT_MC[1][1] = 0.00096073;
	MUON_ID_UNC_HIGHPT_MC[1][2] = 0.000663289;
	MUON_ID_UNC_HIGHPT_MC[1][3] = 0.000529133;
	MUON_ID_UNC_HIGHPT_MC[1][4] = 0.00122616;
	MUON_ID_UNC_HIGHPT_MC[1][5] = 0.000514422;
	MUON_ID_UNC_HIGHPT_MC[1][6] = 0.000789877;
	MUON_ISO_UNC_LOOSE_MC[1][0] = 0.00281129;
	MUON_ISO_UNC_LOOSE_MC[1][1] = 0.00282331;
	MUON_ISO_UNC_LOOSE_MC[1][2] = 0.00249845;
	MUON_ISO_UNC_LOOSE_MC[1][3] = 0.000305684;
	MUON_ISO_UNC_LOOSE_MC[1][4] = 0.000242939;
	MUON_ISO_UNC_LOOSE_MC[1][5] = 0.00028389;
	MUON_ISO_UNC_LOOSE_MC[1][6] = 0.00033901;
	MUON_ISO_UNC_HIGHPT_MC[1][0] = 0.00277066;
	MUON_ISO_UNC_HIGHPT_MC[1][1] = 0.00284822;
	MUON_ISO_UNC_HIGHPT_MC[1][2] = 0.00251193;
	MUON_ISO_UNC_HIGHPT_MC[1][3] = 0.000309191;
	MUON_ISO_UNC_HIGHPT_MC[1][4] = 0.000246851;
	MUON_ISO_UNC_HIGHPT_MC[1][5] = 0.000304037;
	MUON_ISO_UNC_HIGHPT_MC[1][6] = 0.000344613;
	MUON_TRIGGER_UNC_LOOSE_MC[1][0] = 2.54147e-05;
	MUON_TRIGGER_UNC_LOOSE_MC[1][1] = 0.000257801;
	MUON_TRIGGER_UNC_LOOSE_MC[1][2] = 0.00102941;
	MUON_TRIGGER_UNC_LOOSE_MC[1][3] = 0.000762124;
	MUON_TRIGGER_UNC_LOOSE_MC[1][4] = 0.000993233;
	MUON_TRIGGER_UNC_LOOSE_MC[1][5] = 0.00082869;
	MUON_TRIGGER_UNC_LOOSE_MC[1][6] = 0.00160749;
	MUON_TRIGGER_UNC_LOOSE_MC[1][7] = 0.00371748;
	MUON_TRIGGER_UNC_LOOSE_MC[1][8] = 0.0104518;
	MUON_TRIGGER_UNC_HIGHPT_MC[1][0] = 2.59705e-05;
	MUON_TRIGGER_UNC_HIGHPT_MC[1][1] = 0.000336764;
	MUON_TRIGGER_UNC_HIGHPT_MC[1][2] = 0.00104061;
	MUON_TRIGGER_UNC_HIGHPT_MC[1][3] = 0.00142153;
	MUON_TRIGGER_UNC_HIGHPT_MC[1][4] = 0.000862553;
	MUON_TRIGGER_UNC_HIGHPT_MC[1][5] = 0.000715838;
	MUON_TRIGGER_UNC_HIGHPT_MC[1][6] = 0.0016116;
	MUON_TRIGGER_UNC_HIGHPT_MC[1][7] = 0.00429441;
	MUON_TRIGGER_UNC_HIGHPT_MC[1][8] = 0.009853;
	MUON_ID_UNC_LOOSE_DATA[2][0] = 9.47267e-09;
	MUON_ID_UNC_LOOSE_DATA[2][1] = 0.00188671;
	MUON_ID_UNC_LOOSE_DATA[2][2] = 0.000816311;
	MUON_ID_UNC_LOOSE_DATA[2][3] = 0.000284072;
	MUON_ID_UNC_LOOSE_DATA[2][4] = 0.000141374;
	MUON_ID_UNC_LOOSE_DATA[2][5] = 0.00120885;
	MUON_ID_UNC_LOOSE_DATA[2][6] = 0.00112137;
	MUON_ID_UNC_HIGHPT_DATA[2][0] = 4.17611e-05;
	MUON_ID_UNC_HIGHPT_DATA[2][1] = 0.00372222;
	MUON_ID_UNC_HIGHPT_DATA[2][2] = 0.00219508;
	MUON_ID_UNC_HIGHPT_DATA[2][3] = 0.000328664;
	MUON_ID_UNC_HIGHPT_DATA[2][4] = 0.000367534;
	MUON_ID_UNC_HIGHPT_DATA[2][5] = 0.000699183;
	MUON_ID_UNC_HIGHPT_DATA[2][6] = 0.00283879;
	MUON_ISO_UNC_LOOSE_DATA[2][0] = 0.00138848;
	MUON_ISO_UNC_LOOSE_DATA[2][1] = 0.00269604;
	MUON_ISO_UNC_LOOSE_DATA[2][2] = 0.00145499;
	MUON_ISO_UNC_LOOSE_DATA[2][3] = 0.000290477;
	MUON_ISO_UNC_LOOSE_DATA[2][4] = 0.000101933;
	MUON_ISO_UNC_LOOSE_DATA[2][5] = 0.000119841;
	MUON_ISO_UNC_LOOSE_DATA[2][6] = 0.000182267;
	MUON_ISO_UNC_HIGHPT_DATA[2][0] = 0.00132633;
	MUON_ISO_UNC_HIGHPT_DATA[2][1] = 0.00270745;
	MUON_ISO_UNC_HIGHPT_DATA[2][2] = 0.00148272;
	MUON_ISO_UNC_HIGHPT_DATA[2][3] = 0.000296092;
	MUON_ISO_UNC_HIGHPT_DATA[2][4] = 0.000330449;
	MUON_ISO_UNC_HIGHPT_DATA[2][5] = 0.000122285;
	MUON_ISO_UNC_HIGHPT_DATA[2][6] = 0.000175753;
	// MUON_TRIGGER_UNC_LOOSE_DATA[2][0] = 7.29044e-05;
	// MUON_TRIGGER_UNC_LOOSE_DATA[2][1] = 0.000186464;
	// MUON_TRIGGER_UNC_LOOSE_DATA[2][2] = 0.000804265;
	// MUON_TRIGGER_UNC_LOOSE_DATA[2][3] = 0.00111985;
	// MUON_TRIGGER_UNC_LOOSE_DATA[2][4] = 0.00115773;
	// MUON_TRIGGER_UNC_LOOSE_DATA[2][5] = 0.000846352;
	// MUON_TRIGGER_UNC_LOOSE_DATA[2][6] = 0.00110357;
	// MUON_TRIGGER_UNC_LOOSE_DATA[2][7] = 0.00278552;
	// MUON_TRIGGER_UNC_LOOSE_DATA[2][8] = 0.0111735;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[2][0] = 4.45533e-05;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[2][1] = 0.000198784;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[2][2] = 0.000802738;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[2][3] = 0.000732562;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[2][4] = 0.000691983;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[2][5] = 0.000555259;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[2][6] = 0.00105095;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[2][7] = 0.00271737;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[2][8] = 0.0107173;
	MUON_ID_UNC_LOOSE_MC[2][0] = 1.77661e-07;
	MUON_ID_UNC_LOOSE_MC[2][1] = 0.000294365;
	MUON_ID_UNC_LOOSE_MC[2][2] = 0.000153194;
	MUON_ID_UNC_LOOSE_MC[2][3] = 6.33398e-05;
	MUON_ID_UNC_LOOSE_MC[2][4] = 3.59452e-05;
	MUON_ID_UNC_LOOSE_MC[2][5] = 0.000124879;
	MUON_ID_UNC_LOOSE_MC[2][6] = 0.000476615;
	MUON_ID_UNC_HIGHPT_MC[2][0] = 0.000858046;
	MUON_ID_UNC_HIGHPT_MC[2][1] = 0.000359223;
	MUON_ID_UNC_HIGHPT_MC[2][2] = 0.000218803;
	MUON_ID_UNC_HIGHPT_MC[2][3] = 0.000107999;
	MUON_ID_UNC_HIGHPT_MC[2][4] = 9.36772e-05;
	MUON_ID_UNC_HIGHPT_MC[2][5] = 0.000236184;
	MUON_ID_UNC_HIGHPT_MC[2][6] = 0.001035;
	MUON_ISO_UNC_LOOSE_MC[2][0] = 0.00161753;
	MUON_ISO_UNC_LOOSE_MC[2][1] = 0.00315383;
	MUON_ISO_UNC_LOOSE_MC[2][2] = 0.00182294;
	MUON_ISO_UNC_LOOSE_MC[2][3] = 0.000260295;
	MUON_ISO_UNC_LOOSE_MC[2][4] = 0.000316687;
	MUON_ISO_UNC_LOOSE_MC[2][5] = 0.000274337;
	MUON_ISO_UNC_LOOSE_MC[2][6] = 0.000214955;
	MUON_ISO_UNC_HIGHPT_MC[2][0] = 0.00161359;
	MUON_ISO_UNC_HIGHPT_MC[2][1] = 0.00316236;
	MUON_ISO_UNC_HIGHPT_MC[2][2] = 0.00181962;
	MUON_ISO_UNC_HIGHPT_MC[2][3] = 0.000260032;
	MUON_ISO_UNC_HIGHPT_MC[2][4] = 0.000126016;
	MUON_ISO_UNC_HIGHPT_MC[2][5] = 0.000202768;
	MUON_ISO_UNC_HIGHPT_MC[2][6] = 0.000215184;
	MUON_TRIGGER_UNC_LOOSE_MC[2][0] = 0.000128333;
	MUON_TRIGGER_UNC_LOOSE_MC[2][1] = 0.000104633;
	MUON_TRIGGER_UNC_LOOSE_MC[2][2] = 0.000975013;
	MUON_TRIGGER_UNC_LOOSE_MC[2][3] = 0.000853799;
	MUON_TRIGGER_UNC_LOOSE_MC[2][4] = 0.000982831;
	MUON_TRIGGER_UNC_LOOSE_MC[2][5] = 0.000970405;
	MUON_TRIGGER_UNC_LOOSE_MC[2][6] = 0.00175532;
	MUON_TRIGGER_UNC_LOOSE_MC[2][7] = 0.00510913;
	MUON_TRIGGER_UNC_LOOSE_MC[2][8] = 0.0159022;
	MUON_TRIGGER_UNC_HIGHPT_MC[2][0] = 1.97961e-05;
	MUON_TRIGGER_UNC_HIGHPT_MC[2][1] = 0.000382968;
	MUON_TRIGGER_UNC_HIGHPT_MC[2][2] = 0.00101425;
	MUON_TRIGGER_UNC_HIGHPT_MC[2][3] = 0.000798395;
	MUON_TRIGGER_UNC_HIGHPT_MC[2][4] = 0.000921592;
	MUON_TRIGGER_UNC_HIGHPT_MC[2][5] = 0.000933577;
	MUON_TRIGGER_UNC_HIGHPT_MC[2][6] = 0.00175327;
	MUON_TRIGGER_UNC_HIGHPT_MC[2][7] = 0.00454474;
	MUON_TRIGGER_UNC_HIGHPT_MC[2][8] = 0.0164573;
	MUON_ID_UNC_LOOSE_DATA[3][0] = 0.000327127;
	MUON_ID_UNC_LOOSE_DATA[3][1] = 0.00248376;
	MUON_ID_UNC_LOOSE_DATA[3][2] = 0.00106578;
	MUON_ID_UNC_LOOSE_DATA[3][3] = 0.000395573;
	MUON_ID_UNC_LOOSE_DATA[3][4] = 0.000425719;
	MUON_ID_UNC_LOOSE_DATA[3][5] = 0.00166202;
	MUON_ID_UNC_LOOSE_DATA[3][6] = 0.00478026;
	MUON_ID_UNC_HIGHPT_DATA[3][0] = 0.00760157;
	MUON_ID_UNC_HIGHPT_DATA[3][1] = 0.0033534;
	MUON_ID_UNC_HIGHPT_DATA[3][2] = 0.00126964;
	MUON_ID_UNC_HIGHPT_DATA[3][3] = 0.000671425;
	MUON_ID_UNC_HIGHPT_DATA[3][4] = 0.00120291;
	MUON_ID_UNC_HIGHPT_DATA[3][5] = 0.00264803;
	MUON_ID_UNC_HIGHPT_DATA[3][6] = 0.0063471;
	MUON_ISO_UNC_LOOSE_DATA[3][0] = 0.00165665;
	MUON_ISO_UNC_LOOSE_DATA[3][1] = 0.00174742;
	MUON_ISO_UNC_LOOSE_DATA[3][2] = 0.00116005;
	MUON_ISO_UNC_LOOSE_DATA[3][3] = 0.000183353;
	MUON_ISO_UNC_LOOSE_DATA[3][4] = 0.000127166;
	MUON_ISO_UNC_LOOSE_DATA[3][5] = 0.0019833;
	MUON_ISO_UNC_LOOSE_DATA[3][6] = 0.000349701;
	MUON_ISO_UNC_HIGHPT_DATA[3][0] = 0.00170625;
	MUON_ISO_UNC_HIGHPT_DATA[3][1] = 0.00181967;
	MUON_ISO_UNC_HIGHPT_DATA[3][2] = 0.00119516;
	MUON_ISO_UNC_HIGHPT_DATA[3][3] = 0.000190283;
	MUON_ISO_UNC_HIGHPT_DATA[3][4] = 0.000126766;
	MUON_ISO_UNC_HIGHPT_DATA[3][5] = 0.00192279;
	MUON_ISO_UNC_HIGHPT_DATA[3][6] = 0.000346093;
	// MUON_TRIGGER_UNC_LOOSE_DATA[3][0] = 0.000272809;
	// MUON_TRIGGER_UNC_LOOSE_DATA[3][1] = 0.000519751;
	// MUON_TRIGGER_UNC_LOOSE_DATA[3][2] = 0.00183915;
	// MUON_TRIGGER_UNC_LOOSE_DATA[3][3] = 0.00153492;
	// MUON_TRIGGER_UNC_LOOSE_DATA[3][4] = 0.00188493;
	// MUON_TRIGGER_UNC_LOOSE_DATA[3][5] = 0.00166369;
	// MUON_TRIGGER_UNC_LOOSE_DATA[3][6] = 0.00498534;
	// MUON_TRIGGER_UNC_LOOSE_DATA[3][7] = 0.0120593;
	// MUON_TRIGGER_UNC_LOOSE_DATA[3][8] = 0.113365;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[3][0] = 0.000254734;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[3][1] = 0.000498787;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[3][2] = 0.00180377;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[3][3] = 0.00150755;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[3][4] = 0.001639;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[3][5] = 0.00167261;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[3][6] = 0.0041381;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[3][7] = 0.0144408;
	// MUON_TRIGGER_UNC_HIGHPT_DATA[3][8] = 0.134977;
	MUON_ID_UNC_LOOSE_MC[3][0] = 1.70338e-07;
	MUON_ID_UNC_LOOSE_MC[3][1] = 0.0005476;
	MUON_ID_UNC_LOOSE_MC[3][2] = 0.000431454;
	MUON_ID_UNC_LOOSE_MC[3][3] = 0.000149121;
	MUON_ID_UNC_LOOSE_MC[3][4] = 9.88745e-05;
	MUON_ID_UNC_LOOSE_MC[3][5] = 0.000328711;
	MUON_ID_UNC_LOOSE_MC[3][6] = 0.000101169;
	MUON_ID_UNC_HIGHPT_MC[3][0] = 0.00135165;
	MUON_ID_UNC_HIGHPT_MC[3][1] = 0.000754565;
	MUON_ID_UNC_HIGHPT_MC[3][2] = 0.000650332;
	MUON_ID_UNC_HIGHPT_MC[3][3] = 0.000630187;
	MUON_ID_UNC_HIGHPT_MC[3][4] = 0.000336806;
	MUON_ID_UNC_HIGHPT_MC[3][5] = 0.00105996;
	MUON_ID_UNC_HIGHPT_MC[3][6] = 0.001947;
	MUON_ISO_UNC_LOOSE_MC[3][0] = 0.00192606;
	MUON_ISO_UNC_LOOSE_MC[3][1] = 0.00274209;
	MUON_ISO_UNC_LOOSE_MC[3][2] = 0.00173103;
	MUON_ISO_UNC_LOOSE_MC[3][3] = 0.000264847;
	MUON_ISO_UNC_LOOSE_MC[3][4] = 0.000174817;
	MUON_ISO_UNC_LOOSE_MC[3][5] = 0.00034236;
	MUON_ISO_UNC_LOOSE_MC[3][6] = 0.000497012;
	MUON_ISO_UNC_HIGHPT_MC[3][0] = 0.00193697;
	MUON_ISO_UNC_HIGHPT_MC[3][1] = 0.00269202;
	MUON_ISO_UNC_HIGHPT_MC[3][2] = 0.0017377;
	MUON_ISO_UNC_HIGHPT_MC[3][3] = 0.000269272;
	MUON_ISO_UNC_HIGHPT_MC[3][4] = 0.000135044;
	MUON_ISO_UNC_HIGHPT_MC[3][5] = 0.00032112;
	MUON_ISO_UNC_HIGHPT_MC[3][6] = 0.000481024;
	MUON_TRIGGER_UNC_LOOSE_MC[3][0] = 0.000117325;
	MUON_TRIGGER_UNC_LOOSE_MC[3][1] = 0.000365343;
	MUON_TRIGGER_UNC_LOOSE_MC[3][2] = 0.00246411;
	MUON_TRIGGER_UNC_LOOSE_MC[3][3] = 0.00344048;
	MUON_TRIGGER_UNC_LOOSE_MC[3][4] = 0.00612802;
	MUON_TRIGGER_UNC_LOOSE_MC[3][5] = 0.00228855;
	MUON_TRIGGER_UNC_LOOSE_MC[3][6] = 0.00531149;
	MUON_TRIGGER_UNC_LOOSE_MC[3][7] = 0.0190855;
	MUON_TRIGGER_UNC_LOOSE_MC[3][8] = 0.0290383;
	MUON_TRIGGER_UNC_HIGHPT_MC[3][0] = 9.15683e-05;
	MUON_TRIGGER_UNC_HIGHPT_MC[3][1] = 0.000324622;
	MUON_TRIGGER_UNC_HIGHPT_MC[3][2] = 0.0024508;
	MUON_TRIGGER_UNC_HIGHPT_MC[3][3] = 0.00399178;
	MUON_TRIGGER_UNC_HIGHPT_MC[3][4] = 0.00213123;
	MUON_TRIGGER_UNC_HIGHPT_MC[3][5] = 0.00323214;
	MUON_TRIGGER_UNC_HIGHPT_MC[3][6] = 0.00514845;
	MUON_TRIGGER_UNC_HIGHPT_MC[3][7] = 0.0183984;
	MUON_TRIGGER_UNC_HIGHPT_MC[3][8] = 0.029857;

	MUON_TRIGGER_UNC_LOOSE_DATA[0][0] = 0.00025771;
	MUON_TRIGGER_UNC_LOOSE_DATA[0][1] = 0.0112463;
	MUON_TRIGGER_UNC_LOOSE_DATA[0][2] = 0.0545832;
	MUON_TRIGGER_UNC_LOOSE_DATA[0][3] = 0.0312646;
	MUON_TRIGGER_UNC_LOOSE_DATA[0][4] = 0.022554;
	MUON_TRIGGER_UNC_LOOSE_DATA[0][5] = 0.0194172;
	MUON_TRIGGER_UNC_LOOSE_DATA[0][6] = 0.0205655;
	MUON_TRIGGER_UNC_LOOSE_DATA[0][7] = 0.0175683;
	MUON_TRIGGER_UNC_LOOSE_DATA[0][8] = 0.0241179;
	MUON_TRIGGER_UNC_HIGHPT_DATA[0][0] = 0.00024913;
	MUON_TRIGGER_UNC_HIGHPT_DATA[0][1] = 0.0158308;
	MUON_TRIGGER_UNC_HIGHPT_DATA[0][2] = 0.0522997;
	MUON_TRIGGER_UNC_HIGHPT_DATA[0][3] = 0.0280362;
	MUON_TRIGGER_UNC_HIGHPT_DATA[0][4] = 0.019869;
	MUON_TRIGGER_UNC_HIGHPT_DATA[0][5] = 0.0178436;
	MUON_TRIGGER_UNC_HIGHPT_DATA[0][6] = 0.0187292;
	MUON_TRIGGER_UNC_HIGHPT_DATA[0][7] = 0.0171203;
	MUON_TRIGGER_UNC_HIGHPT_DATA[0][8] = 0.0210967;
	MUON_TRIGGER_UNC_LOOSE_DATA[1][0] = 0.000438138;
	MUON_TRIGGER_UNC_LOOSE_DATA[1][1] = 0.00827235;
	MUON_TRIGGER_UNC_LOOSE_DATA[1][2] = 0.0524941;
	MUON_TRIGGER_UNC_LOOSE_DATA[1][3] = 0.0354184;
	MUON_TRIGGER_UNC_LOOSE_DATA[1][4] = 0.0306967;
	MUON_TRIGGER_UNC_LOOSE_DATA[1][5] = 0.0258783;
	MUON_TRIGGER_UNC_LOOSE_DATA[1][6] = 0.0256991;
	MUON_TRIGGER_UNC_LOOSE_DATA[1][7] = 0.0283698;
	MUON_TRIGGER_UNC_LOOSE_DATA[1][8] = 0.0362159;
	MUON_TRIGGER_UNC_HIGHPT_DATA[1][0] = 0.000438513;
	MUON_TRIGGER_UNC_HIGHPT_DATA[1][1] = 0.00856556;
	MUON_TRIGGER_UNC_HIGHPT_DATA[1][2] = 0.0498327;
	MUON_TRIGGER_UNC_HIGHPT_DATA[1][3] = 0.0320156;
	MUON_TRIGGER_UNC_HIGHPT_DATA[1][4] = 0.0275741;
	MUON_TRIGGER_UNC_HIGHPT_DATA[1][5] = 0.0240026;
	MUON_TRIGGER_UNC_HIGHPT_DATA[1][6] = 0.0241727;
	MUON_TRIGGER_UNC_HIGHPT_DATA[1][7] = 0.0254632;
	MUON_TRIGGER_UNC_HIGHPT_DATA[1][8] = 0.0377162;
	MUON_TRIGGER_UNC_LOOSE_DATA[2][0] = 0.0013083;
	MUON_TRIGGER_UNC_LOOSE_DATA[2][1] = 0.0159401;
	MUON_TRIGGER_UNC_LOOSE_DATA[2][2] = 0.0407275;
	MUON_TRIGGER_UNC_LOOSE_DATA[2][3] = 0.0247242;
	MUON_TRIGGER_UNC_LOOSE_DATA[2][4] = 0.018929;
	MUON_TRIGGER_UNC_LOOSE_DATA[2][5] = 0.0150843;
	MUON_TRIGGER_UNC_LOOSE_DATA[2][6] = 0.0158304;
	MUON_TRIGGER_UNC_LOOSE_DATA[2][7] = 0.0230819;
	MUON_TRIGGER_UNC_LOOSE_DATA[2][8] = 0.0237871;
	MUON_TRIGGER_UNC_HIGHPT_DATA[2][0] = 0.00125836;
	MUON_TRIGGER_UNC_HIGHPT_DATA[2][1] = 0.015794;
	MUON_TRIGGER_UNC_HIGHPT_DATA[2][2] = 0.0396372;
	MUON_TRIGGER_UNC_HIGHPT_DATA[2][3] = 0.023034;
	MUON_TRIGGER_UNC_HIGHPT_DATA[2][4] = 0.0169239;
	MUON_TRIGGER_UNC_HIGHPT_DATA[2][5] = 0.014952;
	MUON_TRIGGER_UNC_HIGHPT_DATA[2][6] = 0.0152026;
	MUON_TRIGGER_UNC_HIGHPT_DATA[2][7] = 0.0230045;
	MUON_TRIGGER_UNC_HIGHPT_DATA[2][8] = 0.0190414;
	MUON_TRIGGER_UNC_LOOSE_DATA[3][0] = 0.00767607;
	MUON_TRIGGER_UNC_LOOSE_DATA[3][1] = 0.0289149;
	MUON_TRIGGER_UNC_LOOSE_DATA[3][2] = 0.0584159;
	MUON_TRIGGER_UNC_LOOSE_DATA[3][3] = 0.0571069;
	MUON_TRIGGER_UNC_LOOSE_DATA[3][4] = 0.0421992;
	MUON_TRIGGER_UNC_LOOSE_DATA[3][5] = 0.0302835;
	MUON_TRIGGER_UNC_LOOSE_DATA[3][6] = 0.0232046;
	MUON_TRIGGER_UNC_LOOSE_DATA[3][7] = 0.0506215;
	MUON_TRIGGER_UNC_LOOSE_DATA[3][8] = 0.134946;
	MUON_TRIGGER_UNC_HIGHPT_DATA[3][0] = 0.00790248;
	MUON_TRIGGER_UNC_HIGHPT_DATA[3][1] = 0.0298159;
	MUON_TRIGGER_UNC_HIGHPT_DATA[3][2] = 0.0585101;
	MUON_TRIGGER_UNC_HIGHPT_DATA[3][3] = 0.0575135;
	MUON_TRIGGER_UNC_HIGHPT_DATA[3][4] = 0.0414995;
	MUON_TRIGGER_UNC_HIGHPT_DATA[3][5] = 0.0285501;
	MUON_TRIGGER_UNC_HIGHPT_DATA[3][6] = 0.0199034;
	MUON_TRIGGER_UNC_HIGHPT_DATA[3][7] = 0.0462357;
	MUON_TRIGGER_UNC_HIGHPT_DATA[3][8] = 0.176572;
}

double ScaleFactor::getElectronScaleFactor(const NtupleElectron& electron) {
	return getElectronScaleFactorForID(electron) * getElectronScaleFactorForRECO(electron);
}


double ScaleFactor::getElectronScaleFactorForSystematics(const ElectronSelector* ee) {
	// double eff0_ID = getElectronScaleFactorForID(ee->at(0));
	// double eff0_RECO = getElectronScaleFactorForRECO(ee->at(0));
	double unc0_ID = getElectronScaleFactorUncForID(ee->at(0));
	double unc0_RECO = getElectronScaleFactorUncForRECO(ee->at(0));
	// double eff1_ID = getElectronScaleFactorForID(ee->at(1));
	// double eff1_RECO = getElectronScaleFactorForRECO(ee->at(1));
	double unc1_ID = getElectronScaleFactorUncForID(ee->at(1));
	double unc1_RECO = getElectronScaleFactorUncForRECO(ee->at(1));

	return 1. + sqrt( unc0_ID*unc0_ID + unc1_ID*unc1_ID + unc0_RECO*unc0_RECO + unc1_RECO*unc1_RECO );
}

double ScaleFactor::getElectronScaleFactorForID(const NtupleElectron& electron) {
	int nEtaBins = electronNEtaBinsForID;
	int nPtBins = electronNPtBinsForID;

	for(auto etaIdx=0; etaIdx<nEtaBins; ++etaIdx){
		if( electron.etaSC >= ELECTRON_ETABINS_ID[etaIdx] && electron.etaSC < ELECTRON_ETABINS_ID[etaIdx+1] ) {
			for(auto ptIdx=0; ptIdx<nPtBins; ++ptIdx) {
				if( (electron.pt >= ELECTRON_PTBINS_ID[ptIdx] && electron.pt < ELECTRON_PTBINS_ID[ptIdx+1]) || ptIdx == nPtBins-1 ) return ELECTRON_SF_ID[etaIdx][ptIdx];
			}
		}
	}
	return -999.;
}

double ScaleFactor::getElectronScaleFactorUncForID(const NtupleElectron& electron) {
	int nEtaBins = electronNEtaBinsForID;
	int nPtBins = electronNPtBinsForID;

	for(auto etaIdx=0; etaIdx<nEtaBins; ++etaIdx){
		if( electron.etaSC >= ELECTRON_ETABINS_ID[etaIdx] && electron.etaSC < ELECTRON_ETABINS_ID[etaIdx+1] ) {
			for(auto ptIdx=0; ptIdx<nPtBins; ++ptIdx) {
				if( (electron.pt >= ELECTRON_PTBINS_ID[ptIdx] && electron.pt < ELECTRON_PTBINS_ID[ptIdx+1]) || ptIdx == nPtBins-1 ) return ELECTRON_SF_UNC_ID[etaIdx][ptIdx]/ELECTRON_SF_ID[etaIdx][ptIdx];
			}
		}
	}
	return -999.;
}

double ScaleFactor::getElectronScaleFactorForRECO(const NtupleElectron& electron) {
	int nEtaBins = electronNEtaBinsForRECO;

	for(auto etaIdx=0; etaIdx<nEtaBins; ++etaIdx){
		if( electron.etaSC >= ELECTRON_ETABINS_RECO[etaIdx] && electron.etaSC < ELECTRON_ETABINS_RECO[etaIdx+1] ) {
			return ELECTRON_SF_RECO[etaIdx];
		}
	}
	return -999.;
}

double ScaleFactor::getElectronScaleFactorUncForRECO(const NtupleElectron& electron) {
	int nEtaBins = electronNEtaBinsForRECO;

	for(auto etaIdx=0; etaIdx<nEtaBins; ++etaIdx){
		if( electron.etaSC >= ELECTRON_ETABINS_RECO[etaIdx] && electron.etaSC < ELECTRON_ETABINS_RECO[etaIdx+1] ) {
			return (electron.pt < 20 || electron.pt > 80) ? ELECTRON_SF_UNC_RECO[etaIdx]/ELECTRON_SF_RECO[etaIdx] + 0.01 : ELECTRON_SF_UNC_RECO[etaIdx]/ELECTRON_SF_RECO[etaIdx];
		}
	}
	return -999.;
}

double ScaleFactor::getPhotonScaleFactor(const NtuplePhoton& pho) {
	int nEtaBins = photonNEtaBins;
	int nPtBins = photonNPtBins;

	for(auto etaIdx=0; etaIdx<nEtaBins; ++etaIdx){
		if( pho.etaSC >= PHOTON_ETABINS[etaIdx] && pho.etaSC < PHOTON_ETABINS[etaIdx+1] ) {
			for(auto ptIdx=0; ptIdx<nPtBins; ++ptIdx) {
				if( (pho.pt >= PHOTON_PTBINS[ptIdx] && pho.pt < PHOTON_PTBINS[ptIdx+1]) || ptIdx == nPtBins-1 ) return PHOTON_SF[etaIdx][ptIdx];
			}
		}
	}
	return -999.;
}


double ScaleFactor::getPhotonScaleFactorForSystematics(const NtuplePhoton& pho, const bool& plusSigma) {
	int nEtaBins = photonNEtaBins;
	int nPtBins = photonNPtBins;

	for(auto etaIdx=0; etaIdx<nEtaBins; ++etaIdx){
		if( pho.etaSC >= PHOTON_ETABINS[etaIdx] && pho.etaSC < PHOTON_ETABINS[etaIdx+1] ) {
			for(auto ptIdx=0; ptIdx<nPtBins; ++ptIdx) {
				if( (pho.pt >= PHOTON_PTBINS[ptIdx] && pho.pt < PHOTON_PTBINS[ptIdx+1]) || ptIdx == nPtBins-1 ) {
					if( pho.pt < PHOTON_PTBINS[nPtBins] ) {
						if(plusSigma) return (1. + PHOTON_SF_UNC[etaIdx][ptIdx] / PHOTON_SF[etaIdx][ptIdx]);
						else return (1. - PHOTON_SF_UNC[etaIdx][ptIdx] / PHOTON_SF[etaIdx][ptIdx]);
					}
					else {
						if(plusSigma) return 1.05;
						else return 0.95;
					}
				}
			}
			
		}
	}
	return -999.;
}

double ScaleFactor::getMuonScaleFactors(const MuonSelector* mm) {

	std::vector<double> effsForIDIso1 = getMuonEffsForIDIso(mm->at(0));
	std::vector<double> effsForIDIso2 = getMuonEffsForIDIso(mm->at(1));

	assert( effsForIDIso1.size()==4 && effsForIDIso2.size()==4 );

	double dataForIDIsoLoose1 = effsForIDIso1[0];
	double mcForIDIsoLoose1 = effsForIDIso1[1];
	double dataForIDIsoHighPt1 = effsForIDIso1[2];
	double mcForIDIsoHighPt1 = effsForIDIso1[3];

	double dataForIDIsoLoose2 = effsForIDIso2[0];
	double mcForIDIsoLoose2 = effsForIDIso2[1];
	double dataForIDIsoHighPt2 = effsForIDIso2[2];
	double mcForIDIsoHighPt2 = effsForIDIso2[3];

	// std::cout << "IDIso:" << std::endl;
	// std::cout << dataForIDIsoLoose1 << std::endl;
	// std::cout << mcForIDIsoLoose1 << std::endl;
	// std::cout << dataForIDIsoHighPt1 << std::endl;
	// std::cout << mcForIDIsoHighPt1 << std::endl;
	// std::cout << dataForIDIsoLoose2 << std::endl;
	// std::cout << mcForIDIsoLoose2 << std::endl;
	// std::cout << dataForIDIsoHighPt2 << std::endl;
	// std::cout << mcForIDIsoHighPt2 << std::endl;

	std::vector<double> effsForTrigger1 = getMuonEffsForTrigger(mm->at(0));
	std::vector<double> effsForTrigger2 = getMuonEffsForTrigger(mm->at(1));

	assert( effsForTrigger1.size()==4 && effsForTrigger2.size()==4 );

	double dataForTriggerLoose1 = effsForTrigger1[0];
	double mcForTriggerLoose1 = effsForTrigger1[1];
	double dataForTriggerHighPt1 = effsForTrigger1[2];
	double mcForTriggerHighPt1 = effsForTrigger1[3];

	double dataForTriggerLoose2 = effsForTrigger2[0];
	double mcForTriggerLoose2 = effsForTrigger2[1];
	double dataForTriggerHighPt2 = effsForTrigger2[2];
	double mcForTriggerHighPt2 = effsForTrigger2[3];

	// std::cout << "Trigger:" << std::endl;
	// std::cout << dataForTriggerLoose1 << std::endl;
	// std::cout << mcForTriggerLoose1 << std::endl;
	// std::cout << dataForTriggerHighPt1 << std::endl;
	// std::cout << mcForTriggerHighPt1 << std::endl;
	// std::cout << dataForTriggerLoose2 << std::endl;
	// std::cout << mcForTriggerLoose2 << std::endl;
	// std::cout << dataForTriggerHighPt2 << std::endl;
	// std::cout << mcForTriggerHighPt2 << std::endl;

	double dataForLoose1HighPt2 = dataForIDIsoLoose1 * dataForIDIsoHighPt2 * ( 1. - (1.-dataForTriggerLoose1)*(1.-dataForTriggerHighPt2) );
	double dataForHighPt1Loose2 = dataForIDIsoHighPt1 * dataForIDIsoLoose2 * ( 1. - (1.-dataForTriggerHighPt1)*(1.-dataForTriggerLoose2) );
	double dataForBothHigh = dataForIDIsoHighPt1 * dataForIDIsoHighPt2 * ( 1. - (1.-dataForTriggerHighPt1)*(1.-dataForTriggerHighPt2) );
	double mcForLoose1HighPt2 = mcForIDIsoLoose1 * mcForIDIsoHighPt2 * ( 1. - (1.-mcForTriggerLoose1)*(1.-mcForTriggerHighPt2) );
	double mcForHighPt1Loose2 = mcForIDIsoHighPt1 * mcForIDIsoLoose2 * ( 1. - (1.-mcForTriggerHighPt1)*(1.-mcForTriggerLoose2) );
	double mcForBothHigh = mcForIDIsoHighPt1 * mcForIDIsoHighPt2 * ( 1. - (1.-mcForTriggerHighPt1)*(1.-mcForTriggerHighPt2) );

	// std::cout << (dataForLoose1HighPt2+dataForHighPt1Loose2-dataForBothHigh)/(mcForLoose1HighPt2+mcForHighPt1Loose2-mcForBothHigh) << std::endl;
	return (dataForLoose1HighPt2+dataForHighPt1Loose2-dataForBothHigh)/(mcForLoose1HighPt2+mcForHighPt1Loose2-mcForBothHigh);
}

std::vector<double> ScaleFactor::getMuonEffsForIDIso(const NtupleMuon& mu) {
	int nAbsetaBins = ScaleFactor::muonNAbsetaBins;
	int nPtBins = ScaleFactor::muonNPtBinsForIDIso;

	int absetaIdx = 0, ptIdx = 0;
	for(; absetaIdx<nAbsetaBins; ++absetaIdx){
		if( fabs(mu.eta) >= MUON_ABSETABINS[absetaIdx] && fabs(mu.eta) < MUON_ABSETABINS[absetaIdx+1] ) {
			for(; ptIdx<nPtBins; ++ptIdx) {
				if( (mu.pt >= MUON_PTBINS_IDISO[ptIdx] && mu.pt < MUON_PTBINS_IDISO[ptIdx+1]) || ptIdx == nPtBins-1 ) {
					return { MUON_ID_LOOSE_DATA[absetaIdx][ptIdx]*MUON_ISO_LOOSE_DATA[absetaIdx][ptIdx],
					         MUON_ID_LOOSE_MC[absetaIdx][ptIdx]*MUON_ISO_LOOSE_MC[absetaIdx][ptIdx],
					         MUON_ID_HIGHPT_DATA[absetaIdx][ptIdx]*MUON_ISO_HIGHPT_DATA[absetaIdx][ptIdx],
					         MUON_ID_HIGHPT_MC[absetaIdx][ptIdx]*MUON_ISO_HIGHPT_MC[absetaIdx][ptIdx] };
				}
			}
		}
	}
	return {};
}

std::vector<double> ScaleFactor::getMuonEffsForTrigger(const NtupleMuon& mu) {
	int nAbsetaBins = ScaleFactor::muonNAbsetaBins;
	int nPtBins = ScaleFactor::muonNPtBinsForTrigger;

	int absetaIdx = 0, ptIdx = 0;
	for(; absetaIdx<nAbsetaBins; ++absetaIdx){
		if( fabs(mu.eta) >= MUON_ABSETABINS[absetaIdx] && fabs(mu.eta) < MUON_ABSETABINS[absetaIdx+1] ) {
			for(; ptIdx<nPtBins; ++ptIdx) {
				if( (mu.pt >= MUON_PTBINS_TRIGGER[ptIdx] && mu.pt < MUON_PTBINS_TRIGGER[ptIdx+1]) || ptIdx == nPtBins-1 ) {
					return { MUON_TRIGGER_LOOSE_DATA[absetaIdx][ptIdx], 
						     MUON_TRIGGER_LOOSE_MC[absetaIdx][ptIdx],
						     MUON_TRIGGER_HIGHPT_DATA[absetaIdx][ptIdx], 
						     MUON_TRIGGER_HIGHPT_MC[absetaIdx][ptIdx] };
				}
			}
		}
	}
	return {};
}

double ScaleFactor::getMuonScaleFactorVariations(const MuonSelector* mm, const bool& isIDSF) {

	std::vector<double> effsForIDIso1 = getMuonEffsForIDIso(mm->at(0));
	std::vector<double> effsForIDIso2 = getMuonEffsForIDIso(mm->at(1));

	std::vector<double> uncsForIDIso1 = getMuonEffUncsForIDIso(mm->at(0));
	std::vector<double> uncsForIDIso2 = getMuonEffUncsForIDIso(mm->at(1));

	assert( effsForIDIso1.size()==4 && effsForIDIso2.size()==4 );
	assert( uncsForIDIso1.size()==4 && uncsForIDIso2.size()==4 );

	double effDataForIDIsoLoose1  = effsForIDIso1[0];
	double effMCForIDIsoLoose1    = effsForIDIso1[1];
	double effDataForIDIsoHighPt1 = effsForIDIso1[2];
	double effMCForIDIsoHighPt1   = effsForIDIso1[3];

	double effDataForIDIsoLoose2  = effsForIDIso2[0];
	double effMCForIDIsoLoose2    = effsForIDIso2[1];
	double effDataForIDIsoHighPt2 = effsForIDIso2[2];
	double effMCForIDIsoHighPt2   = effsForIDIso2[3];

	double uncDataForIDIsoLoose1  = uncsForIDIso1[0];
	double uncMCForIDIsoLoose1    = uncsForIDIso1[1];
	double uncDataForIDIsoHighPt1 = uncsForIDIso1[2];
	double uncMCForIDIsoHighPt1   = uncsForIDIso1[3];

	double uncDataForIDIsoLoose2  = uncsForIDIso2[0];
	double uncMCForIDIsoLoose2    = uncsForIDIso2[1];
	double uncDataForIDIsoHighPt2 = uncsForIDIso2[2];
	double uncMCForIDIsoHighPt2   = uncsForIDIso2[3];

	std::vector<double> effsForTrigger1 = getMuonEffsForTrigger(mm->at(0));
	std::vector<double> effsForTrigger2 = getMuonEffsForTrigger(mm->at(1));

	std::vector<double> uncsForTrigger1 = getMuonEffUncsForTrigger(mm->at(0));
	std::vector<double> uncsForTrigger2 = getMuonEffUncsForTrigger(mm->at(1));

	assert( effsForTrigger1.size()==4 && effsForTrigger2.size()==4 );
	assert( uncsForTrigger1.size()==4 && uncsForTrigger2.size()==4 );

	double effDataForTriggerLoose1  = effsForTrigger1[0];
	double effMCForTriggerLoose1    = effsForTrigger1[1];
	double effDataForTriggerHighPt1 = effsForTrigger1[2];
	double effMCForTriggerHighPt1   = effsForTrigger1[3];

	double effDataForTriggerLoose2  = effsForTrigger2[0];
	double effMCForTriggerLoose2    = effsForTrigger2[1];
	double effDataForTriggerHighPt2 = effsForTrigger2[2];
	double effMCForTriggerHighPt2   = effsForTrigger2[3];

	double uncDataForTriggerLoose1  = uncsForTrigger1[0];
	double uncMCForTriggerLoose1    = uncsForTrigger1[1];
	double uncDataForTriggerHighPt1 = uncsForTrigger1[2];
	double uncMCForTriggerHighPt1   = uncsForTrigger1[3];

	double uncDataForTriggerLoose2  = uncsForTrigger2[0];
	double uncMCForTriggerLoose2    = uncsForTrigger2[1];
	double uncDataForTriggerHighPt2 = uncsForTrigger2[2];
	double uncMCForTriggerHighPt2   = uncsForTrigger2[3];

	double effDataForLoose1HighPt2 = effDataForIDIsoLoose1  * effDataForIDIsoHighPt2 * ( 1. - (1.-effDataForTriggerLoose1)  * (1.-effDataForTriggerHighPt2) );
	double effDataForHighPt1Loose2 = effDataForIDIsoHighPt1 * effDataForIDIsoLoose2  * ( 1. - (1.-effDataForTriggerHighPt1) * (1.-effDataForTriggerLoose2) );
	double effDataForBothHigh      = effDataForIDIsoHighPt1 * effDataForIDIsoHighPt2 * ( 1. - (1.-effDataForTriggerHighPt1) * (1.-effDataForTriggerHighPt2) );
	double effMCForLoose1HighPt2   = effMCForIDIsoLoose1    * effMCForIDIsoHighPt2   * ( 1. - (1.-effMCForTriggerLoose1)    * (1.-effMCForTriggerHighPt2) );
	double effMCForHighPt1Loose2   = effMCForIDIsoHighPt1   * effMCForIDIsoLoose2    * ( 1. - (1.-effMCForTriggerHighPt1)   * (1.-effMCForTriggerLoose2) );
	double effMCForBothHigh        = effMCForIDIsoHighPt1   * effMCForIDIsoHighPt2   * ( 1. - (1.-effMCForTriggerHighPt1)   * (1.-effMCForTriggerHighPt2) );

	// double uncDataForLoose1HighPt2 = sqrt( uncDataForIDIsoLoose1*uncDataForIDIsoLoose1   + uncDataForIDIsoHighPt2*uncDataForIDIsoHighPt2 + uncDataForTriggerLoose1*uncDataForTriggerLoose1   + uncDataForTriggerHighPt2*uncDataForTriggerHighPt2 );
	// double uncDataForHighPt1Loose2 = sqrt( uncDataForIDIsoHighPt1*uncDataForIDIsoHighPt1 + uncDataForIDIsoLoose2*uncDataForIDIsoLoose2   + uncDataForTriggerHighPt1*uncDataForTriggerHighPt1 + uncDataForTriggerLoose2*uncDataForTriggerLoose2 );
	// double uncDataForBothHigh      = sqrt( uncDataForIDIsoHighPt1*uncDataForIDIsoHighPt1 + uncDataForIDIsoHighPt2*uncDataForIDIsoHighPt2 + uncDataForTriggerHighPt1*uncDataForTriggerHighPt1 + uncDataForTriggerHighPt2*uncDataForTriggerHighPt2 );
	// double uncMCForLoose1HighPt2   = sqrt( uncMCForIDIsoLoose1*uncMCForIDIsoLoose1       + uncMCForIDIsoHighPt2*uncMCForIDIsoHighPt2     + uncMCForTriggerLoose1*uncMCForTriggerLoose1       + uncMCForTriggerHighPt2*uncMCForTriggerHighPt2 );
	// double uncMCForHighPt1Loose2   = sqrt( uncMCForIDIsoHighPt1*uncMCForIDIsoHighPt1     + uncMCForIDIsoLoose2*uncMCForIDIsoLoose2       + uncMCForTriggerHighPt1*uncMCForTriggerHighPt1     + uncMCForTriggerLoose2*uncMCForTriggerLoose2 );
	// double uncMCForBothHigh        = sqrt( uncMCForIDIsoHighPt1*uncMCForIDIsoHighPt1     + uncMCForIDIsoHighPt2*uncMCForIDIsoHighPt2     + uncMCForTriggerHighPt1*uncMCForTriggerHighPt1     + uncMCForTriggerHighPt2*uncMCForTriggerHighPt2 );

	double uncDataForLoose1HighPt2; 
	double uncDataForHighPt1Loose2;
	double uncDataForBothHigh;
	double uncMCForLoose1HighPt2;
	double uncMCForHighPt1Loose2;
	double uncMCForBothHigh;

	if(isIDSF) {
		uncDataForLoose1HighPt2 = sqrt( uncDataForIDIsoLoose1*uncDataForIDIsoLoose1   + uncDataForIDIsoHighPt2*uncDataForIDIsoHighPt2 );
		uncDataForHighPt1Loose2 = sqrt( uncDataForIDIsoHighPt1*uncDataForIDIsoHighPt1 + uncDataForIDIsoLoose2*uncDataForIDIsoLoose2   );
		uncDataForBothHigh      = sqrt( uncDataForIDIsoHighPt1*uncDataForIDIsoHighPt1 + uncDataForIDIsoHighPt2*uncDataForIDIsoHighPt2 );
		uncMCForLoose1HighPt2   = sqrt( uncMCForIDIsoLoose1*uncMCForIDIsoLoose1       + uncMCForIDIsoHighPt2*uncMCForIDIsoHighPt2     );
		uncMCForHighPt1Loose2   = sqrt( uncMCForIDIsoHighPt1*uncMCForIDIsoHighPt1     + uncMCForIDIsoLoose2*uncMCForIDIsoLoose2       );
		uncMCForBothHigh        = sqrt( uncMCForIDIsoHighPt1*uncMCForIDIsoHighPt1     + uncMCForIDIsoHighPt2*uncMCForIDIsoHighPt2     );
	}
	else {
		uncDataForLoose1HighPt2 = sqrt( effDataForTriggerLoose1*effDataForTriggerLoose1*uncDataForTriggerLoose1*uncDataForTriggerLoose1     + effDataForTriggerHighPt2*effDataForTriggerHighPt2*uncDataForTriggerHighPt2*uncDataForTriggerHighPt2 ) / ( effDataForTriggerLoose1  + effDataForTriggerHighPt2 );
		uncDataForHighPt1Loose2 = sqrt( effDataForTriggerHighPt1*effDataForTriggerHighPt1*uncDataForTriggerHighPt1*uncDataForTriggerHighPt1 + effDataForTriggerLoose2*effDataForTriggerLoose2*uncDataForTriggerLoose2*uncDataForTriggerLoose2     ) / ( effDataForTriggerHighPt1 + effDataForTriggerLoose2  );
		uncDataForBothHigh      = sqrt( effDataForTriggerHighPt1*effDataForTriggerHighPt1*uncDataForTriggerHighPt1*uncDataForTriggerHighPt1 + effDataForTriggerHighPt2*effDataForTriggerHighPt2*uncDataForTriggerHighPt2*uncDataForTriggerHighPt2 ) / ( effDataForTriggerHighPt1 + effDataForTriggerHighPt2 );
		uncMCForLoose1HighPt2   = sqrt( effMCForTriggerLoose1*effMCForTriggerLoose1*uncMCForTriggerLoose1*uncMCForTriggerLoose1             + effMCForTriggerHighPt2*effMCForTriggerHighPt2*uncMCForTriggerHighPt2*uncMCForTriggerHighPt2         ) / ( effMCForTriggerLoose1    + effMCForTriggerHighPt2   );
		uncMCForHighPt1Loose2   = sqrt( effMCForTriggerHighPt1*effMCForTriggerHighPt1*uncMCForTriggerHighPt1*uncMCForTriggerHighPt1         + effMCForTriggerLoose2*effMCForTriggerLoose2*uncMCForTriggerLoose2*uncMCForTriggerLoose2             ) / ( effMCForTriggerHighPt1   + effMCForTriggerLoose2    );
		uncMCForBothHigh        = sqrt( effMCForTriggerHighPt1*effMCForTriggerHighPt1*uncMCForTriggerHighPt1*uncMCForTriggerHighPt1         + effMCForTriggerHighPt2*effMCForTriggerHighPt2*uncMCForTriggerHighPt2*uncMCForTriggerHighPt2         ) / ( effMCForTriggerHighPt1   + effMCForTriggerHighPt2   );
	}

	double dataEff = effDataForLoose1HighPt2 + effDataForHighPt1Loose2 - effDataForBothHigh;
	double mcEff   = effMCForLoose1HighPt2   + effMCForHighPt1Loose2   - effMCForBothHigh;

	double dataUnc = sqrt( (effDataForLoose1HighPt2*uncDataForLoose1HighPt2)*(effDataForLoose1HighPt2*uncDataForLoose1HighPt2) + (effDataForHighPt1Loose2*uncDataForHighPt1Loose2)*(effDataForHighPt1Loose2*uncDataForHighPt1Loose2) + (effDataForBothHigh*uncDataForBothHigh)*(effDataForBothHigh*uncDataForBothHigh) );
	double mcUnc   = sqrt( (effMCForLoose1HighPt2*uncMCForLoose1HighPt2)*(effMCForLoose1HighPt2*uncMCForLoose1HighPt2)         + (effMCForHighPt1Loose2*uncMCForHighPt1Loose2)*(effMCForHighPt1Loose2*uncMCForHighPt1Loose2)         + (effMCForBothHigh*uncMCForBothHigh)*(effMCForBothHigh*uncMCForBothHigh) );

	double unc = sqrt( (dataUnc/dataEff)*(dataUnc/dataEff) + (mcUnc/mcEff)*(mcUnc/mcEff) );
	if(unc > 1.) unc = 1.;
	return 1. + unc;
	// else return 1.05;

}

std::vector<double> ScaleFactor::getMuonEffUncsForIDIso(const NtupleMuon& mu) {
	int nAbsetaBins = ScaleFactor::muonNAbsetaBins;
	int nPtBins = ScaleFactor::muonNPtBinsForIDIso;

	double relUncIDLooseData, relUncIDLooseMC, relUncIsoLooseData, relUncIsoLooseMC;
	double relUncIDHighPtData, relUncIDHighPtMC, relUncIsoHighPtData, relUncIsoHighPtMC;

	int absetaIdx = 0, ptIdx = 0;
	for(; absetaIdx<nAbsetaBins; ++absetaIdx){
		if( fabs(mu.eta) >= MUON_ABSETABINS[absetaIdx] && fabs(mu.eta) < MUON_ABSETABINS[absetaIdx+1] ) {
			for(; ptIdx<nPtBins; ++ptIdx) {
				if( (mu.pt >= MUON_PTBINS_IDISO[ptIdx] && mu.pt < MUON_PTBINS_IDISO[ptIdx+1]) || ptIdx == nPtBins-1 ) {
					relUncIDLooseData   = MUON_ID_UNC_LOOSE_DATA[absetaIdx][ptIdx]   / MUON_ID_LOOSE_DATA[absetaIdx][ptIdx];
					relUncIDLooseMC     = MUON_ID_UNC_LOOSE_MC[absetaIdx][ptIdx]     / MUON_ID_LOOSE_MC[absetaIdx][ptIdx];
					relUncIDHighPtData  = MUON_ID_UNC_HIGHPT_DATA[absetaIdx][ptIdx]  / MUON_ID_HIGHPT_DATA[absetaIdx][ptIdx];
					relUncIDHighPtMC    = MUON_ID_UNC_HIGHPT_MC[absetaIdx][ptIdx]    / MUON_ID_HIGHPT_MC[absetaIdx][ptIdx];
					relUncIsoLooseData  = MUON_ISO_UNC_LOOSE_DATA[absetaIdx][ptIdx]  / MUON_ISO_LOOSE_DATA[absetaIdx][ptIdx];
					relUncIsoLooseMC    = MUON_ISO_UNC_LOOSE_MC[absetaIdx][ptIdx]    / MUON_ISO_LOOSE_MC[absetaIdx][ptIdx];
					relUncIsoHighPtData = MUON_ISO_UNC_HIGHPT_DATA[absetaIdx][ptIdx] / MUON_ISO_HIGHPT_DATA[absetaIdx][ptIdx];
					relUncIsoHighPtMC   = MUON_ISO_UNC_HIGHPT_MC[absetaIdx][ptIdx]   / MUON_ISO_HIGHPT_MC[absetaIdx][ptIdx];
					return { sqrt(relUncIDLooseData*relUncIDLooseData   + relUncIsoLooseData*relUncIsoLooseData),
						 sqrt(relUncIDLooseMC*relUncIDLooseMC       + relUncIsoLooseMC*relUncIsoLooseMC),
						 sqrt(relUncIDHighPtData*relUncIDHighPtData + relUncIsoHighPtData*relUncIsoHighPtData),
						 sqrt(relUncIDHighPtMC*relUncIDHighPtMC     + relUncIsoHighPtMC*relUncIsoHighPtMC) 
						};
				}
			}
		}
	}
	return {};
}

std::vector<double> ScaleFactor::getMuonEffUncsForTrigger(const NtupleMuon& mu) {
	int nAbsetaBins = ScaleFactor::muonNAbsetaBins;
	int nPtBins = ScaleFactor::muonNPtBinsForIDIso;

	double relUncLooseData, relUncLooseMC;
	double relUncHighPtData, relUncHighPtMC;

	int absetaIdx = 0, ptIdx = 0;
	for(; absetaIdx<nAbsetaBins; ++absetaIdx){
		if( fabs(mu.eta) >= MUON_ABSETABINS[absetaIdx] && fabs(mu.eta) < MUON_ABSETABINS[absetaIdx+1] ) {
			for(; ptIdx<nPtBins; ++ptIdx) {
				if( (mu.pt >= MUON_PTBINS_TRIGGER[ptIdx] && mu.pt < MUON_PTBINS_TRIGGER[ptIdx+1]) || ptIdx == nPtBins-1 ) {
					relUncLooseData   = MUON_TRIGGER_UNC_LOOSE_DATA[absetaIdx][ptIdx]   / MUON_TRIGGER_LOOSE_DATA[absetaIdx][ptIdx];
					relUncLooseMC     = MUON_TRIGGER_UNC_LOOSE_MC[absetaIdx][ptIdx]     / MUON_TRIGGER_LOOSE_MC[absetaIdx][ptIdx];
					relUncHighPtData  = MUON_TRIGGER_UNC_HIGHPT_DATA[absetaIdx][ptIdx]  / MUON_TRIGGER_HIGHPT_DATA[absetaIdx][ptIdx];
					relUncHighPtMC    = MUON_TRIGGER_UNC_HIGHPT_MC[absetaIdx][ptIdx]    / MUON_TRIGGER_HIGHPT_MC[absetaIdx][ptIdx];
					return { relUncLooseData, relUncLooseMC, relUncHighPtData, relUncHighPtMC };
				}
			}
		}
	}
	return {};
}
