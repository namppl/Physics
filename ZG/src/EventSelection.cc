#include "Physics/ZG/interface/EventSelection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <assert.h>

EventSelector::EventSelector( const bool& _verbose): verbose{_verbose}, elPtCut{65.}, elExtraPtCut{10.}, elVetoPtCut{10.}, phoPtCut_ee{65.}, muPtCut{26.}, muExtraPtCut{10.}, muVetoPtCut{10.}, phoPtCut_mm{40.}, lowerMassCut{50}, upperMassCut{130}, phoPtOverMZG{40./150.}, multiPho{false}
{}

EventSelector::EventSelector( const bool& _verbose, const double& _elPtCut, const double& _elExtraPtCut, const double& _elVetoPtCut, const double& _phoPtCut_ee, const double& _muPtCut, const double& _muExtraPtCut, const double& _muVetoPtCut, const double& _phoPtCut_mm, const double& _lowerMassCut, const double& _upperMassCut, const double& _phoPtOverMZG, const bool& _multiPho ):
verbose{_verbose}, elPtCut{_elPtCut}, elExtraPtCut{_elExtraPtCut}, elVetoPtCut{_elVetoPtCut}, phoPtCut_ee{_phoPtCut_ee}, muPtCut{_muPtCut}, muExtraPtCut{_muExtraPtCut}, muVetoPtCut{_muVetoPtCut}, phoPtCut_mm{_phoPtCut_mm}, lowerMassCut{_lowerMassCut}, upperMassCut{_upperMassCut}, phoPtOverMZG{_phoPtOverMZG}, multiPho{_multiPho}
{
	std::cout << std::endl;
	std::cout << ":: selection rule ::" << std::endl;
	std::cout << "lower mass threshold: " << lowerMassCut << std::endl; 
	std::cout << "upper mass threshold: " << upperMassCut << std::endl; 
	std::cout << "photon pt over ZG mass: " << phoPtOverMZG << std::endl;
	std::cout << "multiple photon: "<< multiPho << std::endl;
	ee = new ElectronSelector(elPtCut, 2.4, elExtraPtCut, elVetoPtCut);
	mm = new MuonSelector(muPtCut, 2.4, muExtraPtCut, muVetoPtCut);
	gamma_ee = new PhotonSelector(phoPtCut_ee, 2.5);
	gamma_mm = new PhotonSelector(phoPtCut_mm, 2.5);
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
	eventNum = event->event;
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

void EventSelector::selectMoriond17() {
	
	run = event->run;
	lumi = event->lumi;
	eventNum = event->event;
	weight = event->weight;	
	nPU = event->nPU;
	nVert = event->nVertices;

	muon = false;
	electron = false;

	ee->select(event->electrons);
	mm->selectMoriond17(event->muons);
	gamma_ee->select(event->photons);
	gamma_mm->select(event->photons);

	if( ee->passSelection() && mm->nMuons()==0 && gamma_ee->passSelection() ) selectDielGamma();
	else if( mm->passSelection() && ee->nElectrons()==0 && gamma_mm->passSelection() ) selectDimuGamma();

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

	if( el1.charge*el2.charge < 0 && dR1 > 0.4 && dR2 > 0.4 && z_mass > lowerMassCut && z_mass < upperMassCut && pho.pt/boss_mass > phoPtOverMZG ) {
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