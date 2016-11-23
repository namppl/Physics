#include "Physics/ZG/interface/EventSelection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <assert.h>

EventSelector::EventSelector( const bool& _verbose): verbose{_verbose}, dilepton_mass_lower{50}, dilepton_mass_upper{130}, pho_pt_over_mass{40./150.}
{}

EventSelector::EventSelector( const bool& _verbose, const double& mu_pt_cut1, const double& mu_pt_cut2, const double& mu_pt_veto, const double& el_pt_cut1, const double& el_pt_cut2, const double& el_pt_veto, const double& pho_pt_cut, const double& _dilepton_mass_lower, const double& _dilepton_mass_upper, const double& _pho_pt_over_mass ):
verbose{_verbose}, dilepton_mass_lower{_dilepton_mass_lower}, dilepton_mass_upper{_dilepton_mass_upper}, pho_pt_over_mass{_pho_pt_over_mass}
{
	std::cout << std::endl;
	std::cout << ":: selection rule ::" << std::endl;
	std::cout << "lower limit on dilepton mass: " << dilepton_mass_lower << std::endl; 
	std::cout << "upper limit on dilepton mass: " << dilepton_mass_upper << std::endl; 
	std::cout << "photon pt over ZG mass: " << pho_pt_over_mass << std::endl;
	mm = new MuonSelector(mu_pt_cut1, 2.4, mu_pt_cut2, mu_pt_veto);
	ee = new ElectronSelector(el_pt_cut1, 2.4, el_pt_cut2, el_pt_veto);
	gamma = new PhotonSelector(pho_pt_cut, 2.5);
	std::cout << std::endl;
}

EventSelector::~EventSelector() {
	delete mm;
	delete ee;
	delete gamma;
}

void EventSelector::setEvent(NtupleEvent* _event) {
	event = _event;
}

void EventSelector::select() {
	
	run = event->run;
	lumi = event->lumi;
	eventNum = event->event;
	weight = event->weight;	
	nPU = event->nPU;
	nVert = event->nVertices;

	muon = false;
	electron = false;

	mm->select(event->muons);
	ee->select(event->electrons);
	gamma->select(event->photons);

	if( mm->passSelection() && ee->nElectrons()==0 && gamma->passSelection() ) {

		const NtupleMuon& mu1{mm->at(0)};
		const NtupleMuon& mu2{mm->at(1)};
		const NtuplePhoton& pho{gamma->at(0)};

		double z_mass = ( momentum(mu1,muon_mass) + momentum(mu2,muon_mass) ).M();
		double boss_mass = ( momentum(mu1,muon_mass) + momentum(mu2,muon_mass) + momentum(pho,0.) ).M();
		double dR1 = deltaR(mu1.eta,mu1.phi,pho.eta,pho.phi); 
		double dR2 = deltaR(mu2.eta,mu2.phi,pho.eta,pho.phi);

		if( mu1.charge*mu2.charge < 0 && dR1 > 0.4 && dR2 > 0.4 && z_mass > dilepton_mass_lower && z_mass < dilepton_mass_upper && pho.pt/boss_mass > pho_pt_over_mass ) {
			muon = true;
		}

	}
	else if( ee->passSelection() && mm->nMuons()==0 && gamma->passSelection() ) {	

		const NtupleElectron& el1{ee->at(0)};
		const NtupleElectron& el2{ee->at(1)};
		const NtuplePhoton& pho{gamma->at(0)};

		double z_mass = ( momentum(el1,electron_mass) + momentum(el2,electron_mass) ).M();
		double boss_mass = ( momentum(el1,electron_mass) + momentum(el2,electron_mass) + momentum(pho,0.) ).M();
		double dR1 = deltaR(el1.eta,el1.phi,pho.eta,pho.phi); 
		double dR2 = deltaR(el2.eta,el2.phi,pho.eta,pho.phi);

		if( el1.charge*el2.charge < 0 && dR1 > 0.4 && dR2 > 0.4 && z_mass > dilepton_mass_lower && z_mass < dilepton_mass_upper && pho.pt/boss_mass > pho_pt_over_mass ) {
			electron = true;
		}
		
	}

}

bool EventSelector::isTriggered( const TString& triggerName ) {
	for( auto& trigger : event->triggers ) {
	    if( trigger.name.find(triggerName) == 0 ) return trigger.isFired;
	}
    return false;
}

void EventSelector::setGenMass() {
	TLorentzVector mom_system;
	int nGenEl = 0;
	int nGenMu = 0;
	int nGenPho = 0;
	for( auto& par : event->genparticles ) {
		//if( par.fromHardProcessFinalState ) {
		if( par.isHardProcess ) {
			TLorentzVector mom;
			if( fabs(par.id) == 11 ) {
				mom.SetPtEtaPhiM(par.pt, par.eta, par.phi, electron_mass);
				++nGenEl;
			}
			else if( fabs(par.id) == 13 ) {
				mom.SetPtEtaPhiM(par.pt, par.eta, par.phi, muon_mass);
				++nGenMu;
			}
			else if( fabs(par.id) == 22 ) {
				mom.SetPtEtaPhiM(par.pt, par.eta, par.phi, 0.);
				++nGenPho;
			}
			mom_system += mom;
		} 
	}
	assert( (nGenEl==2 && nGenMu==0 && nGenPho==1) || (nGenEl==0 && nGenMu==2 && nGenPho==1) || (nGenEl==0 && nGenMu==0 && nGenPho==1) );
	
	if(nGenEl==2) {
		channel = "ee";
		genMass = mom_system.M();
	}
	else if(nGenMu==2) {
		channel = "mm";
		genMass = mom_system.M();
	}
	else {
		channel = "tt";
		genMass = -999;
	}
}
