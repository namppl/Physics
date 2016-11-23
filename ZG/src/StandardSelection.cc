#include <vector>
#include <iostream>
#include "Physics/ZG/interface/StandardSelection.h"

using namespace std;

Selector::Selector(): ptCut{0}, etaCut{999} {}
Selector::Selector( const double& _ptCut, const double& _etaCut ): ptCut{_ptCut}, etaCut{_etaCut} {}
bool Selector::passSelection() const {
	return pass;
}

LeptonSelector::LeptonSelector(): Selector{}, extraPtCut{0}, vetoPtCut{0} {}
LeptonSelector::LeptonSelector( const double& _ptCut, const double& _etaCut, const double& _extraPtCut, const double& _vetoPtCut ):
Selector( _ptCut, _etaCut ), extraPtCut{_extraPtCut}, vetoPtCut{_vetoPtCut} 
{}

MuonSelector::MuonSelector(): LeptonSelector{} {}
MuonSelector::MuonSelector( const double& _ptCut, const double& _etaCut, const double& _extraPtCut, const double& _vetoPtCut ):
LeptonSelector( _ptCut, _etaCut, _extraPtCut, _vetoPtCut )
{
	std::cout << "muon pt threshold: " << ptCut << " GeV, " << extraPtCut << " GeV (veto " << vetoPtCut <<  " GeV)" << std::endl;
	std::cout << "muon eta restriction: " << etaCut << std::endl;
}

bool MuonSelector::select( const std::vector<NtupleMuon>& candidates ) {

	bool leadingMu = false, tightMu = false;
	selected.clear();

	for( auto& mu : candidates ) {
		if( acceptance(mu,vetoPtCut,etaCut) && looseMuonID(mu,true) && miniIsolation(mu,0.2) ) {
			selected.push_back(&mu);
			if( mu.pt > ptCut ) leadingMu = true;
			if( tightMuonID(mu) ) tightMu= true;
		} 
	} 

	pass = false;
	if( leadingMu && tightMu && selected.size() == 2 ) {
		if( selected[0]->pt > extraPtCut && selected[1]->pt > extraPtCut ) pass = true;
	} 
	return pass;

}

const NtupleMuon& MuonSelector::at(int i) const {
	return *selected[i];
}

int MuonSelector::nMuons() {
	return selected.size();
}

ElectronSelector::ElectronSelector(): LeptonSelector{} {}
ElectronSelector::ElectronSelector( const double& _ptCut, const double& _etaCut, const double& _extraPtCut, const double& _vetoPtCut ):
LeptonSelector( _ptCut, _etaCut, _extraPtCut, _vetoPtCut )
{
	std::cout << "electron pt threshold: " << ptCut << " GeV, " << extraPtCut << " GeV (veto " << vetoPtCut <<  " GeV)" << std::endl;
	std::cout << "electron eta restriction: " << etaCut << std::endl;
}

bool ElectronSelector::select( const std::vector<NtupleElectron>& candidates ) {

	bool leadingEl	= false;
	selected.clear();

	for( auto& el : candidates ) {
		if( acceptance(el,vetoPtCut,etaCut) && excludeECALGap(el) && WPLoose(el,true) && miniIsolation(el,0.1) ) {
			selected.push_back(&el);
			if(el.pt > ptCut) leadingEl = true;
		}
	}

	pass = false;
	if( leadingEl && selected.size() == 2 ) {
		if( selected[0]->pt > extraPtCut && selected[1]->pt > extraPtCut ) pass = true;
	}
	return pass;

}

const NtupleElectron& ElectronSelector::at(int i) const {
	return *selected[i];
}

int ElectronSelector::nElectrons() {
	return selected.size();
}

PhotonSelector::PhotonSelector(): Selector{} {}
PhotonSelector::PhotonSelector( const double& _ptCut, const double& _etaCut ): Selector( _ptCut, _etaCut ) {
	std::cout << "photon pt threshold: " << ptCut << std::endl;
	std::cout << "photon eta restriction: " << etaCut << std::endl;
}


bool PhotonSelector::select( const std::vector<NtuplePhoton>& candidates ) {

	selected.clear();
	
	for( auto& pho : candidates ) {
		if( acceptance(pho,ptCut,etaCut) && excludeECALGap(pho) && WPLoose(pho,true) ) selected.push_back(&pho);
	} 

	pass = ( selected.size() == 1 ) ? true : false;
	return pass;

}

const NtuplePhoton& PhotonSelector::at(int i) const {
	return *selected[i];
}

int PhotonSelector::nPhotons() {
	return selected.size();
}

// template<class T> TLorentzVector momentum(const T& particle, const double& mass){
// 	TLorentzVector mom;
// 	mom.SetPtEtaPhiM(particle.pt,particle.eta,particle.phi,mass);
// 	return mom;
// }

TLorentzVector momentum(const NtupleMuon& mu, const double& mass){
	TLorentzVector mom;
	mom.SetPtEtaPhiM(mu.pt,mu.eta,mu.phi,mass);
	return mom;
}

TLorentzVector momentum(const NtupleElectron& el, const double& mass){
	TLorentzVector mom;
	mom.SetPtEtaPhiM(el.pt,el.eta,el.phi,mass);
	return mom;
}

TLorentzVector momentum(const NtuplePhoton& pho, const double& mass){
	TLorentzVector mom;
	mom.SetPtEtaPhiM(pho.pt,pho.eta,pho.phi,mass);
	return mom;
}

template<class T> bool acceptance(const T& particle, const double& ptCut, const double& etaCut) {
	return ( particle.pt > ptCut && fabs(particle.eta) < etaCut );
}

template<class T> bool excludeECALGap(const T& particle) {
	return ( fabs(particle.eta) < 1.4442 || fabs(particle.eta) > 1.566 );
}

template<class T> bool miniIsolation(const T& particle, const double& iso) {
	return ( particle.miniIso/particle.pt < iso );
}

bool trackerIsolation(const NtupleMuon& mu, const double& iso) {
	return ( mu.isolationR03_sumpt/mu.pt < iso );
}

bool looseMuonID(const NtupleMuon& mu, const bool& isZG) {
    return ( mu.isPFMuon && (mu.isGlobalMuon||mu.isTrackerMuon) && ( !isZG || (fabs(mu.dxyVTX) < 0.2 && fabs(mu.dzVTX) < 0.5) ) );
}

bool tightMuonID(const NtupleMuon& mu) {
	return ( mu.isPFMuon 
		  && mu.isGlobalMuon 
		  && mu.normalizedChi2 < 10 
		  && mu.nValidMuonHits>0 
		  && mu.nMatchedStations>1 
		  && mu.nValidPixelHits>0 
		  && mu.nTrackerLayers>5 
		  && fabs(mu.dxyVTX) < 0.2 
		  && fabs(mu.dzVTX) < 0.5 ) ? true : false;
}

bool highPtMuonID(const NtupleMuon& mu) {
	return ( mu.isGlobalMuon 
		  && mu.isTrackerMuon 
		  && mu.nValidMuonHits>0 
		  && mu.nMatchedStations>1 
		  && mu.nValidPixelHits>0 
		  && mu.nTrackerLayers>5 
		  && fabs(mu.dxyVTX)<0.2 
		  && mu.muonBestTrack_ptError/mu.muonBestTrack_pt<0.3 ) ? true : false;
}

bool WPLoose(const NtupleElectron& el, const bool& isZG) {
    if( fabs(el.eta) <= 1.479 ) { //Barrel
        return ( el.sigmaIetaIeta < 0.0103
              && fabs(el.dEtaIn) < 0.0105
              && fabs(el.dPhiIn) < 0.115
              && el.hOverE < 0.104
              && el.ooEmooP < 0.102
              && fabs(el.d0) < 0.0261
              && fabs(el.dz) < 0.41
              && el.expectedMissingInnerHits <= 2
              && el.passConversionVeto
              && (isZG || el.isoRho < 0.0893) );
    }
    else { //Endcap
        return ( el.sigmaIetaIeta < 0.0301
              && fabs(el.dEtaIn) < 0.00814
              && fabs(el.dPhiIn) < 0.182
              && el.hOverE < 0.0897
              && el.ooEmooP < 0.126
              && fabs(el.d0) < 0.118
              && fabs(el.dz) < 0.822
              && el.expectedMissingInnerHits <= 1
              && el.passConversionVeto 
              && (isZG || el.isoRho < 0.121) );
    }
    return false;
}

bool WPLoose(const NtuplePhoton& pho, const bool& isZG) {
    if( fabs(pho.eta) <= 1.479 ) { // Barrel
        if ( pho.HoverE < 0.05 && pho.Full5x5_SigmaIEtaIEta < 0.0102 && !pho.hasPixelSeed ) {
        	if(isZG) return ( pho.ChIso < 2.5 );
        	else ( pho.ChIso < 3.32 && pho.NhIsoWithEA < 1.92 + 0.014 * pho.pt + 0.000019 * pho.pt * pho.pt && pho.PhIsoWithEA < 0.81 + 0.0053 * pho.pt );
    	}
    }
    else { // Endcap
        if( pho.HoverE < 0.05 && pho.Full5x5_SigmaIEtaIEta < 0.0274 && !pho.hasPixelSeed ){
            if(isZG) return ( pho.ChIso < 2.5 );
            else ( pho.ChIso < 1.97 && pho.NhIsoWithEA < 11.86 + 0.0139 * pho.pt + 0.000025 * pho.pt * pho.pt && pho.PhIsoWithEA < 0.83 + 0.0034 * pho.pt );
        }
    }
    return false;
}
