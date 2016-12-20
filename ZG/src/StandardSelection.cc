#include "Physics/ZG/interface/StandardSelection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <assert.h>
#include <iostream>
#include <vector>

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

bool MuonSelector::selectICHEP16( const std::vector<NtupleMuon>& candidates ) {

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

bool MuonSelector::selectMoriond17( const std::vector<NtupleMuon>& candidates ) {

	selected.clear();
	boostedTrackerIsolation(candidates);
	// for(auto& mu:candidates) {
	// 	if( miniIsolation(mu,0.2) ) selected.push_back(&mu);
	// }
	std::vector<const NtupleMuon*> passed;
	bool leadingMu = false, highPt = false;
	int extraPt = 0;
	for( auto& mu : selected ) {
		if( acceptance(*mu,vetoPtCut,etaCut) && (looseMuonID(*mu,false)||highPtMuonID(*mu)) ) {
		//if( acceptance(*mu,vetoPtCut,etaCut) && (looseMuonID(*mu,true)||highPtMuonID(*mu)) ) {
		// if( acceptance(*mu,vetoPtCut,etaCut) && looseMuonID(*mu,true) ) {
			passed.push_back(mu);
			if(mu->pt > extraPtCut) {
				++extraPt;
				if(mu->pt > ptCut) leadingMu = true;
				if(highPtMuonID(*mu)) highPt = true;
			}
		}
	} 

	selected = passed;
	pass = (selected.size()==2 && extraPt==2 && leadingMu && highPt);

	return pass;

	// selected.clear();

	// for( auto& mu : candidates ) {
	// 	if( acceptance(mu,vetoPtCut,etaCut) && (looseMuonID(mu,true)||highPtMuonID(mu)) ) selected.push_back(&mu);
	// } 

	// pass = (selected.size() > 1);
	// if(pass) pass = boostedTrackerIsolation();
	// if(pass) pass = finalMuonSelection();

	// return pass;

}

void MuonSelector::boostedTrackerIsolation( const std::vector<NtupleMuon>& candidates ) {
	selected.clear();
	for( unsigned i = 0; i < candidates.size(); ++i ) {
		if( trackerIsolation(candidates[i],0.1) ) selected.push_back(&candidates[i]);
		else {
			double iso = candidates[i].isolationR03_sumpt;
			for( unsigned j = 0; j < candidates.size(); ++j ) {
				if(i==j) continue;
				if( deltaR(candidates[i].eta,candidates[i].phi,candidates[j].eta,candidates[j].phi) < 0.3 ) {
					iso -= candidates[j].innerTrack_pt;
				}
			}
			if( iso > -0.1*candidates[i].pt && iso/candidates[i].pt < 0.1 ) selected.push_back(&candidates[i]);
		}
	}
}

// bool MuonSelector::boostedTrackerIsolation() {

// 	std::pair<const NtupleMuon*, const NtupleMuon*> cands{nullptr, nullptr};

// 	std::cout << "### ordering by pt ###" << std::endl;
// 	double leading = 0., trailing = 0.;
// 	for(auto& mu : selected) {
// 		std::cout << mu->pt << std::endl;
// 		if(mu->pt > leading) {
// 			leading = mu->pt;
// 			cands.second = cands.first;
// 			cands.first = mu;
// 		}
// 		else if(mu->pt > trailing) {
// 			trailing = mu->pt;
// 			cands.second = mu;
// 		}
// 	}

// 	assert(cands.first!=nullptr&&cands.second!=nullptr);

// 	std::cout << "selected:" << selected.size() << std::endl;
// 	selected.clear();
// 	selected[0] = cands.first;
// 	selected[1] = cands.second;
// 	std::cout << selected[0]->pt << ", " << selected[1]->pt << std::endl;


// 	if( trackerIsolation(*selected[0],0.1) && trackerIsolation(*selected[1],0.1) ) {
// 		std::cout << "isolated" << std::endl;
// 		return true;
// 	}
// 	else {
// 		if( deltaR(selected[0]->eta,selected[0]->phi,selected[1]->eta,selected[1]->phi) < 0.3 ) {
// 			double iso0 = selected[0]->isolationR03_sumpt - selected[1]->pt;
// 			double iso1 = selected[1]->isolationR03_sumpt - selected[0]->pt;
// 			// if(iso0 < -0.1*selected[1]->pt) iso0 = selected[0]->isolationR03_sumpt;
// 			// if(iso1 < -0.1*selected[0]->pt) iso1 = selected[1]->isolationR03_sumpt;
// 			return ( iso0/selected[0]->pt < 0.1 && iso1/selected[1]->pt < 0.1 );
// 		}
// 		else return false;
// 	}

// }

bool MuonSelector::finalMuonSelection() {

	// std::cout << "### final selection ###" << std::endl;
	bool leadingMu = false, highPt = false;
	int extraPt = 0;
	for(auto& mu : selected) {
		if(mu->pt > extraPtCut) {
			++extraPt;
			if(mu->pt > ptCut) leadingMu = true;
			if(highPtMuonID(*mu)) highPt = true;
		}
	}
	return (extraPt==2 && leadingMu && highPt);

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
		if( acceptance(el,vetoPtCut,etaCut) && excludeECALGap(el) && cutBasedWPLoose(el,true) && miniIsolation(el,0.1) ) {
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
		if( acceptance(pho,ptCut,etaCut) && excludeECALGap(pho) && cutBasedWPLoose(pho,true) ) selected.push_back(&pho);
	} 

	pass = (selected.size()==1);
	return pass;

}

bool PhotonSelector::select( const std::vector<NtuplePhoton>& candidates, const NtupleElectron& el0, const NtupleElectron& el1 ) {

	selected.clear();
	
	for( auto& pho : candidates ) {
		if( acceptance(pho,ptCut,etaCut) && excludeECALGap(pho) && cutBasedWPLoose(pho,true) 
			&& deltaR(pho.eta,pho.phi,el0.eta,el0.phi) > 0.4 
			&& deltaR(pho.eta,pho.phi,el1.eta,el1.phi) > 0.4 ) selected.push_back(&pho);
	} 

	pass = (selected.size()==1);
	return pass;

}

bool PhotonSelector::select( const std::vector<NtuplePhoton>& candidates, const NtupleMuon& mu0, const NtupleMuon& mu1 ) {

	selected.clear();
	
	for( auto& pho : candidates ) {
		if( acceptance(pho,ptCut,etaCut) && excludeECALGap(pho) && cutBasedWPLoose(pho,true) 
			&& deltaR(pho.eta,pho.phi,mu0.eta,mu0.phi) > 0.4 
			&& deltaR(pho.eta,pho.phi,mu1.eta,mu1.phi) > 0.4 ) selected.push_back(&pho);
	} 

	pass = (selected.size()==1);
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
	return ( fabs(particle.etaSC) < 1.4442 || fabs(particle.etaSC) > 1.566 );
}

template<class T> bool miniIsolation(const T& particle, const double& iso) {
	return ( particle.miniIsoAbs/particle.pt < iso );
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
		  && fabs(mu.dzVTX) < 0.5 );
}

bool highPtMuonID(const NtupleMuon& mu) {
	return ( mu.isGlobalMuon 
		  && mu.isTrackerMuon 
		  && mu.nValidMuonHits>0 
		  && mu.nMatchedStations>1 
		  && mu.nValidPixelHits>0 
		  && mu.nTrackerLayers>5 
		  && fabs(mu.dxyVTX)<0.2 
		  && mu.muonBestTrack_ptError/mu.muonBestTrack_pt<0.3 );
}

bool cutBasedWPLoose(const NtupleElectron& el, const bool& isZG) {
    if( fabs(el.eta) <= 1.479 ) { //Barrel
        return ( el.sigmaIetaIeta < 0.011
              && fabs(el.dEtaIn) < 0.00477
              && fabs(el.dPhiIn) < 0.222
              && el.hOverE < 0.298
              && el.ooEmooP < 0.241
              // && fabs(el.d0) < 0.05
              // && fabs(el.dz) < 0.10
              && el.expectedMissingInnerHits <= 2
              && el.passConversionVeto
              && (isZG || el.isoRho < 0.0994) );
    }
    else { //Endcap
        return ( el.sigmaIetaIeta < 0.0314
              && fabs(el.dEtaIn) < 0.00868
              && fabs(el.dPhiIn) < 0.213
              && el.hOverE < 0.101
              && el.ooEmooP < 0.14
              // && fabs(el.d0) < 0.1
              // && fabs(el.dz) < 0.2
              && el.expectedMissingInnerHits <= 1
              && el.passConversionVeto 
              && (isZG || el.isoRho < 0.107) );
    }
    return false;

    // Spring15 //
    // if( fabs(el.eta) <= 1.479 ) { //Barrel
    //     return ( el.sigmaIetaIeta < 0.0103
    //           && fabs(el.dEtaIn) < 0.0105
    //           && fabs(el.dPhiIn) < 0.115
    //           && el.hOverE < 0.104
    //           && el.ooEmooP < 0.102
    //           && fabs(el.d0) < 0.0261
    //           && fabs(el.dz) < 0.41
    //           && el.expectedMissingInnerHits <= 2
    //           && el.passConversionVeto
    //           && (isZG || el.isoRho < 0.0893) );
    // }
    // else { //Endcap
    //     return ( el.sigmaIetaIeta < 0.0301
    //           && fabs(el.dEtaIn) < 0.00814
    //           && fabs(el.dPhiIn) < 0.182
    //           && el.hOverE < 0.0897
    //           && el.ooEmooP < 0.126
    //           && fabs(el.d0) < 0.118
    //           && fabs(el.dz) < 0.822
    //           && el.expectedMissingInnerHits <= 1
    //           && el.passConversionVeto 
    //           && (isZG || el.isoRho < 0.121) );
    // }
    // return false;
}

bool cutBasedWPLoose(const NtuplePhoton& pho, const bool& isZG) {
    if( fabs(pho.eta) <= 1.479 ) { // Barrel
        if ( pho.HoverE < 0.0597 && pho.Full5x5_SigmaIEtaIEta < 0.01031 ) {
        	if(isZG) return ( pho.ChIso < 2.5 && !pho.hasPixelSeed );
        	else return ( pho.ChIso < 1.295 && pho.NhIsoWithEA < 10.910 + 0.0148 * pho.pt + 0.000017 * pho.pt * pho.pt && pho.PhIsoWithEA < 3.630 + 0.0047 * pho.pt );
    	}
    }
    else { // Endcap
        if( pho.HoverE < 0.0481 && pho.Full5x5_SigmaIEtaIEta < 0.03013 ){
            if(isZG) return ( pho.ChIso < 2.5 && !pho.hasPixelSeed );
            else return ( pho.ChIso < 1.011 && pho.NhIsoWithEA < 5.931 + 0.0163 * pho.pt + 0.000014 * pho.pt * pho.pt && pho.PhIsoWithEA < 6.641 + 0.0034 * pho.pt );
        }
    }
    return false;

    // Spring15 //
    // if( fabs(pho.eta) <= 1.479 ) { // Barrel
    //     if ( pho.HoverE < 0.05 && pho.Full5x5_SigmaIEtaIEta < 0.0102 && !pho.hasPixelSeed ) {
    //     	if(isZG) return ( pho.ChIso < 2.5 );
    //     	else return ( pho.ChIso < 3.32 && pho.NhIsoWithEA < 1.92 + 0.014 * pho.pt + 0.000019 * pho.pt * pho.pt && pho.PhIsoWithEA < 0.81 + 0.0053 * pho.pt );
    // 	}
    // }
    // else { // Endcap
    //     if( pho.HoverE < 0.05 && pho.Full5x5_SigmaIEtaIEta < 0.0274 && !pho.hasPixelSeed ){
    //         if(isZG) return ( pho.ChIso < 2.5 );
    //         else return ( pho.ChIso < 1.97 && pho.NhIsoWithEA < 11.86 + 0.0139 * pho.pt + 0.000025 * pho.pt * pho.pt && pho.PhIsoWithEA < 0.83 + 0.0034 * pho.pt );
    //     }
    // }
    // return false;
}
