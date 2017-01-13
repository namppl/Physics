#include "Physics/ZG/interface/ZGTree.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <assert.h>

void TreeFiller::setTree(TTree* _recoTree, TTree* _genTree) {

	recoTree = _recoTree;
	genTree = _genTree;

	recoTree->Branch("run",&run);
	recoTree->Branch("lumi",&lumi);
	recoTree->Branch("event",&eventNum);
	recoTree->Branch("leptType",&leptType);
	recoTree->Branch("weight",&weight);
	recoTree->Branch("nPU",&nPU);
	recoTree->Branch("nVert",&nVert);
	recoTree->Branch("genmass",&genmass);

	recoTree->Branch("lept0_pt",&lept0_pt);
	recoTree->Branch("lept0_eta",&lept0_eta);
	recoTree->Branch("lept0_phi",&lept0_phi);
	recoTree->Branch("lept0_miniRelIso",&lept0_miniRelIso);
	recoTree->Branch("lept0_trkIso",&lept0_trkIso);
	recoTree->Branch("lept0_pdgId",&lept0_pdgId);
	recoTree->Branch("lept0_normalizedChi2",&lept0_normalizedChi2);
	recoTree->Branch("lept0_nValidMuonHits",&lept0_nValidMuonHits);
	recoTree->Branch("lept0_nMatchedStations",&lept0_nMatchedStations);
	recoTree->Branch("lept0_nValidPixelHits",&lept0_nValidPixelHits);
	recoTree->Branch("lept0_nTrackerLayers",&lept0_nTrackerLayers);
	recoTree->Branch("lept0_muonBestTrack_dxyVTX",&lept0_muonBestTrack_dxyVTX);
	recoTree->Branch("lept0_muonBestTrack_dzVTX",&lept0_muonBestTrack_dzVTX);
	recoTree->Branch("lept0_ptError",&lept0_ptError);
	recoTree->Branch("lept0_sigmaIetaIeta",&lept0_sigmaIetaIeta);
	recoTree->Branch("lept0_dEtaIn",&lept0_dEtaIn);
	recoTree->Branch("lept0_dPhiIn",&lept0_dPhiIn);
	recoTree->Branch("lept0_hOverE",&lept0_hOverE);
	recoTree->Branch("lept0_ooEmooP",&lept0_ooEmooP);
	recoTree->Branch("lept0_d0",&lept0_d0);
	recoTree->Branch("lept0_dz",&lept0_dz);
	recoTree->Branch("lept0_expectedMissingInnerHits",&lept0_expectedMissingInnerHits);

	recoTree->Branch("lept1_pt",&lept1_pt);
	recoTree->Branch("lept1_eta",&lept1_eta);
	recoTree->Branch("lept1_phi",&lept1_phi);
	recoTree->Branch("lept1_miniRelIso",&lept1_miniRelIso);
	recoTree->Branch("lept1_trkIso",&lept1_trkIso);
	recoTree->Branch("lept1_pdgId",&lept1_pdgId);
	recoTree->Branch("lept1_normalizedChi2",&lept1_normalizedChi2);
	recoTree->Branch("lept1_nValidMuonHits",&lept1_nValidMuonHits);
	recoTree->Branch("lept1_nMatchedStations",&lept1_nMatchedStations);
	recoTree->Branch("lept1_nValidPixelHits",&lept1_nValidPixelHits);
	recoTree->Branch("lept1_nTrackerLayers",&lept1_nTrackerLayers);
	recoTree->Branch("lept1_muonBestTrack_dxyVTX",&lept1_muonBestTrack_dxyVTX);
	recoTree->Branch("lept1_muonBestTrack_dzVTX",&lept1_muonBestTrack_dzVTX);
	recoTree->Branch("lept1_ptError",&lept1_ptError);
	recoTree->Branch("lept1_sigmaIetaIeta",&lept1_sigmaIetaIeta);
	recoTree->Branch("lept1_dEtaIn",&lept1_dEtaIn);
	recoTree->Branch("lept1_dPhiIn",&lept1_dPhiIn);
	recoTree->Branch("lept1_hOverE",&lept1_hOverE);
	recoTree->Branch("lept1_ooEmooP",&lept1_ooEmooP);
	recoTree->Branch("lept1_d0",&lept1_d0);
	recoTree->Branch("lept1_dz",&lept1_dz);
	recoTree->Branch("lept1_expectedMissingInnerHits",&lept1_expectedMissingInnerHits);

	recoTree->Branch("deltaR_lept",&deltaR_lept);

	recoTree->Branch("gamma_pt",&gamma_pt);
	recoTree->Branch("gamma_eta",&gamma_eta);
	recoTree->Branch("gamma_phi",&gamma_phi);
	recoTree->Branch("gamma_iso",&gamma_iso);
	recoTree->Branch("gamma_HoverE",&gamma_HoverE);
	recoTree->Branch("gamma_Full5x5_SigmaIEtaIEta",&gamma_Full5x5_SigmaIEtaIEta);

	recoTree->Branch("z_pt",&z_pt);
	recoTree->Branch("z_eta",&z_eta);
	recoTree->Branch("z_phi",&z_phi);
	recoTree->Branch("z_mass",&z_mass);

	recoTree->Branch("boss_pt",&boss_pt);
	recoTree->Branch("boss_eta",&boss_eta);
	recoTree->Branch("boss_phi",&boss_phi);
	recoTree->Branch("boss_mass",&boss_mass);

	genTree->Branch("genmass",&genmass);
	genTree->Branch("gen_weight",&gen_weight);
	genTree->Branch("gen_leptType",&gen_leptType);
	genTree->Branch("gen_lept0_pt",&gen_lept0_pt);
	genTree->Branch("gen_lept0_eta",&gen_lept0_eta);
	genTree->Branch("gen_lept0_phi",&gen_lept0_phi);
	genTree->Branch("gen_lept1_pt",&gen_lept1_pt);
	genTree->Branch("gen_lept1_eta",&gen_lept1_eta);
	genTree->Branch("gen_lept1_phi",&gen_lept1_phi);
	genTree->Branch("gen_gamma_pt",&gen_gamma_pt);
	genTree->Branch("gen_gamma_eta",&gen_gamma_eta);
	genTree->Branch("gen_gamma_phi",&gen_gamma_phi);
	genTree->Branch("gen_z_pt",&gen_z_pt);
	genTree->Branch("gen_z_eta",&gen_z_eta);
	genTree->Branch("gen_z_phi",&gen_z_phi);
	genTree->Branch("gen_z_mass",&gen_z_mass);
	genTree->Branch("gen_boss_pt",&gen_boss_pt);
	genTree->Branch("gen_boss_eta",&gen_boss_eta);
	genTree->Branch("gen_boss_phi",&gen_boss_phi);
	genTree->Branch("gen_boss_mass",&gen_boss_mass);

}

void TreeFiller::fillVariables(EventSelector& selector) {

	assert(selector.electron||selector.muon);

	run = selector.run;
	lumi = selector.lumi;
	eventNum = selector.eventNum;
	weight = selector.weight;	
	nPU = selector.nPU;
	nVert = selector.nVert;

	if(selector.electron) {
		leptType = 11;

		bool leading = ( selector.ee->at(0).pt > selector.ee->at(1).pt );
		const NtupleElectron& el1 = ( leading ) ? selector.ee->at(0) : selector.ee->at(1);
		const NtupleElectron& el2 = ( leading ) ? selector.ee->at(1) : selector.ee->at(0);
		const NtuplePhoton& pho = selector.gamma_ee->at(0);
		//std::cout << "pt: " << pho.pt << ", eta: " << pho.eta << ", mva: " << pho.mvaValue << std::endl; 
		
		TLorentzVector z = momentum(el1,electron_mass) + momentum(el2,electron_mass);
		TLorentzVector boss = momentum(el1,electron_mass) + momentum(el2,electron_mass) + momentum(pho,0.);

		lept0_pt  = el1.pt;
		lept0_eta = el1.eta;
		lept0_phi = el1.phi;
		lept0_miniRelIso = el1.miniIsoAbs/el1.pt;
		lept0_pdgId = 11*el1.charge;

		lept0_sigmaIetaIeta = el1.sigmaIetaIeta;
		lept0_dEtaIn = el1.dEtaIn;
		lept0_dPhiIn = el1.dPhiIn;
		lept0_hOverE = el1.hOverE;
		lept0_ooEmooP = el1.ooEmooP;
		lept0_d0 = el1.d0;
		lept0_dz = el1.dz;
		lept0_expectedMissingInnerHits = el1.expectedMissingInnerHits;

		lept1_pt  = el2.pt;
		lept1_eta = el2.eta;
		lept1_phi = el2.phi;
		lept1_miniRelIso = el2.miniIsoAbs/el2.pt;
		lept1_pdgId = 11*el2.charge;

		lept1_sigmaIetaIeta = el2.sigmaIetaIeta;
		lept1_dEtaIn = el2.dEtaIn;
		lept1_dPhiIn = el2.dPhiIn;
		lept1_hOverE = el2.hOverE;
		lept1_ooEmooP = el2.ooEmooP;
		lept1_d0 = el2.d0;
		lept1_dz = el2.dz;
		lept1_expectedMissingInnerHits = el2.expectedMissingInnerHits;

		deltaR_lept = deltaR(el1.eta,el1.phi,el2.eta,el2.phi); 

		gamma_pt = pho.pt;
		gamma_eta = pho.eta;
		gamma_phi = pho.phi;
		gamma_iso = pho.ChIso;
		gamma_HoverE = pho.HoverE;
		gamma_Full5x5_SigmaIEtaIEta = pho.Full5x5_SigmaIEtaIEta;
	
		z_pt = z.Pt();
		z_eta = z.Eta();
		z_phi = z.Phi();
		z_mass = z.M(); 

		boss_pt = boss.Pt();
		boss_eta = boss.Eta();
		boss_phi = boss.Phi();
		boss_mass = boss.M(); 

		recoTree->Fill();
	}
	else {
		leptType = 13;

		bool leading = ( selector.mm->at(0).pt > selector.mm->at(1).pt );
		const NtupleMuon& mu1 = ( leading ) ? selector.mm->at(0) : selector.mm->at(1);
		const NtupleMuon& mu2 = ( leading ) ? selector.mm->at(1) : selector.mm->at(0);
		const NtuplePhoton& pho = selector.gamma_mm->at(0);
		
		TLorentzVector z = momentum(mu1,muon_mass) + momentum(mu2,muon_mass);
		TLorentzVector boss = momentum(mu1,muon_mass) + momentum(mu2,muon_mass) + momentum(pho,0.);

		lept0_pt  = mu1.pt;
		lept0_eta = mu1.eta;
		lept0_phi = mu1.phi;
		lept0_miniRelIso = mu1.miniIsoAbs/mu1.pt;
		lept0_trkIso = mu1.isolationR03_sumpt/mu1.pt;
		lept0_pdgId = 13*mu1.charge;
		lept0_normalizedChi2 = mu1.normalizedChi2;
		lept0_nValidMuonHits = mu1.nValidMuonHits;
		lept0_nMatchedStations = mu1.nMatchedStations;
		lept0_nValidPixelHits = mu1.nValidPixelHits;
		lept0_nTrackerLayers = mu1.nTrackerLayers;
		lept0_muonBestTrack_dxyVTX = mu1.dxyVTX;
		lept0_muonBestTrack_dzVTX = mu1.dzVTX;
		lept0_ptError = mu1.muonBestTrack_ptError;

		lept1_pt  = mu2.pt;
		lept1_eta = mu2.eta;
		lept1_phi = mu2.phi;
		lept1_miniRelIso = mu2.miniIsoAbs/mu2.pt;
		lept1_trkIso = mu2.isolationR03_sumpt/mu2.pt;
		lept1_pdgId = 13*mu2.charge;
		lept1_normalizedChi2 = mu2.normalizedChi2;
		lept1_nValidMuonHits = mu2.nValidMuonHits;
		lept1_nMatchedStations = mu2.nMatchedStations;
		lept1_nValidPixelHits = mu2.nValidPixelHits;
		lept1_nTrackerLayers = mu2.nTrackerLayers;
		lept1_muonBestTrack_dxyVTX = mu2.dxyVTX;
		lept1_muonBestTrack_dzVTX = mu2.dzVTX;
		lept1_ptError = mu2.muonBestTrack_ptError;

		deltaR_lept = deltaR(mu1.eta,mu1.phi,mu2.eta,mu2.phi); 

		gamma_pt = pho.pt;
		gamma_eta = pho.eta;
		gamma_phi = pho.phi;
		gamma_iso = pho.ChIso;
		gamma_HoverE = pho.HoverE;
		gamma_Full5x5_SigmaIEtaIEta = pho.Full5x5_SigmaIEtaIEta;

		z_pt = z.Pt();
		z_eta = z.Eta();
		z_phi = z.Phi();
		z_mass = z.M(); 

		boss_pt = boss.Pt();
		boss_eta = boss.Eta();
		boss_phi = boss.Phi();
		boss_mass = boss.M(); 

		recoTree->Fill();
	}

}

void TreeFiller::fillGenVariables(EventSelector& selector) {
	gen_weight = selector.genWeight;	
	genmass = selector.genMass;
	const NtupleGenParticle& lept0{*selector.genLept[0]};
	const NtupleGenParticle& lept1{*selector.genLept[1]};
	const NtupleGenParticle& gamma{*selector.genPho};

	gen_leptType = fabs(lept0.id);
	//std::cout << lept0.id << ":" << lept1.id << std::endl;
	assert(gen_leptType==fabs(lept1.id));

	gen_lept0_pt = lept0.pt;
	gen_lept0_eta = lept0.eta;
	gen_lept0_phi = lept0.phi;
	gen_lept1_pt = lept1.pt;
	gen_lept1_eta = lept1.eta;
	gen_lept1_phi = lept1.phi;
	gen_gamma_pt = gamma.pt;
	gen_gamma_eta = gamma.eta;
	gen_gamma_phi = gamma.phi;

	TLorentzVector lept0_mom;
	TLorentzVector lept1_mom;
	TLorentzVector gamma_mom;
	lept0_mom.SetPtEtaPhiM(gen_lept0_pt,gen_lept0_eta,gen_lept0_phi,lept0.mass);
	lept0_mom.SetPtEtaPhiM(gen_lept0_pt,gen_lept0_eta,gen_lept0_phi,lept0.mass);
	gamma_mom.SetPtEtaPhiM(gen_gamma_pt,gen_gamma_eta,gen_gamma_phi,0);

	TLorentzVector z = (lept0_mom + lept1_mom);
	TLorentzVector boss = (lept0_mom + lept1_mom + gamma_mom);
	
	gen_z_pt = z.Pt();
	gen_z_eta = z.Eta();
	gen_z_phi = z.Phi();
	gen_z_mass = z.M(); 

	gen_boss_pt = boss.Pt();
	gen_boss_eta = boss.Eta();
	gen_boss_phi = boss.Phi();
	gen_boss_mass = boss.M(); 
	genTree->Fill();
}
