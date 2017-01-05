#ifndef ZGTREE_H
#define ZGTREE_H

#include "TTree.h"
#include "Physics/ZG/interface/EventSelection.h"

class TreeFiller {
public:
	TreeFiller(){}
	void setTree(TTree* _recoTree, TTree* _genTree);
	void fillVariables(EventSelector& selector);
	void fillGenVariables(EventSelector& selector);
	TTree* recoTree;
	TTree* genTree;
	
	int run;
	int lumi;
	UInt_t eventNum;
	int leptType;
	double weight;
	int nPU;
	int nVert;
	double genmass;

	double lept0_pt;
	double lept0_eta;
	double lept0_phi;
	double lept0_miniRelIso;
	double lept0_trkIso;
	double lept0_pdgId;
	double lept0_normalizedChi2;
	double lept0_nValidMuonHits;
	double lept0_nMatchedStations;
	double lept0_nValidPixelHits;
	double lept0_nTrackerLayers;
	double lept0_muonBestTrack_dxyVTX;
	double lept0_muonBestTrack_dzVTX;
	double lept0_ptError;
	double lept0_sigmaIetaIeta;
	double lept0_dEtaIn;
	double lept0_dPhiIn;
	double lept0_hOverE;
	double lept0_ooEmooP;
	double lept0_d0;
	double lept0_dz;
	double lept0_expectedMissingInnerHits;

	double lept1_pt;
	double lept1_eta;
	double lept1_phi;
	double lept1_miniRelIso;
	double lept1_trkIso;
	double lept1_pdgId;
	double lept1_normalizedChi2;
	double lept1_nValidMuonHits;
	double lept1_nMatchedStations;
	double lept1_nValidPixelHits;
	double lept1_nTrackerLayers;
	double lept1_muonBestTrack_dxyVTX;
	double lept1_muonBestTrack_dzVTX;
	double lept1_ptError;
	double lept1_sigmaIetaIeta;
	double lept1_dEtaIn;
	double lept1_dPhiIn;
	double lept1_hOverE;
	double lept1_ooEmooP;
	double lept1_d0;
	double lept1_dz;
	double lept1_expectedMissingInnerHits;

	double deltaR_lept;

	double gamma_pt;
	double gamma_eta;
	double gamma_phi;
	double gamma_iso;
	double gamma_HoverE;
	double gamma_Full5x5_SigmaIEtaIEta;

	double z_pt;
	double z_eta;
	double z_phi;
	double z_mass;

	double boss_pt;
	double boss_eta;
	double boss_phi;
	double boss_mass;

	double gen_weight;
	int gen_leptType;

	double gen_lept0_pt;
	double gen_lept0_eta;
	double gen_lept0_phi;
	double gen_lept1_pt;
	double gen_lept1_eta;
	double gen_lept1_phi;

	double gen_gamma_pt;
	double gen_gamma_eta;
	double gen_gamma_phi;

	double gen_z_pt;
	double gen_z_eta;
	double gen_z_phi;
	double gen_z_mass;

	double gen_boss_pt;
	double gen_boss_eta;
	double gen_boss_phi;
	double gen_boss_mass;
};

#endif
