#ifndef EVENT_SELECTION_H
#define EVENT_SELECTION_H
#include "Physics/ZG/interface/StandardSelection.h"

class ScaleFactor {
public:
	ScaleFactor();
	~ScaleFactor();

	const static int electronNEtaBinsForID = 6;
	const static int electronNPtBinsForID = 8;

	double ELECTRON_ETABINS_ID[electronNEtaBinsForID+1] = {-2.5, -1.566, -1.444, 0, 1.444, 1.566, 2.5};
	double ELECTRON_PTBINS_ID[electronNPtBinsForID+1]   = {10, 20, 30, 40, 50, 80, 120, 200, 500};
	double ELECTRON_SF_ID[electronNEtaBinsForID][electronNPtBinsForID];
	double ELECTRON_SF_UNC_ID[electronNEtaBinsForID][electronNPtBinsForID];

	const static int electronNEtaBinsForRECO = 30;

	double ELECTRON_ETABINS_RECO[electronNEtaBinsForRECO+1] = {-2.5, -2.45, -2.4, -2.3, -2.2, -2, -1.8, -1.63, -1.566, -1.444, -1.2, -1, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 1, 1.2, 1.444, 1.566, 1.63, 1.8, 2, 2.2, 2.3, 2.4, 2.45, 2.5};
	double ELECTRON_SF_RECO[electronNEtaBinsForRECO];
	double ELECTRON_SF_UNC_RECO[electronNEtaBinsForRECO];
	
	const static int photonNEtaBins = 10;
	const static int photonNPtBins = 5;

	double PHOTON_ETABINS[photonNEtaBins+1] = {-2.5, -2, -1.566, -1.444, -0.8, 0, 0.8, 1.444, 1.566, 2, 2.5};
	double PHOTON_PTBINS[photonNPtBins+1] = {20, 35, 50, 90, 150, 500};
	double PHOTON_SF[photonNEtaBins][photonNPtBins];
	double PHOTON_SF_UNC[photonNEtaBins][photonNPtBins];

	const static int muonNAbsetaBins = 4;
	const static int muonNPtBinsForIDIso = 7;
	const static int muonNPtBinsForTrigger = 9;

	double MUON_ABSETABINS[muonNAbsetaBins+1] = {0., 0.9, 1.2, 2.1, 2.4};
	double MUON_PTBINS_IDISO[muonNPtBinsForIDIso+1] = {10, 20, 25, 30, 40, 50, 60, 120};
	double MUON_PTBINS_TRIGGER[muonNPtBinsForTrigger+1] = {10, 40, 50, 52, 55, 60, 80, 120, 200, 500};

	double MUON_ID_LOOSE_DATA[muonNAbsetaBins][muonNPtBinsForIDIso];
	double MUON_ID_LOOSE_MC[muonNAbsetaBins][muonNPtBinsForIDIso];
	double MUON_ID_HIGHPT_DATA[muonNAbsetaBins][muonNPtBinsForIDIso];
	double MUON_ID_HIGHPT_MC[muonNAbsetaBins][muonNPtBinsForIDIso];

	double MUON_ISO_LOOSE_DATA[muonNAbsetaBins][muonNPtBinsForIDIso];
	double MUON_ISO_LOOSE_MC[muonNAbsetaBins][muonNPtBinsForIDIso];
	double MUON_ISO_HIGHPT_DATA[muonNAbsetaBins][muonNPtBinsForIDIso];
	double MUON_ISO_HIGHPT_MC[muonNAbsetaBins][muonNPtBinsForIDIso];

	double MUON_TRIGGER_LOOSE_DATA[muonNAbsetaBins][muonNPtBinsForTrigger];
	double MUON_TRIGGER_LOOSE_MC[muonNAbsetaBins][muonNPtBinsForTrigger];
	double MUON_TRIGGER_HIGHPT_DATA[muonNAbsetaBins][muonNPtBinsForTrigger];
	double MUON_TRIGGER_HIGHPT_MC[muonNAbsetaBins][muonNPtBinsForTrigger];

	double MUON_ID_UNC_LOOSE_DATA[muonNAbsetaBins][muonNPtBinsForIDIso];
	double MUON_ID_UNC_LOOSE_MC[muonNAbsetaBins][muonNPtBinsForIDIso];
	double MUON_ID_UNC_HIGHPT_DATA[muonNAbsetaBins][muonNPtBinsForIDIso];
	double MUON_ID_UNC_HIGHPT_MC[muonNAbsetaBins][muonNPtBinsForIDIso];

	double MUON_ISO_UNC_LOOSE_DATA[muonNAbsetaBins][muonNPtBinsForIDIso];
	double MUON_ISO_UNC_LOOSE_MC[muonNAbsetaBins][muonNPtBinsForIDIso];
	double MUON_ISO_UNC_HIGHPT_DATA[muonNAbsetaBins][muonNPtBinsForIDIso];
	double MUON_ISO_UNC_HIGHPT_MC[muonNAbsetaBins][muonNPtBinsForIDIso];

	double MUON_TRIGGER_UNC_LOOSE_DATA[muonNAbsetaBins][muonNPtBinsForTrigger];
	double MUON_TRIGGER_UNC_LOOSE_MC[muonNAbsetaBins][muonNPtBinsForTrigger];
	double MUON_TRIGGER_UNC_HIGHPT_DATA[muonNAbsetaBins][muonNPtBinsForTrigger];
	double MUON_TRIGGER_UNC_HIGHPT_MC[muonNAbsetaBins][muonNPtBinsForTrigger];

	double getElectronScaleFactor(const NtupleElectron& electron);
	double getElectronScaleFactorForSystematics(const ElectronSelector* ee);
	double getElectronScaleFactorForID(const NtupleElectron& electron);
	double getElectronScaleFactorForRECO(const NtupleElectron& electron);
	double getElectronScaleFactorUncForID(const NtupleElectron& electron);
	double getElectronScaleFactorUncForRECO(const NtupleElectron& electron);
	double getPhotonScaleFactor(const NtuplePhoton& pho);
	double getPhotonScaleFactorForSystematics(const NtuplePhoton& pho, const bool& plusSigma);
	double getMuonScaleFactors(const MuonSelector* mm);
	double getMuonScaleFactorVariations(const MuonSelector* mm, const bool& isIDSF);
	std::vector<double> getMuonEffsForIDIso(const NtupleMuon& mu);
	std::vector<double> getMuonEffsForTrigger(const NtupleMuon& mu);
	std::vector<double> getMuonEffUncsForIDIso(const NtupleMuon& mu);
	std::vector<double> getMuonEffUncsForTrigger(const NtupleMuon& mu);

};

class EventSelector {
public:
	EventSelector( const bool& _verbose );
	EventSelector( const bool& _verbose, const double& _elPtCut, const double& _elExtraPtCut, const double& _elVetoPtCut, const double& _phoPtCut_ee, const double& _muPtCut, const double& _muExtraPtCut, const double& _muVetoPtCut, const double& _phoPtCut_mm, const double& _lowerMassCut, const double& _upperMassCut, const double& _phoPtOverMZG, const bool& _multiPho, const double& _norm );
	~EventSelector();
	void setEvent( NtupleEvent* _event );
	void selectICHEP16();
	void selectMoriond17(const int& option);
	void selectDimuGamma();
	void selectDielGamma();
	bool isTriggered( const TString& triggerName );
	void setGenMass();
	bool matchedAcceptance();
	void applyScaleFactors();
	void fillEventVariables();
	int run;
	int lumi;
	unsigned long long eventN;
	UInt_t eventNum;
	double aMCNLO;
	double weight;
	double genWeight;
	int nPU;
	int nVert;
	double genMass;
	bool Flag_badMuons;
	bool Flag_duplicatedMuons;
	bool Flag_noBadMuons;
	bool Flag_dupECALClusters;
	bool Flag_noHitsNotReplaced;
	std::vector<NtupleGenParticle*> genLept;
	NtupleGenParticle* genPho;
	TString channel;
	NtupleEvent* event;
	ElectronSelector* ee;
	MuonSelector* mm;
	PhotonSelector* gamma_ee;
	PhotonSelector* gamma_mm;
	ScaleFactor* scaleFactors;
	bool muon;
	bool electron;
	bool accepted;
	const bool verbose;
	const double elPtCut;
	const double elExtraPtCut;
	const double elVetoPtCut;
	const double phoPtCut_ee;
	const double muPtCut;
	const double muExtraPtCut;
	const double muVetoPtCut;
	const double phoPtCut_mm;
	const double lowerMassCut;
	const double upperMassCut;
	const double phoPtOverMZG;
	const bool multiPho;
	const double norm; // normalzed to 1/fb

};

double PUReweight(const int& nPU);

#endif
