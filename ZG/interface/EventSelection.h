#ifndef EVENT_SELECTION_H
#define EVENT_SELECTION_H
#include "Physics/ZG/interface/StandardSelection.h"

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
	std::vector<NtupleGenParticle*> genLept;
	NtupleGenParticle* genPho;
	TString channel;
	NtupleEvent* event;
	ElectronSelector* ee;
	MuonSelector* mm;
	PhotonSelector* gamma_ee;
	PhotonSelector* gamma_mm;
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
	const double norm;

};

double PUReweight(const int& nPU);

#endif
