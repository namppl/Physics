#ifndef EVENT_SELECTION_H
#define EVENT_SELECTION_H
#include "Physics/ZG/interface/StandardSelection.h"

class EventSelector {
public:
	EventSelector( const bool& _verbose );
	EventSelector( const bool& _verbose, const double& mu_pt_cut1, const double& mu_pt_cut2, const double& mu_pt_veto, const double& el_pt_cut1, const double& el_pt_cut2, const double& el_pt_veto, const double& pho_pt_cut, const double& _dilepton_mass_lower, const double& _dilepton_mass_upper, const double& _pho_pt_over_mass );
	~EventSelector();
	void setEvent( NtupleEvent* _event );
	void select();
	bool isTriggered( const TString& triggerName );
	void setGenMass();
	int run;
	int lumi;
	UInt_t eventNum;
	double weight;
	double genWeight;
	int nPU;
	int nVert;
	double genMass;
	std::vector<NtupleGenParticle*> genLept;
	NtupleGenParticle* genPho;
	TString channel;
	const bool verbose;
	const double dilepton_mass_lower;
	const double dilepton_mass_upper;
	const double pho_pt_over_mass;
	NtupleEvent* event;
	MuonSelector* mm;
	ElectronSelector* ee;
	PhotonSelector* gamma;
	bool muon;
	bool electron;

};

#endif
