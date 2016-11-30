#ifndef StandardSelection_H
#define StandardSelection_H

#include "TLorentzVector.h"
#include "Physics/NtupleMaker/interface/NtupleClasses.h"

// template<class T> TLorentzVector momentum(const T& particle, const double& mass);
TLorentzVector momentum(const NtupleMuon& mu, const double& mass);
TLorentzVector momentum(const NtupleElectron& el, const double& mass);
TLorentzVector momentum(const NtuplePhoton& pho, const double& mass);
template<class T> bool acceptance(const T& particle, const double& ptCut, const double& etaCut);
template<class T> bool excludeECALGap(const T& particle);
template<class T> bool miniIsolation(const T& particle, const double& iso);
bool trackerIsolation(const NtupleMuon& mu, const double& iso);
bool looseMuonID(const NtupleMuon& mu, const bool& isZG);
bool tightMuonID(const NtupleMuon& mu);
bool highPtMuonID(const NtupleMuon& mu);
bool WPLoose(const NtupleElectron& el, const bool& isZG);
bool WPLoose(const NtuplePhoton& pho, const bool& isZG);

class Selector {
protected:
	const double ptCut;
	const double etaCut;
	bool pass;
public:
	Selector();
	Selector(const double& _ptCut, const double& _etaCut);
	bool passSelection() const;
};

class LeptonSelector : public Selector {
protected:
	const double extraPtCut;
	const double vetoPtCut;
public:
	LeptonSelector();
	LeptonSelector(const double& _ptCut, const double& _etaCut, const double& _extraPtCut, const double& _vetoPtCut);
};

class MuonSelector : public LeptonSelector {
protected:
	std::vector<const NtupleMuon*> selected;
public:
	MuonSelector();
	MuonSelector(const double& _ptCut, const double& _etaCut, const double& _extraPtCut, const double& _vetoPtCut);
	bool select(const std::vector<NtupleMuon>& candidates);
	const NtupleMuon& at(int i) const;
	int nMuons();
};

class ElectronSelector : public LeptonSelector { 
protected:
	std::vector<const NtupleElectron*> selected; 
public:
	ElectronSelector();
	ElectronSelector(const double& _ptCut, const double& _etaCut, const double& _extraPtCut, const double& _vetoPtCut);
	bool select(const std::vector<NtupleElectron>& candidates);
	const NtupleElectron& at(int i) const;
	int nElectrons();
};

class PhotonSelector : public Selector {
protected:
	std::vector<const NtuplePhoton*> selected;
public:
	PhotonSelector();
	PhotonSelector(const double& _ptCut, const double& _etaCut);
	bool select(const std::vector<NtuplePhoton>& candidates);
	const NtuplePhoton& at(int i) const;
	int nPhotons();
};

#endif