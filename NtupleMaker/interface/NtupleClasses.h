#ifndef NtupleMaker_H
#define NtupleMaker_H

#include <vector>
#include <string>

const double muon_mass     = 0.1056583715;
const double electron_mass = 0.000510998;

class NtupleTrigger {

public:
    bool isFired;
    std::string name;

    NtupleTrigger(){};
    virtual ~NtupleTrigger(){};

};

class NtupleTriggerObject {

public:
    double eta;
    double phi;
    std::string name;

    NtupleTriggerObject(){};
    virtual ~NtupleTriggerObject(){};

};

class NtupleMuon {

public:
    bool isGlobalMuon;
    bool isStandAloneMuon;
    bool isTrackerMuon;
    bool isPFMuon;
    double px;
    double py;
    double pz;
    double pt;
    double phi;
    double eta;
    int charge;
    int nChambers;
    int stationMask;
    int nMatchedStations;
    double normalizedChi2;
    int nValidHits;
    int nValidTrackerHits;
    int nValidPixelHits;
    int nTrackerLayers;
    int nValidMuonHits;
    double qoverp;
    double theta;
    double lambda;
    double dxy;
    double d0;
    double dsz;
    double dz;
    double dxyBS;
    double dzBS;
    double dszBS;
    double dxyVTX;
    double dzVTX;
    double dszVTX;
    double vx;
    double vy;
    double vz;

    // muon tracks
    double muonBestTrack_px;
    double muonBestTrack_py;
    double muonBestTrack_pz;
    double muonBestTrack_pt;
    double muonBestTrack_ptError;
    int muonBestTrack_nValidPixelHits;
    int muonBestTrack_nTrackerLayers;
    int muonBestTrack_nValidMuonHits;
    double muonBestTrack_dxyVTX;
    double muonBestTrack_dzVTX;
    double muonBestTrack_dszVTX;
    double tunePMuonBestTrack_px;
    double tunePMuonBestTrack_py;
    double tunePMuonBestTrack_pz;
    double tunePMuonBestTrack_pt;
    double tunePMuonBestTrack_ptError;
    int tunePMuonBestTrack_nValidPixelHits;
    int tunePMuonBestTrack_nTrackerLayers;
    int tunePMuonBestTrack_nValidMuonHits;
    double tunePMuonBestTrack_dxyVTX;
    double tunePMuonBestTrack_dzVTX;
    double tunePMuonBestTrack_dszVTX;
    double innerTrack_px;
    double innerTrack_py;
    double innerTrack_pz;
    double innerTrack_pt;
    double innerTrack_ptError;
    int innerTrack_nValidPixelHits;
    int innerTrack_nTrackerLayers;
    double innerTrack_dxyVTX;
    double innerTrack_dzVTX;
    double innerTrack_dszVTX;
    double outerTrack_px;
    double outerTrack_py;
    double outerTrack_pz;
    double outerTrack_pt;
    double outerTrack_ptError;
    int outerTrack_nValidMuonHits;


    // muon isolation
    double miniIsoRel;
    double miniIsoAbs;
    double miniIsoCh;
    double miniIsoPh;
    double miniIsoNh;
    double miniIsoRho;
    double miniIsoPt;
    double isolationR03_sumpt;
    double isolationR03_hadEt;
    double isolationR03_emEt;
    double isolationR05_sumpt;
    double isolationR05_hadEt;
    double isolationR05_emEt;
    double PfChargedHadronIsoR05;
    double PfNeutralHadronIsoR05;
    double PfGammaIsoR05;
    double PfChargedHadronIsoR04;
    double PfNeutralHadronIsoR04;
    double PfGammaIsoR04;
    double PfPUPtR04;
    double PfChargedHadronIsoR03;
    double PfNeutralHadronIsoR03;
    double PfGammaIsoR03;
    double PfPUPtR03;

    NtupleMuon(){};
    virtual ~NtupleMuon(){};

};
   
class NtupleElectron{

public:
    double pt;
    double eta;
    double rap;
    double phi;
    double E;
    int charge;
    double enSC;
    double preEnSC;
    double rawEnSC;
    double etSC;
    double etaSC;
    double phiSC;
    double sigmaIetaIeta;
    double E1x5;
    double E2x5;
    double E5x5;
    double hOverE;
    double etaScWidth;
    double phiScWidth;
    double r9;
    double dEtaIn;
    double dPhiIn;
    double ooEmooP;
    double isoChargedHadrons;
    double isoNeutralHadrons;
    double isoPhotons;
    double isoChargedFromPU;
    double isoDeltaBeta;
    double isoRho;
    double d0;
    double dz;
    double expectedMissingInnerHits;
    bool passConversionVeto;
    double brem;
    bool passVetoId;
    bool passLooseId;
    bool passMediumId;
    bool passTightId;
    double eleInBarrel;
    double eleInEndcap;
    bool eleEcalDrivenSeed;

    double dr03ecalRecHitSumEt;
    double dr03hcalDepth1TowerSumEt;
    double dr03hcalDepth1TowerSumEtBc;
    double dr03hcalDepth2TowerSumEt;
    double dr03hcalDepth2TowerSumEtBc;
    double dr03tkSumPt;

    double dr04ecalRecHitSumEt;
    double dr04hcalDepth1TowerSumEt;
    double dr04hcalDepth1TowerSumEtBc;
    double dr04hcalDepth2TowerSumEt;
    double dr04hcalDepth2TowerSumEtBc;
    double dr04tkSumPt;

    double mvaValue;
    int mvaCategory;

    double miniIsoRel;
    double miniIsoAbs;
    double miniIsoPt;

    NtupleElectron(){};
    virtual ~NtupleElectron(){};

};


class NtupleDimuon {

public:
    int X;
    int Y;
    double openingAngle;
    double vertexFitChi2;
    double vertexFitNdof;
    double vertexFitChi2Ndof;
    double vertexFitProb;

    NtupleDimuon(){};
    virtual ~NtupleDimuon(){};

};

class NtupleDielectron {

public:
    int X;
    int Y;
    double openingAngle;
    double vertexFitChi2;
    double vertexFitNdof;
    double vertexFitChi2Ndof;
    double vertexFitProb;

    NtupleDielectron(){};
    virtual ~NtupleDielectron(){};

};

class NtupleEmu {
    
public:
    int X;
    int Y;
    double openingAngle;
    double vertexFitChi2;
    double vertexFitNdof;
    double vertexFitChi2Ndof;
    double vertexFitProb;

    NtupleEmu(){};
    virtual ~NtupleEmu(){};

};


class NtupleGenParticle {

public:
    double px;
    double py;
    double pz;
    double pt;
    double phi;
    double eta;
    double energy;
    double mass;
    int mother;
    int charge;
    int status;
    int id;
    bool fromHardProcessFinalState;
    bool fromHardProcessDecayed;
    bool fromHardProcessBeforeFSR;
    bool isHardProcess;
    bool isLastCopy;
    bool isLastCopyBeforeFSR;
    bool isPromptDecayed;
    bool isPromptFinalState;
    bool isDirectHardProcessTauDecayProductFinalState;
    bool isDirectPromptTauDecayProductFinalState;

    NtupleGenParticle(){};
    virtual ~NtupleGenParticle(){};

};

class NtuplePhoton {

public:
    bool hasPixelSeed;
    bool passElectronVeto;
    double pt;
    double eta;
    double phi;
    //double et;
    double etaSC;
    double phiSC;
    double HoverE;
    double Full5x5_SigmaIEtaIEta;
	double r9;
    bool passMediumId;
    double mvaValue;
    double mvaCategory;
    double ChIso;
    double NhIso;
    double PhIso;
    double ChIsoWithEA;
    double NhIsoWithEA;
    double PhIsoWithEA;

    NtuplePhoton(){};
    virtual ~NtuplePhoton(){};

};


class NtupleJet {

public:
    double pt;
    double eta;
    double phi;
    int charge;
    int flavor;
    double bTag;
    double CHfrac;
    double NHfrac;
    double NHEMfrac;
    double CHEMfrac;
    int CHmulti;
    int NHmulti; 

    NtupleJet(){};
    virtual ~NtupleJet(){};

};

class NtupleMET {

public:
    double pt;
    double phi;
    double px;
    double py;
    double sumEt;
    double Type1_pt;
    double Type1_phi;
    double Type1_px;
    double Type1_py;
    double Type1_sumEt;

    NtupleMET(){};
    virtual ~NtupleMET(){};

};


class NtupleEvent {

public:
    int run;
    unsigned long long event;
    int lumi;
    int nVertices;
    int nPU;    
    int nPUin;    
    int nMuons;
    int nElectrons;    
    int nPhotons;
    int nJets;    
    double rho;        
    double weight;
    int nGenParticles;

    std::vector<NtupleTrigger> triggers;
    std::vector<NtupleTriggerObject> triggerobjects;
    std::vector<NtupleMuon> muons;
    std::vector<NtupleElectron> electrons;
    std::vector<NtupleDimuon> dimuons;
    std::vector<NtupleDielectron> dielectrons;
    std::vector<NtupleEmu> emus;
    std::vector<NtuplePhoton> photons;
    std::vector<NtupleJet> jets;
    std::vector<NtupleMET> MET;
    std::vector<NtupleGenParticle> genparticles;

    NtupleEvent(){};
    virtual ~NtupleEvent(){};

};


#endif
