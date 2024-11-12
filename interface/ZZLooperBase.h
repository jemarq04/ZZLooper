#ifndef ZZLOOPERBASE_H
#define ZZLOOPERBASE_H

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "Math/Vector4D.h"
#include "Math/Boost.h"

//NOTE: Only accepts 2022 skimmed ntuples
class ZZLooperBase{
  public:
    ZZLooperBase(const char *name, const char *channel, const char *filename);
    ~ZZLooperBase();

    void SetMC(bool isMC=true){_isMC = isMC;}
    void SetLumi(float lumi=1.){_lumi = lumi;}
    void SetXsec(float xsec=1.){_xsec = xsec;}
		void SetKfac(float kfac=1.){_kfac = kfac;}

    virtual void Loop() = 0;
  protected:
    void Init();

    void SetTitles(TH1 *hist, std::string x, std::string y="Number of Entries") const;
    double GetDeltaPhi(double phi1, double phi2) const;
    double GetPolCosTheta(const ROOT::Math::PtEtaPhiEVector &lep, const ROOT::Math::PtEtaPhiEVector &z) const;

    void MakeIndexFile(std::string filename) const;

    const double Z_MASS = 91.187; // GeV, given by PDG 2020

    TFile *_infile;
    TTree *_ntuple;
		Float_t _lumi=1.;
    Float_t _xsec=1., _kfac=1.;
    std::string _name;
    std::string _channel;

    bool _hasInit = false;
    bool _isMC = false;

    //Shared branches
    Float_t summedWeights;
    //Float_t Energy;
    //Float_t EnergyNoFSR;
    //Float_t Eta;
    //Float_t EtaNoFSR;
    //Float_t Flag_BadChargedCandidateFilterPass;
    //Float_t Flag_BadPFMuonFilterPass;
    Float_t Mass;
    //Float_t MassNoFSR;
    //Float_t MtToMET;
    //Float_t Phi;
    //Float_t PhiNoFSR;
    //Float_t Pt;
    //Float_t PtNoFSR;
    //Float_t deltaEtajj;
    //Float_t deltaEtajj_jerDown;
    //Float_t deltaEtajj_jerUp;
    //Float_t deltaEtajj_jesDown;
    //Float_t deltaEtajj_jesUp;
    //Float_t deltaPhiTojj;
    //Float_t deltaPhiTojj_jerDown;
    //Float_t deltaPhiTojj_jerUp;
    //Float_t deltaPhiTojj_jesDown;
    //Float_t deltaPhiTojj_jesUp;
    //Float_t etajj;
    //Float_t etajj_jerDown;
    //Float_t etajj_jerUp;
    //Float_t etajj_jesDown;
    //Float_t etajj_jesUp;
    Float_t genWeight;
    //Float_t jetPUSFmulfac;
    //Float_t maxScaleWeight;
    //Float_t minScaleWeight;
    //Float_t mjj;
    //Float_t mjj_jerDown;
    //Float_t mjj_jerUp;
    //Float_t mjj_jesDown;
    //Float_t mjj_jesUp;
    //Float_t nTruePU;
    //Float_t originalXWGTUP;
    //Float_t phijj;
    //Float_t phijj_jerDown;
    //Float_t phijj_jerUp;
    //Float_t phijj_jesDown;
    //Float_t phijj_jesUp;
    //Float_t ptjj;
    //Float_t ptjj_jerDown;
    //Float_t ptjj_jerUp;
    //Float_t ptjj_jesDown;
    //Float_t ptjj_jesUp;
    //Float_t pvRho;
    //Float_t pvZ;
    //Float_t pvndof;
    //Float_t rapidity;
    //Float_t type1_pfMETEt;
    //Float_t type1_pfMETPhi;
    //Float_t zeppenfeld;
    //Float_t zeppenfeld_jerDown;
    //Float_t zeppenfeld_jerUp;
    //Float_t zeppenfeld_jesDown;
    //Float_t zeppenfeld_jesUp;
    //Float_t zeppenfeldj3;
    //Float_t zeppenfeldj3_jerDown;
    //Float_t zeppenfeldj3_jerUp;
    //Float_t zeppenfeldj3_jesDown;
    //Float_t zeppenfeldj3_jesUp;
    //Bool_t pvIdFake;
    //Bool_t pvIsValid;
    //Int_t Charge;
    //Int_t PdgId;
    UInt_t lumi;
    //UInt_t nJets;
    //UInt_t nJets_jerDown;
    //UInt_t nJets_jerUp;
    //UInt_t nJets_jesDown;
    //UInt_t nJets_jesUp;
    //UInt_t nZZLooseElec;
    //UInt_t nZZLooseMu;
    //UInt_t nZZTightElec;
    //UInt_t nZZTightIsoElec;
    //UInt_t nZZTightIsoMu;
    //UInt_t nZZTightMu;
    UInt_t nvtx;
    UInt_t run;
    ULong64_t evt;
    //std::vector<float> *jetDeepCSV = 0;
    //std::vector<float> *jetEta = 0;
    //std::vector<float> *jetEta_jerDown = 0;
    //std::vector<float> *jetEta_jerUp = 0;
    //std::vector<float> *jetEta_jesDown = 0;
    //std::vector<float> *jetEta_jesUp = 0;
    //std::vector<float> *jetPhi = 0;
    //std::vector<float> *jetPt = 0;
    //std::vector<float> *jetPt_jerDown = 0;
    //std::vector<float> *jetPt_jerUp = 0;
    //std::vector<float> *jetPt_jesDown = 0;
    //std::vector<float> *jetPt_jesUp = 0;
    //std::vector<float> *jetRapidity = 0;
    //std::vector<float> *jetRapidity_jerDown = 0;
    //std::vector<float> *jetRapidity_jerUp = 0;
    //std::vector<float> *jetRapidity_jesDown = 0;
    //std::vector<float> *jetRapidity_jesUp = 0;
    //std::vector<float> *pdfWeights = 0;
    //std::vector<float> *scaleWeightIDs = 0;
    //std::vector<float> *scaleWeights = 0;
    //std::vector<int> *isGenJetMatched = 0;
    //std::vector<int> *jetHadronFlavor = 0;
    //std::vector<int> *jetPUID = 0;
    //std::vector<int> *jetPUID_jerDown = 0;
    //std::vector<int> *jetPUID_jerUp = 0;
    //std::vector<int> *jetPUID_jesDown = 0;
    //std::vector<int> *jetPUID_jesUp = 0;
    //Bool_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVLPass;
    //Bool_t HLT_DoubleEle25_CaloIdL_MWPass;
    //Bool_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8Pass;
    //Bool_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVLPass;
    //Bool_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass;
    //Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass;
    //Bool_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass;
    //Bool_t HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZPass;
    //Bool_t HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZPass;
    //Bool_t HLT_TripleMu_10_5_5_DZPass;
    //Bool_t HLT_TripleMu_12_10_5Pass;
    //Bool_t HLT_Ele30_WPTight_GsfPass;
    //Bool_t HLT_Ele32_WPTight_GsfPass;
    //Bool_t HLT_IsoMu24Pass;
    //Bool_t Flag_goodVerticesPass;
    //Bool_t Flag_globalSuperTightHalo2016FilterPass;
    //Bool_t Flag_EcalDeadCellTriggerPrimitiveFilterPass;
    //Float_t Flag_BadPFMuonFilterPass;
    //Bool_t Flag_BadPFMuonDzFilterPass;
    //Bool_t Flag_hfNoisyHitsFilterPass;
    //Bool_t Flag_eeBadScFilterPass;

    //eeee
    //Float_t e1_e3_DR;
    //Float_t e1_e3_Mass;
    //Float_t e1_e3_MassNoFSR;
    //Float_t e1_e4_DR;
    //Float_t e1_e4_Mass;
    //Float_t e1_e4_MassNoFSR;
    //Float_t e2_e3_DR;
    //Float_t e2_e3_Mass;
    //Float_t e2_e3_MassNoFSR;
    //Float_t e2_e4_DR;
    //Float_t e2_e4_Mass;
    //Float_t e2_e4_MassNoFSR;
    //Bool_t e1_e3_SS;
    //Bool_t e1_e4_SS;
    //Bool_t e2_e3_SS;
    //Bool_t e2_e4_SS;
    //Float_t e3_e4_DR;
    //Float_t e3_e4_Energy;
    //Float_t e3_e4_EnergyNoFSR;
    //Float_t e3_e4_Eta;
    //Float_t e3_e4_EtaNoFSR;
    Float_t e3_e4_Mass;
    //Float_t e3_e4_MassNoFSR;
    //Float_t e3_e4_Mt;
    //Float_t e3_e4_MtToMET;
    //Float_t e3_e4_Phi;
    //Float_t e3_e4_PhiNoFSR;
    //Float_t e3_e4_Pt;
    //Float_t e3_e4_PtNoFSR;
    //Bool_t e3_e4_SS;
    //Int_t e3_e4_Charge;
    //Int_t e3_e4_PdgId;
    //Float_t e3EffectiveArea;
    Float_t e3Energy;
    Float_t e3Eta;
    //Float_t e3FSRDREt;
    //Float_t e3FSREt;
    //Float_t e3FSREta;
    //Float_t e3FSRPhi;
    //Float_t e3GenEta;
    //Float_t e3GenPhi;
    //Float_t e3GenPt;
    //Float_t e3IP2D;
    //Float_t e3IP2DUncertainty;
    //Float_t e3IP3D;
    //Float_t e3IP3DUncertainty;
    //Float_t e3Loose;
    //Float_t e3Mass;
    //Float_t e3Medium;
    //Float_t e3MtToMET;
    //Float_t e3PFChargedIso;
    //Float_t e3PFNeutralIso;
    //Float_t e3PFPUIso;
    //Float_t e3PFPhotonIso;
    //Float_t e3PVDXY;
    //Float_t e3PVDZ;
    Float_t e3Phi;
    Float_t e3Pt;
    //Float_t e3RelPFIsoRho;
    //Float_t e3Rho;
    //Float_t e3SCEnergy;
    //Float_t e3SCEta;
    //Float_t e3SCPhi;
    //Float_t e3SCRawEnergy;
    //Float_t e3SIP2D;
    //Float_t e3SIP3D;
    //Float_t e3ScaleGainDn;
    //Float_t e3ScaleGainUp;
    //Float_t e3ScaleStatDn;
    //Float_t e3ScaleStatUp;
    //Float_t e3ScaleSystDn;
    //Float_t e3ScaleSystUp;
    //Float_t e3ScaleTotDn;
    //Float_t e3ScaleTotUp;
    //Float_t e3SigmaPhiDn;
    //Float_t e3SigmaPhiUp;
    //Float_t e3SigmaRhoDn;
    //Float_t e3SigmaRhoUp;
    //Float_t e3SigmaTotDn;
    //Float_t e3SigmaTotUp;
    //Float_t e3Tight;
    //Float_t e3UnCorrPt;
    //Float_t e3Veto;
    //Float_t e3ZZIso;
    //Bool_t e3IsEB;
    //Bool_t e3IsGap;
    //Bool_t e3ZZIsoPass;
    //Bool_t e3ZZLooseID;
    //Bool_t e3ZZLooseIDNoVtx;
    //Bool_t e3ZZTightID;
    //Bool_t e3ZZTightIDNoVtx;
    //Int_t e3Charge;
    //Int_t e3GenCharge;
    //Int_t e3GenPdgId;
    Int_t e3PdgId;
    //UInt_t e3MissingHits;
    //Float_t e4EffectiveArea;
    Float_t e4Energy;
    Float_t e4Eta;
    //Float_t e4FSRDREt;
    //Float_t e4FSREt;
    //Float_t e4FSREta;
    //Float_t e4FSRPhi;
    //Float_t e4GenEta;
    //Float_t e4GenPhi;
    //Float_t e4GenPt;
    //Float_t e4IP2D;
    //Float_t e4IP2DUncertainty;
    //Float_t e4IP3D;
    //Float_t e4IP3DUncertainty;
    //Float_t e4Loose;
    //Float_t e4Mass;
    //Float_t e4Medium;
    //Float_t e4MtToMET;
    //Float_t e4PFChargedIso;
    //Float_t e4PFNeutralIso;
    //Float_t e4PFPUIso;
    //Float_t e4PFPhotonIso;
    //Float_t e4PVDXY;
    //Float_t e4PVDZ;
    Float_t e4Phi;
    Float_t e4Pt;
    //Float_t e4RelPFIsoRho;
    //Float_t e4Rho;
    //Float_t e4SCEnergy;
    //Float_t e4SCEta;
    //Float_t e4SCPhi;
    //Float_t e4SCRawEnergy;
    //Float_t e4SIP2D;
    //Float_t e4SIP3D;
    //Float_t e4ScaleGainDn;
    //Float_t e4ScaleGainUp;
    //Float_t e4ScaleStatDn;
    //Float_t e4ScaleStatUp;
    //Float_t e4ScaleSystDn;
    //Float_t e4ScaleSystUp;
    //Float_t e4ScaleTotDn;
    //Float_t e4ScaleTotUp;
    //Float_t e4SigmaPhiDn;
    //Float_t e4SigmaPhiUp;
    //Float_t e4SigmaRhoDn;
    //Float_t e4SigmaRhoUp;
    //Float_t e4SigmaTotDn;
    //Float_t e4SigmaTotUp;
    //Float_t e4Tight;
    //Float_t e4UnCorrPt;
    //Float_t e4Veto;
    //Float_t e4ZZIso;
    //Bool_t e4IsEB;
    //Bool_t e4IsGap;
    //Bool_t e4ZZIsoPass;
    //Bool_t e4ZZLooseID;
    //Bool_t e4ZZLooseIDNoVtx;
    //Bool_t e4ZZTightID;
    //Bool_t e4ZZTightIDNoVtx;
    //Int_t e4Charge;
    //Int_t e4GenCharge;
    //Int_t e4GenPdgId;
    Int_t e4PdgId;
    //UInt_t e4MissingHits;

    //eemm
    //Float_t e1_m1_DR;
    //Float_t e1_m1_Mass;
    //Float_t e1_m1_MassNoFSR;
    //Float_t e1_m2_DR;
    //Float_t e1_m2_Mass;
    //Float_t e1_m2_MassNoFSR;
    //Float_t e2_m1_DR;
    //Float_t e2_m1_Mass;
    //Float_t e2_m1_MassNoFSR;
    //Float_t e2_m2_DR;
    //Float_t e2_m2_Mass;
    //Float_t e2_m2_MassNoFSR;
    //Bool_t e1_m1_SS;
    //Bool_t e1_m2_SS;
    //Bool_t e2_m1_SS;
    //Bool_t e2_m2_SS;

    //mmmm
    //Float_t m1_m3_DR;
    //Float_t m1_m3_Mass;
    //Float_t m1_m3_MassNoFSR;
    //Float_t m1_m4_DR;
    //Float_t m1_m4_Mass;
    //Float_t m1_m4_MassNoFSR;
    //Float_t m2_m3_DR;
    //Float_t m2_m3_Mass;
    //Float_t m2_m3_MassNoFSR;
    //Float_t m2_m4_DR;
    //Float_t m2_m4_Mass;
    //Float_t m2_m4_MassNoFSR;
    //Bool_t m1_m3_SS;
    //Bool_t m1_m4_SS;
    //Bool_t m2_m3_SS;
    //Bool_t m2_m4_SS;
    //Float_t m3_m4_DR;
    //Float_t m3_m4_Energy;
    //Float_t m3_m4_EnergyNoFSR;
    //Float_t m3_m4_Eta;
    //Float_t m3_m4_EtaNoFSR;
    Float_t m3_m4_Mass;
    //Float_t m3_m4_MassNoFSR;
    //Float_t m3_m4_Mt;
    //Float_t m3_m4_MtToMET;
    //Float_t m3_m4_Phi;
    //Float_t m3_m4_PhiNoFSR;
    //Float_t m3_m4_Pt;
    //Float_t m3_m4_PtNoFSR;
    //Bool_t m3_m4_SS;
    //Int_t m3_m4_Charge;
    //Int_t m3_m4_PdgId;
    //Float_t m3EffScaleFactor;
    //Float_t m3EffScaleFactorError;
    Float_t m3Energy;
    Float_t m3Eta;
    //Float_t m3FSRDREt;
    //Float_t m3FSREt;
    //Float_t m3FSREta;
    //Float_t m3FSRPhi;
    //Float_t m3GenEta;
    //Float_t m3GenPhi;
    //Float_t m3GenPt;
    //Float_t m3IP2D;
    //Float_t m3IP2DUncertainty;
    //Float_t m3IP3D;
    //Float_t m3IP3DUncertainty;
    //Float_t m3Mass;
    //Float_t m3MtToMET;
    //Float_t m3PFChargedIso;
    //Float_t m3PFNeutralIso;
    //Float_t m3PFPUIso;
    //Float_t m3PFPhotonIso;
    //Float_t m3PVDXY;
    //Float_t m3PVDZ;
    Float_t m3Phi;
    Float_t m3Pt;
    //Float_t m3PtErr;
    //Float_t m3PtErrUncorrected;
    //Float_t m3PtUncorrected;
    //Float_t m3RelPFIsoDBR03;
    //Float_t m3RelPFIsoDBR04;
    //Float_t m3RochesterScaleUncDn;
    //Float_t m3RochesterScaleUncUp;
    //Float_t m3RochesterSmearUncDn;
    //Float_t m3RochesterSmearUncUp;
    //Float_t m3SIP2D;
    //Float_t m3SIP3D;
    //Float_t m3TrackIso;
    //Float_t m3ZZIso;
    //Bool_t m3HighPtID;
    //Bool_t m3HighPtIDNoVtx;
    //Bool_t m3IsGlobal;
    //Bool_t m3IsLoose;
    //Bool_t m3IsMedium;
    //Bool_t m3IsPFMuon;
    //Bool_t m3IsTracker;
    //Bool_t m3PASTightID;
    //Bool_t m3PASTightIDNoVtx;
    //Bool_t m3PFID;
    //Bool_t m3PFIDNoVtx;
    //Bool_t m3ZZIsoPass;
    //Bool_t m3ZZLooseID;
    //Bool_t m3ZZLooseIDNoVtx;
    //Bool_t m3ZZTightID;
    //Bool_t m3ZZTightIDNoVtx;
    //Int_t m3Charge;
    //Int_t m3GenCharge;
    //Int_t m3GenPdgId;
    Int_t m3PdgId;
    //UInt_t m3BestTrackType;
    //UInt_t m3CutBasedLoose;
    //UInt_t m3CutBasedMedium;
    //UInt_t m3CutBasedTight;
    //UInt_t m3IsTight;
    //UInt_t m3MatchedStations;
    //UInt_t m3NoOfMatches;
    //UInt_t m3PFIsoLoose;
    //UInt_t m3PFIsoMedium;
    //UInt_t m3PFIsoTight;
    //UInt_t m3PFIsoVTight;
    //Float_t m4EffScaleFactor;
    //Float_t m4EffScaleFactorError;
    Float_t m4Energy;
    Float_t m4Eta;
    //Float_t m4FSRDREt;
    //Float_t m4FSREt;
    //Float_t m4FSREta;
    //Float_t m4FSRPhi;
    //Float_t m4GenEta;
    //Float_t m4GenPhi;
    //Float_t m4GenPt;
    //Float_t m4IP2D;
    //Float_t m4IP2DUncertainty;
    //Float_t m4IP3D;
    //Float_t m4IP3DUncertainty;
    //Float_t m4Mass;
    //Float_t m4MtToMET;
    //Float_t m4PFChargedIso;
    //Float_t m4PFNeutralIso;
    //Float_t m4PFPUIso;
    //Float_t m4PFPhotonIso;
    //Float_t m4PVDXY;
    //Float_t m4PVDZ;
    Float_t m4Phi;
    Float_t m4Pt;
    //Float_t m4PtErr;
    //Float_t m4PtErrUncorrected;
    //Float_t m4PtUncorrected;
    //Float_t m4RelPFIsoDBR03;
    //Float_t m4RelPFIsoDBR04;
    //Float_t m4RochesterScaleUncDn;
    //Float_t m4RochesterScaleUncUp;
    //Float_t m4RochesterSmearUncDn;
    //Float_t m4RochesterSmearUncUp;
    //Float_t m4SIP2D;
    //Float_t m4SIP3D;
    //Float_t m4TrackIso;
    //Float_t m4ZZIso;
    //Bool_t m4HighPtID;
    //Bool_t m4HighPtIDNoVtx;
    //Bool_t m4IsGlobal;
    //Bool_t m4IsLoose;
    //Bool_t m4IsMedium;
    //Bool_t m4IsPFMuon;
    //Bool_t m4IsTracker;
    //Bool_t m4PASTightID;
    //Bool_t m4PASTightIDNoVtx;
    //Bool_t m4PFID;
    //Bool_t m4PFIDNoVtx;
    //Bool_t m4ZZIsoPass;
    //Bool_t m4ZZLooseID;
    //Bool_t m4ZZLooseIDNoVtx;
    //Bool_t m4ZZTightID;
    //Bool_t m4ZZTightIDNoVtx;
    //Int_t m4Charge;
    //Int_t m4GenCharge;
    //Int_t m4GenPdgId;
    Int_t m4PdgId;
    //UInt_t m4BestTrackType;
    //UInt_t m4CutBasedLoose;
    //UInt_t m4CutBasedMedium;
    //UInt_t m4CutBasedTight;
    //UInt_t m4IsTight;
    //UInt_t m4MatchedStations;
    //UInt_t m4NoOfMatches;
    //UInt_t m4PFIsoLoose;
    //UInt_t m4PFIsoMedium;
    //UInt_t m4PFIsoTight;
    //UInt_t m4PFIsoVTight;

    //eeee or eemm
    //Float_t e1_e2_DR;
    //Float_t e1_e2_Energy;
    //Float_t e1_e2_EnergyNoFSR;
    //Float_t e1_e2_Eta;
    //Float_t e1_e2_EtaNoFSR;
    Float_t e1_e2_Mass;
    //Float_t e1_e2_MassNoFSR;
    //Float_t e1_e2_Mt;
    //Float_t e1_e2_MtToMET;
    //Float_t e1_e2_Phi;
    //Float_t e1_e2_PhiNoFSR;
    //Float_t e1_e2_Pt;
    //Float_t e1_e2_PtNoFSR;
    //Bool_t e1_e2_SS;
    //Int_t e1_e2_Charge;
    //Int_t e1_e2_PdgId;
    //Float_t e1EffectiveArea;
    Float_t e1Energy;
    Float_t e1Eta;
    //Float_t e1FSRDREt;
    //Float_t e1FSREt;
    //Float_t e1FSREta;
    //Float_t e1FSRPhi;
    //Float_t e1GenEta;
    //Float_t e1GenPhi;
    //Float_t e1GenPt;
    //Float_t e1IP2D;
    //Float_t e1IP2DUncertainty;
    //Float_t e1IP3D;
    //Float_t e1IP3DUncertainty;
    //Float_t e1Loose;
    //Float_t e1Mass;
    //Float_t e1Medium;
    //Float_t e1MtToMET;
    //Float_t e1PFChargedIso;
    //Float_t e1PFNeutralIso;
    //Float_t e1PFPUIso;
    //Float_t e1PFPhotonIso;
    //Float_t e1PVDXY;
    //Float_t e1PVDZ;
    Float_t e1Phi;
    Float_t e1Pt;
    //Float_t e1RelPFIsoRho;
    //Float_t e1Rho;
    //Float_t e1SCEnergy;
    //Float_t e1SCEta;
    //Float_t e1SCPhi;
    //Float_t e1SCRawEnergy;
    //Float_t e1SIP2D;
    //Float_t e1SIP3D;
    //Float_t e1ScaleGainDn;
    //Float_t e1ScaleGainUp;
    //Float_t e1ScaleStatDn;
    //Float_t e1ScaleStatUp;
    //Float_t e1ScaleSystDn;
    //Float_t e1ScaleSystUp;
    //Float_t e1ScaleTotDn;
    //Float_t e1ScaleTotUp;
    //Float_t e1SigmaPhiDn;
    //Float_t e1SigmaPhiUp;
    //Float_t e1SigmaRhoDn;
    //Float_t e1SigmaRhoUp;
    //Float_t e1SigmaTotDn;
    //Float_t e1SigmaTotUp;
    //Float_t e1Tight;
    //Float_t e1UnCorrPt;
    //Float_t e1Veto;
    //Float_t e1ZZIso;
    //Bool_t e1IsEB;
    //Bool_t e1IsGap;
    //Bool_t e1ZZIsoPass;
    //Bool_t e1ZZLooseID;
    //Bool_t e1ZZLooseIDNoVtx;
    //Bool_t e1ZZTightID;
    //Bool_t e1ZZTightIDNoVtx;
    //Int_t e1Charge;
    //Int_t e1GenCharge;
    //Int_t e1GenPdgId;
    Int_t e1PdgId;
    //UInt_t e1MissingHits;
    //Float_t e2EffectiveArea;
    Float_t e2Energy;
    Float_t e2Eta;
    //Float_t e2FSRDREt;
    //Float_t e2FSREt;
    //Float_t e2FSREta;
    //Float_t e2FSRPhi;
    //Float_t e2GenEta;
    //Float_t e2GenPhi;
    //Float_t e2GenPt;
    //Float_t e2IP2D;
    //Float_t e2IP2DUncertainty;
    //Float_t e2IP3D;
    //Float_t e2IP3DUncertainty;
    //Float_t e2Loose;
    //Float_t e2Mass;
    //Float_t e2Medium;
    //Float_t e2MtToMET;
    //Float_t e2PFChargedIso;
    //Float_t e2PFNeutralIso;
    //Float_t e2PFPUIso;
    //Float_t e2PFPhotonIso;
    //Float_t e2PVDXY;
    //Float_t e2PVDZ;
    Float_t e2Phi;
    Float_t e2Pt;
    //Float_t e2RelPFIsoRho;
    //Float_t e2Rho;
    //Float_t e2SCEnergy;
    //Float_t e2SCEta;
    //Float_t e2SCPhi;
    //Float_t e2SCRawEnergy;
    //Float_t e2SIP2D;
    //Float_t e2SIP3D;
    //Float_t e2ScaleGainDn;
    //Float_t e2ScaleGainUp;
    //Float_t e2ScaleStatDn;
    //Float_t e2ScaleStatUp;
    //Float_t e2ScaleSystDn;
    //Float_t e2ScaleSystUp;
    //Float_t e2ScaleTotDn;
    //Float_t e2ScaleTotUp;
    //Float_t e2SigmaPhiDn;
    //Float_t e2SigmaPhiUp;
    //Float_t e2SigmaRhoDn;
    //Float_t e2SigmaRhoUp;
    //Float_t e2SigmaTotDn;
    //Float_t e2SigmaTotUp;
    //Float_t e2Tight;
    //Float_t e2UnCorrPt;
    //Float_t e2Veto;
    //Float_t e2ZZIso;
    //Bool_t e2IsEB;
    //Bool_t e2IsGap;
    //Bool_t e2ZZIsoPass;
    //Bool_t e2ZZLooseID;
    //Bool_t e2ZZLooseIDNoVtx;
    //Bool_t e2ZZTightID;
    //Bool_t e2ZZTightIDNoVtx;
    //Int_t e2Charge;
    //Int_t e2GenCharge;
    //Int_t e2GenPdgId;
    Int_t e2PdgId;
    //UInt_t e2MissingHits;

    //mmmm or eemm
    //Float_t m1_m2_DR;
    //Float_t m1_m2_Energy;
    //Float_t m1_m2_EnergyNoFSR;
    //Float_t m1_m2_Eta;
    //Float_t m1_m2_EtaNoFSR;
    Float_t m1_m2_Mass;
    //Float_t m1_m2_MassNoFSR;
    //Float_t m1_m2_Mt;
    //Float_t m1_m2_MtToMET;
    //Float_t m1_m2_Phi;
    //Float_t m1_m2_PhiNoFSR;
    //Float_t m1_m2_Pt;
    //Float_t m1_m2_PtNoFSR;
    //Bool_t m1_m2_SS;
    //Int_t m1_m2_Charge;
    //Int_t m1_m2_PdgId;
    //Float_t m1EffScaleFactor;
    //Float_t m1EffScaleFactorError;
    Float_t m1Energy;
    Float_t m1Eta;
    //Float_t m1FSRDREt;
    //Float_t m1FSREt;
    //Float_t m1FSREta;
    //Float_t m1FSRPhi;
    //Float_t m1GenEta;
    //Float_t m1GenPhi;
    //Float_t m1GenPt;
    //Float_t m1IP2D;
    //Float_t m1IP2DUncertainty;
    //Float_t m1IP3D;
    //Float_t m1IP3DUncertainty;
    //Float_t m1Mass;
    //Float_t m1MtToMET;
    //Float_t m1PFChargedIso;
    //Float_t m1PFNeutralIso;
    //Float_t m1PFPUIso;
    //Float_t m1PFPhotonIso;
    //Float_t m1PVDXY;
    //Float_t m1PVDZ;
    Float_t m1Phi;
    Float_t m1Pt;
    //Float_t m1PtErr;
    //Float_t m1PtErrUncorrected;
    //Float_t m1PtUncorrected;
    //Float_t m1RelPFIsoDBR03;
    //Float_t m1RelPFIsoDBR04;
    //Float_t m1RochesterScaleUncDn;
    //Float_t m1RochesterScaleUncUp;
    //Float_t m1RochesterSmearUncDn;
    //Float_t m1RochesterSmearUncUp;
    //Float_t m1SIP2D;
    //Float_t m1SIP3D;
    //Float_t m1TrackIso;
    //Float_t m1ZZIso;
    //Bool_t m1HighPtID;
    //Bool_t m1HighPtIDNoVtx;
    //Bool_t m1IsGlobal;
    //Bool_t m1IsLoose;
    //Bool_t m1IsMedium;
    //Bool_t m1IsPFMuon;
    //Bool_t m1IsTracker;
    //Bool_t m1PASTightID;
    //Bool_t m1PASTightIDNoVtx;
    //Bool_t m1PFID;
    //Bool_t m1PFIDNoVtx;
    //Bool_t m1ZZIsoPass;
    //Bool_t m1ZZLooseID;
    //Bool_t m1ZZLooseIDNoVtx;
    //Bool_t m1ZZTightID;
    //Bool_t m1ZZTightIDNoVtx;
    //Int_t m1Charge;
    //Int_t m1GenCharge;
    //Int_t m1GenPdgId;
    Int_t m1PdgId;
    //UInt_t m1BestTrackType;
    //UInt_t m1CutBasedLoose;
    //UInt_t m1CutBasedMedium;
    //UInt_t m1CutBasedTight;
    //UInt_t m1IsTight;
    //UInt_t m1MatchedStations;
    //UInt_t m1NoOfMatches;
    //UInt_t m1PFIsoLoose;
    //UInt_t m1PFIsoMedium;
    //UInt_t m1PFIsoTight;
    //UInt_t m1PFIsoVTight;
    //Float_t m2EffScaleFactor;
    //Float_t m2EffScaleFactorError;
    Float_t m2Energy;
    Float_t m2Eta;
    //Float_t m2FSRDREt;
    //Float_t m2FSREt;
    //Float_t m2FSREta;
    //Float_t m2FSRPhi;
    //Float_t m2GenEta;
    //Float_t m2GenPhi;
    //Float_t m2GenPt;
    //Float_t m2IP2D;
    //Float_t m2IP2DUncertainty;
    //Float_t m2IP3D;
    //Float_t m2IP3DUncertainty;
    //Float_t m2Mass;
    //Float_t m2MtToMET;
    //Float_t m2PFChargedIso;
    //Float_t m2PFNeutralIso;
    //Float_t m2PFPUIso;
    //Float_t m2PFPhotonIso;
    //Float_t m2PVDXY;
    //Float_t m2PVDZ;
    Float_t m2Phi;
    Float_t m2Pt;
    //Float_t m2PtErr;
    //Float_t m2PtErrUncorrected;
    //Float_t m2PtUncorrected;
    //Float_t m2RelPFIsoDBR03;
    //Float_t m2RelPFIsoDBR04;
    //Float_t m2RochesterScaleUncDn;
    //Float_t m2RochesterScaleUncUp;
    //Float_t m2RochesterSmearUncDn;
    //Float_t m2RochesterSmearUncUp;
    //Float_t m2SIP2D;
    //Float_t m2SIP3D;
    //Float_t m2TrackIso;
    //Float_t m2ZZIso;
    //Bool_t m2HighPtID;
    //Bool_t m2HighPtIDNoVtx;
    //Bool_t m2IsGlobal;
    //Bool_t m2IsLoose;
    //Bool_t m2IsMedium;
    //Bool_t m2IsPFMuon;
    //Bool_t m2IsTracker;
    //Bool_t m2PASTightID;
    //Bool_t m2PASTightIDNoVtx;
    //Bool_t m2PFID;
    //Bool_t m2PFIDNoVtx;
    //Bool_t m2ZZIsoPass;
    //Bool_t m2ZZLooseID;
    //Bool_t m2ZZLooseIDNoVtx;
    //Bool_t m2ZZTightID;
    //Bool_t m2ZZTightIDNoVtx;
    //Int_t m2Charge;
    //Int_t m2GenCharge;
    //Int_t m2GenPdgId;
    Int_t m2PdgId;
    //UInt_t m2BestTrackType;
    //UInt_t m2CutBasedLoose;
    //UInt_t m2CutBasedMedium;
    //UInt_t m2CutBasedTight;
    //UInt_t m2IsTight;
    //UInt_t m2MatchedStations;
    //UInt_t m2NoOfMatches;
    //UInt_t m2PFIsoLoose;
    //UInt_t m2PFIsoMedium;
    //UInt_t m2PFIsoTight;
    //UInt_t m2PFIsoVTight;
};

ZZLooperBase::ZZLooperBase(const char *name, const char *channel, const char *filename)
  : _name(name), _channel(channel){
  _infile = TFile::Open(filename);
  _ntuple = (TTree*)_infile->Get((_channel + "/ntuple").c_str());
}

ZZLooperBase::~ZZLooperBase(){_infile->Close();}

void ZZLooperBase::Init(){
  if (_hasInit) return;
  _hasInit = true;
  //Shared branches
  if (_isMC) _ntuple->SetBranchAddress("summedWeights", &summedWeights);
  //_ntuple->SetBranchAddress("Energy", &Energy);
  //_ntuple->SetBranchAddress("EnergyNoFSR", &EnergyNoFSR);
  //_ntuple->SetBranchAddress("Eta", &Eta);
  //_ntuple->SetBranchAddress("EtaNoFSR", &EtaNoFSR);
  //_ntuple->SetBranchAddress("Flag_BadChargedCandidateFilterPass", &Flag_BadChargedCandidateFilterPass);
  //_ntuple->SetBranchAddress("Flag_BadPFMuonFilterPass", &Flag_BadPFMuonFilterPass);
  _ntuple->SetBranchAddress("Mass", &Mass);
  //_ntuple->SetBranchAddress("MassNoFSR", &MassNoFSR);
  //_ntuple->SetBranchAddress("MtToMET", &MtToMET);
  //_ntuple->SetBranchAddress("Phi", &Phi);
  //_ntuple->SetBranchAddress("PhiNoFSR", &PhiNoFSR);
  //_ntuple->SetBranchAddress("Pt", &Pt);
  //_ntuple->SetBranchAddress("PtNoFSR", &PtNoFSR);
  //_ntuple->SetBranchAddress("deltaEtajj", &deltaEtajj);
  if (_isMC){
    //_ntuple->SetBranchAddress("deltaEtajj_jerDown", &deltaEtajj_jerDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("deltaEtajj_jerUp", &deltaEtajj_jerUp);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("deltaEtajj_jesDown", &deltaEtajj_jesDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("deltaEtajj_jesUp", &deltaEtajj_jesUp);
  }
  //_ntuple->SetBranchAddress("deltaPhiTojj", &deltaPhiTojj);
  if (_isMC){
    //_ntuple->SetBranchAddress("deltaPhiTojj_jerDown", &deltaPhiTojj_jerDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("deltaPhiTojj_jerUp", &deltaPhiTojj_jerUp);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("deltaPhiTojj_jesDown", &deltaPhiTojj_jesDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("deltaPhiTojj_jesUp", &deltaPhiTojj_jesUp);
  }
  //_ntuple->SetBranchAddress("etajj", &etajj);
  if (_isMC){
    //_ntuple->SetBranchAddress("etajj_jerDown", &etajj_jerDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("etajj_jerUp", &etajj_jerUp);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("etajj_jesDown", &etajj_jesDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("etajj_jesUp", &etajj_jesUp);
  }
  if (_isMC){
    _ntuple->SetBranchAddress("genWeight", &genWeight);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("jetPUSFmulfac", &jetPUSFmulfac);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("maxScaleWeight", &maxScaleWeight);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("minScaleWeight", &minScaleWeight);
  }
  //_ntuple->SetBranchAddress("mjj", &mjj);
  if (_isMC){
    //_ntuple->SetBranchAddress("mjj_jerDown", &mjj_jerDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("mjj_jerUp", &mjj_jerUp);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("mjj_jesDown", &mjj_jesDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("mjj_jesUp", &mjj_jesUp);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("nTruePU", &nTruePU);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("originalXWGTUP", &originalXWGTUP);
  }
  //_ntuple->SetBranchAddress("phijj", &phijj);
  if (_isMC){
    //_ntuple->SetBranchAddress("phijj_jerDown", &phijj_jerDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("phijj_jerUp", &phijj_jerUp);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("phijj_jesDown", &phijj_jesDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("phijj_jesUp", &phijj_jesUp);
  }
  //_ntuple->SetBranchAddress("ptjj", &ptjj);
  if (_isMC){
    //_ntuple->SetBranchAddress("ptjj_jerDown", &ptjj_jerDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("ptjj_jerUp", &ptjj_jerUp);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("ptjj_jesDown", &ptjj_jesDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("ptjj_jesUp", &ptjj_jesUp);
  }
  //_ntuple->SetBranchAddress("pvRho", &pvRho);
  //_ntuple->SetBranchAddress("pvZ", &pvZ);
  //_ntuple->SetBranchAddress("pvndof", &pvndof);
  //_ntuple->SetBranchAddress("rapidity", &rapidity);
  //_ntuple->SetBranchAddress("type1_pfMETEt", &type1_pfMETEt);
  //_ntuple->SetBranchAddress("type1_pfMETPhi", &type1_pfMETPhi);
  //_ntuple->SetBranchAddress("zeppenfeld", &zeppenfeld);
  if (_isMC){
    //_ntuple->SetBranchAddress("zeppenfeld_jerDown", &zeppenfeld_jerDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("zeppenfeld_jerUp", &zeppenfeld_jerUp);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("zeppenfeld_jesDown", &zeppenfeld_jesDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("zeppenfeld_jesUp", &zeppenfeld_jesUp);
  }
  //_ntuple->SetBranchAddress("zeppenfeldj3", &zeppenfeldj3);
  if (_isMC){
    //_ntuple->SetBranchAddress("zeppenfeldj3_jerDown", &zeppenfeldj3_jerDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("zeppenfeldj3_jerUp", &zeppenfeldj3_jerUp);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("zeppenfeldj3_jesDown", &zeppenfeldj3_jesDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("zeppenfeldj3_jesUp", &zeppenfeldj3_jesUp);
  }
  //_ntuple->SetBranchAddress("pvIdFake", &pvIdFake);
  //_ntuple->SetBranchAddress("pvIsValid", &pvIsValid);
  //_ntuple->SetBranchAddress("Charge", &Charge);
  //_ntuple->SetBranchAddress("PdgId", &PdgId);
  _ntuple->SetBranchAddress("lumi", &lumi);
  //_ntuple->SetBranchAddress("nJets", &nJets);
  if (_isMC){
    //_ntuple->SetBranchAddress("nJets_jerDown", &nJets_jerDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("nJets_jerUp", &nJets_jerUp);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("nJets_jesDown", &nJets_jesDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("nJets_jesUp", &nJets_jesUp);
  }
  //_ntuple->SetBranchAddress("nZZLooseElec", &nZZLooseElec);
  //_ntuple->SetBranchAddress("nZZLooseMu", &nZZLooseMu);
  //_ntuple->SetBranchAddress("nZZTightElec", &nZZTightElec);
  //_ntuple->SetBranchAddress("nZZTightIsoElec", &nZZTightIsoElec);
  //_ntuple->SetBranchAddress("nZZTightIsoMu", &nZZTightIsoMu);
  //_ntuple->SetBranchAddress("nZZTightMu", &nZZTightMu);
  _ntuple->SetBranchAddress("nvtx", &nvtx);
  _ntuple->SetBranchAddress("run", &run);
  _ntuple->SetBranchAddress("evt", &evt);
  //_ntuple->SetBranchAddress("jetDeepCSV", &jetDeepCSV);
  //_ntuple->SetBranchAddress("jetEta", &jetEta);
  if (_isMC){
    //_ntuple->SetBranchAddress("jetEta_jerDown", &jetEta_jerDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("jetEta_jerUp", &jetEta_jerUp);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("jetEta_jesDown", &jetEta_jesDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("jetEta_jesUp", &jetEta_jesUp);
  }
  //_ntuple->SetBranchAddress("jetPhi", &jetPhi);
  //_ntuple->SetBranchAddress("jetPt", &jetPt);
  if (_isMC){
    //_ntuple->SetBranchAddress("jetPt_jerDown", &jetPt_jerDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("jetPt_jerUp", &jetPt_jerUp);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("jetPt_jesDown", &jetPt_jesDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("jetPt_jesUp", &jetPt_jesUp);
  }
  //_ntuple->SetBranchAddress("jetRapidity", &jetRapidity);
  if (_isMC){
    //_ntuple->SetBranchAddress("jetRapidity_jerDown", &jetRapidity_jerDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("jetRapidity_jerUp", &jetRapidity_jerUp);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("jetRapidity_jesDown", &jetRapidity_jesDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("jetRapidity_jesUp", &jetRapidity_jesUp);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("pdfWeights", &pdfWeights);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("scaleWeightIDs", &scaleWeightIDs);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("scaleWeights", &scaleWeights);
  }
  //_ntuple->SetBranchAddress("isGenJetMatched", &isGenJetMatched);
  //_ntuple->SetBranchAddress("jetHadronFlavor", &jetHadronFlavor);
  //_ntuple->SetBranchAddress("jetPUID", &jetPUID);
  if (_isMC){
    //_ntuple->SetBranchAddress("jetPUID_jerDown", &jetPUID_jerDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("jetPUID_jerUp", &jetPUID_jerUp);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("jetPUID_jesDown", &jetPUID_jesDown);
  }
  if (_isMC){
    //_ntuple->SetBranchAddress("jetPUID_jesUp", &jetPUID_jesUp);
  }
  //_ntuple->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVLPass", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVLPass);
  //_ntuple->SetBranchAddress("HLT_DoubleEle25_CaloIdL_MWPass", &HLT_DoubleEle25_CaloIdL_MWPass);
  //_ntuple->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8Pass", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8Pass);
  //_ntuple->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVLPass", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVLPass);
  //_ntuple->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass);
  //_ntuple->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass);
  //_ntuple->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass);
  //_ntuple->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZPass", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZPass);
  //_ntuple->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZPass", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZPass);
  //_ntuple->SetBranchAddress("HLT_TripleMu_10_5_5_DZPass", &HLT_TripleMu_10_5_5_DZPass);
  //_ntuple->SetBranchAddress("HLT_TripleMu_12_10_5Pass", &HLT_TripleMu_12_10_5Pass);
  //_ntuple->SetBranchAddress("HLT_Ele30_WPTight_GsfPass", &HLT_Ele30_WPTight_GsfPass);
  //_ntuple->SetBranchAddress("HLT_Ele32_WPTight_GsfPass", &HLT_Ele32_WPTight_GsfPass);
  //_ntuple->SetBranchAddress("HLT_IsoMu24Pass", &HLT_IsoMu24Pass);
  //_ntuple->SetBranchAddress("Flag_goodVerticesPass", &Flag_goodVerticesPass);
  //_ntuple->SetBranchAddress("Flag_globalSuperTightHalo2016FilterPass", &Flag_globalSuperTightHalo2016FilterPass);
  //_ntuple->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilterPass", &Flag_EcalDeadCellTriggerPrimitiveFilterPass);
  //_ntuple->SetBranchAddress("Flag_BadPFMuonFilterPass", &Flag_BadPFMuonFilterPass);
  //_ntuple->SetBranchAddress("Flag_BadPFMuonDzFilterPass", &Flag_BadPFMuonDzFilterPass);
  //_ntuple->SetBranchAddress("Flag_hfNoisyHitsFilterPass", &Flag_hfNoisyHitsFilterPass);
  //_ntuple->SetBranchAddress("Flag_eeBadScFilterPass", &Flag_eeBadScFilterPass);

  if (_channel == "eeee"){
    //_ntuple->SetBranchAddress("e1_e3_DR", &e1_e3_DR);
    //_ntuple->SetBranchAddress("e1_e3_Mass", &e1_e3_Mass);
    //_ntuple->SetBranchAddress("e1_e3_MassNoFSR", &e1_e3_MassNoFSR);
    //_ntuple->SetBranchAddress("e1_e4_DR", &e1_e4_DR);
    //_ntuple->SetBranchAddress("e1_e4_Mass", &e1_e4_Mass);
    //_ntuple->SetBranchAddress("e1_e4_MassNoFSR", &e1_e4_MassNoFSR);
    //_ntuple->SetBranchAddress("e2_e3_DR", &e2_e3_DR);
    //_ntuple->SetBranchAddress("e2_e3_Mass", &e2_e3_Mass);
    //_ntuple->SetBranchAddress("e2_e3_MassNoFSR", &e2_e3_MassNoFSR);
    //_ntuple->SetBranchAddress("e2_e4_DR", &e2_e4_DR);
    //_ntuple->SetBranchAddress("e2_e4_Mass", &e2_e4_Mass);
    //_ntuple->SetBranchAddress("e2_e4_MassNoFSR", &e2_e4_MassNoFSR);
    //_ntuple->SetBranchAddress("e1_e3_SS", &e1_e3_SS);
    //_ntuple->SetBranchAddress("e1_e4_SS", &e1_e4_SS);
    //_ntuple->SetBranchAddress("e2_e3_SS", &e2_e3_SS);
    //_ntuple->SetBranchAddress("e2_e4_SS", &e2_e4_SS);
    //_ntuple->SetBranchAddress("e3_e4_DR", &e3_e4_DR);
    //_ntuple->SetBranchAddress("e3_e4_Energy", &e3_e4_Energy);
    //_ntuple->SetBranchAddress("e3_e4_EnergyNoFSR", &e3_e4_EnergyNoFSR);
    //_ntuple->SetBranchAddress("e3_e4_Eta", &e3_e4_Eta);
    //_ntuple->SetBranchAddress("e3_e4_EtaNoFSR", &e3_e4_EtaNoFSR);
    _ntuple->SetBranchAddress("e3_e4_Mass", &e3_e4_Mass);
    //_ntuple->SetBranchAddress("e3_e4_MassNoFSR", &e3_e4_MassNoFSR);
    //_ntuple->SetBranchAddress("e3_e4_Mt", &e3_e4_Mt);
    //_ntuple->SetBranchAddress("e3_e4_MtToMET", &e3_e4_MtToMET);
    //_ntuple->SetBranchAddress("e3_e4_Phi", &e3_e4_Phi);
    //_ntuple->SetBranchAddress("e3_e4_PhiNoFSR", &e3_e4_PhiNoFSR);
    //_ntuple->SetBranchAddress("e3_e4_Pt", &e3_e4_Pt);
    //_ntuple->SetBranchAddress("e3_e4_PtNoFSR", &e3_e4_PtNoFSR);
    //_ntuple->SetBranchAddress("e3_e4_SS", &e3_e4_SS);
    //_ntuple->SetBranchAddress("e3_e4_Charge", &e3_e4_Charge);
    //_ntuple->SetBranchAddress("e3_e4_PdgId", &e3_e4_PdgId);
    //_ntuple->SetBranchAddress("e3EffectiveArea", &e3EffectiveArea);
    _ntuple->SetBranchAddress("e3Energy", &e3Energy);
    _ntuple->SetBranchAddress("e3Eta", &e3Eta);
    //_ntuple->SetBranchAddress("e3FSRDREt", &e3FSRDREt);
    //_ntuple->SetBranchAddress("e3FSREt", &e3FSREt);
    //_ntuple->SetBranchAddress("e3FSREta", &e3FSREta);
    //_ntuple->SetBranchAddress("e3FSRPhi", &e3FSRPhi);
    if (_isMC){
      //_ntuple->SetBranchAddress("e3GenEta", &e3GenEta);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("e3GenPhi", &e3GenPhi);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("e3GenPt", &e3GenPt);
    }
    //_ntuple->SetBranchAddress("e3IP2D", &e3IP2D);
    //_ntuple->SetBranchAddress("e3IP2DUncertainty", &e3IP2DUncertainty);
    //_ntuple->SetBranchAddress("e3IP3D", &e3IP3D);
    //_ntuple->SetBranchAddress("e3IP3DUncertainty", &e3IP3DUncertainty);
    //_ntuple->SetBranchAddress("e3Loose", &e3Loose);
    //_ntuple->SetBranchAddress("e3Mass", &e3Mass);
    //_ntuple->SetBranchAddress("e3Medium", &e3Medium);
    //_ntuple->SetBranchAddress("e3MtToMET", &e3MtToMET);
    //_ntuple->SetBranchAddress("e3PFChargedIso", &e3PFChargedIso);
    //_ntuple->SetBranchAddress("e3PFNeutralIso", &e3PFNeutralIso);
    //_ntuple->SetBranchAddress("e3PFPUIso", &e3PFPUIso);
    //_ntuple->SetBranchAddress("e3PFPhotonIso", &e3PFPhotonIso);
    //_ntuple->SetBranchAddress("e3PVDXY", &e3PVDXY);
    //_ntuple->SetBranchAddress("e3PVDZ", &e3PVDZ);
    _ntuple->SetBranchAddress("e3Phi", &e3Phi);
    _ntuple->SetBranchAddress("e3Pt", &e3Pt);
    //_ntuple->SetBranchAddress("e3RelPFIsoRho", &e3RelPFIsoRho);
    //_ntuple->SetBranchAddress("e3Rho", &e3Rho);
    //_ntuple->SetBranchAddress("e3SCEnergy", &e3SCEnergy);
    //_ntuple->SetBranchAddress("e3SCEta", &e3SCEta);
    //_ntuple->SetBranchAddress("e3SCPhi", &e3SCPhi);
    //_ntuple->SetBranchAddress("e3SCRawEnergy", &e3SCRawEnergy);
    //_ntuple->SetBranchAddress("e3SIP2D", &e3SIP2D);
    //_ntuple->SetBranchAddress("e3SIP3D", &e3SIP3D);
    //_ntuple->SetBranchAddress("e3ScaleGainDn", &e3ScaleGainDn);
    //_ntuple->SetBranchAddress("e3ScaleGainUp", &e3ScaleGainUp);
    //_ntuple->SetBranchAddress("e3ScaleStatDn", &e3ScaleStatDn);
    //_ntuple->SetBranchAddress("e3ScaleStatUp", &e3ScaleStatUp);
    //_ntuple->SetBranchAddress("e3ScaleSystDn", &e3ScaleSystDn);
    //_ntuple->SetBranchAddress("e3ScaleSystUp", &e3ScaleSystUp);
    //_ntuple->SetBranchAddress("e3ScaleTotDn", &e3ScaleTotDn);
    //_ntuple->SetBranchAddress("e3ScaleTotUp", &e3ScaleTotUp);
    //_ntuple->SetBranchAddress("e3SigmaPhiDn", &e3SigmaPhiDn);
    //_ntuple->SetBranchAddress("e3SigmaPhiUp", &e3SigmaPhiUp);
    //_ntuple->SetBranchAddress("e3SigmaRhoDn", &e3SigmaRhoDn);
    //_ntuple->SetBranchAddress("e3SigmaRhoUp", &e3SigmaRhoUp);
    //_ntuple->SetBranchAddress("e3SigmaTotDn", &e3SigmaTotDn);
    //_ntuple->SetBranchAddress("e3SigmaTotUp", &e3SigmaTotUp);
    //_ntuple->SetBranchAddress("e3Tight", &e3Tight);
    //_ntuple->SetBranchAddress("e3UnCorrPt", &e3UnCorrPt);
    //_ntuple->SetBranchAddress("e3Veto", &e3Veto);
    //_ntuple->SetBranchAddress("e3ZZIso", &e3ZZIso);
    //_ntuple->SetBranchAddress("e3IsEB", &e3IsEB);
    //_ntuple->SetBranchAddress("e3IsGap", &e3IsGap);
    //_ntuple->SetBranchAddress("e3ZZIsoPass", &e3ZZIsoPass);
    //_ntuple->SetBranchAddress("e3ZZLooseID", &e3ZZLooseID);
    //_ntuple->SetBranchAddress("e3ZZLooseIDNoVtx", &e3ZZLooseIDNoVtx);
    //_ntuple->SetBranchAddress("e3ZZTightID", &e3ZZTightID);
    //_ntuple->SetBranchAddress("e3ZZTightIDNoVtx", &e3ZZTightIDNoVtx);
    //_ntuple->SetBranchAddress("e3Charge", &e3Charge);
    if (_isMC){
      //_ntuple->SetBranchAddress("e3GenCharge", &e3GenCharge);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("e3GenPdgId", &e3GenPdgId);
    }
    _ntuple->SetBranchAddress("e3PdgId", &e3PdgId);
    //_ntuple->SetBranchAddress("e3MissingHits", &e3MissingHits);
    //_ntuple->SetBranchAddress("e4EffectiveArea", &e4EffectiveArea);
    _ntuple->SetBranchAddress("e4Energy", &e4Energy);
    _ntuple->SetBranchAddress("e4Eta", &e4Eta);
    //_ntuple->SetBranchAddress("e4FSRDREt", &e4FSRDREt);
    //_ntuple->SetBranchAddress("e4FSREt", &e4FSREt);
    //_ntuple->SetBranchAddress("e4FSREta", &e4FSREta);
    //_ntuple->SetBranchAddress("e4FSRPhi", &e4FSRPhi);
    if (_isMC){
      //_ntuple->SetBranchAddress("e4GenEta", &e4GenEta);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("e4GenPhi", &e4GenPhi);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("e4GenPt", &e4GenPt);
    }
    //_ntuple->SetBranchAddress("e4IP2D", &e4IP2D);
    //_ntuple->SetBranchAddress("e4IP2DUncertainty", &e4IP2DUncertainty);
    //_ntuple->SetBranchAddress("e4IP3D", &e4IP3D);
    //_ntuple->SetBranchAddress("e4IP3DUncertainty", &e4IP3DUncertainty);
    //_ntuple->SetBranchAddress("e4Loose", &e4Loose);
    //_ntuple->SetBranchAddress("e4Mass", &e4Mass);
    //_ntuple->SetBranchAddress("e4Medium", &e4Medium);
    //_ntuple->SetBranchAddress("e4MtToMET", &e4MtToMET);
    //_ntuple->SetBranchAddress("e4PFChargedIso", &e4PFChargedIso);
    //_ntuple->SetBranchAddress("e4PFNeutralIso", &e4PFNeutralIso);
    //_ntuple->SetBranchAddress("e4PFPUIso", &e4PFPUIso);
    //_ntuple->SetBranchAddress("e4PFPhotonIso", &e4PFPhotonIso);
    //_ntuple->SetBranchAddress("e4PVDXY", &e4PVDXY);
    //_ntuple->SetBranchAddress("e4PVDZ", &e4PVDZ);
    _ntuple->SetBranchAddress("e4Phi", &e4Phi);
    _ntuple->SetBranchAddress("e4Pt", &e4Pt);
    //_ntuple->SetBranchAddress("e4RelPFIsoRho", &e4RelPFIsoRho);
    //_ntuple->SetBranchAddress("e4Rho", &e4Rho);
    //_ntuple->SetBranchAddress("e4SCEnergy", &e4SCEnergy);
    //_ntuple->SetBranchAddress("e4SCEta", &e4SCEta);
    //_ntuple->SetBranchAddress("e4SCPhi", &e4SCPhi);
    //_ntuple->SetBranchAddress("e4SCRawEnergy", &e4SCRawEnergy);
    //_ntuple->SetBranchAddress("e4SIP2D", &e4SIP2D);
    //_ntuple->SetBranchAddress("e4SIP3D", &e4SIP3D);
    //_ntuple->SetBranchAddress("e4ScaleGainDn", &e4ScaleGainDn);
    //_ntuple->SetBranchAddress("e4ScaleGainUp", &e4ScaleGainUp);
    //_ntuple->SetBranchAddress("e4ScaleStatDn", &e4ScaleStatDn);
    //_ntuple->SetBranchAddress("e4ScaleStatUp", &e4ScaleStatUp);
    //_ntuple->SetBranchAddress("e4ScaleSystDn", &e4ScaleSystDn);
    //_ntuple->SetBranchAddress("e4ScaleSystUp", &e4ScaleSystUp);
    //_ntuple->SetBranchAddress("e4ScaleTotDn", &e4ScaleTotDn);
    //_ntuple->SetBranchAddress("e4ScaleTotUp", &e4ScaleTotUp);
    //_ntuple->SetBranchAddress("e4SigmaPhiDn", &e4SigmaPhiDn);
    //_ntuple->SetBranchAddress("e4SigmaPhiUp", &e4SigmaPhiUp);
    //_ntuple->SetBranchAddress("e4SigmaRhoDn", &e4SigmaRhoDn);
    //_ntuple->SetBranchAddress("e4SigmaRhoUp", &e4SigmaRhoUp);
    //_ntuple->SetBranchAddress("e4SigmaTotDn", &e4SigmaTotDn);
    //_ntuple->SetBranchAddress("e4SigmaTotUp", &e4SigmaTotUp);
    //_ntuple->SetBranchAddress("e4Tight", &e4Tight);
    //_ntuple->SetBranchAddress("e4UnCorrPt", &e4UnCorrPt);
    //_ntuple->SetBranchAddress("e4Veto", &e4Veto);
    //_ntuple->SetBranchAddress("e4ZZIso", &e4ZZIso);
    //_ntuple->SetBranchAddress("e4IsEB", &e4IsEB);
    //_ntuple->SetBranchAddress("e4IsGap", &e4IsGap);
    //_ntuple->SetBranchAddress("e4ZZIsoPass", &e4ZZIsoPass);
    //_ntuple->SetBranchAddress("e4ZZLooseID", &e4ZZLooseID);
    //_ntuple->SetBranchAddress("e4ZZLooseIDNoVtx", &e4ZZLooseIDNoVtx);
    //_ntuple->SetBranchAddress("e4ZZTightID", &e4ZZTightID);
    //_ntuple->SetBranchAddress("e4ZZTightIDNoVtx", &e4ZZTightIDNoVtx);
    //_ntuple->SetBranchAddress("e4Charge", &e4Charge);
    if (_isMC){
      //_ntuple->SetBranchAddress("e4GenCharge", &e4GenCharge);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("e4GenPdgId", &e4GenPdgId);
    }
    _ntuple->SetBranchAddress("e4PdgId", &e4PdgId);
    //_ntuple->SetBranchAddress("e4MissingHits", &e4MissingHits);
  }

  else if (_channel == "eemm"){
    //_ntuple->SetBranchAddress("e1_m1_DR", &e1_m1_DR);
    //_ntuple->SetBranchAddress("e1_m1_Mass", &e1_m1_Mass);
    //_ntuple->SetBranchAddress("e1_m1_MassNoFSR", &e1_m1_MassNoFSR);
    //_ntuple->SetBranchAddress("e1_m2_DR", &e1_m2_DR);
    //_ntuple->SetBranchAddress("e1_m2_Mass", &e1_m2_Mass);
    //_ntuple->SetBranchAddress("e1_m2_MassNoFSR", &e1_m2_MassNoFSR);
    //_ntuple->SetBranchAddress("e2_m1_DR", &e2_m1_DR);
    //_ntuple->SetBranchAddress("e2_m1_Mass", &e2_m1_Mass);
    //_ntuple->SetBranchAddress("e2_m1_MassNoFSR", &e2_m1_MassNoFSR);
    //_ntuple->SetBranchAddress("e2_m2_DR", &e2_m2_DR);
    //_ntuple->SetBranchAddress("e2_m2_Mass", &e2_m2_Mass);
    //_ntuple->SetBranchAddress("e2_m2_MassNoFSR", &e2_m2_MassNoFSR);
    //_ntuple->SetBranchAddress("e1_m1_SS", &e1_m1_SS);
    //_ntuple->SetBranchAddress("e1_m2_SS", &e1_m2_SS);
    //_ntuple->SetBranchAddress("e2_m1_SS", &e2_m1_SS);
    //_ntuple->SetBranchAddress("e2_m2_SS", &e2_m2_SS);
  }

  else if (_channel == "mmmm"){
    //_ntuple->SetBranchAddress("m1_m3_DR", &m1_m3_DR);
    //_ntuple->SetBranchAddress("m1_m3_Mass", &m1_m3_Mass);
    //_ntuple->SetBranchAddress("m1_m3_MassNoFSR", &m1_m3_MassNoFSR);
    //_ntuple->SetBranchAddress("m1_m4_DR", &m1_m4_DR);
    //_ntuple->SetBranchAddress("m1_m4_Mass", &m1_m4_Mass);
    //_ntuple->SetBranchAddress("m1_m4_MassNoFSR", &m1_m4_MassNoFSR);
    //_ntuple->SetBranchAddress("m2_m3_DR", &m2_m3_DR);
    //_ntuple->SetBranchAddress("m2_m3_Mass", &m2_m3_Mass);
    //_ntuple->SetBranchAddress("m2_m3_MassNoFSR", &m2_m3_MassNoFSR);
    //_ntuple->SetBranchAddress("m2_m4_DR", &m2_m4_DR);
    //_ntuple->SetBranchAddress("m2_m4_Mass", &m2_m4_Mass);
    //_ntuple->SetBranchAddress("m2_m4_MassNoFSR", &m2_m4_MassNoFSR);
    //_ntuple->SetBranchAddress("m1_m3_SS", &m1_m3_SS);
    //_ntuple->SetBranchAddress("m1_m4_SS", &m1_m4_SS);
    //_ntuple->SetBranchAddress("m2_m3_SS", &m2_m3_SS);
    //_ntuple->SetBranchAddress("m2_m4_SS", &m2_m4_SS);
    //_ntuple->SetBranchAddress("m3_m4_DR", &m3_m4_DR);
    //_ntuple->SetBranchAddress("m3_m4_Energy", &m3_m4_Energy);
    //_ntuple->SetBranchAddress("m3_m4_EnergyNoFSR", &m3_m4_EnergyNoFSR);
    //_ntuple->SetBranchAddress("m3_m4_Eta", &m3_m4_Eta);
    //_ntuple->SetBranchAddress("m3_m4_EtaNoFSR", &m3_m4_EtaNoFSR);
    _ntuple->SetBranchAddress("m3_m4_Mass", &m3_m4_Mass);
    //_ntuple->SetBranchAddress("m3_m4_MassNoFSR", &m3_m4_MassNoFSR);
    //_ntuple->SetBranchAddress("m3_m4_Mt", &m3_m4_Mt);
    //_ntuple->SetBranchAddress("m3_m4_MtToMET", &m3_m4_MtToMET);
    //_ntuple->SetBranchAddress("m3_m4_Phi", &m3_m4_Phi);
    //_ntuple->SetBranchAddress("m3_m4_PhiNoFSR", &m3_m4_PhiNoFSR);
    //_ntuple->SetBranchAddress("m3_m4_Pt", &m3_m4_Pt);
    //_ntuple->SetBranchAddress("m3_m4_PtNoFSR", &m3_m4_PtNoFSR);
    //_ntuple->SetBranchAddress("m3_m4_SS", &m3_m4_SS);
    //_ntuple->SetBranchAddress("m3_m4_Charge", &m3_m4_Charge);
    //_ntuple->SetBranchAddress("m3_m4_PdgId", &m3_m4_PdgId);
    //_ntuple->SetBranchAddress("m3EffScaleFactor", &m3EffScaleFactor);
    //_ntuple->SetBranchAddress("m3EffScaleFactorError", &m3EffScaleFactorError);
    _ntuple->SetBranchAddress("m3Energy", &m3Energy);
    _ntuple->SetBranchAddress("m3Eta", &m3Eta);
    //_ntuple->SetBranchAddress("m3FSRDREt", &m3FSRDREt);
    //_ntuple->SetBranchAddress("m3FSREt", &m3FSREt);
    //_ntuple->SetBranchAddress("m3FSREta", &m3FSREta);
    //_ntuple->SetBranchAddress("m3FSRPhi", &m3FSRPhi);
    if (_isMC){
      //_ntuple->SetBranchAddress("m3GenEta", &m3GenEta);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("m3GenPhi", &m3GenPhi);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("m3GenPt", &m3GenPt);
    }
    //_ntuple->SetBranchAddress("m3IP2D", &m3IP2D);
    //_ntuple->SetBranchAddress("m3IP2DUncertainty", &m3IP2DUncertainty);
    //_ntuple->SetBranchAddress("m3IP3D", &m3IP3D);
    //_ntuple->SetBranchAddress("m3IP3DUncertainty", &m3IP3DUncertainty);
    //_ntuple->SetBranchAddress("m3Mass", &m3Mass);
    //_ntuple->SetBranchAddress("m3MtToMET", &m3MtToMET);
    //_ntuple->SetBranchAddress("m3PFChargedIso", &m3PFChargedIso);
    //_ntuple->SetBranchAddress("m3PFNeutralIso", &m3PFNeutralIso);
    //_ntuple->SetBranchAddress("m3PFPUIso", &m3PFPUIso);
    //_ntuple->SetBranchAddress("m3PFPhotonIso", &m3PFPhotonIso);
    //_ntuple->SetBranchAddress("m3PVDXY", &m3PVDXY);
    //_ntuple->SetBranchAddress("m3PVDZ", &m3PVDZ);
    _ntuple->SetBranchAddress("m3Phi", &m3Phi);
    _ntuple->SetBranchAddress("m3Pt", &m3Pt);
    //_ntuple->SetBranchAddress("m3PtErr", &m3PtErr);
    //_ntuple->SetBranchAddress("m3PtErrUncorrected", &m3PtErrUncorrected);
    //_ntuple->SetBranchAddress("m3PtUncorrected", &m3PtUncorrected);
    //_ntuple->SetBranchAddress("m3RelPFIsoDBR03", &m3RelPFIsoDBR03);
    //_ntuple->SetBranchAddress("m3RelPFIsoDBR04", &m3RelPFIsoDBR04);
    //_ntuple->SetBranchAddress("m3RochesterScaleUncDn", &m3RochesterScaleUncDn);
    //_ntuple->SetBranchAddress("m3RochesterScaleUncUp", &m3RochesterScaleUncUp);
    //_ntuple->SetBranchAddress("m3RochesterSmearUncDn", &m3RochesterSmearUncDn);
    //_ntuple->SetBranchAddress("m3RochesterSmearUncUp", &m3RochesterSmearUncUp);
    //_ntuple->SetBranchAddress("m3SIP2D", &m3SIP2D);
    //_ntuple->SetBranchAddress("m3SIP3D", &m3SIP3D);
    //_ntuple->SetBranchAddress("m3TrackIso", &m3TrackIso);
    //_ntuple->SetBranchAddress("m3ZZIso", &m3ZZIso);
    //_ntuple->SetBranchAddress("m3HighPtID", &m3HighPtID);
    //_ntuple->SetBranchAddress("m3HighPtIDNoVtx", &m3HighPtIDNoVtx);
    //_ntuple->SetBranchAddress("m3IsGlobal", &m3IsGlobal);
    //_ntuple->SetBranchAddress("m3IsLoose", &m3IsLoose);
    //_ntuple->SetBranchAddress("m3IsMedium", &m3IsMedium);
    //_ntuple->SetBranchAddress("m3IsPFMuon", &m3IsPFMuon);
    //_ntuple->SetBranchAddress("m3IsTracker", &m3IsTracker);
    //_ntuple->SetBranchAddress("m3PASTightID", &m3PASTightID);
    //_ntuple->SetBranchAddress("m3PASTightIDNoVtx", &m3PASTightIDNoVtx);
    //_ntuple->SetBranchAddress("m3PFID", &m3PFID);
    //_ntuple->SetBranchAddress("m3PFIDNoVtx", &m3PFIDNoVtx);
    //_ntuple->SetBranchAddress("m3ZZIsoPass", &m3ZZIsoPass);
    //_ntuple->SetBranchAddress("m3ZZLooseID", &m3ZZLooseID);
    //_ntuple->SetBranchAddress("m3ZZLooseIDNoVtx", &m3ZZLooseIDNoVtx);
    //_ntuple->SetBranchAddress("m3ZZTightID", &m3ZZTightID);
    //_ntuple->SetBranchAddress("m3ZZTightIDNoVtx", &m3ZZTightIDNoVtx);
    //_ntuple->SetBranchAddress("m3Charge", &m3Charge);
    if (_isMC){
      //_ntuple->SetBranchAddress("m3GenCharge", &m3GenCharge);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("m3GenPdgId", &m3GenPdgId);
    }
    _ntuple->SetBranchAddress("m3PdgId", &m3PdgId);
    //_ntuple->SetBranchAddress("m3BestTrackType", &m3BestTrackType);
    //_ntuple->SetBranchAddress("m3CutBasedLoose", &m3CutBasedLoose);
    //_ntuple->SetBranchAddress("m3CutBasedMedium", &m3CutBasedMedium);
    //_ntuple->SetBranchAddress("m3CutBasedTight", &m3CutBasedTight);
    //_ntuple->SetBranchAddress("m3IsTight", &m3IsTight);
    //_ntuple->SetBranchAddress("m3MatchedStations", &m3MatchedStations);
    //_ntuple->SetBranchAddress("m3NoOfMatches", &m3NoOfMatches);
    //_ntuple->SetBranchAddress("m3PFIsoLoose", &m3PFIsoLoose);
    //_ntuple->SetBranchAddress("m3PFIsoMedium", &m3PFIsoMedium);
    //_ntuple->SetBranchAddress("m3PFIsoTight", &m3PFIsoTight);
    //_ntuple->SetBranchAddress("m3PFIsoVTight", &m3PFIsoVTight);
    //_ntuple->SetBranchAddress("m4EffScaleFactor", &m4EffScaleFactor);
    //_ntuple->SetBranchAddress("m4EffScaleFactorError", &m4EffScaleFactorError);
    _ntuple->SetBranchAddress("m4Energy", &m4Energy);
    _ntuple->SetBranchAddress("m4Eta", &m4Eta);
    //_ntuple->SetBranchAddress("m4FSRDREt", &m4FSRDREt);
    //_ntuple->SetBranchAddress("m4FSREt", &m4FSREt);
    //_ntuple->SetBranchAddress("m4FSREta", &m4FSREta);
    //_ntuple->SetBranchAddress("m4FSRPhi", &m4FSRPhi);
    if (_isMC){
      //_ntuple->SetBranchAddress("m4GenEta", &m4GenEta);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("m4GenPhi", &m4GenPhi);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("m4GenPt", &m4GenPt);
    }
    //_ntuple->SetBranchAddress("m4IP2D", &m4IP2D);
    //_ntuple->SetBranchAddress("m4IP2DUncertainty", &m4IP2DUncertainty);
    //_ntuple->SetBranchAddress("m4IP3D", &m4IP3D);
    //_ntuple->SetBranchAddress("m4IP3DUncertainty", &m4IP3DUncertainty);
    //_ntuple->SetBranchAddress("m4Mass", &m4Mass);
    //_ntuple->SetBranchAddress("m4MtToMET", &m4MtToMET);
    //_ntuple->SetBranchAddress("m4PFChargedIso", &m4PFChargedIso);
    //_ntuple->SetBranchAddress("m4PFNeutralIso", &m4PFNeutralIso);
    //_ntuple->SetBranchAddress("m4PFPUIso", &m4PFPUIso);
    //_ntuple->SetBranchAddress("m4PFPhotonIso", &m4PFPhotonIso);
    //_ntuple->SetBranchAddress("m4PVDXY", &m4PVDXY);
    //_ntuple->SetBranchAddress("m4PVDZ", &m4PVDZ);
    _ntuple->SetBranchAddress("m4Phi", &m4Phi);
    _ntuple->SetBranchAddress("m4Pt", &m4Pt);
    //_ntuple->SetBranchAddress("m4PtErr", &m4PtErr);
    //_ntuple->SetBranchAddress("m4PtErrUncorrected", &m4PtErrUncorrected);
    //_ntuple->SetBranchAddress("m4PtUncorrected", &m4PtUncorrected);
    //_ntuple->SetBranchAddress("m4RelPFIsoDBR03", &m4RelPFIsoDBR03);
    //_ntuple->SetBranchAddress("m4RelPFIsoDBR04", &m4RelPFIsoDBR04);
    //_ntuple->SetBranchAddress("m4RochesterScaleUncDn", &m4RochesterScaleUncDn);
    //_ntuple->SetBranchAddress("m4RochesterScaleUncUp", &m4RochesterScaleUncUp);
    //_ntuple->SetBranchAddress("m4RochesterSmearUncDn", &m4RochesterSmearUncDn);
    //_ntuple->SetBranchAddress("m4RochesterSmearUncUp", &m4RochesterSmearUncUp);
    //_ntuple->SetBranchAddress("m4SIP2D", &m4SIP2D);
    //_ntuple->SetBranchAddress("m4SIP3D", &m4SIP3D);
    //_ntuple->SetBranchAddress("m4TrackIso", &m4TrackIso);
    //_ntuple->SetBranchAddress("m4ZZIso", &m4ZZIso);
    //_ntuple->SetBranchAddress("m4HighPtID", &m4HighPtID);
    //_ntuple->SetBranchAddress("m4HighPtIDNoVtx", &m4HighPtIDNoVtx);
    //_ntuple->SetBranchAddress("m4IsGlobal", &m4IsGlobal);
    //_ntuple->SetBranchAddress("m4IsLoose", &m4IsLoose);
    //_ntuple->SetBranchAddress("m4IsMedium", &m4IsMedium);
    //_ntuple->SetBranchAddress("m4IsPFMuon", &m4IsPFMuon);
    //_ntuple->SetBranchAddress("m4IsTracker", &m4IsTracker);
    //_ntuple->SetBranchAddress("m4PASTightID", &m4PASTightID);
    //_ntuple->SetBranchAddress("m4PASTightIDNoVtx", &m4PASTightIDNoVtx);
    //_ntuple->SetBranchAddress("m4PFID", &m4PFID);
    //_ntuple->SetBranchAddress("m4PFIDNoVtx", &m4PFIDNoVtx);
    //_ntuple->SetBranchAddress("m4ZZIsoPass", &m4ZZIsoPass);
    //_ntuple->SetBranchAddress("m4ZZLooseID", &m4ZZLooseID);
    //_ntuple->SetBranchAddress("m4ZZLooseIDNoVtx", &m4ZZLooseIDNoVtx);
    //_ntuple->SetBranchAddress("m4ZZTightID", &m4ZZTightID);
    //_ntuple->SetBranchAddress("m4ZZTightIDNoVtx", &m4ZZTightIDNoVtx);
    //_ntuple->SetBranchAddress("m4Charge", &m4Charge);
    if (_isMC){
      //_ntuple->SetBranchAddress("m4GenCharge", &m4GenCharge);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("m4GenPdgId", &m4GenPdgId);
    }
    _ntuple->SetBranchAddress("m4PdgId", &m4PdgId);
    //_ntuple->SetBranchAddress("m4BestTrackType", &m4BestTrackType);
    //_ntuple->SetBranchAddress("m4CutBasedLoose", &m4CutBasedLoose);
    //_ntuple->SetBranchAddress("m4CutBasedMedium", &m4CutBasedMedium);
    //_ntuple->SetBranchAddress("m4CutBasedTight", &m4CutBasedTight);
    //_ntuple->SetBranchAddress("m4IsTight", &m4IsTight);
    //_ntuple->SetBranchAddress("m4MatchedStations", &m4MatchedStations);
    //_ntuple->SetBranchAddress("m4NoOfMatches", &m4NoOfMatches);
    //_ntuple->SetBranchAddress("m4PFIsoLoose", &m4PFIsoLoose);
    //_ntuple->SetBranchAddress("m4PFIsoMedium", &m4PFIsoMedium);
    //_ntuple->SetBranchAddress("m4PFIsoTight", &m4PFIsoTight);
    //_ntuple->SetBranchAddress("m4PFIsoVTight", &m4PFIsoVTight);
  }

  if (_channel == "eeee" || _channel == "eemm"){
    //_ntuple->SetBranchAddress("e1_e2_DR", &e1_e2_DR);
    //_ntuple->SetBranchAddress("e1_e2_Energy", &e1_e2_Energy);
    //_ntuple->SetBranchAddress("e1_e2_EnergyNoFSR", &e1_e2_EnergyNoFSR);
    //_ntuple->SetBranchAddress("e1_e2_Eta", &e1_e2_Eta);
    //_ntuple->SetBranchAddress("e1_e2_EtaNoFSR", &e1_e2_EtaNoFSR);
    _ntuple->SetBranchAddress("e1_e2_Mass", &e1_e2_Mass);
    //_ntuple->SetBranchAddress("e1_e2_MassNoFSR", &e1_e2_MassNoFSR);
    //_ntuple->SetBranchAddress("e1_e2_Mt", &e1_e2_Mt);
    //_ntuple->SetBranchAddress("e1_e2_MtToMET", &e1_e2_MtToMET);
    //_ntuple->SetBranchAddress("e1_e2_Phi", &e1_e2_Phi);
    //_ntuple->SetBranchAddress("e1_e2_PhiNoFSR", &e1_e2_PhiNoFSR);
    //_ntuple->SetBranchAddress("e1_e2_Pt", &e1_e2_Pt);
    //_ntuple->SetBranchAddress("e1_e2_PtNoFSR", &e1_e2_PtNoFSR);
    //_ntuple->SetBranchAddress("e1_e2_SS", &e1_e2_SS);
    //_ntuple->SetBranchAddress("e1_e2_Charge", &e1_e2_Charge);
    //_ntuple->SetBranchAddress("e1_e2_PdgId", &e1_e2_PdgId);
    //_ntuple->SetBranchAddress("e1EffectiveArea", &e1EffectiveArea);
    _ntuple->SetBranchAddress("e1Energy", &e1Energy);
    _ntuple->SetBranchAddress("e1Eta", &e1Eta);
    //_ntuple->SetBranchAddress("e1FSRDREt", &e1FSRDREt);
    //_ntuple->SetBranchAddress("e1FSREt", &e1FSREt);
    //_ntuple->SetBranchAddress("e1FSREta", &e1FSREta);
    //_ntuple->SetBranchAddress("e1FSRPhi", &e1FSRPhi);
    if (_isMC){
      //_ntuple->SetBranchAddress("e1GenEta", &e1GenEta);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("e1GenPhi", &e1GenPhi);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("e1GenPt", &e1GenPt);
    }
    //_ntuple->SetBranchAddress("e1IP2D", &e1IP2D);
    //_ntuple->SetBranchAddress("e1IP2DUncertainty", &e1IP2DUncertainty);
    //_ntuple->SetBranchAddress("e1IP3D", &e1IP3D);
    //_ntuple->SetBranchAddress("e1IP3DUncertainty", &e1IP3DUncertainty);
    //_ntuple->SetBranchAddress("e1Loose", &e1Loose);
    //_ntuple->SetBranchAddress("e1Mass", &e1Mass);
    //_ntuple->SetBranchAddress("e1Medium", &e1Medium);
    //_ntuple->SetBranchAddress("e1MtToMET", &e1MtToMET);
    //_ntuple->SetBranchAddress("e1PFChargedIso", &e1PFChargedIso);
    //_ntuple->SetBranchAddress("e1PFNeutralIso", &e1PFNeutralIso);
    //_ntuple->SetBranchAddress("e1PFPUIso", &e1PFPUIso);
    //_ntuple->SetBranchAddress("e1PFPhotonIso", &e1PFPhotonIso);
    //_ntuple->SetBranchAddress("e1PVDXY", &e1PVDXY);
    //_ntuple->SetBranchAddress("e1PVDZ", &e1PVDZ);
    _ntuple->SetBranchAddress("e1Phi", &e1Phi);
    _ntuple->SetBranchAddress("e1Pt", &e1Pt);
    //_ntuple->SetBranchAddress("e1RelPFIsoRho", &e1RelPFIsoRho);
    //_ntuple->SetBranchAddress("e1Rho", &e1Rho);
    //_ntuple->SetBranchAddress("e1SCEnergy", &e1SCEnergy);
    //_ntuple->SetBranchAddress("e1SCEta", &e1SCEta);
    //_ntuple->SetBranchAddress("e1SCPhi", &e1SCPhi);
    //_ntuple->SetBranchAddress("e1SCRawEnergy", &e1SCRawEnergy);
    //_ntuple->SetBranchAddress("e1SIP2D", &e1SIP2D);
    //_ntuple->SetBranchAddress("e1SIP3D", &e1SIP3D);
    //_ntuple->SetBranchAddress("e1ScaleGainDn", &e1ScaleGainDn);
    //_ntuple->SetBranchAddress("e1ScaleGainUp", &e1ScaleGainUp);
    //_ntuple->SetBranchAddress("e1ScaleStatDn", &e1ScaleStatDn);
    //_ntuple->SetBranchAddress("e1ScaleStatUp", &e1ScaleStatUp);
    //_ntuple->SetBranchAddress("e1ScaleSystDn", &e1ScaleSystDn);
    //_ntuple->SetBranchAddress("e1ScaleSystUp", &e1ScaleSystUp);
    //_ntuple->SetBranchAddress("e1ScaleTotDn", &e1ScaleTotDn);
    //_ntuple->SetBranchAddress("e1ScaleTotUp", &e1ScaleTotUp);
    //_ntuple->SetBranchAddress("e1SigmaPhiDn", &e1SigmaPhiDn);
    //_ntuple->SetBranchAddress("e1SigmaPhiUp", &e1SigmaPhiUp);
    //_ntuple->SetBranchAddress("e1SigmaRhoDn", &e1SigmaRhoDn);
    //_ntuple->SetBranchAddress("e1SigmaRhoUp", &e1SigmaRhoUp);
    //_ntuple->SetBranchAddress("e1SigmaTotDn", &e1SigmaTotDn);
    //_ntuple->SetBranchAddress("e1SigmaTotUp", &e1SigmaTotUp);
    //_ntuple->SetBranchAddress("e1Tight", &e1Tight);
    //_ntuple->SetBranchAddress("e1UnCorrPt", &e1UnCorrPt);
    //_ntuple->SetBranchAddress("e1Veto", &e1Veto);
    //_ntuple->SetBranchAddress("e1ZZIso", &e1ZZIso);
    //_ntuple->SetBranchAddress("e1IsEB", &e1IsEB);
    //_ntuple->SetBranchAddress("e1IsGap", &e1IsGap);
    //_ntuple->SetBranchAddress("e1ZZIsoPass", &e1ZZIsoPass);
    //_ntuple->SetBranchAddress("e1ZZLooseID", &e1ZZLooseID);
    //_ntuple->SetBranchAddress("e1ZZLooseIDNoVtx", &e1ZZLooseIDNoVtx);
    //_ntuple->SetBranchAddress("e1ZZTightID", &e1ZZTightID);
    //_ntuple->SetBranchAddress("e1ZZTightIDNoVtx", &e1ZZTightIDNoVtx);
    //_ntuple->SetBranchAddress("e1Charge", &e1Charge);
    if (_isMC){
      //_ntuple->SetBranchAddress("e1GenCharge", &e1GenCharge);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("e1GenPdgId", &e1GenPdgId);
    }
    _ntuple->SetBranchAddress("e1PdgId", &e1PdgId);
    //_ntuple->SetBranchAddress("e1MissingHits", &e1MissingHits);
    //_ntuple->SetBranchAddress("e2EffectiveArea", &e2EffectiveArea);
    _ntuple->SetBranchAddress("e2Energy", &e2Energy);
    _ntuple->SetBranchAddress("e2Eta", &e2Eta);
    //_ntuple->SetBranchAddress("e2FSRDREt", &e2FSRDREt);
    //_ntuple->SetBranchAddress("e2FSREt", &e2FSREt);
    //_ntuple->SetBranchAddress("e2FSREta", &e2FSREta);
    //_ntuple->SetBranchAddress("e2FSRPhi", &e2FSRPhi);
    if (_isMC){
      //_ntuple->SetBranchAddress("e2GenEta", &e2GenEta);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("e2GenPhi", &e2GenPhi);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("e2GenPt", &e2GenPt);
    }
    //_ntuple->SetBranchAddress("e2IP2D", &e2IP2D);
    //_ntuple->SetBranchAddress("e2IP2DUncertainty", &e2IP2DUncertainty);
    //_ntuple->SetBranchAddress("e2IP3D", &e2IP3D);
    //_ntuple->SetBranchAddress("e2IP3DUncertainty", &e2IP3DUncertainty);
    //_ntuple->SetBranchAddress("e2Loose", &e2Loose);
    //_ntuple->SetBranchAddress("e2Mass", &e2Mass);
    //_ntuple->SetBranchAddress("e2Medium", &e2Medium);
    //_ntuple->SetBranchAddress("e2MtToMET", &e2MtToMET);
    //_ntuple->SetBranchAddress("e2PFChargedIso", &e2PFChargedIso);
    //_ntuple->SetBranchAddress("e2PFNeutralIso", &e2PFNeutralIso);
    //_ntuple->SetBranchAddress("e2PFPUIso", &e2PFPUIso);
    //_ntuple->SetBranchAddress("e2PFPhotonIso", &e2PFPhotonIso);
    //_ntuple->SetBranchAddress("e2PVDXY", &e2PVDXY);
    //_ntuple->SetBranchAddress("e2PVDZ", &e2PVDZ);
    _ntuple->SetBranchAddress("e2Phi", &e2Phi);
    _ntuple->SetBranchAddress("e2Pt", &e2Pt);
    //_ntuple->SetBranchAddress("e2RelPFIsoRho", &e2RelPFIsoRho);
    //_ntuple->SetBranchAddress("e2Rho", &e2Rho);
    //_ntuple->SetBranchAddress("e2SCEnergy", &e2SCEnergy);
    //_ntuple->SetBranchAddress("e2SCEta", &e2SCEta);
    //_ntuple->SetBranchAddress("e2SCPhi", &e2SCPhi);
    //_ntuple->SetBranchAddress("e2SCRawEnergy", &e2SCRawEnergy);
    //_ntuple->SetBranchAddress("e2SIP2D", &e2SIP2D);
    //_ntuple->SetBranchAddress("e2SIP3D", &e2SIP3D);
    //_ntuple->SetBranchAddress("e2ScaleGainDn", &e2ScaleGainDn);
    //_ntuple->SetBranchAddress("e2ScaleGainUp", &e2ScaleGainUp);
    //_ntuple->SetBranchAddress("e2ScaleStatDn", &e2ScaleStatDn);
    //_ntuple->SetBranchAddress("e2ScaleStatUp", &e2ScaleStatUp);
    //_ntuple->SetBranchAddress("e2ScaleSystDn", &e2ScaleSystDn);
    //_ntuple->SetBranchAddress("e2ScaleSystUp", &e2ScaleSystUp);
    //_ntuple->SetBranchAddress("e2ScaleTotDn", &e2ScaleTotDn);
    //_ntuple->SetBranchAddress("e2ScaleTotUp", &e2ScaleTotUp);
    //_ntuple->SetBranchAddress("e2SigmaPhiDn", &e2SigmaPhiDn);
    //_ntuple->SetBranchAddress("e2SigmaPhiUp", &e2SigmaPhiUp);
    //_ntuple->SetBranchAddress("e2SigmaRhoDn", &e2SigmaRhoDn);
    //_ntuple->SetBranchAddress("e2SigmaRhoUp", &e2SigmaRhoUp);
    //_ntuple->SetBranchAddress("e2SigmaTotDn", &e2SigmaTotDn);
    //_ntuple->SetBranchAddress("e2SigmaTotUp", &e2SigmaTotUp);
    //_ntuple->SetBranchAddress("e2Tight", &e2Tight);
    //_ntuple->SetBranchAddress("e2UnCorrPt", &e2UnCorrPt);
    //_ntuple->SetBranchAddress("e2Veto", &e2Veto);
    //_ntuple->SetBranchAddress("e2ZZIso", &e2ZZIso);
    //_ntuple->SetBranchAddress("e2IsEB", &e2IsEB);
    //_ntuple->SetBranchAddress("e2IsGap", &e2IsGap);
    //_ntuple->SetBranchAddress("e2ZZIsoPass", &e2ZZIsoPass);
    //_ntuple->SetBranchAddress("e2ZZLooseID", &e2ZZLooseID);
    //_ntuple->SetBranchAddress("e2ZZLooseIDNoVtx", &e2ZZLooseIDNoVtx);
    //_ntuple->SetBranchAddress("e2ZZTightID", &e2ZZTightID);
    //_ntuple->SetBranchAddress("e2ZZTightIDNoVtx", &e2ZZTightIDNoVtx);
    //_ntuple->SetBranchAddress("e2Charge", &e2Charge);
    if (_isMC){
      //_ntuple->SetBranchAddress("e2GenCharge", &e2GenCharge);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("e2GenPdgId", &e2GenPdgId);
    }
    _ntuple->SetBranchAddress("e2PdgId", &e2PdgId);
    //_ntuple->SetBranchAddress("e2MissingHits", &e2MissingHits);
  }

  if (_channel == "mmmm" || _channel == "eemm"){
    //_ntuple->SetBranchAddress("m1_m2_DR", &m1_m2_DR);
    //_ntuple->SetBranchAddress("m1_m2_Energy", &m1_m2_Energy);
    //_ntuple->SetBranchAddress("m1_m2_EnergyNoFSR", &m1_m2_EnergyNoFSR);
    //_ntuple->SetBranchAddress("m1_m2_Eta", &m1_m2_Eta);
    //_ntuple->SetBranchAddress("m1_m2_EtaNoFSR", &m1_m2_EtaNoFSR);
    _ntuple->SetBranchAddress("m1_m2_Mass", &m1_m2_Mass);
    //_ntuple->SetBranchAddress("m1_m2_MassNoFSR", &m1_m2_MassNoFSR);
    //_ntuple->SetBranchAddress("m1_m2_Mt", &m1_m2_Mt);
    //_ntuple->SetBranchAddress("m1_m2_MtToMET", &m1_m2_MtToMET);
    //_ntuple->SetBranchAddress("m1_m2_Phi", &m1_m2_Phi);
    //_ntuple->SetBranchAddress("m1_m2_PhiNoFSR", &m1_m2_PhiNoFSR);
    //_ntuple->SetBranchAddress("m1_m2_Pt", &m1_m2_Pt);
    //_ntuple->SetBranchAddress("m1_m2_PtNoFSR", &m1_m2_PtNoFSR);
    //_ntuple->SetBranchAddress("m1_m2_SS", &m1_m2_SS);
    //_ntuple->SetBranchAddress("m1_m2_Charge", &m1_m2_Charge);
    //_ntuple->SetBranchAddress("m1_m2_PdgId", &m1_m2_PdgId);
    //_ntuple->SetBranchAddress("m1EffScaleFactor", &m1EffScaleFactor);
    //_ntuple->SetBranchAddress("m1EffScaleFactorError", &m1EffScaleFactorError);
    _ntuple->SetBranchAddress("m1Energy", &m1Energy);
    _ntuple->SetBranchAddress("m1Eta", &m1Eta);
    //_ntuple->SetBranchAddress("m1FSRDREt", &m1FSRDREt);
    //_ntuple->SetBranchAddress("m1FSREt", &m1FSREt);
    //_ntuple->SetBranchAddress("m1FSREta", &m1FSREta);
    //_ntuple->SetBranchAddress("m1FSRPhi", &m1FSRPhi);
    if (_isMC){
      //_ntuple->SetBranchAddress("m1GenEta", &m1GenEta);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("m1GenPhi", &m1GenPhi);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("m1GenPt", &m1GenPt);
    }
    //_ntuple->SetBranchAddress("m1IP2D", &m1IP2D);
    //_ntuple->SetBranchAddress("m1IP2DUncertainty", &m1IP2DUncertainty);
    //_ntuple->SetBranchAddress("m1IP3D", &m1IP3D);
    //_ntuple->SetBranchAddress("m1IP3DUncertainty", &m1IP3DUncertainty);
    //_ntuple->SetBranchAddress("m1Mass", &m1Mass);
    //_ntuple->SetBranchAddress("m1MtToMET", &m1MtToMET);
    //_ntuple->SetBranchAddress("m1PFChargedIso", &m1PFChargedIso);
    //_ntuple->SetBranchAddress("m1PFNeutralIso", &m1PFNeutralIso);
    //_ntuple->SetBranchAddress("m1PFPUIso", &m1PFPUIso);
    //_ntuple->SetBranchAddress("m1PFPhotonIso", &m1PFPhotonIso);
    //_ntuple->SetBranchAddress("m1PVDXY", &m1PVDXY);
    //_ntuple->SetBranchAddress("m1PVDZ", &m1PVDZ);
    _ntuple->SetBranchAddress("m1Phi", &m1Phi);
    _ntuple->SetBranchAddress("m1Pt", &m1Pt);
    //_ntuple->SetBranchAddress("m1PtErr", &m1PtErr);
    //_ntuple->SetBranchAddress("m1PtErrUncorrected", &m1PtErrUncorrected);
    //_ntuple->SetBranchAddress("m1PtUncorrected", &m1PtUncorrected);
    //_ntuple->SetBranchAddress("m1RelPFIsoDBR03", &m1RelPFIsoDBR03);
    //_ntuple->SetBranchAddress("m1RelPFIsoDBR04", &m1RelPFIsoDBR04);
    //_ntuple->SetBranchAddress("m1RochesterScaleUncDn", &m1RochesterScaleUncDn);
    //_ntuple->SetBranchAddress("m1RochesterScaleUncUp", &m1RochesterScaleUncUp);
    //_ntuple->SetBranchAddress("m1RochesterSmearUncDn", &m1RochesterSmearUncDn);
    //_ntuple->SetBranchAddress("m1RochesterSmearUncUp", &m1RochesterSmearUncUp);
    //_ntuple->SetBranchAddress("m1SIP2D", &m1SIP2D);
    //_ntuple->SetBranchAddress("m1SIP3D", &m1SIP3D);
    //_ntuple->SetBranchAddress("m1TrackIso", &m1TrackIso);
    //_ntuple->SetBranchAddress("m1ZZIso", &m1ZZIso);
    //_ntuple->SetBranchAddress("m1HighPtID", &m1HighPtID);
    //_ntuple->SetBranchAddress("m1HighPtIDNoVtx", &m1HighPtIDNoVtx);
    //_ntuple->SetBranchAddress("m1IsGlobal", &m1IsGlobal);
    //_ntuple->SetBranchAddress("m1IsLoose", &m1IsLoose);
    //_ntuple->SetBranchAddress("m1IsMedium", &m1IsMedium);
    //_ntuple->SetBranchAddress("m1IsPFMuon", &m1IsPFMuon);
    //_ntuple->SetBranchAddress("m1IsTracker", &m1IsTracker);
    //_ntuple->SetBranchAddress("m1PASTightID", &m1PASTightID);
    //_ntuple->SetBranchAddress("m1PASTightIDNoVtx", &m1PASTightIDNoVtx);
    //_ntuple->SetBranchAddress("m1PFID", &m1PFID);
    //_ntuple->SetBranchAddress("m1PFIDNoVtx", &m1PFIDNoVtx);
    //_ntuple->SetBranchAddress("m1ZZIsoPass", &m1ZZIsoPass);
    //_ntuple->SetBranchAddress("m1ZZLooseID", &m1ZZLooseID);
    //_ntuple->SetBranchAddress("m1ZZLooseIDNoVtx", &m1ZZLooseIDNoVtx);
    //_ntuple->SetBranchAddress("m1ZZTightID", &m1ZZTightID);
    //_ntuple->SetBranchAddress("m1ZZTightIDNoVtx", &m1ZZTightIDNoVtx);
    //_ntuple->SetBranchAddress("m1Charge", &m1Charge);
    if (_isMC){
      //_ntuple->SetBranchAddress("m1GenCharge", &m1GenCharge);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("m1GenPdgId", &m1GenPdgId);
    }
    _ntuple->SetBranchAddress("m1PdgId", &m1PdgId);
    //_ntuple->SetBranchAddress("m1BestTrackType", &m1BestTrackType);
    //_ntuple->SetBranchAddress("m1CutBasedLoose", &m1CutBasedLoose);
    //_ntuple->SetBranchAddress("m1CutBasedMedium", &m1CutBasedMedium);
    //_ntuple->SetBranchAddress("m1CutBasedTight", &m1CutBasedTight);
    //_ntuple->SetBranchAddress("m1IsTight", &m1IsTight);
    //_ntuple->SetBranchAddress("m1MatchedStations", &m1MatchedStations);
    //_ntuple->SetBranchAddress("m1NoOfMatches", &m1NoOfMatches);
    //_ntuple->SetBranchAddress("m1PFIsoLoose", &m1PFIsoLoose);
    //_ntuple->SetBranchAddress("m1PFIsoMedium", &m1PFIsoMedium);
    //_ntuple->SetBranchAddress("m1PFIsoTight", &m1PFIsoTight);
    //_ntuple->SetBranchAddress("m1PFIsoVTight", &m1PFIsoVTight);
    //_ntuple->SetBranchAddress("m2EffScaleFactor", &m2EffScaleFactor);
    //_ntuple->SetBranchAddress("m2EffScaleFactorError", &m2EffScaleFactorError);
    _ntuple->SetBranchAddress("m2Energy", &m2Energy);
    _ntuple->SetBranchAddress("m2Eta", &m2Eta);
    //_ntuple->SetBranchAddress("m2FSRDREt", &m2FSRDREt);
    //_ntuple->SetBranchAddress("m2FSREt", &m2FSREt);
    //_ntuple->SetBranchAddress("m2FSREta", &m2FSREta);
    //_ntuple->SetBranchAddress("m2FSRPhi", &m2FSRPhi);
    if (_isMC){
      //_ntuple->SetBranchAddress("m2GenEta", &m2GenEta);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("m2GenPhi", &m2GenPhi);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("m2GenPt", &m2GenPt);
    }
    //_ntuple->SetBranchAddress("m2IP2D", &m2IP2D);
    //_ntuple->SetBranchAddress("m2IP2DUncertainty", &m2IP2DUncertainty);
    //_ntuple->SetBranchAddress("m2IP3D", &m2IP3D);
    //_ntuple->SetBranchAddress("m2IP3DUncertainty", &m2IP3DUncertainty);
    //_ntuple->SetBranchAddress("m2Mass", &m2Mass);
    //_ntuple->SetBranchAddress("m2MtToMET", &m2MtToMET);
    //_ntuple->SetBranchAddress("m2PFChargedIso", &m2PFChargedIso);
    //_ntuple->SetBranchAddress("m2PFNeutralIso", &m2PFNeutralIso);
    //_ntuple->SetBranchAddress("m2PFPUIso", &m2PFPUIso);
    //_ntuple->SetBranchAddress("m2PFPhotonIso", &m2PFPhotonIso);
    //_ntuple->SetBranchAddress("m2PVDXY", &m2PVDXY);
    //_ntuple->SetBranchAddress("m2PVDZ", &m2PVDZ);
    _ntuple->SetBranchAddress("m2Phi", &m2Phi);
    _ntuple->SetBranchAddress("m2Pt", &m2Pt);
    //_ntuple->SetBranchAddress("m2PtErr", &m2PtErr);
    //_ntuple->SetBranchAddress("m2PtErrUncorrected", &m2PtErrUncorrected);
    //_ntuple->SetBranchAddress("m2PtUncorrected", &m2PtUncorrected);
    //_ntuple->SetBranchAddress("m2RelPFIsoDBR03", &m2RelPFIsoDBR03);
    //_ntuple->SetBranchAddress("m2RelPFIsoDBR04", &m2RelPFIsoDBR04);
    //_ntuple->SetBranchAddress("m2RochesterScaleUncDn", &m2RochesterScaleUncDn);
    //_ntuple->SetBranchAddress("m2RochesterScaleUncUp", &m2RochesterScaleUncUp);
    //_ntuple->SetBranchAddress("m2RochesterSmearUncDn", &m2RochesterSmearUncDn);
    //_ntuple->SetBranchAddress("m2RochesterSmearUncUp", &m2RochesterSmearUncUp);
    //_ntuple->SetBranchAddress("m2SIP2D", &m2SIP2D);
    //_ntuple->SetBranchAddress("m2SIP3D", &m2SIP3D);
    //_ntuple->SetBranchAddress("m2TrackIso", &m2TrackIso);
    //_ntuple->SetBranchAddress("m2ZZIso", &m2ZZIso);
    //_ntuple->SetBranchAddress("m2HighPtID", &m2HighPtID);
    //_ntuple->SetBranchAddress("m2HighPtIDNoVtx", &m2HighPtIDNoVtx);
    //_ntuple->SetBranchAddress("m2IsGlobal", &m2IsGlobal);
    //_ntuple->SetBranchAddress("m2IsLoose", &m2IsLoose);
    //_ntuple->SetBranchAddress("m2IsMedium", &m2IsMedium);
    //_ntuple->SetBranchAddress("m2IsPFMuon", &m2IsPFMuon);
    //_ntuple->SetBranchAddress("m2IsTracker", &m2IsTracker);
    //_ntuple->SetBranchAddress("m2PASTightID", &m2PASTightID);
    //_ntuple->SetBranchAddress("m2PASTightIDNoVtx", &m2PASTightIDNoVtx);
    //_ntuple->SetBranchAddress("m2PFID", &m2PFID);
    //_ntuple->SetBranchAddress("m2PFIDNoVtx", &m2PFIDNoVtx);
    //_ntuple->SetBranchAddress("m2ZZIsoPass", &m2ZZIsoPass);
    //_ntuple->SetBranchAddress("m2ZZLooseID", &m2ZZLooseID);
    //_ntuple->SetBranchAddress("m2ZZLooseIDNoVtx", &m2ZZLooseIDNoVtx);
    //_ntuple->SetBranchAddress("m2ZZTightID", &m2ZZTightID);
    //_ntuple->SetBranchAddress("m2ZZTightIDNoVtx", &m2ZZTightIDNoVtx);
    //_ntuple->SetBranchAddress("m2Charge", &m2Charge);
    if (_isMC){
      //_ntuple->SetBranchAddress("m2GenCharge", &m2GenCharge);
    }
    if (_isMC){
      //_ntuple->SetBranchAddress("m2GenPdgId", &m2GenPdgId);
    }
    _ntuple->SetBranchAddress("m2PdgId", &m2PdgId);
    //_ntuple->SetBranchAddress("m2BestTrackType", &m2BestTrackType);
    //_ntuple->SetBranchAddress("m2CutBasedLoose", &m2CutBasedLoose);
    //_ntuple->SetBranchAddress("m2CutBasedMedium", &m2CutBasedMedium);
    //_ntuple->SetBranchAddress("m2CutBasedTight", &m2CutBasedTight);
    //_ntuple->SetBranchAddress("m2IsTight", &m2IsTight);
    //_ntuple->SetBranchAddress("m2MatchedStations", &m2MatchedStations);
    //_ntuple->SetBranchAddress("m2NoOfMatches", &m2NoOfMatches);
    //_ntuple->SetBranchAddress("m2PFIsoLoose", &m2PFIsoLoose);
    //_ntuple->SetBranchAddress("m2PFIsoMedium", &m2PFIsoMedium);
    //_ntuple->SetBranchAddress("m2PFIsoTight", &m2PFIsoTight);
    //_ntuple->SetBranchAddress("m2PFIsoVTight", &m2PFIsoVTight);
  }
}

void ZZLooperBase::SetTitles(TH1 *hist, std::string x, std::string y) const{
  if (hist->InheritsFrom("TH3")) return;
  hist->GetXaxis()->SetTitle(x.c_str());
  hist->GetYaxis()->SetTitle(y.c_str());
  //hist->GetYaxis()->SetTitleOffset(1.4);
}

double ZZLooperBase::GetDeltaPhi(double phi1, double phi2) const{
  double deltaphi = phi1 - phi2;
  while (deltaphi < 0) deltaphi += (2 * M_PI);
  while (deltaphi > M_PI) deltaphi -= (2*M_PI);
  return deltaphi;
}

double ZZLooperBase::GetPolCosTheta(const ROOT::Math::PtEtaPhiEVector &lp, const ROOT::Math::PtEtaPhiEVector &ln) const{
  ROOT::Math::PtEtaPhiEVector z = lp+ln;
  ROOT::Math::Boost boost(z.BoostToCM());
  ROOT::Math::PtEtaPhiEVector lep_boost = boost(ln);
  return lep_boost.Vect().Dot(z.Vect()) / std::sqrt( lep_boost.Vect().Mag2() * z.Vect().Mag2() );
}

void ZZLooperBase::MakeIndexFile(std::string filename) const{
  try{
    std::ifstream fin("inputs/index.php");
    std::ofstream fout((filename + "/index.php").c_str());
    fout << fin.rdbuf();
  }
  catch (...){
    std::cerr << "Error making index file " << filename << std::endl;
  }
}

#endif//ZZLOOPERBASE_H
