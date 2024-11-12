#ifndef COMPLOOPERBASE_H
#define COMPLOOPERBASE_H

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
class CompLooperBase{
  public:
    CompLooperBase(const char *name, const char *channel, const char *filename1, const char *filename2);
    ~CompLooperBase();

    void SetLumi(float lumi=1.){_lumi = lumi;}
    void SetMC1(bool isMC=true){_isT1MC = isMC;}
    void SetMC2(bool isMC=true){_isT2MC = isMC;}
    void SetXsec1(float xsec=1.){_xsec1 = xsec;}
    void SetXsec2(float xsec=1.){_xsec2 = xsec;}
		void SetKfac1(float kfac=1.){_kfac1 = kfac;}
		void SetKfac2(float kfac=1.){_kfac2 = kfac;}
    void SetMC(bool isMC1=true, bool isMC2=true){SetMC1(isMC1); SetMC2(isMC2);}
    void SetXsec(float xsec1=1., float xsec2=1.){SetXsec1(xsec1); SetXsec2(xsec2);}
		void SetKfac(float kfac1=1., float kfac2=1.){SetKfac1(kfac1); SetKfac2(kfac2);}

    virtual void Loop() = 0;
  protected:
    void Init();

    void SetTitles(TH1 *hist, std::string x, std::string y="Number of Entries") const;
    double GetDeltaPhi(double phi1, double phi2) const;
    double GetPolCosTheta(const ROOT::Math::PtEtaPhiEVector &lep, const ROOT::Math::PtEtaPhiEVector &z) const;

    void MakeIndexFile(std::string filename) const;

    const double Z_MASS = 91.187; // GeV, given by PDG 2020

    TFile *_infile1, *_infile2;
    TTree *_ntuple1, *_ntuple2;
		Float_t _lumi=1.;
    Float_t _xsec1=1., _xsec2=1.;
		Float_t _kfac1=1., _kfac2=1.;
    std::string _name;
    std::string _channel;

    bool _hasInit = false;
    bool _isT1MC = false, _isT2MC = false;

    //Shared branches
    Float_t summedWeights1;
    Float_t summedWeights2;
    //Float_t Energy1;
    //Float_t Energy2;
    //Float_t EnergyNoFSR1;
    //Float_t EnergyNoFSR2;
    //Float_t Eta1;
    //Float_t Eta2;
    //Float_t EtaNoFSR1;
    //Float_t EtaNoFSR2;
    //Float_t Flag_BadChargedCandidateFilterPass1;
    //Float_t Flag_BadChargedCandidateFilterPass2;
    //Float_t Flag_BadPFMuonFilterPass1;
    //Float_t Flag_BadPFMuonFilterPass2;
    Float_t Mass1;
    Float_t Mass2;
    //Float_t MassNoFSR1;
    //Float_t MassNoFSR2;
    //Float_t MtToMET1;
    //Float_t MtToMET2;
    //Float_t Phi1;
    //Float_t Phi2;
    //Float_t PhiNoFSR1;
    //Float_t PhiNoFSR2;
    //Float_t Pt1;
    //Float_t Pt2;
    //Float_t PtNoFSR1;
    //Float_t PtNoFSR2;
    //Float_t deltaEtajj1;
    //Float_t deltaEtajj2;
    //Float_t deltaEtajj_jerDown1;
    //Float_t deltaEtajj_jerDown2;
    //Float_t deltaEtajj_jerUp1;
    //Float_t deltaEtajj_jerUp2;
    //Float_t deltaEtajj_jesDown1;
    //Float_t deltaEtajj_jesDown2;
    //Float_t deltaEtajj_jesUp1;
    //Float_t deltaEtajj_jesUp2;
    //Float_t deltaPhiTojj1;
    //Float_t deltaPhiTojj2;
    //Float_t deltaPhiTojj_jerDown1;
    //Float_t deltaPhiTojj_jerDown2;
    //Float_t deltaPhiTojj_jerUp1;
    //Float_t deltaPhiTojj_jerUp2;
    //Float_t deltaPhiTojj_jesDown1;
    //Float_t deltaPhiTojj_jesDown2;
    //Float_t deltaPhiTojj_jesUp1;
    //Float_t deltaPhiTojj_jesUp2;
    //Float_t etajj1;
    //Float_t etajj2;
    //Float_t etajj_jerDown1;
    //Float_t etajj_jerDown2;
    //Float_t etajj_jerUp1;
    //Float_t etajj_jerUp2;
    //Float_t etajj_jesDown1;
    //Float_t etajj_jesDown2;
    //Float_t etajj_jesUp1;
    //Float_t etajj_jesUp2;
    Float_t genWeight1;
    Float_t genWeight2;
    //Float_t jetPUSFmulfac1;
    //Float_t jetPUSFmulfac2;
    //Float_t maxScaleWeight1;
    //Float_t maxScaleWeight2;
    //Float_t minScaleWeight1;
    //Float_t minScaleWeight2;
    //Float_t mjj1;
    //Float_t mjj2;
    //Float_t mjj_jerDown1;
    //Float_t mjj_jerDown2;
    //Float_t mjj_jerUp1;
    //Float_t mjj_jerUp2;
    //Float_t mjj_jesDown1;
    //Float_t mjj_jesDown2;
    //Float_t mjj_jesUp1;
    //Float_t mjj_jesUp2;
    //Float_t nTruePU1;
    //Float_t nTruePU2;
    //Float_t originalXWGTUP1;
    //Float_t originalXWGTUP2;
    //Float_t phijj1;
    //Float_t phijj2;
    //Float_t phijj_jerDown1;
    //Float_t phijj_jerDown2;
    //Float_t phijj_jerUp1;
    //Float_t phijj_jerUp2;
    //Float_t phijj_jesDown1;
    //Float_t phijj_jesDown2;
    //Float_t phijj_jesUp1;
    //Float_t phijj_jesUp2;
    //Float_t ptjj1;
    //Float_t ptjj2;
    //Float_t ptjj_jerDown1;
    //Float_t ptjj_jerDown2;
    //Float_t ptjj_jerUp1;
    //Float_t ptjj_jerUp2;
    //Float_t ptjj_jesDown1;
    //Float_t ptjj_jesDown2;
    //Float_t ptjj_jesUp1;
    //Float_t ptjj_jesUp2;
    //Float_t pvRho1;
    //Float_t pvRho2;
    //Float_t pvZ1;
    //Float_t pvZ2;
    //Float_t pvndof1;
    //Float_t pvndof2;
    //Float_t rapidity1;
    //Float_t rapidity2;
    //Float_t type1_pfMETEt1;
    //Float_t type1_pfMETEt2;
    //Float_t type1_pfMETPhi1;
    //Float_t type1_pfMETPhi2;
    //Float_t zeppenfeld1;
    //Float_t zeppenfeld2;
    //Float_t zeppenfeld_jerDown1;
    //Float_t zeppenfeld_jerDown2;
    //Float_t zeppenfeld_jerUp1;
    //Float_t zeppenfeld_jerUp2;
    //Float_t zeppenfeld_jesDown1;
    //Float_t zeppenfeld_jesDown2;
    //Float_t zeppenfeld_jesUp1;
    //Float_t zeppenfeld_jesUp2;
    //Float_t zeppenfeldj31;
    //Float_t zeppenfeldj32;
    //Float_t zeppenfeldj3_jerDown1;
    //Float_t zeppenfeldj3_jerDown2;
    //Float_t zeppenfeldj3_jerUp1;
    //Float_t zeppenfeldj3_jerUp2;
    //Float_t zeppenfeldj3_jesDown1;
    //Float_t zeppenfeldj3_jesDown2;
    //Float_t zeppenfeldj3_jesUp1;
    //Float_t zeppenfeldj3_jesUp2;
    //Bool_t pvIdFake1;
    //Bool_t pvIdFake2;
    //Bool_t pvIsValid1;
    //Bool_t pvIsValid2;
    //Int_t Charge1;
    //Int_t Charge2;
    //Int_t PdgId1;
    //Int_t PdgId2;
    UInt_t lumi1;
    UInt_t lumi2;
    //UInt_t nJets1;
    //UInt_t nJets2;
    //UInt_t nJets_jerDown1;
    //UInt_t nJets_jerDown2;
    //UInt_t nJets_jerUp1;
    //UInt_t nJets_jerUp2;
    //UInt_t nJets_jesDown1;
    //UInt_t nJets_jesDown2;
    //UInt_t nJets_jesUp1;
    //UInt_t nJets_jesUp2;
    //UInt_t nZZLooseElec1;
    //UInt_t nZZLooseElec2;
    //UInt_t nZZLooseMu1;
    //UInt_t nZZLooseMu2;
    //UInt_t nZZTightElec1;
    //UInt_t nZZTightElec2;
    //UInt_t nZZTightIsoElec1;
    //UInt_t nZZTightIsoElec2;
    //UInt_t nZZTightIsoMu1;
    //UInt_t nZZTightIsoMu2;
    //UInt_t nZZTightMu1;
    //UInt_t nZZTightMu2;
    UInt_t nvtx1;
    UInt_t nvtx2;
    UInt_t run1;
    UInt_t run2;
    ULong64_t evt1;
    ULong64_t evt2;
    //std::vector<float> *jetDeepCSV1 = 0;
    //std::vector<float> *jetDeepCSV2 = 0;
    //std::vector<float> *jetEta1 = 0;
    //std::vector<float> *jetEta2 = 0;
    //std::vector<float> *jetEta_jerDown1 = 0;
    //std::vector<float> *jetEta_jerDown2 = 0;
    //std::vector<float> *jetEta_jerUp1 = 0;
    //std::vector<float> *jetEta_jerUp2 = 0;
    //std::vector<float> *jetEta_jesDown1 = 0;
    //std::vector<float> *jetEta_jesDown2 = 0;
    //std::vector<float> *jetEta_jesUp1 = 0;
    //std::vector<float> *jetEta_jesUp2 = 0;
    //std::vector<float> *jetPhi1 = 0;
    //std::vector<float> *jetPhi2 = 0;
    //std::vector<float> *jetPt1 = 0;
    //std::vector<float> *jetPt2 = 0;
    //std::vector<float> *jetPt_jerDown1 = 0;
    //std::vector<float> *jetPt_jerDown2 = 0;
    //std::vector<float> *jetPt_jerUp1 = 0;
    //std::vector<float> *jetPt_jerUp2 = 0;
    //std::vector<float> *jetPt_jesDown1 = 0;
    //std::vector<float> *jetPt_jesDown2 = 0;
    //std::vector<float> *jetPt_jesUp1 = 0;
    //std::vector<float> *jetPt_jesUp2 = 0;
    //std::vector<float> *jetRapidity1 = 0;
    //std::vector<float> *jetRapidity2 = 0;
    //std::vector<float> *jetRapidity_jerDown1 = 0;
    //std::vector<float> *jetRapidity_jerDown2 = 0;
    //std::vector<float> *jetRapidity_jerUp1 = 0;
    //std::vector<float> *jetRapidity_jerUp2 = 0;
    //std::vector<float> *jetRapidity_jesDown1 = 0;
    //std::vector<float> *jetRapidity_jesDown2 = 0;
    //std::vector<float> *jetRapidity_jesUp1 = 0;
    //std::vector<float> *jetRapidity_jesUp2 = 0;
    //std::vector<float> *pdfWeights1 = 0;
    //std::vector<float> *pdfWeights2 = 0;
    //std::vector<float> *scaleWeightIDs1 = 0;
    //std::vector<float> *scaleWeightIDs2 = 0;
    //std::vector<float> *scaleWeights1 = 0;
    //std::vector<float> *scaleWeights2 = 0;
    //std::vector<int> *isGenJetMatched1 = 0;
    //std::vector<int> *isGenJetMatched2 = 0;
    //std::vector<int> *jetHadronFlavor1 = 0;
    //std::vector<int> *jetHadronFlavor2 = 0;
    //std::vector<int> *jetPUID1 = 0;
    //std::vector<int> *jetPUID2 = 0;
    //std::vector<int> *jetPUID_jerDown1 = 0;
    //std::vector<int> *jetPUID_jerDown2 = 0;
    //std::vector<int> *jetPUID_jerUp1 = 0;
    //std::vector<int> *jetPUID_jerUp2 = 0;
    //std::vector<int> *jetPUID_jesDown1 = 0;
    //std::vector<int> *jetPUID_jesDown2 = 0;
    //std::vector<int> *jetPUID_jesUp1 = 0;
    //std::vector<int> *jetPUID_jesUp2 = 0;
    //Bool_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVLPass1;
    //Bool_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVLPass2;
    //Bool_t HLT_DoubleEle25_CaloIdL_MWPass1;
    //Bool_t HLT_DoubleEle25_CaloIdL_MWPass2;
    //Bool_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8Pass1;
    //Bool_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8Pass2;
    //Bool_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVLPass1;
    //Bool_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVLPass2;
    //Bool_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass1;
    //Bool_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass2;
    //Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass1;
    //Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass2;
    //Bool_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass1;
    //Bool_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass2;
    //Bool_t HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZPass1;
    //Bool_t HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZPass2;
    //Bool_t HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZPass1;
    //Bool_t HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZPass2;
    //Bool_t HLT_TripleMu_10_5_5_DZPass1;
    //Bool_t HLT_TripleMu_10_5_5_DZPass2;
    //Bool_t HLT_TripleMu_12_10_5Pass1;
    //Bool_t HLT_TripleMu_12_10_5Pass2;
    //Bool_t HLT_Ele30_WPTight_GsfPass1;
    //Bool_t HLT_Ele30_WPTight_GsfPass2;
    //Bool_t HLT_Ele32_WPTight_GsfPass1;
    //Bool_t HLT_Ele32_WPTight_GsfPass2;
    //Bool_t HLT_IsoMu24Pass1;
    //Bool_t HLT_IsoMu24Pass2;
    //Bool_t Flag_goodVerticesPass1;
    //Bool_t Flag_goodVerticesPass2;
    //Bool_t Flag_globalSuperTightHalo2016FilterPass1;
    //Bool_t Flag_globalSuperTightHalo2016FilterPass2;
    //Bool_t Flag_EcalDeadCellTriggerPrimitiveFilterPass1;
    //Bool_t Flag_EcalDeadCellTriggerPrimitiveFilterPass2;
    //Float_t Flag_BadPFMuonFilterPass1;
    //Float_t Flag_BadPFMuonFilterPass2;
    //Bool_t Flag_BadPFMuonDzFilterPass1;
    //Bool_t Flag_BadPFMuonDzFilterPass2;
    //Bool_t Flag_hfNoisyHitsFilterPass1;
    //Bool_t Flag_hfNoisyHitsFilterPass2;
    //Bool_t Flag_eeBadScFilterPass1;
    //Bool_t Flag_eeBadScFilterPass2;

    //eeee
    //Float_t e1_e3_DR1;
    //Float_t e1_e3_DR2;
    //Float_t e1_e3_Mass1;
    //Float_t e1_e3_Mass2;
    //Float_t e1_e3_MassNoFSR1;
    //Float_t e1_e3_MassNoFSR2;
    //Float_t e1_e4_DR1;
    //Float_t e1_e4_DR2;
    //Float_t e1_e4_Mass1;
    //Float_t e1_e4_Mass2;
    //Float_t e1_e4_MassNoFSR1;
    //Float_t e1_e4_MassNoFSR2;
    //Float_t e2_e3_DR1;
    //Float_t e2_e3_DR2;
    //Float_t e2_e3_Mass1;
    //Float_t e2_e3_Mass2;
    //Float_t e2_e3_MassNoFSR1;
    //Float_t e2_e3_MassNoFSR2;
    //Float_t e2_e4_DR1;
    //Float_t e2_e4_DR2;
    //Float_t e2_e4_Mass1;
    //Float_t e2_e4_Mass2;
    //Float_t e2_e4_MassNoFSR1;
    //Float_t e2_e4_MassNoFSR2;
    //Bool_t e1_e3_SS1;
    //Bool_t e1_e3_SS2;
    //Bool_t e1_e4_SS1;
    //Bool_t e1_e4_SS2;
    //Bool_t e2_e3_SS1;
    //Bool_t e2_e3_SS2;
    //Bool_t e2_e4_SS1;
    //Bool_t e2_e4_SS2;
    //Float_t e3_e4_DR1;
    //Float_t e3_e4_DR2;
    //Float_t e3_e4_Energy1;
    //Float_t e3_e4_Energy2;
    //Float_t e3_e4_EnergyNoFSR1;
    //Float_t e3_e4_EnergyNoFSR2;
    //Float_t e3_e4_Eta1;
    //Float_t e3_e4_Eta2;
    //Float_t e3_e4_EtaNoFSR1;
    //Float_t e3_e4_EtaNoFSR2;
    Float_t e3_e4_Mass1;
    Float_t e3_e4_Mass2;
    //Float_t e3_e4_MassNoFSR1;
    //Float_t e3_e4_MassNoFSR2;
    //Float_t e3_e4_Mt1;
    //Float_t e3_e4_Mt2;
    //Float_t e3_e4_MtToMET1;
    //Float_t e3_e4_MtToMET2;
    //Float_t e3_e4_Phi1;
    //Float_t e3_e4_Phi2;
    //Float_t e3_e4_PhiNoFSR1;
    //Float_t e3_e4_PhiNoFSR2;
    //Float_t e3_e4_Pt1;
    //Float_t e3_e4_Pt2;
    //Float_t e3_e4_PtNoFSR1;
    //Float_t e3_e4_PtNoFSR2;
    //Bool_t e3_e4_SS1;
    //Bool_t e3_e4_SS2;
    //Int_t e3_e4_Charge1;
    //Int_t e3_e4_Charge2;
    //Int_t e3_e4_PdgId1;
    //Int_t e3_e4_PdgId2;
    //Float_t e3EffectiveArea1;
    //Float_t e3EffectiveArea2;
    Float_t e3Energy1;
    Float_t e3Energy2;
    Float_t e3Eta1;
    Float_t e3Eta2;
    //Float_t e3FSRDREt1;
    //Float_t e3FSRDREt2;
    //Float_t e3FSREt1;
    //Float_t e3FSREt2;
    //Float_t e3FSREta1;
    //Float_t e3FSREta2;
    //Float_t e3FSRPhi1;
    //Float_t e3FSRPhi2;
    //Float_t e3GenEta1;
    //Float_t e3GenEta2;
    //Float_t e3GenPhi1;
    //Float_t e3GenPhi2;
    //Float_t e3GenPt1;
    //Float_t e3GenPt2;
    //Float_t e3IP2D1;
    //Float_t e3IP2D2;
    //Float_t e3IP2DUncertainty1;
    //Float_t e3IP2DUncertainty2;
    //Float_t e3IP3D1;
    //Float_t e3IP3D2;
    //Float_t e3IP3DUncertainty1;
    //Float_t e3IP3DUncertainty2;
    //Float_t e3Loose1;
    //Float_t e3Loose2;
    //Float_t e3Mass1;
    //Float_t e3Mass2;
    //Float_t e3Medium1;
    //Float_t e3Medium2;
    //Float_t e3MtToMET1;
    //Float_t e3MtToMET2;
    //Float_t e3PFChargedIso1;
    //Float_t e3PFChargedIso2;
    //Float_t e3PFNeutralIso1;
    //Float_t e3PFNeutralIso2;
    //Float_t e3PFPUIso1;
    //Float_t e3PFPUIso2;
    //Float_t e3PFPhotonIso1;
    //Float_t e3PFPhotonIso2;
    //Float_t e3PVDXY1;
    //Float_t e3PVDXY2;
    //Float_t e3PVDZ1;
    //Float_t e3PVDZ2;
    Float_t e3Phi1;
    Float_t e3Phi2;
    Float_t e3Pt1;
    Float_t e3Pt2;
    //Float_t e3RelPFIsoRho1;
    //Float_t e3RelPFIsoRho2;
    //Float_t e3Rho1;
    //Float_t e3Rho2;
    //Float_t e3SCEnergy1;
    //Float_t e3SCEnergy2;
    //Float_t e3SCEta1;
    //Float_t e3SCEta2;
    //Float_t e3SCPhi1;
    //Float_t e3SCPhi2;
    //Float_t e3SCRawEnergy1;
    //Float_t e3SCRawEnergy2;
    //Float_t e3SIP2D1;
    //Float_t e3SIP2D2;
    //Float_t e3SIP3D1;
    //Float_t e3SIP3D2;
    //Float_t e3ScaleGainDn1;
    //Float_t e3ScaleGainDn2;
    //Float_t e3ScaleGainUp1;
    //Float_t e3ScaleGainUp2;
    //Float_t e3ScaleStatDn1;
    //Float_t e3ScaleStatDn2;
    //Float_t e3ScaleStatUp1;
    //Float_t e3ScaleStatUp2;
    //Float_t e3ScaleSystDn1;
    //Float_t e3ScaleSystDn2;
    //Float_t e3ScaleSystUp1;
    //Float_t e3ScaleSystUp2;
    //Float_t e3ScaleTotDn1;
    //Float_t e3ScaleTotDn2;
    //Float_t e3ScaleTotUp1;
    //Float_t e3ScaleTotUp2;
    //Float_t e3SigmaPhiDn1;
    //Float_t e3SigmaPhiDn2;
    //Float_t e3SigmaPhiUp1;
    //Float_t e3SigmaPhiUp2;
    //Float_t e3SigmaRhoDn1;
    //Float_t e3SigmaRhoDn2;
    //Float_t e3SigmaRhoUp1;
    //Float_t e3SigmaRhoUp2;
    //Float_t e3SigmaTotDn1;
    //Float_t e3SigmaTotDn2;
    //Float_t e3SigmaTotUp1;
    //Float_t e3SigmaTotUp2;
    //Float_t e3Tight1;
    //Float_t e3Tight2;
    //Float_t e3UnCorrPt1;
    //Float_t e3UnCorrPt2;
    //Float_t e3Veto1;
    //Float_t e3Veto2;
    //Float_t e3ZZIso1;
    //Float_t e3ZZIso2;
    //Bool_t e3IsEB1;
    //Bool_t e3IsEB2;
    //Bool_t e3IsGap1;
    //Bool_t e3IsGap2;
    //Bool_t e3ZZIsoPass1;
    //Bool_t e3ZZIsoPass2;
    //Bool_t e3ZZLooseID1;
    //Bool_t e3ZZLooseID2;
    //Bool_t e3ZZLooseIDNoVtx1;
    //Bool_t e3ZZLooseIDNoVtx2;
    //Bool_t e3ZZTightID1;
    //Bool_t e3ZZTightID2;
    //Bool_t e3ZZTightIDNoVtx1;
    //Bool_t e3ZZTightIDNoVtx2;
    //Int_t e3Charge1;
    //Int_t e3Charge2;
    //Int_t e3GenCharge1;
    //Int_t e3GenCharge2;
    //Int_t e3GenPdgId1;
    //Int_t e3GenPdgId2;
    Int_t e3PdgId1;
    Int_t e3PdgId2;
    //UInt_t e3MissingHits1;
    //UInt_t e3MissingHits2;
    //Float_t e4EffectiveArea1;
    //Float_t e4EffectiveArea2;
    Float_t e4Energy1;
    Float_t e4Energy2;
    Float_t e4Eta1;
    Float_t e4Eta2;
    //Float_t e4FSRDREt1;
    //Float_t e4FSRDREt2;
    //Float_t e4FSREt1;
    //Float_t e4FSREt2;
    //Float_t e4FSREta1;
    //Float_t e4FSREta2;
    //Float_t e4FSRPhi1;
    //Float_t e4FSRPhi2;
    //Float_t e4GenEta1;
    //Float_t e4GenEta2;
    //Float_t e4GenPhi1;
    //Float_t e4GenPhi2;
    //Float_t e4GenPt1;
    //Float_t e4GenPt2;
    //Float_t e4IP2D1;
    //Float_t e4IP2D2;
    //Float_t e4IP2DUncertainty1;
    //Float_t e4IP2DUncertainty2;
    //Float_t e4IP3D1;
    //Float_t e4IP3D2;
    //Float_t e4IP3DUncertainty1;
    //Float_t e4IP3DUncertainty2;
    //Float_t e4Loose1;
    //Float_t e4Loose2;
    //Float_t e4Mass1;
    //Float_t e4Mass2;
    //Float_t e4Medium1;
    //Float_t e4Medium2;
    //Float_t e4MtToMET1;
    //Float_t e4MtToMET2;
    //Float_t e4PFChargedIso1;
    //Float_t e4PFChargedIso2;
    //Float_t e4PFNeutralIso1;
    //Float_t e4PFNeutralIso2;
    //Float_t e4PFPUIso1;
    //Float_t e4PFPUIso2;
    //Float_t e4PFPhotonIso1;
    //Float_t e4PFPhotonIso2;
    //Float_t e4PVDXY1;
    //Float_t e4PVDXY2;
    //Float_t e4PVDZ1;
    //Float_t e4PVDZ2;
    Float_t e4Phi1;
    Float_t e4Phi2;
    Float_t e4Pt1;
    Float_t e4Pt2;
    //Float_t e4RelPFIsoRho1;
    //Float_t e4RelPFIsoRho2;
    //Float_t e4Rho1;
    //Float_t e4Rho2;
    //Float_t e4SCEnergy1;
    //Float_t e4SCEnergy2;
    //Float_t e4SCEta1;
    //Float_t e4SCEta2;
    //Float_t e4SCPhi1;
    //Float_t e4SCPhi2;
    //Float_t e4SCRawEnergy1;
    //Float_t e4SCRawEnergy2;
    //Float_t e4SIP2D1;
    //Float_t e4SIP2D2;
    //Float_t e4SIP3D1;
    //Float_t e4SIP3D2;
    //Float_t e4ScaleGainDn1;
    //Float_t e4ScaleGainDn2;
    //Float_t e4ScaleGainUp1;
    //Float_t e4ScaleGainUp2;
    //Float_t e4ScaleStatDn1;
    //Float_t e4ScaleStatDn2;
    //Float_t e4ScaleStatUp1;
    //Float_t e4ScaleStatUp2;
    //Float_t e4ScaleSystDn1;
    //Float_t e4ScaleSystDn2;
    //Float_t e4ScaleSystUp1;
    //Float_t e4ScaleSystUp2;
    //Float_t e4ScaleTotDn1;
    //Float_t e4ScaleTotDn2;
    //Float_t e4ScaleTotUp1;
    //Float_t e4ScaleTotUp2;
    //Float_t e4SigmaPhiDn1;
    //Float_t e4SigmaPhiDn2;
    //Float_t e4SigmaPhiUp1;
    //Float_t e4SigmaPhiUp2;
    //Float_t e4SigmaRhoDn1;
    //Float_t e4SigmaRhoDn2;
    //Float_t e4SigmaRhoUp1;
    //Float_t e4SigmaRhoUp2;
    //Float_t e4SigmaTotDn1;
    //Float_t e4SigmaTotDn2;
    //Float_t e4SigmaTotUp1;
    //Float_t e4SigmaTotUp2;
    //Float_t e4Tight1;
    //Float_t e4Tight2;
    //Float_t e4UnCorrPt1;
    //Float_t e4UnCorrPt2;
    //Float_t e4Veto1;
    //Float_t e4Veto2;
    //Float_t e4ZZIso1;
    //Float_t e4ZZIso2;
    //Bool_t e4IsEB1;
    //Bool_t e4IsEB2;
    //Bool_t e4IsGap1;
    //Bool_t e4IsGap2;
    //Bool_t e4ZZIsoPass1;
    //Bool_t e4ZZIsoPass2;
    //Bool_t e4ZZLooseID1;
    //Bool_t e4ZZLooseID2;
    //Bool_t e4ZZLooseIDNoVtx1;
    //Bool_t e4ZZLooseIDNoVtx2;
    //Bool_t e4ZZTightID1;
    //Bool_t e4ZZTightID2;
    //Bool_t e4ZZTightIDNoVtx1;
    //Bool_t e4ZZTightIDNoVtx2;
    //Int_t e4Charge1;
    //Int_t e4Charge2;
    //Int_t e4GenCharge1;
    //Int_t e4GenCharge2;
    //Int_t e4GenPdgId1;
    //Int_t e4GenPdgId2;
    Int_t e4PdgId1;
    Int_t e4PdgId2;
    //UInt_t e4MissingHits1;
    //UInt_t e4MissingHits2;

    //eemm
    //Float_t e1_m1_DR1;
    //Float_t e1_m1_DR2;
    //Float_t e1_m1_Mass1;
    //Float_t e1_m1_Mass2;
    //Float_t e1_m1_MassNoFSR1;
    //Float_t e1_m1_MassNoFSR2;
    //Float_t e1_m2_DR1;
    //Float_t e1_m2_DR2;
    //Float_t e1_m2_Mass1;
    //Float_t e1_m2_Mass2;
    //Float_t e1_m2_MassNoFSR1;
    //Float_t e1_m2_MassNoFSR2;
    //Float_t e2_m1_DR1;
    //Float_t e2_m1_DR2;
    //Float_t e2_m1_Mass1;
    //Float_t e2_m1_Mass2;
    //Float_t e2_m1_MassNoFSR1;
    //Float_t e2_m1_MassNoFSR2;
    //Float_t e2_m2_DR1;
    //Float_t e2_m2_DR2;
    //Float_t e2_m2_Mass1;
    //Float_t e2_m2_Mass2;
    //Float_t e2_m2_MassNoFSR1;
    //Float_t e2_m2_MassNoFSR2;
    //Bool_t e1_m1_SS1;
    //Bool_t e1_m1_SS2;
    //Bool_t e1_m2_SS1;
    //Bool_t e1_m2_SS2;
    //Bool_t e2_m1_SS1;
    //Bool_t e2_m1_SS2;
    //Bool_t e2_m2_SS1;
    //Bool_t e2_m2_SS2;

    //mmmm
    //Float_t m1_m3_DR1;
    //Float_t m1_m3_DR2;
    //Float_t m1_m3_Mass1;
    //Float_t m1_m3_Mass2;
    //Float_t m1_m3_MassNoFSR1;
    //Float_t m1_m3_MassNoFSR2;
    //Float_t m1_m4_DR1;
    //Float_t m1_m4_DR2;
    //Float_t m1_m4_Mass1;
    //Float_t m1_m4_Mass2;
    //Float_t m1_m4_MassNoFSR1;
    //Float_t m1_m4_MassNoFSR2;
    //Float_t m2_m3_DR1;
    //Float_t m2_m3_DR2;
    //Float_t m2_m3_Mass1;
    //Float_t m2_m3_Mass2;
    //Float_t m2_m3_MassNoFSR1;
    //Float_t m2_m3_MassNoFSR2;
    //Float_t m2_m4_DR1;
    //Float_t m2_m4_DR2;
    //Float_t m2_m4_Mass1;
    //Float_t m2_m4_Mass2;
    //Float_t m2_m4_MassNoFSR1;
    //Float_t m2_m4_MassNoFSR2;
    //Bool_t m1_m3_SS1;
    //Bool_t m1_m3_SS2;
    //Bool_t m1_m4_SS1;
    //Bool_t m1_m4_SS2;
    //Bool_t m2_m3_SS1;
    //Bool_t m2_m3_SS2;
    //Bool_t m2_m4_SS1;
    //Bool_t m2_m4_SS2;
    //Float_t m3_m4_DR1;
    //Float_t m3_m4_DR2;
    //Float_t m3_m4_Energy1;
    //Float_t m3_m4_Energy2;
    //Float_t m3_m4_EnergyNoFSR1;
    //Float_t m3_m4_EnergyNoFSR2;
    //Float_t m3_m4_Eta1;
    //Float_t m3_m4_Eta2;
    //Float_t m3_m4_EtaNoFSR1;
    //Float_t m3_m4_EtaNoFSR2;
    Float_t m3_m4_Mass1;
    Float_t m3_m4_Mass2;
    //Float_t m3_m4_MassNoFSR1;
    //Float_t m3_m4_MassNoFSR2;
    //Float_t m3_m4_Mt1;
    //Float_t m3_m4_Mt2;
    //Float_t m3_m4_MtToMET1;
    //Float_t m3_m4_MtToMET2;
    //Float_t m3_m4_Phi1;
    //Float_t m3_m4_Phi2;
    //Float_t m3_m4_PhiNoFSR1;
    //Float_t m3_m4_PhiNoFSR2;
    //Float_t m3_m4_Pt1;
    //Float_t m3_m4_Pt2;
    //Float_t m3_m4_PtNoFSR1;
    //Float_t m3_m4_PtNoFSR2;
    //Bool_t m3_m4_SS1;
    //Bool_t m3_m4_SS2;
    //Int_t m3_m4_Charge1;
    //Int_t m3_m4_Charge2;
    //Int_t m3_m4_PdgId1;
    //Int_t m3_m4_PdgId2;
    //Float_t m3EffScaleFactor1;
    //Float_t m3EffScaleFactor2;
    //Float_t m3EffScaleFactorError1;
    //Float_t m3EffScaleFactorError2;
    Float_t m3Energy1;
    Float_t m3Energy2;
    Float_t m3Eta1;
    Float_t m3Eta2;
    //Float_t m3FSRDREt1;
    //Float_t m3FSRDREt2;
    //Float_t m3FSREt1;
    //Float_t m3FSREt2;
    //Float_t m3FSREta1;
    //Float_t m3FSREta2;
    //Float_t m3FSRPhi1;
    //Float_t m3FSRPhi2;
    //Float_t m3GenEta1;
    //Float_t m3GenEta2;
    //Float_t m3GenPhi1;
    //Float_t m3GenPhi2;
    //Float_t m3GenPt1;
    //Float_t m3GenPt2;
    //Float_t m3IP2D1;
    //Float_t m3IP2D2;
    //Float_t m3IP2DUncertainty1;
    //Float_t m3IP2DUncertainty2;
    //Float_t m3IP3D1;
    //Float_t m3IP3D2;
    //Float_t m3IP3DUncertainty1;
    //Float_t m3IP3DUncertainty2;
    //Float_t m3Mass1;
    //Float_t m3Mass2;
    //Float_t m3MtToMET1;
    //Float_t m3MtToMET2;
    //Float_t m3PFChargedIso1;
    //Float_t m3PFChargedIso2;
    //Float_t m3PFNeutralIso1;
    //Float_t m3PFNeutralIso2;
    //Float_t m3PFPUIso1;
    //Float_t m3PFPUIso2;
    //Float_t m3PFPhotonIso1;
    //Float_t m3PFPhotonIso2;
    //Float_t m3PVDXY1;
    //Float_t m3PVDXY2;
    //Float_t m3PVDZ1;
    //Float_t m3PVDZ2;
    Float_t m3Phi1;
    Float_t m3Phi2;
    Float_t m3Pt1;
    Float_t m3Pt2;
    //Float_t m3PtErr1;
    //Float_t m3PtErr2;
    //Float_t m3PtErrUncorrected1;
    //Float_t m3PtErrUncorrected2;
    //Float_t m3PtUncorrected1;
    //Float_t m3PtUncorrected2;
    //Float_t m3RelPFIsoDBR031;
    //Float_t m3RelPFIsoDBR032;
    //Float_t m3RelPFIsoDBR041;
    //Float_t m3RelPFIsoDBR042;
    //Float_t m3RochesterScaleUncDn1;
    //Float_t m3RochesterScaleUncDn2;
    //Float_t m3RochesterScaleUncUp1;
    //Float_t m3RochesterScaleUncUp2;
    //Float_t m3RochesterSmearUncDn1;
    //Float_t m3RochesterSmearUncDn2;
    //Float_t m3RochesterSmearUncUp1;
    //Float_t m3RochesterSmearUncUp2;
    //Float_t m3SIP2D1;
    //Float_t m3SIP2D2;
    //Float_t m3SIP3D1;
    //Float_t m3SIP3D2;
    //Float_t m3TrackIso1;
    //Float_t m3TrackIso2;
    //Float_t m3ZZIso1;
    //Float_t m3ZZIso2;
    //Bool_t m3HighPtID1;
    //Bool_t m3HighPtID2;
    //Bool_t m3HighPtIDNoVtx1;
    //Bool_t m3HighPtIDNoVtx2;
    //Bool_t m3IsGlobal1;
    //Bool_t m3IsGlobal2;
    //Bool_t m3IsLoose1;
    //Bool_t m3IsLoose2;
    //Bool_t m3IsMedium1;
    //Bool_t m3IsMedium2;
    //Bool_t m3IsPFMuon1;
    //Bool_t m3IsPFMuon2;
    //Bool_t m3IsTracker1;
    //Bool_t m3IsTracker2;
    //Bool_t m3PASTightID1;
    //Bool_t m3PASTightID2;
    //Bool_t m3PASTightIDNoVtx1;
    //Bool_t m3PASTightIDNoVtx2;
    //Bool_t m3PFID1;
    //Bool_t m3PFID2;
    //Bool_t m3PFIDNoVtx1;
    //Bool_t m3PFIDNoVtx2;
    //Bool_t m3ZZIsoPass1;
    //Bool_t m3ZZIsoPass2;
    //Bool_t m3ZZLooseID1;
    //Bool_t m3ZZLooseID2;
    //Bool_t m3ZZLooseIDNoVtx1;
    //Bool_t m3ZZLooseIDNoVtx2;
    //Bool_t m3ZZTightID1;
    //Bool_t m3ZZTightID2;
    //Bool_t m3ZZTightIDNoVtx1;
    //Bool_t m3ZZTightIDNoVtx2;
    //Int_t m3Charge1;
    //Int_t m3Charge2;
    //Int_t m3GenCharge1;
    //Int_t m3GenCharge2;
    //Int_t m3GenPdgId1;
    //Int_t m3GenPdgId2;
    Int_t m3PdgId1;
    Int_t m3PdgId2;
    //UInt_t m3BestTrackType1;
    //UInt_t m3BestTrackType2;
    //UInt_t m3CutBasedLoose1;
    //UInt_t m3CutBasedLoose2;
    //UInt_t m3CutBasedMedium1;
    //UInt_t m3CutBasedMedium2;
    //UInt_t m3CutBasedTight1;
    //UInt_t m3CutBasedTight2;
    //UInt_t m3IsTight1;
    //UInt_t m3IsTight2;
    //UInt_t m3MatchedStations1;
    //UInt_t m3MatchedStations2;
    //UInt_t m3NoOfMatches1;
    //UInt_t m3NoOfMatches2;
    //UInt_t m3PFIsoLoose1;
    //UInt_t m3PFIsoLoose2;
    //UInt_t m3PFIsoMedium1;
    //UInt_t m3PFIsoMedium2;
    //UInt_t m3PFIsoTight1;
    //UInt_t m3PFIsoTight2;
    //UInt_t m3PFIsoVTight1;
    //UInt_t m3PFIsoVTight2;
    //Float_t m4EffScaleFactor1;
    //Float_t m4EffScaleFactor2;
    //Float_t m4EffScaleFactorError1;
    //Float_t m4EffScaleFactorError2;
    Float_t m4Energy1;
    Float_t m4Energy2;
    Float_t m4Eta1;
    Float_t m4Eta2;
    //Float_t m4FSRDREt1;
    //Float_t m4FSRDREt2;
    //Float_t m4FSREt1;
    //Float_t m4FSREt2;
    //Float_t m4FSREta1;
    //Float_t m4FSREta2;
    //Float_t m4FSRPhi1;
    //Float_t m4FSRPhi2;
    //Float_t m4GenEta1;
    //Float_t m4GenEta2;
    //Float_t m4GenPhi1;
    //Float_t m4GenPhi2;
    //Float_t m4GenPt1;
    //Float_t m4GenPt2;
    //Float_t m4IP2D1;
    //Float_t m4IP2D2;
    //Float_t m4IP2DUncertainty1;
    //Float_t m4IP2DUncertainty2;
    //Float_t m4IP3D1;
    //Float_t m4IP3D2;
    //Float_t m4IP3DUncertainty1;
    //Float_t m4IP3DUncertainty2;
    //Float_t m4Mass1;
    //Float_t m4Mass2;
    //Float_t m4MtToMET1;
    //Float_t m4MtToMET2;
    //Float_t m4PFChargedIso1;
    //Float_t m4PFChargedIso2;
    //Float_t m4PFNeutralIso1;
    //Float_t m4PFNeutralIso2;
    //Float_t m4PFPUIso1;
    //Float_t m4PFPUIso2;
    //Float_t m4PFPhotonIso1;
    //Float_t m4PFPhotonIso2;
    //Float_t m4PVDXY1;
    //Float_t m4PVDXY2;
    //Float_t m4PVDZ1;
    //Float_t m4PVDZ2;
    Float_t m4Phi1;
    Float_t m4Phi2;
    Float_t m4Pt1;
    Float_t m4Pt2;
    //Float_t m4PtErr1;
    //Float_t m4PtErr2;
    //Float_t m4PtErrUncorrected1;
    //Float_t m4PtErrUncorrected2;
    //Float_t m4PtUncorrected1;
    //Float_t m4PtUncorrected2;
    //Float_t m4RelPFIsoDBR031;
    //Float_t m4RelPFIsoDBR032;
    //Float_t m4RelPFIsoDBR041;
    //Float_t m4RelPFIsoDBR042;
    //Float_t m4RochesterScaleUncDn1;
    //Float_t m4RochesterScaleUncDn2;
    //Float_t m4RochesterScaleUncUp1;
    //Float_t m4RochesterScaleUncUp2;
    //Float_t m4RochesterSmearUncDn1;
    //Float_t m4RochesterSmearUncDn2;
    //Float_t m4RochesterSmearUncUp1;
    //Float_t m4RochesterSmearUncUp2;
    //Float_t m4SIP2D1;
    //Float_t m4SIP2D2;
    //Float_t m4SIP3D1;
    //Float_t m4SIP3D2;
    //Float_t m4TrackIso1;
    //Float_t m4TrackIso2;
    //Float_t m4ZZIso1;
    //Float_t m4ZZIso2;
    //Bool_t m4HighPtID1;
    //Bool_t m4HighPtID2;
    //Bool_t m4HighPtIDNoVtx1;
    //Bool_t m4HighPtIDNoVtx2;
    //Bool_t m4IsGlobal1;
    //Bool_t m4IsGlobal2;
    //Bool_t m4IsLoose1;
    //Bool_t m4IsLoose2;
    //Bool_t m4IsMedium1;
    //Bool_t m4IsMedium2;
    //Bool_t m4IsPFMuon1;
    //Bool_t m4IsPFMuon2;
    //Bool_t m4IsTracker1;
    //Bool_t m4IsTracker2;
    //Bool_t m4PASTightID1;
    //Bool_t m4PASTightID2;
    //Bool_t m4PASTightIDNoVtx1;
    //Bool_t m4PASTightIDNoVtx2;
    //Bool_t m4PFID1;
    //Bool_t m4PFID2;
    //Bool_t m4PFIDNoVtx1;
    //Bool_t m4PFIDNoVtx2;
    //Bool_t m4ZZIsoPass1;
    //Bool_t m4ZZIsoPass2;
    //Bool_t m4ZZLooseID1;
    //Bool_t m4ZZLooseID2;
    //Bool_t m4ZZLooseIDNoVtx1;
    //Bool_t m4ZZLooseIDNoVtx2;
    //Bool_t m4ZZTightID1;
    //Bool_t m4ZZTightID2;
    //Bool_t m4ZZTightIDNoVtx1;
    //Bool_t m4ZZTightIDNoVtx2;
    //Int_t m4Charge1;
    //Int_t m4Charge2;
    //Int_t m4GenCharge1;
    //Int_t m4GenCharge2;
    //Int_t m4GenPdgId1;
    //Int_t m4GenPdgId2;
    Int_t m4PdgId1;
    Int_t m4PdgId2;
    //UInt_t m4BestTrackType1;
    //UInt_t m4BestTrackType2;
    //UInt_t m4CutBasedLoose1;
    //UInt_t m4CutBasedLoose2;
    //UInt_t m4CutBasedMedium1;
    //UInt_t m4CutBasedMedium2;
    //UInt_t m4CutBasedTight1;
    //UInt_t m4CutBasedTight2;
    //UInt_t m4IsTight1;
    //UInt_t m4IsTight2;
    //UInt_t m4MatchedStations1;
    //UInt_t m4MatchedStations2;
    //UInt_t m4NoOfMatches1;
    //UInt_t m4NoOfMatches2;
    //UInt_t m4PFIsoLoose1;
    //UInt_t m4PFIsoLoose2;
    //UInt_t m4PFIsoMedium1;
    //UInt_t m4PFIsoMedium2;
    //UInt_t m4PFIsoTight1;
    //UInt_t m4PFIsoTight2;
    //UInt_t m4PFIsoVTight1;
    //UInt_t m4PFIsoVTight2;

    //eeee or eemm
    //Float_t e1_e2_DR1;
    //Float_t e1_e2_DR2;
    //Float_t e1_e2_Energy1;
    //Float_t e1_e2_Energy2;
    //Float_t e1_e2_EnergyNoFSR1;
    //Float_t e1_e2_EnergyNoFSR2;
    //Float_t e1_e2_Eta1;
    //Float_t e1_e2_Eta2;
    //Float_t e1_e2_EtaNoFSR1;
    //Float_t e1_e2_EtaNoFSR2;
    Float_t e1_e2_Mass1;
    Float_t e1_e2_Mass2;
    //Float_t e1_e2_MassNoFSR1;
    //Float_t e1_e2_MassNoFSR2;
    //Float_t e1_e2_Mt1;
    //Float_t e1_e2_Mt2;
    //Float_t e1_e2_MtToMET1;
    //Float_t e1_e2_MtToMET2;
    //Float_t e1_e2_Phi1;
    //Float_t e1_e2_Phi2;
    //Float_t e1_e2_PhiNoFSR1;
    //Float_t e1_e2_PhiNoFSR2;
    //Float_t e1_e2_Pt1;
    //Float_t e1_e2_Pt2;
    //Float_t e1_e2_PtNoFSR1;
    //Float_t e1_e2_PtNoFSR2;
    //Bool_t e1_e2_SS1;
    //Bool_t e1_e2_SS2;
    //Int_t e1_e2_Charge1;
    //Int_t e1_e2_Charge2;
    //Int_t e1_e2_PdgId1;
    //Int_t e1_e2_PdgId2;
    //Float_t e1EffectiveArea1;
    //Float_t e1EffectiveArea2;
    Float_t e1Energy1;
    Float_t e1Energy2;
    Float_t e1Eta1;
    Float_t e1Eta2;
    //Float_t e1FSRDREt1;
    //Float_t e1FSRDREt2;
    //Float_t e1FSREt1;
    //Float_t e1FSREt2;
    //Float_t e1FSREta1;
    //Float_t e1FSREta2;
    //Float_t e1FSRPhi1;
    //Float_t e1FSRPhi2;
    //Float_t e1GenEta1;
    //Float_t e1GenEta2;
    //Float_t e1GenPhi1;
    //Float_t e1GenPhi2;
    //Float_t e1GenPt1;
    //Float_t e1GenPt2;
    //Float_t e1IP2D1;
    //Float_t e1IP2D2;
    //Float_t e1IP2DUncertainty1;
    //Float_t e1IP2DUncertainty2;
    //Float_t e1IP3D1;
    //Float_t e1IP3D2;
    //Float_t e1IP3DUncertainty1;
    //Float_t e1IP3DUncertainty2;
    //Float_t e1Loose1;
    //Float_t e1Loose2;
    //Float_t e1Mass1;
    //Float_t e1Mass2;
    //Float_t e1Medium1;
    //Float_t e1Medium2;
    //Float_t e1MtToMET1;
    //Float_t e1MtToMET2;
    //Float_t e1PFChargedIso1;
    //Float_t e1PFChargedIso2;
    //Float_t e1PFNeutralIso1;
    //Float_t e1PFNeutralIso2;
    //Float_t e1PFPUIso1;
    //Float_t e1PFPUIso2;
    //Float_t e1PFPhotonIso1;
    //Float_t e1PFPhotonIso2;
    //Float_t e1PVDXY1;
    //Float_t e1PVDXY2;
    //Float_t e1PVDZ1;
    //Float_t e1PVDZ2;
    Float_t e1Phi1;
    Float_t e1Phi2;
    Float_t e1Pt1;
    Float_t e1Pt2;
    //Float_t e1RelPFIsoRho1;
    //Float_t e1RelPFIsoRho2;
    //Float_t e1Rho1;
    //Float_t e1Rho2;
    //Float_t e1SCEnergy1;
    //Float_t e1SCEnergy2;
    //Float_t e1SCEta1;
    //Float_t e1SCEta2;
    //Float_t e1SCPhi1;
    //Float_t e1SCPhi2;
    //Float_t e1SCRawEnergy1;
    //Float_t e1SCRawEnergy2;
    //Float_t e1SIP2D1;
    //Float_t e1SIP2D2;
    //Float_t e1SIP3D1;
    //Float_t e1SIP3D2;
    //Float_t e1ScaleGainDn1;
    //Float_t e1ScaleGainDn2;
    //Float_t e1ScaleGainUp1;
    //Float_t e1ScaleGainUp2;
    //Float_t e1ScaleStatDn1;
    //Float_t e1ScaleStatDn2;
    //Float_t e1ScaleStatUp1;
    //Float_t e1ScaleStatUp2;
    //Float_t e1ScaleSystDn1;
    //Float_t e1ScaleSystDn2;
    //Float_t e1ScaleSystUp1;
    //Float_t e1ScaleSystUp2;
    //Float_t e1ScaleTotDn1;
    //Float_t e1ScaleTotDn2;
    //Float_t e1ScaleTotUp1;
    //Float_t e1ScaleTotUp2;
    //Float_t e1SigmaPhiDn1;
    //Float_t e1SigmaPhiDn2;
    //Float_t e1SigmaPhiUp1;
    //Float_t e1SigmaPhiUp2;
    //Float_t e1SigmaRhoDn1;
    //Float_t e1SigmaRhoDn2;
    //Float_t e1SigmaRhoUp1;
    //Float_t e1SigmaRhoUp2;
    //Float_t e1SigmaTotDn1;
    //Float_t e1SigmaTotDn2;
    //Float_t e1SigmaTotUp1;
    //Float_t e1SigmaTotUp2;
    //Float_t e1Tight1;
    //Float_t e1Tight2;
    //Float_t e1UnCorrPt1;
    //Float_t e1UnCorrPt2;
    //Float_t e1Veto1;
    //Float_t e1Veto2;
    //Float_t e1ZZIso1;
    //Float_t e1ZZIso2;
    //Bool_t e1IsEB1;
    //Bool_t e1IsEB2;
    //Bool_t e1IsGap1;
    //Bool_t e1IsGap2;
    //Bool_t e1ZZIsoPass1;
    //Bool_t e1ZZIsoPass2;
    //Bool_t e1ZZLooseID1;
    //Bool_t e1ZZLooseID2;
    //Bool_t e1ZZLooseIDNoVtx1;
    //Bool_t e1ZZLooseIDNoVtx2;
    //Bool_t e1ZZTightID1;
    //Bool_t e1ZZTightID2;
    //Bool_t e1ZZTightIDNoVtx1;
    //Bool_t e1ZZTightIDNoVtx2;
    //Int_t e1Charge1;
    //Int_t e1Charge2;
    //Int_t e1GenCharge1;
    //Int_t e1GenCharge2;
    //Int_t e1GenPdgId1;
    //Int_t e1GenPdgId2;
    Int_t e1PdgId1;
    Int_t e1PdgId2;
    //UInt_t e1MissingHits1;
    //UInt_t e1MissingHits2;
    //Float_t e2EffectiveArea1;
    //Float_t e2EffectiveArea2;
    Float_t e2Energy1;
    Float_t e2Energy2;
    Float_t e2Eta1;
    Float_t e2Eta2;
    //Float_t e2FSRDREt1;
    //Float_t e2FSRDREt2;
    //Float_t e2FSREt1;
    //Float_t e2FSREt2;
    //Float_t e2FSREta1;
    //Float_t e2FSREta2;
    //Float_t e2FSRPhi1;
    //Float_t e2FSRPhi2;
    //Float_t e2GenEta1;
    //Float_t e2GenEta2;
    //Float_t e2GenPhi1;
    //Float_t e2GenPhi2;
    //Float_t e2GenPt1;
    //Float_t e2GenPt2;
    //Float_t e2IP2D1;
    //Float_t e2IP2D2;
    //Float_t e2IP2DUncertainty1;
    //Float_t e2IP2DUncertainty2;
    //Float_t e2IP3D1;
    //Float_t e2IP3D2;
    //Float_t e2IP3DUncertainty1;
    //Float_t e2IP3DUncertainty2;
    //Float_t e2Loose1;
    //Float_t e2Loose2;
    //Float_t e2Mass1;
    //Float_t e2Mass2;
    //Float_t e2Medium1;
    //Float_t e2Medium2;
    //Float_t e2MtToMET1;
    //Float_t e2MtToMET2;
    //Float_t e2PFChargedIso1;
    //Float_t e2PFChargedIso2;
    //Float_t e2PFNeutralIso1;
    //Float_t e2PFNeutralIso2;
    //Float_t e2PFPUIso1;
    //Float_t e2PFPUIso2;
    //Float_t e2PFPhotonIso1;
    //Float_t e2PFPhotonIso2;
    //Float_t e2PVDXY1;
    //Float_t e2PVDXY2;
    //Float_t e2PVDZ1;
    //Float_t e2PVDZ2;
    Float_t e2Phi1;
    Float_t e2Phi2;
    Float_t e2Pt1;
    Float_t e2Pt2;
    //Float_t e2RelPFIsoRho1;
    //Float_t e2RelPFIsoRho2;
    //Float_t e2Rho1;
    //Float_t e2Rho2;
    //Float_t e2SCEnergy1;
    //Float_t e2SCEnergy2;
    //Float_t e2SCEta1;
    //Float_t e2SCEta2;
    //Float_t e2SCPhi1;
    //Float_t e2SCPhi2;
    //Float_t e2SCRawEnergy1;
    //Float_t e2SCRawEnergy2;
    //Float_t e2SIP2D1;
    //Float_t e2SIP2D2;
    //Float_t e2SIP3D1;
    //Float_t e2SIP3D2;
    //Float_t e2ScaleGainDn1;
    //Float_t e2ScaleGainDn2;
    //Float_t e2ScaleGainUp1;
    //Float_t e2ScaleGainUp2;
    //Float_t e2ScaleStatDn1;
    //Float_t e2ScaleStatDn2;
    //Float_t e2ScaleStatUp1;
    //Float_t e2ScaleStatUp2;
    //Float_t e2ScaleSystDn1;
    //Float_t e2ScaleSystDn2;
    //Float_t e2ScaleSystUp1;
    //Float_t e2ScaleSystUp2;
    //Float_t e2ScaleTotDn1;
    //Float_t e2ScaleTotDn2;
    //Float_t e2ScaleTotUp1;
    //Float_t e2ScaleTotUp2;
    //Float_t e2SigmaPhiDn1;
    //Float_t e2SigmaPhiDn2;
    //Float_t e2SigmaPhiUp1;
    //Float_t e2SigmaPhiUp2;
    //Float_t e2SigmaRhoDn1;
    //Float_t e2SigmaRhoDn2;
    //Float_t e2SigmaRhoUp1;
    //Float_t e2SigmaRhoUp2;
    //Float_t e2SigmaTotDn1;
    //Float_t e2SigmaTotDn2;
    //Float_t e2SigmaTotUp1;
    //Float_t e2SigmaTotUp2;
    //Float_t e2Tight1;
    //Float_t e2Tight2;
    //Float_t e2UnCorrPt1;
    //Float_t e2UnCorrPt2;
    //Float_t e2Veto1;
    //Float_t e2Veto2;
    //Float_t e2ZZIso1;
    //Float_t e2ZZIso2;
    //Bool_t e2IsEB1;
    //Bool_t e2IsEB2;
    //Bool_t e2IsGap1;
    //Bool_t e2IsGap2;
    //Bool_t e2ZZIsoPass1;
    //Bool_t e2ZZIsoPass2;
    //Bool_t e2ZZLooseID1;
    //Bool_t e2ZZLooseID2;
    //Bool_t e2ZZLooseIDNoVtx1;
    //Bool_t e2ZZLooseIDNoVtx2;
    //Bool_t e2ZZTightID1;
    //Bool_t e2ZZTightID2;
    //Bool_t e2ZZTightIDNoVtx1;
    //Bool_t e2ZZTightIDNoVtx2;
    //Int_t e2Charge1;
    //Int_t e2Charge2;
    //Int_t e2GenCharge1;
    //Int_t e2GenCharge2;
    //Int_t e2GenPdgId1;
    //Int_t e2GenPdgId2;
    Int_t e2PdgId1;
    Int_t e2PdgId2;
    //UInt_t e2MissingHits1;
    //UInt_t e2MissingHits2;

    //mmmm or eemm
    //Float_t m1_m2_DR1;
    //Float_t m1_m2_DR2;
    //Float_t m1_m2_Energy1;
    //Float_t m1_m2_Energy2;
    //Float_t m1_m2_EnergyNoFSR1;
    //Float_t m1_m2_EnergyNoFSR2;
    //Float_t m1_m2_Eta1;
    //Float_t m1_m2_Eta2;
    //Float_t m1_m2_EtaNoFSR1;
    //Float_t m1_m2_EtaNoFSR2;
    Float_t m1_m2_Mass1;
    Float_t m1_m2_Mass2;
    //Float_t m1_m2_MassNoFSR1;
    //Float_t m1_m2_MassNoFSR2;
    //Float_t m1_m2_Mt1;
    //Float_t m1_m2_Mt2;
    //Float_t m1_m2_MtToMET1;
    //Float_t m1_m2_MtToMET2;
    //Float_t m1_m2_Phi1;
    //Float_t m1_m2_Phi2;
    //Float_t m1_m2_PhiNoFSR1;
    //Float_t m1_m2_PhiNoFSR2;
    //Float_t m1_m2_Pt1;
    //Float_t m1_m2_Pt2;
    //Float_t m1_m2_PtNoFSR1;
    //Float_t m1_m2_PtNoFSR2;
    //Bool_t m1_m2_SS1;
    //Bool_t m1_m2_SS2;
    //Int_t m1_m2_Charge1;
    //Int_t m1_m2_Charge2;
    //Int_t m1_m2_PdgId1;
    //Int_t m1_m2_PdgId2;
    //Float_t m1EffScaleFactor1;
    //Float_t m1EffScaleFactor2;
    //Float_t m1EffScaleFactorError1;
    //Float_t m1EffScaleFactorError2;
    Float_t m1Energy1;
    Float_t m1Energy2;
    Float_t m1Eta1;
    Float_t m1Eta2;
    //Float_t m1FSRDREt1;
    //Float_t m1FSRDREt2;
    //Float_t m1FSREt1;
    //Float_t m1FSREt2;
    //Float_t m1FSREta1;
    //Float_t m1FSREta2;
    //Float_t m1FSRPhi1;
    //Float_t m1FSRPhi2;
    //Float_t m1GenEta1;
    //Float_t m1GenEta2;
    //Float_t m1GenPhi1;
    //Float_t m1GenPhi2;
    //Float_t m1GenPt1;
    //Float_t m1GenPt2;
    //Float_t m1IP2D1;
    //Float_t m1IP2D2;
    //Float_t m1IP2DUncertainty1;
    //Float_t m1IP2DUncertainty2;
    //Float_t m1IP3D1;
    //Float_t m1IP3D2;
    //Float_t m1IP3DUncertainty1;
    //Float_t m1IP3DUncertainty2;
    //Float_t m1Mass1;
    //Float_t m1Mass2;
    //Float_t m1MtToMET1;
    //Float_t m1MtToMET2;
    //Float_t m1PFChargedIso1;
    //Float_t m1PFChargedIso2;
    //Float_t m1PFNeutralIso1;
    //Float_t m1PFNeutralIso2;
    //Float_t m1PFPUIso1;
    //Float_t m1PFPUIso2;
    //Float_t m1PFPhotonIso1;
    //Float_t m1PFPhotonIso2;
    //Float_t m1PVDXY1;
    //Float_t m1PVDXY2;
    //Float_t m1PVDZ1;
    //Float_t m1PVDZ2;
    Float_t m1Phi1;
    Float_t m1Phi2;
    Float_t m1Pt1;
    Float_t m1Pt2;
    //Float_t m1PtErr1;
    //Float_t m1PtErr2;
    //Float_t m1PtErrUncorrected1;
    //Float_t m1PtErrUncorrected2;
    //Float_t m1PtUncorrected1;
    //Float_t m1PtUncorrected2;
    //Float_t m1RelPFIsoDBR031;
    //Float_t m1RelPFIsoDBR032;
    //Float_t m1RelPFIsoDBR041;
    //Float_t m1RelPFIsoDBR042;
    //Float_t m1RochesterScaleUncDn1;
    //Float_t m1RochesterScaleUncDn2;
    //Float_t m1RochesterScaleUncUp1;
    //Float_t m1RochesterScaleUncUp2;
    //Float_t m1RochesterSmearUncDn1;
    //Float_t m1RochesterSmearUncDn2;
    //Float_t m1RochesterSmearUncUp1;
    //Float_t m1RochesterSmearUncUp2;
    //Float_t m1SIP2D1;
    //Float_t m1SIP2D2;
    //Float_t m1SIP3D1;
    //Float_t m1SIP3D2;
    //Float_t m1TrackIso1;
    //Float_t m1TrackIso2;
    //Float_t m1ZZIso1;
    //Float_t m1ZZIso2;
    //Bool_t m1HighPtID1;
    //Bool_t m1HighPtID2;
    //Bool_t m1HighPtIDNoVtx1;
    //Bool_t m1HighPtIDNoVtx2;
    //Bool_t m1IsGlobal1;
    //Bool_t m1IsGlobal2;
    //Bool_t m1IsLoose1;
    //Bool_t m1IsLoose2;
    //Bool_t m1IsMedium1;
    //Bool_t m1IsMedium2;
    //Bool_t m1IsPFMuon1;
    //Bool_t m1IsPFMuon2;
    //Bool_t m1IsTracker1;
    //Bool_t m1IsTracker2;
    //Bool_t m1PASTightID1;
    //Bool_t m1PASTightID2;
    //Bool_t m1PASTightIDNoVtx1;
    //Bool_t m1PASTightIDNoVtx2;
    //Bool_t m1PFID1;
    //Bool_t m1PFID2;
    //Bool_t m1PFIDNoVtx1;
    //Bool_t m1PFIDNoVtx2;
    //Bool_t m1ZZIsoPass1;
    //Bool_t m1ZZIsoPass2;
    //Bool_t m1ZZLooseID1;
    //Bool_t m1ZZLooseID2;
    //Bool_t m1ZZLooseIDNoVtx1;
    //Bool_t m1ZZLooseIDNoVtx2;
    //Bool_t m1ZZTightID1;
    //Bool_t m1ZZTightID2;
    //Bool_t m1ZZTightIDNoVtx1;
    //Bool_t m1ZZTightIDNoVtx2;
    //Int_t m1Charge1;
    //Int_t m1Charge2;
    //Int_t m1GenCharge1;
    //Int_t m1GenCharge2;
    //Int_t m1GenPdgId1;
    //Int_t m1GenPdgId2;
    Int_t m1PdgId1;
    Int_t m1PdgId2;
    //UInt_t m1BestTrackType1;
    //UInt_t m1BestTrackType2;
    //UInt_t m1CutBasedLoose1;
    //UInt_t m1CutBasedLoose2;
    //UInt_t m1CutBasedMedium1;
    //UInt_t m1CutBasedMedium2;
    //UInt_t m1CutBasedTight1;
    //UInt_t m1CutBasedTight2;
    //UInt_t m1IsTight1;
    //UInt_t m1IsTight2;
    //UInt_t m1MatchedStations1;
    //UInt_t m1MatchedStations2;
    //UInt_t m1NoOfMatches1;
    //UInt_t m1NoOfMatches2;
    //UInt_t m1PFIsoLoose1;
    //UInt_t m1PFIsoLoose2;
    //UInt_t m1PFIsoMedium1;
    //UInt_t m1PFIsoMedium2;
    //UInt_t m1PFIsoTight1;
    //UInt_t m1PFIsoTight2;
    //UInt_t m1PFIsoVTight1;
    //UInt_t m1PFIsoVTight2;
    //Float_t m2EffScaleFactor1;
    //Float_t m2EffScaleFactor2;
    //Float_t m2EffScaleFactorError1;
    //Float_t m2EffScaleFactorError2;
    Float_t m2Energy1;
    Float_t m2Energy2;
    Float_t m2Eta1;
    Float_t m2Eta2;
    //Float_t m2FSRDREt1;
    //Float_t m2FSRDREt2;
    //Float_t m2FSREt1;
    //Float_t m2FSREt2;
    //Float_t m2FSREta1;
    //Float_t m2FSREta2;
    //Float_t m2FSRPhi1;
    //Float_t m2FSRPhi2;
    //Float_t m2GenEta1;
    //Float_t m2GenEta2;
    //Float_t m2GenPhi1;
    //Float_t m2GenPhi2;
    //Float_t m2GenPt1;
    //Float_t m2GenPt2;
    //Float_t m2IP2D1;
    //Float_t m2IP2D2;
    //Float_t m2IP2DUncertainty1;
    //Float_t m2IP2DUncertainty2;
    //Float_t m2IP3D1;
    //Float_t m2IP3D2;
    //Float_t m2IP3DUncertainty1;
    //Float_t m2IP3DUncertainty2;
    //Float_t m2Mass1;
    //Float_t m2Mass2;
    //Float_t m2MtToMET1;
    //Float_t m2MtToMET2;
    //Float_t m2PFChargedIso1;
    //Float_t m2PFChargedIso2;
    //Float_t m2PFNeutralIso1;
    //Float_t m2PFNeutralIso2;
    //Float_t m2PFPUIso1;
    //Float_t m2PFPUIso2;
    //Float_t m2PFPhotonIso1;
    //Float_t m2PFPhotonIso2;
    //Float_t m2PVDXY1;
    //Float_t m2PVDXY2;
    //Float_t m2PVDZ1;
    //Float_t m2PVDZ2;
    Float_t m2Phi1;
    Float_t m2Phi2;
    Float_t m2Pt1;
    Float_t m2Pt2;
    //Float_t m2PtErr1;
    //Float_t m2PtErr2;
    //Float_t m2PtErrUncorrected1;
    //Float_t m2PtErrUncorrected2;
    //Float_t m2PtUncorrected1;
    //Float_t m2PtUncorrected2;
    //Float_t m2RelPFIsoDBR031;
    //Float_t m2RelPFIsoDBR032;
    //Float_t m2RelPFIsoDBR041;
    //Float_t m2RelPFIsoDBR042;
    //Float_t m2RochesterScaleUncDn1;
    //Float_t m2RochesterScaleUncDn2;
    //Float_t m2RochesterScaleUncUp1;
    //Float_t m2RochesterScaleUncUp2;
    //Float_t m2RochesterSmearUncDn1;
    //Float_t m2RochesterSmearUncDn2;
    //Float_t m2RochesterSmearUncUp1;
    //Float_t m2RochesterSmearUncUp2;
    //Float_t m2SIP2D1;
    //Float_t m2SIP2D2;
    //Float_t m2SIP3D1;
    //Float_t m2SIP3D2;
    //Float_t m2TrackIso1;
    //Float_t m2TrackIso2;
    //Float_t m2ZZIso1;
    //Float_t m2ZZIso2;
    //Bool_t m2HighPtID1;
    //Bool_t m2HighPtID2;
    //Bool_t m2HighPtIDNoVtx1;
    //Bool_t m2HighPtIDNoVtx2;
    //Bool_t m2IsGlobal1;
    //Bool_t m2IsGlobal2;
    //Bool_t m2IsLoose1;
    //Bool_t m2IsLoose2;
    //Bool_t m2IsMedium1;
    //Bool_t m2IsMedium2;
    //Bool_t m2IsPFMuon1;
    //Bool_t m2IsPFMuon2;
    //Bool_t m2IsTracker1;
    //Bool_t m2IsTracker2;
    //Bool_t m2PASTightID1;
    //Bool_t m2PASTightID2;
    //Bool_t m2PASTightIDNoVtx1;
    //Bool_t m2PASTightIDNoVtx2;
    //Bool_t m2PFID1;
    //Bool_t m2PFID2;
    //Bool_t m2PFIDNoVtx1;
    //Bool_t m2PFIDNoVtx2;
    //Bool_t m2ZZIsoPass1;
    //Bool_t m2ZZIsoPass2;
    //Bool_t m2ZZLooseID1;
    //Bool_t m2ZZLooseID2;
    //Bool_t m2ZZLooseIDNoVtx1;
    //Bool_t m2ZZLooseIDNoVtx2;
    //Bool_t m2ZZTightID1;
    //Bool_t m2ZZTightID2;
    //Bool_t m2ZZTightIDNoVtx1;
    //Bool_t m2ZZTightIDNoVtx2;
    //Int_t m2Charge1;
    //Int_t m2Charge2;
    //Int_t m2GenCharge1;
    //Int_t m2GenCharge2;
    //Int_t m2GenPdgId1;
    //Int_t m2GenPdgId2;
    Int_t m2PdgId1;
    Int_t m2PdgId2;
    //UInt_t m2BestTrackType1;
    //UInt_t m2BestTrackType2;
    //UInt_t m2CutBasedLoose1;
    //UInt_t m2CutBasedLoose2;
    //UInt_t m2CutBasedMedium1;
    //UInt_t m2CutBasedMedium2;
    //UInt_t m2CutBasedTight1;
    //UInt_t m2CutBasedTight2;
    //UInt_t m2IsTight1;
    //UInt_t m2IsTight2;
    //UInt_t m2MatchedStations1;
    //UInt_t m2MatchedStations2;
    //UInt_t m2NoOfMatches1;
    //UInt_t m2NoOfMatches2;
    //UInt_t m2PFIsoLoose1;
    //UInt_t m2PFIsoLoose2;
    //UInt_t m2PFIsoMedium1;
    //UInt_t m2PFIsoMedium2;
    //UInt_t m2PFIsoTight1;
    //UInt_t m2PFIsoTight2;
    //UInt_t m2PFIsoVTight1;
    //UInt_t m2PFIsoVTight2;
};

CompLooperBase::CompLooperBase(const char *name, const char *channel, const char *filename1, const char *filename2)
  : _name(name), _channel(channel){
  _infile1 = TFile::Open(filename1);
  if (_infile1 != nullptr)
		_ntuple1 = (TTree*)_infile1->Get((_channel + "/ntuple").c_str());
  _infile2 = TFile::Open(filename2);
  if (_infile2 != nullptr)
		_ntuple2 = (TTree*)_infile2->Get((_channel + "/ntuple").c_str());
}

CompLooperBase::~CompLooperBase(){
  if (_infile1 != nullptr)
		_infile1->Close();
  if (_infile2 != nullptr)
		_infile2->Close();
}

void CompLooperBase::Init(){
  if (_hasInit) return;
  _hasInit = true;
  //Shared branches
  if (_isT1MC) _ntuple1->SetBranchAddress("summedWeights", &summedWeights1);
  if (_isT2MC) _ntuple2->SetBranchAddress("summedWeights", &summedWeights2);
  //_ntuple1->SetBranchAddress("Energy", &Energy1);
  //_ntuple2->SetBranchAddress("Energy", &Energy2);
  //_ntuple1->SetBranchAddress("EnergyNoFSR", &EnergyNoFSR1);
  //_ntuple2->SetBranchAddress("EnergyNoFSR", &EnergyNoFSR2);
  //_ntuple1->SetBranchAddress("Eta", &Eta1);
  //_ntuple2->SetBranchAddress("Eta", &Eta2);
  //_ntuple1->SetBranchAddress("EtaNoFSR", &EtaNoFSR1);
  //_ntuple2->SetBranchAddress("EtaNoFSR", &EtaNoFSR2);
  //_ntuple1->SetBranchAddress("Flag_BadChargedCandidateFilterPass", &Flag_BadChargedCandidateFilterPass1);
  //_ntuple2->SetBranchAddress("Flag_BadChargedCandidateFilterPass", &Flag_BadChargedCandidateFilterPass2);
  //_ntuple1->SetBranchAddress("Flag_BadPFMuonFilterPass", &Flag_BadPFMuonFilterPass1);
  //_ntuple2->SetBranchAddress("Flag_BadPFMuonFilterPass", &Flag_BadPFMuonFilterPass2);
  _ntuple1->SetBranchAddress("Mass", &Mass1);
  _ntuple2->SetBranchAddress("Mass", &Mass2);
  //_ntuple1->SetBranchAddress("MassNoFSR", &MassNoFSR1);
  //_ntuple2->SetBranchAddress("MassNoFSR", &MassNoFSR2);
  //_ntuple1->SetBranchAddress("MtToMET", &MtToMET1);
  //_ntuple2->SetBranchAddress("MtToMET", &MtToMET2);
  //_ntuple1->SetBranchAddress("Phi", &Phi1);
  //_ntuple2->SetBranchAddress("Phi", &Phi2);
  //_ntuple1->SetBranchAddress("PhiNoFSR", &PhiNoFSR1);
  //_ntuple2->SetBranchAddress("PhiNoFSR", &PhiNoFSR2);
  //_ntuple1->SetBranchAddress("Pt", &Pt1);
  //_ntuple2->SetBranchAddress("Pt", &Pt2);
  //_ntuple1->SetBranchAddress("PtNoFSR", &PtNoFSR1);
  //_ntuple2->SetBranchAddress("PtNoFSR", &PtNoFSR2);
  //_ntuple1->SetBranchAddress("deltaEtajj", &deltaEtajj1);
  //_ntuple2->SetBranchAddress("deltaEtajj", &deltaEtajj2);
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("deltaEtajj_jerDown", &deltaEtajj_jerDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("deltaEtajj_jerDown", &deltaEtajj_jerDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("deltaEtajj_jerUp", &deltaEtajj_jerUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("deltaEtajj_jerUp", &deltaEtajj_jerUp2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("deltaEtajj_jesDown", &deltaEtajj_jesDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("deltaEtajj_jesDown", &deltaEtajj_jesDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("deltaEtajj_jesUp", &deltaEtajj_jesUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("deltaEtajj_jesUp", &deltaEtajj_jesUp2);
  }
  //_ntuple1->SetBranchAddress("deltaPhiTojj", &deltaPhiTojj1);
  //_ntuple2->SetBranchAddress("deltaPhiTojj", &deltaPhiTojj2);
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("deltaPhiTojj_jerDown", &deltaPhiTojj_jerDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("deltaPhiTojj_jerDown", &deltaPhiTojj_jerDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("deltaPhiTojj_jerUp", &deltaPhiTojj_jerUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("deltaPhiTojj_jerUp", &deltaPhiTojj_jerUp2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("deltaPhiTojj_jesDown", &deltaPhiTojj_jesDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("deltaPhiTojj_jesDown", &deltaPhiTojj_jesDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("deltaPhiTojj_jesUp", &deltaPhiTojj_jesUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("deltaPhiTojj_jesUp", &deltaPhiTojj_jesUp2);
  }
  //_ntuple1->SetBranchAddress("etajj", &etajj1);
  //_ntuple2->SetBranchAddress("etajj", &etajj2);
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("etajj_jerDown", &etajj_jerDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("etajj_jerDown", &etajj_jerDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("etajj_jerUp", &etajj_jerUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("etajj_jerUp", &etajj_jerUp2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("etajj_jesDown", &etajj_jesDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("etajj_jesDown", &etajj_jesDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("etajj_jesUp", &etajj_jesUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("etajj_jesUp", &etajj_jesUp2);
  }
  if (_isT1MC){
    _ntuple1->SetBranchAddress("genWeight", &genWeight1);
  }
  if (_isT2MC){
    _ntuple2->SetBranchAddress("genWeight", &genWeight2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("jetPUSFmulfac", &jetPUSFmulfac1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("jetPUSFmulfac", &jetPUSFmulfac2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("maxScaleWeight", &maxScaleWeight1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("maxScaleWeight", &maxScaleWeight2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("minScaleWeight", &minScaleWeight1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("minScaleWeight", &minScaleWeight2);
  }
  //_ntuple1->SetBranchAddress("mjj", &mjj1);
  //_ntuple2->SetBranchAddress("mjj", &mjj2);
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("mjj_jerDown", &mjj_jerDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("mjj_jerDown", &mjj_jerDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("mjj_jerUp", &mjj_jerUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("mjj_jerUp", &mjj_jerUp2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("mjj_jesDown", &mjj_jesDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("mjj_jesDown", &mjj_jesDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("mjj_jesUp", &mjj_jesUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("mjj_jesUp", &mjj_jesUp2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("nTruePU", &nTruePU1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("nTruePU", &nTruePU2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("originalXWGTUP", &originalXWGTUP1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("originalXWGTUP", &originalXWGTUP2);
  }
  //_ntuple1->SetBranchAddress("phijj", &phijj1);
  //_ntuple2->SetBranchAddress("phijj", &phijj2);
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("phijj_jerDown", &phijj_jerDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("phijj_jerDown", &phijj_jerDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("phijj_jerUp", &phijj_jerUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("phijj_jerUp", &phijj_jerUp2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("phijj_jesDown", &phijj_jesDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("phijj_jesDown", &phijj_jesDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("phijj_jesUp", &phijj_jesUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("phijj_jesUp", &phijj_jesUp2);
  }
  //_ntuple1->SetBranchAddress("ptjj", &ptjj1);
  //_ntuple2->SetBranchAddress("ptjj", &ptjj2);
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("ptjj_jerDown", &ptjj_jerDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("ptjj_jerDown", &ptjj_jerDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("ptjj_jerUp", &ptjj_jerUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("ptjj_jerUp", &ptjj_jerUp2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("ptjj_jesDown", &ptjj_jesDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("ptjj_jesDown", &ptjj_jesDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("ptjj_jesUp", &ptjj_jesUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("ptjj_jesUp", &ptjj_jesUp2);
  }
  //_ntuple1->SetBranchAddress("pvRho", &pvRho1);
  //_ntuple2->SetBranchAddress("pvRho", &pvRho2);
  //_ntuple1->SetBranchAddress("pvZ", &pvZ1);
  //_ntuple2->SetBranchAddress("pvZ", &pvZ2);
  //_ntuple1->SetBranchAddress("pvndof", &pvndof1);
  //_ntuple2->SetBranchAddress("pvndof", &pvndof2);
  //_ntuple1->SetBranchAddress("rapidity", &rapidity1);
  //_ntuple2->SetBranchAddress("rapidity", &rapidity2);
  //_ntuple1->SetBranchAddress("type1_pfMETEt", &type1_pfMETEt1);
  //_ntuple2->SetBranchAddress("type1_pfMETEt", &type1_pfMETEt2);
  //_ntuple1->SetBranchAddress("type1_pfMETPhi", &type1_pfMETPhi1);
  //_ntuple2->SetBranchAddress("type1_pfMETPhi", &type1_pfMETPhi2);
  //_ntuple1->SetBranchAddress("zeppenfeld", &zeppenfeld1);
  //_ntuple2->SetBranchAddress("zeppenfeld", &zeppenfeld2);
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("zeppenfeld_jerDown", &zeppenfeld_jerDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("zeppenfeld_jerDown", &zeppenfeld_jerDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("zeppenfeld_jerUp", &zeppenfeld_jerUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("zeppenfeld_jerUp", &zeppenfeld_jerUp2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("zeppenfeld_jesDown", &zeppenfeld_jesDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("zeppenfeld_jesDown", &zeppenfeld_jesDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("zeppenfeld_jesUp", &zeppenfeld_jesUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("zeppenfeld_jesUp", &zeppenfeld_jesUp2);
  }
  //_ntuple1->SetBranchAddress("zeppenfeldj3", &zeppenfeldj31);
  //_ntuple2->SetBranchAddress("zeppenfeldj3", &zeppenfeldj32);
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("zeppenfeldj3_jerDown", &zeppenfeldj3_jerDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("zeppenfeldj3_jerDown", &zeppenfeldj3_jerDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("zeppenfeldj3_jerUp", &zeppenfeldj3_jerUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("zeppenfeldj3_jerUp", &zeppenfeldj3_jerUp2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("zeppenfeldj3_jesDown", &zeppenfeldj3_jesDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("zeppenfeldj3_jesDown", &zeppenfeldj3_jesDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("zeppenfeldj3_jesUp", &zeppenfeldj3_jesUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("zeppenfeldj3_jesUp", &zeppenfeldj3_jesUp2);
  }
  //_ntuple1->SetBranchAddress("pvIdFake", &pvIdFake1);
  //_ntuple2->SetBranchAddress("pvIdFake", &pvIdFake2);
  //_ntuple1->SetBranchAddress("pvIsValid", &pvIsValid1);
  //_ntuple2->SetBranchAddress("pvIsValid", &pvIsValid2);
  //_ntuple1->SetBranchAddress("Charge", &Charge1);
  //_ntuple2->SetBranchAddress("Charge", &Charge2);
  //_ntuple1->SetBranchAddress("PdgId", &PdgId1);
  //_ntuple2->SetBranchAddress("PdgId", &PdgId2);
  _ntuple1->SetBranchAddress("lumi", &lumi1);
  _ntuple2->SetBranchAddress("lumi", &lumi2);
  //_ntuple1->SetBranchAddress("nJets", &nJets1);
  //_ntuple2->SetBranchAddress("nJets", &nJets2);
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("nJets_jerDown", &nJets_jerDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("nJets_jerDown", &nJets_jerDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("nJets_jerUp", &nJets_jerUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("nJets_jerUp", &nJets_jerUp2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("nJets_jesDown", &nJets_jesDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("nJets_jesDown", &nJets_jesDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("nJets_jesUp", &nJets_jesUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("nJets_jesUp", &nJets_jesUp2);
  }
  //_ntuple1->SetBranchAddress("nZZLooseElec", &nZZLooseElec1);
  //_ntuple2->SetBranchAddress("nZZLooseElec", &nZZLooseElec2);
  //_ntuple1->SetBranchAddress("nZZLooseMu", &nZZLooseMu1);
  //_ntuple2->SetBranchAddress("nZZLooseMu", &nZZLooseMu2);
  //_ntuple1->SetBranchAddress("nZZTightElec", &nZZTightElec1);
  //_ntuple2->SetBranchAddress("nZZTightElec", &nZZTightElec2);
  //_ntuple1->SetBranchAddress("nZZTightIsoElec", &nZZTightIsoElec1);
  //_ntuple2->SetBranchAddress("nZZTightIsoElec", &nZZTightIsoElec2);
  //_ntuple1->SetBranchAddress("nZZTightIsoMu", &nZZTightIsoMu1);
  //_ntuple2->SetBranchAddress("nZZTightIsoMu", &nZZTightIsoMu2);
  //_ntuple1->SetBranchAddress("nZZTightMu", &nZZTightMu1);
  //_ntuple2->SetBranchAddress("nZZTightMu", &nZZTightMu2);
  _ntuple1->SetBranchAddress("nvtx", &nvtx1);
  _ntuple2->SetBranchAddress("nvtx", &nvtx2);
  _ntuple1->SetBranchAddress("run", &run1);
  _ntuple2->SetBranchAddress("run", &run2);
  _ntuple1->SetBranchAddress("evt", &evt1);
  _ntuple2->SetBranchAddress("evt", &evt2);
  //_ntuple1->SetBranchAddress("jetDeepCSV", &jetDeepCSV1);
  //_ntuple2->SetBranchAddress("jetDeepCSV", &jetDeepCSV2);
  //_ntuple1->SetBranchAddress("jetEta", &jetEta1);
  //_ntuple2->SetBranchAddress("jetEta", &jetEta2);
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("jetEta_jerDown", &jetEta_jerDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("jetEta_jerDown", &jetEta_jerDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("jetEta_jerUp", &jetEta_jerUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("jetEta_jerUp", &jetEta_jerUp2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("jetEta_jesDown", &jetEta_jesDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("jetEta_jesDown", &jetEta_jesDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("jetEta_jesUp", &jetEta_jesUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("jetEta_jesUp", &jetEta_jesUp2);
  }
  //_ntuple1->SetBranchAddress("jetPhi", &jetPhi1);
  //_ntuple2->SetBranchAddress("jetPhi", &jetPhi2);
  //_ntuple1->SetBranchAddress("jetPt", &jetPt1);
  //_ntuple2->SetBranchAddress("jetPt", &jetPt2);
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("jetPt_jerDown", &jetPt_jerDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("jetPt_jerDown", &jetPt_jerDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("jetPt_jerUp", &jetPt_jerUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("jetPt_jerUp", &jetPt_jerUp2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("jetPt_jesDown", &jetPt_jesDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("jetPt_jesDown", &jetPt_jesDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("jetPt_jesUp", &jetPt_jesUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("jetPt_jesUp", &jetPt_jesUp2);
  }
  //_ntuple1->SetBranchAddress("jetRapidity", &jetRapidity1);
  //_ntuple2->SetBranchAddress("jetRapidity", &jetRapidity2);
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("jetRapidity_jerDown", &jetRapidity_jerDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("jetRapidity_jerDown", &jetRapidity_jerDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("jetRapidity_jerUp", &jetRapidity_jerUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("jetRapidity_jerUp", &jetRapidity_jerUp2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("jetRapidity_jesDown", &jetRapidity_jesDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("jetRapidity_jesDown", &jetRapidity_jesDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("jetRapidity_jesUp", &jetRapidity_jesUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("jetRapidity_jesUp", &jetRapidity_jesUp2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("pdfWeights", &pdfWeights1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("pdfWeights", &pdfWeights2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("scaleWeightIDs", &scaleWeightIDs1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("scaleWeightIDs", &scaleWeightIDs2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("scaleWeights", &scaleWeights1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("scaleWeights", &scaleWeights2);
  }
  //_ntuple1->SetBranchAddress("isGenJetMatched", &isGenJetMatched1);
  //_ntuple2->SetBranchAddress("isGenJetMatched", &isGenJetMatched2);
  //_ntuple1->SetBranchAddress("jetHadronFlavor", &jetHadronFlavor1);
  //_ntuple2->SetBranchAddress("jetHadronFlavor", &jetHadronFlavor2);
  //_ntuple1->SetBranchAddress("jetPUID", &jetPUID1);
  //_ntuple2->SetBranchAddress("jetPUID", &jetPUID2);
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("jetPUID_jerDown", &jetPUID_jerDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("jetPUID_jerDown", &jetPUID_jerDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("jetPUID_jerUp", &jetPUID_jerUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("jetPUID_jerUp", &jetPUID_jerUp2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("jetPUID_jesDown", &jetPUID_jesDown1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("jetPUID_jesDown", &jetPUID_jesDown2);
  }
  if (_isT1MC){
    //_ntuple1->SetBranchAddress("jetPUID_jesUp", &jetPUID_jesUp1);
  }
  if (_isT2MC){
    //_ntuple2->SetBranchAddress("jetPUID_jesUp", &jetPUID_jesUp2);
  }
  //_ntuple1->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVLPass", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVLPass1);
  //_ntuple2->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVLPass", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVLPass2);
  //_ntuple1->SetBranchAddress("HLT_DoubleEle25_CaloIdL_MWPass", &HLT_DoubleEle25_CaloIdL_MWPass1);
  //_ntuple2->SetBranchAddress("HLT_DoubleEle25_CaloIdL_MWPass", &HLT_DoubleEle25_CaloIdL_MWPass2);
  //_ntuple1->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8Pass", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8Pass1);
  //_ntuple2->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8Pass", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8Pass2);
  //_ntuple1->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVLPass", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVLPass1);
  //_ntuple2->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVLPass", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVLPass2);
  //_ntuple1->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass1);
  //_ntuple2->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass2);
  //_ntuple1->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass1);
  //_ntuple2->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass2);
  //_ntuple1->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass1);
  //_ntuple2->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass2);
  //_ntuple1->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZPass", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZPass1);
  //_ntuple2->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZPass", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZPass2);
  //_ntuple1->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZPass", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZPass1);
  //_ntuple2->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZPass", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZPass2);
  //_ntuple1->SetBranchAddress("HLT_TripleMu_10_5_5_DZPass", &HLT_TripleMu_10_5_5_DZPass1);
  //_ntuple2->SetBranchAddress("HLT_TripleMu_10_5_5_DZPass", &HLT_TripleMu_10_5_5_DZPass2);
  //_ntuple1->SetBranchAddress("HLT_TripleMu_12_10_5Pass", &HLT_TripleMu_12_10_5Pass1);
  //_ntuple2->SetBranchAddress("HLT_TripleMu_12_10_5Pass", &HLT_TripleMu_12_10_5Pass2);
  //_ntuple1->SetBranchAddress("HLT_Ele30_WPTight_GsfPass", &HLT_Ele30_WPTight_GsfPass1);
  //_ntuple2->SetBranchAddress("HLT_Ele30_WPTight_GsfPass", &HLT_Ele30_WPTight_GsfPass2);
  //_ntuple1->SetBranchAddress("HLT_Ele32_WPTight_GsfPass", &HLT_Ele32_WPTight_GsfPass1);
  //_ntuple2->SetBranchAddress("HLT_Ele32_WPTight_GsfPass", &HLT_Ele32_WPTight_GsfPass2);
  //_ntuple1->SetBranchAddress("HLT_IsoMu24Pass", &HLT_IsoMu24Pass1);
  //_ntuple2->SetBranchAddress("HLT_IsoMu24Pass", &HLT_IsoMu24Pass2);
  //_ntuple1->SetBranchAddress("Flag_goodVerticesPass", &Flag_goodVerticesPass1);
  //_ntuple2->SetBranchAddress("Flag_goodVerticesPass", &Flag_goodVerticesPass2);
  //_ntuple1->SetBranchAddress("Flag_globalSuperTightHalo2016FilterPass", &Flag_globalSuperTightHalo2016FilterPass1);
  //_ntuple2->SetBranchAddress("Flag_globalSuperTightHalo2016FilterPass", &Flag_globalSuperTightHalo2016FilterPass2);
  //_ntuple1->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilterPass", &Flag_EcalDeadCellTriggerPrimitiveFilterPass1);
  //_ntuple2->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilterPass", &Flag_EcalDeadCellTriggerPrimitiveFilterPass2);
  //_ntuple1->SetBranchAddress("Flag_BadPFMuonFilterPass", &Flag_BadPFMuonFilterPass1);
  //_ntuple2->SetBranchAddress("Flag_BadPFMuonFilterPass", &Flag_BadPFMuonFilterPass2);
  //_ntuple1->SetBranchAddress("Flag_BadPFMuonDzFilterPass", &Flag_BadPFMuonDzFilterPass1);
  //_ntuple2->SetBranchAddress("Flag_BadPFMuonDzFilterPass", &Flag_BadPFMuonDzFilterPass2);
  //_ntuple1->SetBranchAddress("Flag_hfNoisyHitsFilterPass", &Flag_hfNoisyHitsFilterPass1);
  //_ntuple2->SetBranchAddress("Flag_hfNoisyHitsFilterPass", &Flag_hfNoisyHitsFilterPass2);
  //_ntuple1->SetBranchAddress("Flag_eeBadScFilterPass", &Flag_eeBadScFilterPass1);
  //_ntuple2->SetBranchAddress("Flag_eeBadScFilterPass", &Flag_eeBadScFilterPass2);

  if (_channel == "eeee"){
    //_ntuple1->SetBranchAddress("e1_e3_DR", &e1_e3_DR1);
    //_ntuple2->SetBranchAddress("e1_e3_DR", &e1_e3_DR2);
    //_ntuple1->SetBranchAddress("e1_e3_Mass", &e1_e3_Mass1);
    //_ntuple2->SetBranchAddress("e1_e3_Mass", &e1_e3_Mass2);
    //_ntuple1->SetBranchAddress("e1_e3_MassNoFSR", &e1_e3_MassNoFSR1);
    //_ntuple2->SetBranchAddress("e1_e3_MassNoFSR", &e1_e3_MassNoFSR2);
    //_ntuple1->SetBranchAddress("e1_e4_DR", &e1_e4_DR1);
    //_ntuple2->SetBranchAddress("e1_e4_DR", &e1_e4_DR2);
    //_ntuple1->SetBranchAddress("e1_e4_Mass", &e1_e4_Mass1);
    //_ntuple2->SetBranchAddress("e1_e4_Mass", &e1_e4_Mass2);
    //_ntuple1->SetBranchAddress("e1_e4_MassNoFSR", &e1_e4_MassNoFSR1);
    //_ntuple2->SetBranchAddress("e1_e4_MassNoFSR", &e1_e4_MassNoFSR2);
    //_ntuple1->SetBranchAddress("e2_e3_DR", &e2_e3_DR1);
    //_ntuple2->SetBranchAddress("e2_e3_DR", &e2_e3_DR2);
    //_ntuple1->SetBranchAddress("e2_e3_Mass", &e2_e3_Mass1);
    //_ntuple2->SetBranchAddress("e2_e3_Mass", &e2_e3_Mass2);
    //_ntuple1->SetBranchAddress("e2_e3_MassNoFSR", &e2_e3_MassNoFSR1);
    //_ntuple2->SetBranchAddress("e2_e3_MassNoFSR", &e2_e3_MassNoFSR2);
    //_ntuple1->SetBranchAddress("e2_e4_DR", &e2_e4_DR1);
    //_ntuple2->SetBranchAddress("e2_e4_DR", &e2_e4_DR2);
    //_ntuple1->SetBranchAddress("e2_e4_Mass", &e2_e4_Mass1);
    //_ntuple2->SetBranchAddress("e2_e4_Mass", &e2_e4_Mass2);
    //_ntuple1->SetBranchAddress("e2_e4_MassNoFSR", &e2_e4_MassNoFSR1);
    //_ntuple2->SetBranchAddress("e2_e4_MassNoFSR", &e2_e4_MassNoFSR2);
    //_ntuple1->SetBranchAddress("e1_e3_SS", &e1_e3_SS1);
    //_ntuple2->SetBranchAddress("e1_e3_SS", &e1_e3_SS2);
    //_ntuple1->SetBranchAddress("e1_e4_SS", &e1_e4_SS1);
    //_ntuple2->SetBranchAddress("e1_e4_SS", &e1_e4_SS2);
    //_ntuple1->SetBranchAddress("e2_e3_SS", &e2_e3_SS1);
    //_ntuple2->SetBranchAddress("e2_e3_SS", &e2_e3_SS2);
    //_ntuple1->SetBranchAddress("e2_e4_SS", &e2_e4_SS1);
    //_ntuple2->SetBranchAddress("e2_e4_SS", &e2_e4_SS2);
    //_ntuple1->SetBranchAddress("e3_e4_DR", &e3_e4_DR1);
    //_ntuple2->SetBranchAddress("e3_e4_DR", &e3_e4_DR2);
    //_ntuple1->SetBranchAddress("e3_e4_Energy", &e3_e4_Energy1);
    //_ntuple2->SetBranchAddress("e3_e4_Energy", &e3_e4_Energy2);
    //_ntuple1->SetBranchAddress("e3_e4_EnergyNoFSR", &e3_e4_EnergyNoFSR1);
    //_ntuple2->SetBranchAddress("e3_e4_EnergyNoFSR", &e3_e4_EnergyNoFSR2);
    //_ntuple1->SetBranchAddress("e3_e4_Eta", &e3_e4_Eta1);
    //_ntuple2->SetBranchAddress("e3_e4_Eta", &e3_e4_Eta2);
    //_ntuple1->SetBranchAddress("e3_e4_EtaNoFSR", &e3_e4_EtaNoFSR1);
    //_ntuple2->SetBranchAddress("e3_e4_EtaNoFSR", &e3_e4_EtaNoFSR2);
    _ntuple1->SetBranchAddress("e3_e4_Mass", &e3_e4_Mass1);
    _ntuple2->SetBranchAddress("e3_e4_Mass", &e3_e4_Mass2);
    //_ntuple1->SetBranchAddress("e3_e4_MassNoFSR", &e3_e4_MassNoFSR1);
    //_ntuple2->SetBranchAddress("e3_e4_MassNoFSR", &e3_e4_MassNoFSR2);
    //_ntuple1->SetBranchAddress("e3_e4_Mt", &e3_e4_Mt1);
    //_ntuple2->SetBranchAddress("e3_e4_Mt", &e3_e4_Mt2);
    //_ntuple1->SetBranchAddress("e3_e4_MtToMET", &e3_e4_MtToMET1);
    //_ntuple2->SetBranchAddress("e3_e4_MtToMET", &e3_e4_MtToMET2);
    //_ntuple1->SetBranchAddress("e3_e4_Phi", &e3_e4_Phi1);
    //_ntuple2->SetBranchAddress("e3_e4_Phi", &e3_e4_Phi2);
    //_ntuple1->SetBranchAddress("e3_e4_PhiNoFSR", &e3_e4_PhiNoFSR1);
    //_ntuple2->SetBranchAddress("e3_e4_PhiNoFSR", &e3_e4_PhiNoFSR2);
    //_ntuple1->SetBranchAddress("e3_e4_Pt", &e3_e4_Pt1);
    //_ntuple2->SetBranchAddress("e3_e4_Pt", &e3_e4_Pt2);
    //_ntuple1->SetBranchAddress("e3_e4_PtNoFSR", &e3_e4_PtNoFSR1);
    //_ntuple2->SetBranchAddress("e3_e4_PtNoFSR", &e3_e4_PtNoFSR2);
    //_ntuple1->SetBranchAddress("e3_e4_SS", &e3_e4_SS1);
    //_ntuple2->SetBranchAddress("e3_e4_SS", &e3_e4_SS2);
    //_ntuple1->SetBranchAddress("e3_e4_Charge", &e3_e4_Charge1);
    //_ntuple2->SetBranchAddress("e3_e4_Charge", &e3_e4_Charge2);
    //_ntuple1->SetBranchAddress("e3_e4_PdgId", &e3_e4_PdgId1);
    //_ntuple2->SetBranchAddress("e3_e4_PdgId", &e3_e4_PdgId2);
    //_ntuple1->SetBranchAddress("e3EffectiveArea", &e3EffectiveArea1);
    //_ntuple2->SetBranchAddress("e3EffectiveArea", &e3EffectiveArea2);
    _ntuple1->SetBranchAddress("e3Energy", &e3Energy1);
    _ntuple2->SetBranchAddress("e3Energy", &e3Energy2);
    _ntuple1->SetBranchAddress("e3Eta", &e3Eta1);
    _ntuple2->SetBranchAddress("e3Eta", &e3Eta2);
    //_ntuple1->SetBranchAddress("e3FSRDREt", &e3FSRDREt1);
    //_ntuple2->SetBranchAddress("e3FSRDREt", &e3FSRDREt2);
    //_ntuple1->SetBranchAddress("e3FSREt", &e3FSREt1);
    //_ntuple2->SetBranchAddress("e3FSREt", &e3FSREt2);
    //_ntuple1->SetBranchAddress("e3FSREta", &e3FSREta1);
    //_ntuple2->SetBranchAddress("e3FSREta", &e3FSREta2);
    //_ntuple1->SetBranchAddress("e3FSRPhi", &e3FSRPhi1);
    //_ntuple2->SetBranchAddress("e3FSRPhi", &e3FSRPhi2);
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e3GenEta", &e3GenEta1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e3GenEta", &e3GenEta2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e3GenPhi", &e3GenPhi1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e3GenPhi", &e3GenPhi2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e3GenPt", &e3GenPt1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e3GenPt", &e3GenPt2);
    }
    //_ntuple1->SetBranchAddress("e3IP2D", &e3IP2D1);
    //_ntuple2->SetBranchAddress("e3IP2D", &e3IP2D2);
    //_ntuple1->SetBranchAddress("e3IP2DUncertainty", &e3IP2DUncertainty1);
    //_ntuple2->SetBranchAddress("e3IP2DUncertainty", &e3IP2DUncertainty2);
    //_ntuple1->SetBranchAddress("e3IP3D", &e3IP3D1);
    //_ntuple2->SetBranchAddress("e3IP3D", &e3IP3D2);
    //_ntuple1->SetBranchAddress("e3IP3DUncertainty", &e3IP3DUncertainty1);
    //_ntuple2->SetBranchAddress("e3IP3DUncertainty", &e3IP3DUncertainty2);
    //_ntuple1->SetBranchAddress("e3Loose", &e3Loose1);
    //_ntuple2->SetBranchAddress("e3Loose", &e3Loose2);
    //_ntuple1->SetBranchAddress("e3Mass", &e3Mass1);
    //_ntuple2->SetBranchAddress("e3Mass", &e3Mass2);
    //_ntuple1->SetBranchAddress("e3Medium", &e3Medium1);
    //_ntuple2->SetBranchAddress("e3Medium", &e3Medium2);
    //_ntuple1->SetBranchAddress("e3MtToMET", &e3MtToMET1);
    //_ntuple2->SetBranchAddress("e3MtToMET", &e3MtToMET2);
    //_ntuple1->SetBranchAddress("e3PFChargedIso", &e3PFChargedIso1);
    //_ntuple2->SetBranchAddress("e3PFChargedIso", &e3PFChargedIso2);
    //_ntuple1->SetBranchAddress("e3PFNeutralIso", &e3PFNeutralIso1);
    //_ntuple2->SetBranchAddress("e3PFNeutralIso", &e3PFNeutralIso2);
    //_ntuple1->SetBranchAddress("e3PFPUIso", &e3PFPUIso1);
    //_ntuple2->SetBranchAddress("e3PFPUIso", &e3PFPUIso2);
    //_ntuple1->SetBranchAddress("e3PFPhotonIso", &e3PFPhotonIso1);
    //_ntuple2->SetBranchAddress("e3PFPhotonIso", &e3PFPhotonIso2);
    //_ntuple1->SetBranchAddress("e3PVDXY", &e3PVDXY1);
    //_ntuple2->SetBranchAddress("e3PVDXY", &e3PVDXY2);
    //_ntuple1->SetBranchAddress("e3PVDZ", &e3PVDZ1);
    //_ntuple2->SetBranchAddress("e3PVDZ", &e3PVDZ2);
    _ntuple1->SetBranchAddress("e3Phi", &e3Phi1);
    _ntuple2->SetBranchAddress("e3Phi", &e3Phi2);
    _ntuple1->SetBranchAddress("e3Pt", &e3Pt1);
    _ntuple2->SetBranchAddress("e3Pt", &e3Pt2);
    //_ntuple1->SetBranchAddress("e3RelPFIsoRho", &e3RelPFIsoRho1);
    //_ntuple2->SetBranchAddress("e3RelPFIsoRho", &e3RelPFIsoRho2);
    //_ntuple1->SetBranchAddress("e3Rho", &e3Rho1);
    //_ntuple2->SetBranchAddress("e3Rho", &e3Rho2);
    //_ntuple1->SetBranchAddress("e3SCEnergy", &e3SCEnergy1);
    //_ntuple2->SetBranchAddress("e3SCEnergy", &e3SCEnergy2);
    //_ntuple1->SetBranchAddress("e3SCEta", &e3SCEta1);
    //_ntuple2->SetBranchAddress("e3SCEta", &e3SCEta2);
    //_ntuple1->SetBranchAddress("e3SCPhi", &e3SCPhi1);
    //_ntuple2->SetBranchAddress("e3SCPhi", &e3SCPhi2);
    //_ntuple1->SetBranchAddress("e3SCRawEnergy", &e3SCRawEnergy1);
    //_ntuple2->SetBranchAddress("e3SCRawEnergy", &e3SCRawEnergy2);
    //_ntuple1->SetBranchAddress("e3SIP2D", &e3SIP2D1);
    //_ntuple2->SetBranchAddress("e3SIP2D", &e3SIP2D2);
    //_ntuple1->SetBranchAddress("e3SIP3D", &e3SIP3D1);
    //_ntuple2->SetBranchAddress("e3SIP3D", &e3SIP3D2);
    //_ntuple1->SetBranchAddress("e3ScaleGainDn", &e3ScaleGainDn1);
    //_ntuple2->SetBranchAddress("e3ScaleGainDn", &e3ScaleGainDn2);
    //_ntuple1->SetBranchAddress("e3ScaleGainUp", &e3ScaleGainUp1);
    //_ntuple2->SetBranchAddress("e3ScaleGainUp", &e3ScaleGainUp2);
    //_ntuple1->SetBranchAddress("e3ScaleStatDn", &e3ScaleStatDn1);
    //_ntuple2->SetBranchAddress("e3ScaleStatDn", &e3ScaleStatDn2);
    //_ntuple1->SetBranchAddress("e3ScaleStatUp", &e3ScaleStatUp1);
    //_ntuple2->SetBranchAddress("e3ScaleStatUp", &e3ScaleStatUp2);
    //_ntuple1->SetBranchAddress("e3ScaleSystDn", &e3ScaleSystDn1);
    //_ntuple2->SetBranchAddress("e3ScaleSystDn", &e3ScaleSystDn2);
    //_ntuple1->SetBranchAddress("e3ScaleSystUp", &e3ScaleSystUp1);
    //_ntuple2->SetBranchAddress("e3ScaleSystUp", &e3ScaleSystUp2);
    //_ntuple1->SetBranchAddress("e3ScaleTotDn", &e3ScaleTotDn1);
    //_ntuple2->SetBranchAddress("e3ScaleTotDn", &e3ScaleTotDn2);
    //_ntuple1->SetBranchAddress("e3ScaleTotUp", &e3ScaleTotUp1);
    //_ntuple2->SetBranchAddress("e3ScaleTotUp", &e3ScaleTotUp2);
    //_ntuple1->SetBranchAddress("e3SigmaPhiDn", &e3SigmaPhiDn1);
    //_ntuple2->SetBranchAddress("e3SigmaPhiDn", &e3SigmaPhiDn2);
    //_ntuple1->SetBranchAddress("e3SigmaPhiUp", &e3SigmaPhiUp1);
    //_ntuple2->SetBranchAddress("e3SigmaPhiUp", &e3SigmaPhiUp2);
    //_ntuple1->SetBranchAddress("e3SigmaRhoDn", &e3SigmaRhoDn1);
    //_ntuple2->SetBranchAddress("e3SigmaRhoDn", &e3SigmaRhoDn2);
    //_ntuple1->SetBranchAddress("e3SigmaRhoUp", &e3SigmaRhoUp1);
    //_ntuple2->SetBranchAddress("e3SigmaRhoUp", &e3SigmaRhoUp2);
    //_ntuple1->SetBranchAddress("e3SigmaTotDn", &e3SigmaTotDn1);
    //_ntuple2->SetBranchAddress("e3SigmaTotDn", &e3SigmaTotDn2);
    //_ntuple1->SetBranchAddress("e3SigmaTotUp", &e3SigmaTotUp1);
    //_ntuple2->SetBranchAddress("e3SigmaTotUp", &e3SigmaTotUp2);
    //_ntuple1->SetBranchAddress("e3Tight", &e3Tight1);
    //_ntuple2->SetBranchAddress("e3Tight", &e3Tight2);
    //_ntuple1->SetBranchAddress("e3UnCorrPt", &e3UnCorrPt1);
    //_ntuple2->SetBranchAddress("e3UnCorrPt", &e3UnCorrPt2);
    //_ntuple1->SetBranchAddress("e3Veto", &e3Veto1);
    //_ntuple2->SetBranchAddress("e3Veto", &e3Veto2);
    //_ntuple1->SetBranchAddress("e3ZZIso", &e3ZZIso1);
    //_ntuple2->SetBranchAddress("e3ZZIso", &e3ZZIso2);
    //_ntuple1->SetBranchAddress("e3IsEB", &e3IsEB1);
    //_ntuple2->SetBranchAddress("e3IsEB", &e3IsEB2);
    //_ntuple1->SetBranchAddress("e3IsGap", &e3IsGap1);
    //_ntuple2->SetBranchAddress("e3IsGap", &e3IsGap2);
    //_ntuple1->SetBranchAddress("e3ZZIsoPass", &e3ZZIsoPass1);
    //_ntuple2->SetBranchAddress("e3ZZIsoPass", &e3ZZIsoPass2);
    //_ntuple1->SetBranchAddress("e3ZZLooseID", &e3ZZLooseID1);
    //_ntuple2->SetBranchAddress("e3ZZLooseID", &e3ZZLooseID2);
    //_ntuple1->SetBranchAddress("e3ZZLooseIDNoVtx", &e3ZZLooseIDNoVtx1);
    //_ntuple2->SetBranchAddress("e3ZZLooseIDNoVtx", &e3ZZLooseIDNoVtx2);
    //_ntuple1->SetBranchAddress("e3ZZTightID", &e3ZZTightID1);
    //_ntuple2->SetBranchAddress("e3ZZTightID", &e3ZZTightID2);
    //_ntuple1->SetBranchAddress("e3ZZTightIDNoVtx", &e3ZZTightIDNoVtx1);
    //_ntuple2->SetBranchAddress("e3ZZTightIDNoVtx", &e3ZZTightIDNoVtx2);
    //_ntuple1->SetBranchAddress("e3Charge", &e3Charge1);
    //_ntuple2->SetBranchAddress("e3Charge", &e3Charge2);
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e3GenCharge", &e3GenCharge1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e3GenCharge", &e3GenCharge2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e3GenPdgId", &e3GenPdgId1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e3GenPdgId", &e3GenPdgId2);
    }
    _ntuple1->SetBranchAddress("e3PdgId", &e3PdgId1);
    _ntuple2->SetBranchAddress("e3PdgId", &e3PdgId2);
    //_ntuple1->SetBranchAddress("e3MissingHits", &e3MissingHits1);
    //_ntuple2->SetBranchAddress("e3MissingHits", &e3MissingHits2);
    //_ntuple1->SetBranchAddress("e4EffectiveArea", &e4EffectiveArea1);
    //_ntuple2->SetBranchAddress("e4EffectiveArea", &e4EffectiveArea2);
    _ntuple1->SetBranchAddress("e4Energy", &e4Energy1);
    _ntuple2->SetBranchAddress("e4Energy", &e4Energy2);
    _ntuple1->SetBranchAddress("e4Eta", &e4Eta1);
    _ntuple2->SetBranchAddress("e4Eta", &e4Eta2);
    //_ntuple1->SetBranchAddress("e4FSRDREt", &e4FSRDREt1);
    //_ntuple2->SetBranchAddress("e4FSRDREt", &e4FSRDREt2);
    //_ntuple1->SetBranchAddress("e4FSREt", &e4FSREt1);
    //_ntuple2->SetBranchAddress("e4FSREt", &e4FSREt2);
    //_ntuple1->SetBranchAddress("e4FSREta", &e4FSREta1);
    //_ntuple2->SetBranchAddress("e4FSREta", &e4FSREta2);
    //_ntuple1->SetBranchAddress("e4FSRPhi", &e4FSRPhi1);
    //_ntuple2->SetBranchAddress("e4FSRPhi", &e4FSRPhi2);
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e4GenEta", &e4GenEta1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e4GenEta", &e4GenEta2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e4GenPhi", &e4GenPhi1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e4GenPhi", &e4GenPhi2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e4GenPt", &e4GenPt1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e4GenPt", &e4GenPt2);
    }
    //_ntuple1->SetBranchAddress("e4IP2D", &e4IP2D1);
    //_ntuple2->SetBranchAddress("e4IP2D", &e4IP2D2);
    //_ntuple1->SetBranchAddress("e4IP2DUncertainty", &e4IP2DUncertainty1);
    //_ntuple2->SetBranchAddress("e4IP2DUncertainty", &e4IP2DUncertainty2);
    //_ntuple1->SetBranchAddress("e4IP3D", &e4IP3D1);
    //_ntuple2->SetBranchAddress("e4IP3D", &e4IP3D2);
    //_ntuple1->SetBranchAddress("e4IP3DUncertainty", &e4IP3DUncertainty1);
    //_ntuple2->SetBranchAddress("e4IP3DUncertainty", &e4IP3DUncertainty2);
    //_ntuple1->SetBranchAddress("e4Loose", &e4Loose1);
    //_ntuple2->SetBranchAddress("e4Loose", &e4Loose2);
    //_ntuple1->SetBranchAddress("e4Mass", &e4Mass1);
    //_ntuple2->SetBranchAddress("e4Mass", &e4Mass2);
    //_ntuple1->SetBranchAddress("e4Medium", &e4Medium1);
    //_ntuple2->SetBranchAddress("e4Medium", &e4Medium2);
    //_ntuple1->SetBranchAddress("e4MtToMET", &e4MtToMET1);
    //_ntuple2->SetBranchAddress("e4MtToMET", &e4MtToMET2);
    //_ntuple1->SetBranchAddress("e4PFChargedIso", &e4PFChargedIso1);
    //_ntuple2->SetBranchAddress("e4PFChargedIso", &e4PFChargedIso2);
    //_ntuple1->SetBranchAddress("e4PFNeutralIso", &e4PFNeutralIso1);
    //_ntuple2->SetBranchAddress("e4PFNeutralIso", &e4PFNeutralIso2);
    //_ntuple1->SetBranchAddress("e4PFPUIso", &e4PFPUIso1);
    //_ntuple2->SetBranchAddress("e4PFPUIso", &e4PFPUIso2);
    //_ntuple1->SetBranchAddress("e4PFPhotonIso", &e4PFPhotonIso1);
    //_ntuple2->SetBranchAddress("e4PFPhotonIso", &e4PFPhotonIso2);
    //_ntuple1->SetBranchAddress("e4PVDXY", &e4PVDXY1);
    //_ntuple2->SetBranchAddress("e4PVDXY", &e4PVDXY2);
    //_ntuple1->SetBranchAddress("e4PVDZ", &e4PVDZ1);
    //_ntuple2->SetBranchAddress("e4PVDZ", &e4PVDZ2);
    _ntuple1->SetBranchAddress("e4Phi", &e4Phi1);
    _ntuple2->SetBranchAddress("e4Phi", &e4Phi2);
    _ntuple1->SetBranchAddress("e4Pt", &e4Pt1);
    _ntuple2->SetBranchAddress("e4Pt", &e4Pt2);
    //_ntuple1->SetBranchAddress("e4RelPFIsoRho", &e4RelPFIsoRho1);
    //_ntuple2->SetBranchAddress("e4RelPFIsoRho", &e4RelPFIsoRho2);
    //_ntuple1->SetBranchAddress("e4Rho", &e4Rho1);
    //_ntuple2->SetBranchAddress("e4Rho", &e4Rho2);
    //_ntuple1->SetBranchAddress("e4SCEnergy", &e4SCEnergy1);
    //_ntuple2->SetBranchAddress("e4SCEnergy", &e4SCEnergy2);
    //_ntuple1->SetBranchAddress("e4SCEta", &e4SCEta1);
    //_ntuple2->SetBranchAddress("e4SCEta", &e4SCEta2);
    //_ntuple1->SetBranchAddress("e4SCPhi", &e4SCPhi1);
    //_ntuple2->SetBranchAddress("e4SCPhi", &e4SCPhi2);
    //_ntuple1->SetBranchAddress("e4SCRawEnergy", &e4SCRawEnergy1);
    //_ntuple2->SetBranchAddress("e4SCRawEnergy", &e4SCRawEnergy2);
    //_ntuple1->SetBranchAddress("e4SIP2D", &e4SIP2D1);
    //_ntuple2->SetBranchAddress("e4SIP2D", &e4SIP2D2);
    //_ntuple1->SetBranchAddress("e4SIP3D", &e4SIP3D1);
    //_ntuple2->SetBranchAddress("e4SIP3D", &e4SIP3D2);
    //_ntuple1->SetBranchAddress("e4ScaleGainDn", &e4ScaleGainDn1);
    //_ntuple2->SetBranchAddress("e4ScaleGainDn", &e4ScaleGainDn2);
    //_ntuple1->SetBranchAddress("e4ScaleGainUp", &e4ScaleGainUp1);
    //_ntuple2->SetBranchAddress("e4ScaleGainUp", &e4ScaleGainUp2);
    //_ntuple1->SetBranchAddress("e4ScaleStatDn", &e4ScaleStatDn1);
    //_ntuple2->SetBranchAddress("e4ScaleStatDn", &e4ScaleStatDn2);
    //_ntuple1->SetBranchAddress("e4ScaleStatUp", &e4ScaleStatUp1);
    //_ntuple2->SetBranchAddress("e4ScaleStatUp", &e4ScaleStatUp2);
    //_ntuple1->SetBranchAddress("e4ScaleSystDn", &e4ScaleSystDn1);
    //_ntuple2->SetBranchAddress("e4ScaleSystDn", &e4ScaleSystDn2);
    //_ntuple1->SetBranchAddress("e4ScaleSystUp", &e4ScaleSystUp1);
    //_ntuple2->SetBranchAddress("e4ScaleSystUp", &e4ScaleSystUp2);
    //_ntuple1->SetBranchAddress("e4ScaleTotDn", &e4ScaleTotDn1);
    //_ntuple2->SetBranchAddress("e4ScaleTotDn", &e4ScaleTotDn2);
    //_ntuple1->SetBranchAddress("e4ScaleTotUp", &e4ScaleTotUp1);
    //_ntuple2->SetBranchAddress("e4ScaleTotUp", &e4ScaleTotUp2);
    //_ntuple1->SetBranchAddress("e4SigmaPhiDn", &e4SigmaPhiDn1);
    //_ntuple2->SetBranchAddress("e4SigmaPhiDn", &e4SigmaPhiDn2);
    //_ntuple1->SetBranchAddress("e4SigmaPhiUp", &e4SigmaPhiUp1);
    //_ntuple2->SetBranchAddress("e4SigmaPhiUp", &e4SigmaPhiUp2);
    //_ntuple1->SetBranchAddress("e4SigmaRhoDn", &e4SigmaRhoDn1);
    //_ntuple2->SetBranchAddress("e4SigmaRhoDn", &e4SigmaRhoDn2);
    //_ntuple1->SetBranchAddress("e4SigmaRhoUp", &e4SigmaRhoUp1);
    //_ntuple2->SetBranchAddress("e4SigmaRhoUp", &e4SigmaRhoUp2);
    //_ntuple1->SetBranchAddress("e4SigmaTotDn", &e4SigmaTotDn1);
    //_ntuple2->SetBranchAddress("e4SigmaTotDn", &e4SigmaTotDn2);
    //_ntuple1->SetBranchAddress("e4SigmaTotUp", &e4SigmaTotUp1);
    //_ntuple2->SetBranchAddress("e4SigmaTotUp", &e4SigmaTotUp2);
    //_ntuple1->SetBranchAddress("e4Tight", &e4Tight1);
    //_ntuple2->SetBranchAddress("e4Tight", &e4Tight2);
    //_ntuple1->SetBranchAddress("e4UnCorrPt", &e4UnCorrPt1);
    //_ntuple2->SetBranchAddress("e4UnCorrPt", &e4UnCorrPt2);
    //_ntuple1->SetBranchAddress("e4Veto", &e4Veto1);
    //_ntuple2->SetBranchAddress("e4Veto", &e4Veto2);
    //_ntuple1->SetBranchAddress("e4ZZIso", &e4ZZIso1);
    //_ntuple2->SetBranchAddress("e4ZZIso", &e4ZZIso2);
    //_ntuple1->SetBranchAddress("e4IsEB", &e4IsEB1);
    //_ntuple2->SetBranchAddress("e4IsEB", &e4IsEB2);
    //_ntuple1->SetBranchAddress("e4IsGap", &e4IsGap1);
    //_ntuple2->SetBranchAddress("e4IsGap", &e4IsGap2);
    //_ntuple1->SetBranchAddress("e4ZZIsoPass", &e4ZZIsoPass1);
    //_ntuple2->SetBranchAddress("e4ZZIsoPass", &e4ZZIsoPass2);
    //_ntuple1->SetBranchAddress("e4ZZLooseID", &e4ZZLooseID1);
    //_ntuple2->SetBranchAddress("e4ZZLooseID", &e4ZZLooseID2);
    //_ntuple1->SetBranchAddress("e4ZZLooseIDNoVtx", &e4ZZLooseIDNoVtx1);
    //_ntuple2->SetBranchAddress("e4ZZLooseIDNoVtx", &e4ZZLooseIDNoVtx2);
    //_ntuple1->SetBranchAddress("e4ZZTightID", &e4ZZTightID1);
    //_ntuple2->SetBranchAddress("e4ZZTightID", &e4ZZTightID2);
    //_ntuple1->SetBranchAddress("e4ZZTightIDNoVtx", &e4ZZTightIDNoVtx1);
    //_ntuple2->SetBranchAddress("e4ZZTightIDNoVtx", &e4ZZTightIDNoVtx2);
    //_ntuple1->SetBranchAddress("e4Charge", &e4Charge1);
    //_ntuple2->SetBranchAddress("e4Charge", &e4Charge2);
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e4GenCharge", &e4GenCharge1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e4GenCharge", &e4GenCharge2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e4GenPdgId", &e4GenPdgId1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e4GenPdgId", &e4GenPdgId2);
    }
    _ntuple1->SetBranchAddress("e4PdgId", &e4PdgId1);
    _ntuple2->SetBranchAddress("e4PdgId", &e4PdgId2);
    //_ntuple1->SetBranchAddress("e4MissingHits", &e4MissingHits1);
    //_ntuple2->SetBranchAddress("e4MissingHits", &e4MissingHits2);
  }

  else if (_channel == "eemm"){
    //_ntuple1->SetBranchAddress("e1_m1_DR", &e1_m1_DR1);
    //_ntuple2->SetBranchAddress("e1_m1_DR", &e1_m1_DR2);
    //_ntuple1->SetBranchAddress("e1_m1_Mass", &e1_m1_Mass1);
    //_ntuple2->SetBranchAddress("e1_m1_Mass", &e1_m1_Mass2);
    //_ntuple1->SetBranchAddress("e1_m1_MassNoFSR", &e1_m1_MassNoFSR1);
    //_ntuple2->SetBranchAddress("e1_m1_MassNoFSR", &e1_m1_MassNoFSR2);
    //_ntuple1->SetBranchAddress("e1_m2_DR", &e1_m2_DR1);
    //_ntuple2->SetBranchAddress("e1_m2_DR", &e1_m2_DR2);
    //_ntuple1->SetBranchAddress("e1_m2_Mass", &e1_m2_Mass1);
    //_ntuple2->SetBranchAddress("e1_m2_Mass", &e1_m2_Mass2);
    //_ntuple1->SetBranchAddress("e1_m2_MassNoFSR", &e1_m2_MassNoFSR1);
    //_ntuple2->SetBranchAddress("e1_m2_MassNoFSR", &e1_m2_MassNoFSR2);
    //_ntuple1->SetBranchAddress("e2_m1_DR", &e2_m1_DR1);
    //_ntuple2->SetBranchAddress("e2_m1_DR", &e2_m1_DR2);
    //_ntuple1->SetBranchAddress("e2_m1_Mass", &e2_m1_Mass1);
    //_ntuple2->SetBranchAddress("e2_m1_Mass", &e2_m1_Mass2);
    //_ntuple1->SetBranchAddress("e2_m1_MassNoFSR", &e2_m1_MassNoFSR1);
    //_ntuple2->SetBranchAddress("e2_m1_MassNoFSR", &e2_m1_MassNoFSR2);
    //_ntuple1->SetBranchAddress("e2_m2_DR", &e2_m2_DR1);
    //_ntuple2->SetBranchAddress("e2_m2_DR", &e2_m2_DR2);
    //_ntuple1->SetBranchAddress("e2_m2_Mass", &e2_m2_Mass1);
    //_ntuple2->SetBranchAddress("e2_m2_Mass", &e2_m2_Mass2);
    //_ntuple1->SetBranchAddress("e2_m2_MassNoFSR", &e2_m2_MassNoFSR1);
    //_ntuple2->SetBranchAddress("e2_m2_MassNoFSR", &e2_m2_MassNoFSR2);
    //_ntuple1->SetBranchAddress("e1_m1_SS", &e1_m1_SS1);
    //_ntuple2->SetBranchAddress("e1_m1_SS", &e1_m1_SS2);
    //_ntuple1->SetBranchAddress("e1_m2_SS", &e1_m2_SS1);
    //_ntuple2->SetBranchAddress("e1_m2_SS", &e1_m2_SS2);
    //_ntuple1->SetBranchAddress("e2_m1_SS", &e2_m1_SS1);
    //_ntuple2->SetBranchAddress("e2_m1_SS", &e2_m1_SS2);
    //_ntuple1->SetBranchAddress("e2_m2_SS", &e2_m2_SS1);
    //_ntuple2->SetBranchAddress("e2_m2_SS", &e2_m2_SS2);
  }

  else if (_channel == "mmmm"){
    //_ntuple1->SetBranchAddress("m1_m3_DR", &m1_m3_DR1);
    //_ntuple2->SetBranchAddress("m1_m3_DR", &m1_m3_DR2);
    //_ntuple1->SetBranchAddress("m1_m3_Mass", &m1_m3_Mass1);
    //_ntuple2->SetBranchAddress("m1_m3_Mass", &m1_m3_Mass2);
    //_ntuple1->SetBranchAddress("m1_m3_MassNoFSR", &m1_m3_MassNoFSR1);
    //_ntuple2->SetBranchAddress("m1_m3_MassNoFSR", &m1_m3_MassNoFSR2);
    //_ntuple1->SetBranchAddress("m1_m4_DR", &m1_m4_DR1);
    //_ntuple2->SetBranchAddress("m1_m4_DR", &m1_m4_DR2);
    //_ntuple1->SetBranchAddress("m1_m4_Mass", &m1_m4_Mass1);
    //_ntuple2->SetBranchAddress("m1_m4_Mass", &m1_m4_Mass2);
    //_ntuple1->SetBranchAddress("m1_m4_MassNoFSR", &m1_m4_MassNoFSR1);
    //_ntuple2->SetBranchAddress("m1_m4_MassNoFSR", &m1_m4_MassNoFSR2);
    //_ntuple1->SetBranchAddress("m2_m3_DR", &m2_m3_DR1);
    //_ntuple2->SetBranchAddress("m2_m3_DR", &m2_m3_DR2);
    //_ntuple1->SetBranchAddress("m2_m3_Mass", &m2_m3_Mass1);
    //_ntuple2->SetBranchAddress("m2_m3_Mass", &m2_m3_Mass2);
    //_ntuple1->SetBranchAddress("m2_m3_MassNoFSR", &m2_m3_MassNoFSR1);
    //_ntuple2->SetBranchAddress("m2_m3_MassNoFSR", &m2_m3_MassNoFSR2);
    //_ntuple1->SetBranchAddress("m2_m4_DR", &m2_m4_DR1);
    //_ntuple2->SetBranchAddress("m2_m4_DR", &m2_m4_DR2);
    //_ntuple1->SetBranchAddress("m2_m4_Mass", &m2_m4_Mass1);
    //_ntuple2->SetBranchAddress("m2_m4_Mass", &m2_m4_Mass2);
    //_ntuple1->SetBranchAddress("m2_m4_MassNoFSR", &m2_m4_MassNoFSR1);
    //_ntuple2->SetBranchAddress("m2_m4_MassNoFSR", &m2_m4_MassNoFSR2);
    //_ntuple1->SetBranchAddress("m1_m3_SS", &m1_m3_SS1);
    //_ntuple2->SetBranchAddress("m1_m3_SS", &m1_m3_SS2);
    //_ntuple1->SetBranchAddress("m1_m4_SS", &m1_m4_SS1);
    //_ntuple2->SetBranchAddress("m1_m4_SS", &m1_m4_SS2);
    //_ntuple1->SetBranchAddress("m2_m3_SS", &m2_m3_SS1);
    //_ntuple2->SetBranchAddress("m2_m3_SS", &m2_m3_SS2);
    //_ntuple1->SetBranchAddress("m2_m4_SS", &m2_m4_SS1);
    //_ntuple2->SetBranchAddress("m2_m4_SS", &m2_m4_SS2);
    //_ntuple1->SetBranchAddress("m3_m4_DR", &m3_m4_DR1);
    //_ntuple2->SetBranchAddress("m3_m4_DR", &m3_m4_DR2);
    //_ntuple1->SetBranchAddress("m3_m4_Energy", &m3_m4_Energy1);
    //_ntuple2->SetBranchAddress("m3_m4_Energy", &m3_m4_Energy2);
    //_ntuple1->SetBranchAddress("m3_m4_EnergyNoFSR", &m3_m4_EnergyNoFSR1);
    //_ntuple2->SetBranchAddress("m3_m4_EnergyNoFSR", &m3_m4_EnergyNoFSR2);
    //_ntuple1->SetBranchAddress("m3_m4_Eta", &m3_m4_Eta1);
    //_ntuple2->SetBranchAddress("m3_m4_Eta", &m3_m4_Eta2);
    //_ntuple1->SetBranchAddress("m3_m4_EtaNoFSR", &m3_m4_EtaNoFSR1);
    //_ntuple2->SetBranchAddress("m3_m4_EtaNoFSR", &m3_m4_EtaNoFSR2);
    _ntuple1->SetBranchAddress("m3_m4_Mass", &m3_m4_Mass1);
    _ntuple2->SetBranchAddress("m3_m4_Mass", &m3_m4_Mass2);
    //_ntuple1->SetBranchAddress("m3_m4_MassNoFSR", &m3_m4_MassNoFSR1);
    //_ntuple2->SetBranchAddress("m3_m4_MassNoFSR", &m3_m4_MassNoFSR2);
    //_ntuple1->SetBranchAddress("m3_m4_Mt", &m3_m4_Mt1);
    //_ntuple2->SetBranchAddress("m3_m4_Mt", &m3_m4_Mt2);
    //_ntuple1->SetBranchAddress("m3_m4_MtToMET", &m3_m4_MtToMET1);
    //_ntuple2->SetBranchAddress("m3_m4_MtToMET", &m3_m4_MtToMET2);
    //_ntuple1->SetBranchAddress("m3_m4_Phi", &m3_m4_Phi1);
    //_ntuple2->SetBranchAddress("m3_m4_Phi", &m3_m4_Phi2);
    //_ntuple1->SetBranchAddress("m3_m4_PhiNoFSR", &m3_m4_PhiNoFSR1);
    //_ntuple2->SetBranchAddress("m3_m4_PhiNoFSR", &m3_m4_PhiNoFSR2);
    //_ntuple1->SetBranchAddress("m3_m4_Pt", &m3_m4_Pt1);
    //_ntuple2->SetBranchAddress("m3_m4_Pt", &m3_m4_Pt2);
    //_ntuple1->SetBranchAddress("m3_m4_PtNoFSR", &m3_m4_PtNoFSR1);
    //_ntuple2->SetBranchAddress("m3_m4_PtNoFSR", &m3_m4_PtNoFSR2);
    //_ntuple1->SetBranchAddress("m3_m4_SS", &m3_m4_SS1);
    //_ntuple2->SetBranchAddress("m3_m4_SS", &m3_m4_SS2);
    //_ntuple1->SetBranchAddress("m3_m4_Charge", &m3_m4_Charge1);
    //_ntuple2->SetBranchAddress("m3_m4_Charge", &m3_m4_Charge2);
    //_ntuple1->SetBranchAddress("m3_m4_PdgId", &m3_m4_PdgId1);
    //_ntuple2->SetBranchAddress("m3_m4_PdgId", &m3_m4_PdgId2);
    //_ntuple1->SetBranchAddress("m3EffScaleFactor", &m3EffScaleFactor1);
    //_ntuple2->SetBranchAddress("m3EffScaleFactor", &m3EffScaleFactor2);
    //_ntuple1->SetBranchAddress("m3EffScaleFactorError", &m3EffScaleFactorError1);
    //_ntuple2->SetBranchAddress("m3EffScaleFactorError", &m3EffScaleFactorError2);
    _ntuple1->SetBranchAddress("m3Energy", &m3Energy1);
    _ntuple2->SetBranchAddress("m3Energy", &m3Energy2);
    _ntuple1->SetBranchAddress("m3Eta", &m3Eta1);
    _ntuple2->SetBranchAddress("m3Eta", &m3Eta2);
    //_ntuple1->SetBranchAddress("m3FSRDREt", &m3FSRDREt1);
    //_ntuple2->SetBranchAddress("m3FSRDREt", &m3FSRDREt2);
    //_ntuple1->SetBranchAddress("m3FSREt", &m3FSREt1);
    //_ntuple2->SetBranchAddress("m3FSREt", &m3FSREt2);
    //_ntuple1->SetBranchAddress("m3FSREta", &m3FSREta1);
    //_ntuple2->SetBranchAddress("m3FSREta", &m3FSREta2);
    //_ntuple1->SetBranchAddress("m3FSRPhi", &m3FSRPhi1);
    //_ntuple2->SetBranchAddress("m3FSRPhi", &m3FSRPhi2);
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m3GenEta", &m3GenEta1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m3GenEta", &m3GenEta2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m3GenPhi", &m3GenPhi1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m3GenPhi", &m3GenPhi2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m3GenPt", &m3GenPt1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m3GenPt", &m3GenPt2);
    }
    //_ntuple1->SetBranchAddress("m3IP2D", &m3IP2D1);
    //_ntuple2->SetBranchAddress("m3IP2D", &m3IP2D2);
    //_ntuple1->SetBranchAddress("m3IP2DUncertainty", &m3IP2DUncertainty1);
    //_ntuple2->SetBranchAddress("m3IP2DUncertainty", &m3IP2DUncertainty2);
    //_ntuple1->SetBranchAddress("m3IP3D", &m3IP3D1);
    //_ntuple2->SetBranchAddress("m3IP3D", &m3IP3D2);
    //_ntuple1->SetBranchAddress("m3IP3DUncertainty", &m3IP3DUncertainty1);
    //_ntuple2->SetBranchAddress("m3IP3DUncertainty", &m3IP3DUncertainty2);
    //_ntuple1->SetBranchAddress("m3Mass", &m3Mass1);
    //_ntuple2->SetBranchAddress("m3Mass", &m3Mass2);
    //_ntuple1->SetBranchAddress("m3MtToMET", &m3MtToMET1);
    //_ntuple2->SetBranchAddress("m3MtToMET", &m3MtToMET2);
    //_ntuple1->SetBranchAddress("m3PFChargedIso", &m3PFChargedIso1);
    //_ntuple2->SetBranchAddress("m3PFChargedIso", &m3PFChargedIso2);
    //_ntuple1->SetBranchAddress("m3PFNeutralIso", &m3PFNeutralIso1);
    //_ntuple2->SetBranchAddress("m3PFNeutralIso", &m3PFNeutralIso2);
    //_ntuple1->SetBranchAddress("m3PFPUIso", &m3PFPUIso1);
    //_ntuple2->SetBranchAddress("m3PFPUIso", &m3PFPUIso2);
    //_ntuple1->SetBranchAddress("m3PFPhotonIso", &m3PFPhotonIso1);
    //_ntuple2->SetBranchAddress("m3PFPhotonIso", &m3PFPhotonIso2);
    //_ntuple1->SetBranchAddress("m3PVDXY", &m3PVDXY1);
    //_ntuple2->SetBranchAddress("m3PVDXY", &m3PVDXY2);
    //_ntuple1->SetBranchAddress("m3PVDZ", &m3PVDZ1);
    //_ntuple2->SetBranchAddress("m3PVDZ", &m3PVDZ2);
    _ntuple1->SetBranchAddress("m3Phi", &m3Phi1);
    _ntuple2->SetBranchAddress("m3Phi", &m3Phi2);
    _ntuple1->SetBranchAddress("m3Pt", &m3Pt1);
    _ntuple2->SetBranchAddress("m3Pt", &m3Pt2);
    //_ntuple1->SetBranchAddress("m3PtErr", &m3PtErr1);
    //_ntuple2->SetBranchAddress("m3PtErr", &m3PtErr2);
    //_ntuple1->SetBranchAddress("m3PtErrUncorrected", &m3PtErrUncorrected1);
    //_ntuple2->SetBranchAddress("m3PtErrUncorrected", &m3PtErrUncorrected2);
    //_ntuple1->SetBranchAddress("m3PtUncorrected", &m3PtUncorrected1);
    //_ntuple2->SetBranchAddress("m3PtUncorrected", &m3PtUncorrected2);
    //_ntuple1->SetBranchAddress("m3RelPFIsoDBR03", &m3RelPFIsoDBR031);
    //_ntuple2->SetBranchAddress("m3RelPFIsoDBR03", &m3RelPFIsoDBR032);
    //_ntuple1->SetBranchAddress("m3RelPFIsoDBR04", &m3RelPFIsoDBR041);
    //_ntuple2->SetBranchAddress("m3RelPFIsoDBR04", &m3RelPFIsoDBR042);
    //_ntuple1->SetBranchAddress("m3RochesterScaleUncDn", &m3RochesterScaleUncDn1);
    //_ntuple2->SetBranchAddress("m3RochesterScaleUncDn", &m3RochesterScaleUncDn2);
    //_ntuple1->SetBranchAddress("m3RochesterScaleUncUp", &m3RochesterScaleUncUp1);
    //_ntuple2->SetBranchAddress("m3RochesterScaleUncUp", &m3RochesterScaleUncUp2);
    //_ntuple1->SetBranchAddress("m3RochesterSmearUncDn", &m3RochesterSmearUncDn1);
    //_ntuple2->SetBranchAddress("m3RochesterSmearUncDn", &m3RochesterSmearUncDn2);
    //_ntuple1->SetBranchAddress("m3RochesterSmearUncUp", &m3RochesterSmearUncUp1);
    //_ntuple2->SetBranchAddress("m3RochesterSmearUncUp", &m3RochesterSmearUncUp2);
    //_ntuple1->SetBranchAddress("m3SIP2D", &m3SIP2D1);
    //_ntuple2->SetBranchAddress("m3SIP2D", &m3SIP2D2);
    //_ntuple1->SetBranchAddress("m3SIP3D", &m3SIP3D1);
    //_ntuple2->SetBranchAddress("m3SIP3D", &m3SIP3D2);
    //_ntuple1->SetBranchAddress("m3TrackIso", &m3TrackIso1);
    //_ntuple2->SetBranchAddress("m3TrackIso", &m3TrackIso2);
    //_ntuple1->SetBranchAddress("m3ZZIso", &m3ZZIso1);
    //_ntuple2->SetBranchAddress("m3ZZIso", &m3ZZIso2);
    //_ntuple1->SetBranchAddress("m3HighPtID", &m3HighPtID1);
    //_ntuple2->SetBranchAddress("m3HighPtID", &m3HighPtID2);
    //_ntuple1->SetBranchAddress("m3HighPtIDNoVtx", &m3HighPtIDNoVtx1);
    //_ntuple2->SetBranchAddress("m3HighPtIDNoVtx", &m3HighPtIDNoVtx2);
    //_ntuple1->SetBranchAddress("m3IsGlobal", &m3IsGlobal1);
    //_ntuple2->SetBranchAddress("m3IsGlobal", &m3IsGlobal2);
    //_ntuple1->SetBranchAddress("m3IsLoose", &m3IsLoose1);
    //_ntuple2->SetBranchAddress("m3IsLoose", &m3IsLoose2);
    //_ntuple1->SetBranchAddress("m3IsMedium", &m3IsMedium1);
    //_ntuple2->SetBranchAddress("m3IsMedium", &m3IsMedium2);
    //_ntuple1->SetBranchAddress("m3IsPFMuon", &m3IsPFMuon1);
    //_ntuple2->SetBranchAddress("m3IsPFMuon", &m3IsPFMuon2);
    //_ntuple1->SetBranchAddress("m3IsTracker", &m3IsTracker1);
    //_ntuple2->SetBranchAddress("m3IsTracker", &m3IsTracker2);
    //_ntuple1->SetBranchAddress("m3PASTightID", &m3PASTightID1);
    //_ntuple2->SetBranchAddress("m3PASTightID", &m3PASTightID2);
    //_ntuple1->SetBranchAddress("m3PASTightIDNoVtx", &m3PASTightIDNoVtx1);
    //_ntuple2->SetBranchAddress("m3PASTightIDNoVtx", &m3PASTightIDNoVtx2);
    //_ntuple1->SetBranchAddress("m3PFID", &m3PFID1);
    //_ntuple2->SetBranchAddress("m3PFID", &m3PFID2);
    //_ntuple1->SetBranchAddress("m3PFIDNoVtx", &m3PFIDNoVtx1);
    //_ntuple2->SetBranchAddress("m3PFIDNoVtx", &m3PFIDNoVtx2);
    //_ntuple1->SetBranchAddress("m3ZZIsoPass", &m3ZZIsoPass1);
    //_ntuple2->SetBranchAddress("m3ZZIsoPass", &m3ZZIsoPass2);
    //_ntuple1->SetBranchAddress("m3ZZLooseID", &m3ZZLooseID1);
    //_ntuple2->SetBranchAddress("m3ZZLooseID", &m3ZZLooseID2);
    //_ntuple1->SetBranchAddress("m3ZZLooseIDNoVtx", &m3ZZLooseIDNoVtx1);
    //_ntuple2->SetBranchAddress("m3ZZLooseIDNoVtx", &m3ZZLooseIDNoVtx2);
    //_ntuple1->SetBranchAddress("m3ZZTightID", &m3ZZTightID1);
    //_ntuple2->SetBranchAddress("m3ZZTightID", &m3ZZTightID2);
    //_ntuple1->SetBranchAddress("m3ZZTightIDNoVtx", &m3ZZTightIDNoVtx1);
    //_ntuple2->SetBranchAddress("m3ZZTightIDNoVtx", &m3ZZTightIDNoVtx2);
    //_ntuple1->SetBranchAddress("m3Charge", &m3Charge1);
    //_ntuple2->SetBranchAddress("m3Charge", &m3Charge2);
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m3GenCharge", &m3GenCharge1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m3GenCharge", &m3GenCharge2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m3GenPdgId", &m3GenPdgId1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m3GenPdgId", &m3GenPdgId2);
    }
    _ntuple1->SetBranchAddress("m3PdgId", &m3PdgId1);
    _ntuple2->SetBranchAddress("m3PdgId", &m3PdgId2);
    //_ntuple1->SetBranchAddress("m3BestTrackType", &m3BestTrackType1);
    //_ntuple2->SetBranchAddress("m3BestTrackType", &m3BestTrackType2);
    //_ntuple1->SetBranchAddress("m3CutBasedLoose", &m3CutBasedLoose1);
    //_ntuple2->SetBranchAddress("m3CutBasedLoose", &m3CutBasedLoose2);
    //_ntuple1->SetBranchAddress("m3CutBasedMedium", &m3CutBasedMedium1);
    //_ntuple2->SetBranchAddress("m3CutBasedMedium", &m3CutBasedMedium2);
    //_ntuple1->SetBranchAddress("m3CutBasedTight", &m3CutBasedTight1);
    //_ntuple2->SetBranchAddress("m3CutBasedTight", &m3CutBasedTight2);
    //_ntuple1->SetBranchAddress("m3IsTight", &m3IsTight1);
    //_ntuple2->SetBranchAddress("m3IsTight", &m3IsTight2);
    //_ntuple1->SetBranchAddress("m3MatchedStations", &m3MatchedStations1);
    //_ntuple2->SetBranchAddress("m3MatchedStations", &m3MatchedStations2);
    //_ntuple1->SetBranchAddress("m3NoOfMatches", &m3NoOfMatches1);
    //_ntuple2->SetBranchAddress("m3NoOfMatches", &m3NoOfMatches2);
    //_ntuple1->SetBranchAddress("m3PFIsoLoose", &m3PFIsoLoose1);
    //_ntuple2->SetBranchAddress("m3PFIsoLoose", &m3PFIsoLoose2);
    //_ntuple1->SetBranchAddress("m3PFIsoMedium", &m3PFIsoMedium1);
    //_ntuple2->SetBranchAddress("m3PFIsoMedium", &m3PFIsoMedium2);
    //_ntuple1->SetBranchAddress("m3PFIsoTight", &m3PFIsoTight1);
    //_ntuple2->SetBranchAddress("m3PFIsoTight", &m3PFIsoTight2);
    //_ntuple1->SetBranchAddress("m3PFIsoVTight", &m3PFIsoVTight1);
    //_ntuple2->SetBranchAddress("m3PFIsoVTight", &m3PFIsoVTight2);
    //_ntuple1->SetBranchAddress("m4EffScaleFactor", &m4EffScaleFactor1);
    //_ntuple2->SetBranchAddress("m4EffScaleFactor", &m4EffScaleFactor2);
    //_ntuple1->SetBranchAddress("m4EffScaleFactorError", &m4EffScaleFactorError1);
    //_ntuple2->SetBranchAddress("m4EffScaleFactorError", &m4EffScaleFactorError2);
    _ntuple1->SetBranchAddress("m4Energy", &m4Energy1);
    _ntuple2->SetBranchAddress("m4Energy", &m4Energy2);
    _ntuple1->SetBranchAddress("m4Eta", &m4Eta1);
    _ntuple2->SetBranchAddress("m4Eta", &m4Eta2);
    //_ntuple1->SetBranchAddress("m4FSRDREt", &m4FSRDREt1);
    //_ntuple2->SetBranchAddress("m4FSRDREt", &m4FSRDREt2);
    //_ntuple1->SetBranchAddress("m4FSREt", &m4FSREt1);
    //_ntuple2->SetBranchAddress("m4FSREt", &m4FSREt2);
    //_ntuple1->SetBranchAddress("m4FSREta", &m4FSREta1);
    //_ntuple2->SetBranchAddress("m4FSREta", &m4FSREta2);
    //_ntuple1->SetBranchAddress("m4FSRPhi", &m4FSRPhi1);
    //_ntuple2->SetBranchAddress("m4FSRPhi", &m4FSRPhi2);
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m4GenEta", &m4GenEta1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m4GenEta", &m4GenEta2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m4GenPhi", &m4GenPhi1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m4GenPhi", &m4GenPhi2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m4GenPt", &m4GenPt1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m4GenPt", &m4GenPt2);
    }
    //_ntuple1->SetBranchAddress("m4IP2D", &m4IP2D1);
    //_ntuple2->SetBranchAddress("m4IP2D", &m4IP2D2);
    //_ntuple1->SetBranchAddress("m4IP2DUncertainty", &m4IP2DUncertainty1);
    //_ntuple2->SetBranchAddress("m4IP2DUncertainty", &m4IP2DUncertainty2);
    //_ntuple1->SetBranchAddress("m4IP3D", &m4IP3D1);
    //_ntuple2->SetBranchAddress("m4IP3D", &m4IP3D2);
    //_ntuple1->SetBranchAddress("m4IP3DUncertainty", &m4IP3DUncertainty1);
    //_ntuple2->SetBranchAddress("m4IP3DUncertainty", &m4IP3DUncertainty2);
    //_ntuple1->SetBranchAddress("m4Mass", &m4Mass1);
    //_ntuple2->SetBranchAddress("m4Mass", &m4Mass2);
    //_ntuple1->SetBranchAddress("m4MtToMET", &m4MtToMET1);
    //_ntuple2->SetBranchAddress("m4MtToMET", &m4MtToMET2);
    //_ntuple1->SetBranchAddress("m4PFChargedIso", &m4PFChargedIso1);
    //_ntuple2->SetBranchAddress("m4PFChargedIso", &m4PFChargedIso2);
    //_ntuple1->SetBranchAddress("m4PFNeutralIso", &m4PFNeutralIso1);
    //_ntuple2->SetBranchAddress("m4PFNeutralIso", &m4PFNeutralIso2);
    //_ntuple1->SetBranchAddress("m4PFPUIso", &m4PFPUIso1);
    //_ntuple2->SetBranchAddress("m4PFPUIso", &m4PFPUIso2);
    //_ntuple1->SetBranchAddress("m4PFPhotonIso", &m4PFPhotonIso1);
    //_ntuple2->SetBranchAddress("m4PFPhotonIso", &m4PFPhotonIso2);
    //_ntuple1->SetBranchAddress("m4PVDXY", &m4PVDXY1);
    //_ntuple2->SetBranchAddress("m4PVDXY", &m4PVDXY2);
    //_ntuple1->SetBranchAddress("m4PVDZ", &m4PVDZ1);
    //_ntuple2->SetBranchAddress("m4PVDZ", &m4PVDZ2);
    _ntuple1->SetBranchAddress("m4Phi", &m4Phi1);
    _ntuple2->SetBranchAddress("m4Phi", &m4Phi2);
    _ntuple1->SetBranchAddress("m4Pt", &m4Pt1);
    _ntuple2->SetBranchAddress("m4Pt", &m4Pt2);
    //_ntuple1->SetBranchAddress("m4PtErr", &m4PtErr1);
    //_ntuple2->SetBranchAddress("m4PtErr", &m4PtErr2);
    //_ntuple1->SetBranchAddress("m4PtErrUncorrected", &m4PtErrUncorrected1);
    //_ntuple2->SetBranchAddress("m4PtErrUncorrected", &m4PtErrUncorrected2);
    //_ntuple1->SetBranchAddress("m4PtUncorrected", &m4PtUncorrected1);
    //_ntuple2->SetBranchAddress("m4PtUncorrected", &m4PtUncorrected2);
    //_ntuple1->SetBranchAddress("m4RelPFIsoDBR03", &m4RelPFIsoDBR031);
    //_ntuple2->SetBranchAddress("m4RelPFIsoDBR03", &m4RelPFIsoDBR032);
    //_ntuple1->SetBranchAddress("m4RelPFIsoDBR04", &m4RelPFIsoDBR041);
    //_ntuple2->SetBranchAddress("m4RelPFIsoDBR04", &m4RelPFIsoDBR042);
    //_ntuple1->SetBranchAddress("m4RochesterScaleUncDn", &m4RochesterScaleUncDn1);
    //_ntuple2->SetBranchAddress("m4RochesterScaleUncDn", &m4RochesterScaleUncDn2);
    //_ntuple1->SetBranchAddress("m4RochesterScaleUncUp", &m4RochesterScaleUncUp1);
    //_ntuple2->SetBranchAddress("m4RochesterScaleUncUp", &m4RochesterScaleUncUp2);
    //_ntuple1->SetBranchAddress("m4RochesterSmearUncDn", &m4RochesterSmearUncDn1);
    //_ntuple2->SetBranchAddress("m4RochesterSmearUncDn", &m4RochesterSmearUncDn2);
    //_ntuple1->SetBranchAddress("m4RochesterSmearUncUp", &m4RochesterSmearUncUp1);
    //_ntuple2->SetBranchAddress("m4RochesterSmearUncUp", &m4RochesterSmearUncUp2);
    //_ntuple1->SetBranchAddress("m4SIP2D", &m4SIP2D1);
    //_ntuple2->SetBranchAddress("m4SIP2D", &m4SIP2D2);
    //_ntuple1->SetBranchAddress("m4SIP3D", &m4SIP3D1);
    //_ntuple2->SetBranchAddress("m4SIP3D", &m4SIP3D2);
    //_ntuple1->SetBranchAddress("m4TrackIso", &m4TrackIso1);
    //_ntuple2->SetBranchAddress("m4TrackIso", &m4TrackIso2);
    //_ntuple1->SetBranchAddress("m4ZZIso", &m4ZZIso1);
    //_ntuple2->SetBranchAddress("m4ZZIso", &m4ZZIso2);
    //_ntuple1->SetBranchAddress("m4HighPtID", &m4HighPtID1);
    //_ntuple2->SetBranchAddress("m4HighPtID", &m4HighPtID2);
    //_ntuple1->SetBranchAddress("m4HighPtIDNoVtx", &m4HighPtIDNoVtx1);
    //_ntuple2->SetBranchAddress("m4HighPtIDNoVtx", &m4HighPtIDNoVtx2);
    //_ntuple1->SetBranchAddress("m4IsGlobal", &m4IsGlobal1);
    //_ntuple2->SetBranchAddress("m4IsGlobal", &m4IsGlobal2);
    //_ntuple1->SetBranchAddress("m4IsLoose", &m4IsLoose1);
    //_ntuple2->SetBranchAddress("m4IsLoose", &m4IsLoose2);
    //_ntuple1->SetBranchAddress("m4IsMedium", &m4IsMedium1);
    //_ntuple2->SetBranchAddress("m4IsMedium", &m4IsMedium2);
    //_ntuple1->SetBranchAddress("m4IsPFMuon", &m4IsPFMuon1);
    //_ntuple2->SetBranchAddress("m4IsPFMuon", &m4IsPFMuon2);
    //_ntuple1->SetBranchAddress("m4IsTracker", &m4IsTracker1);
    //_ntuple2->SetBranchAddress("m4IsTracker", &m4IsTracker2);
    //_ntuple1->SetBranchAddress("m4PASTightID", &m4PASTightID1);
    //_ntuple2->SetBranchAddress("m4PASTightID", &m4PASTightID2);
    //_ntuple1->SetBranchAddress("m4PASTightIDNoVtx", &m4PASTightIDNoVtx1);
    //_ntuple2->SetBranchAddress("m4PASTightIDNoVtx", &m4PASTightIDNoVtx2);
    //_ntuple1->SetBranchAddress("m4PFID", &m4PFID1);
    //_ntuple2->SetBranchAddress("m4PFID", &m4PFID2);
    //_ntuple1->SetBranchAddress("m4PFIDNoVtx", &m4PFIDNoVtx1);
    //_ntuple2->SetBranchAddress("m4PFIDNoVtx", &m4PFIDNoVtx2);
    //_ntuple1->SetBranchAddress("m4ZZIsoPass", &m4ZZIsoPass1);
    //_ntuple2->SetBranchAddress("m4ZZIsoPass", &m4ZZIsoPass2);
    //_ntuple1->SetBranchAddress("m4ZZLooseID", &m4ZZLooseID1);
    //_ntuple2->SetBranchAddress("m4ZZLooseID", &m4ZZLooseID2);
    //_ntuple1->SetBranchAddress("m4ZZLooseIDNoVtx", &m4ZZLooseIDNoVtx1);
    //_ntuple2->SetBranchAddress("m4ZZLooseIDNoVtx", &m4ZZLooseIDNoVtx2);
    //_ntuple1->SetBranchAddress("m4ZZTightID", &m4ZZTightID1);
    //_ntuple2->SetBranchAddress("m4ZZTightID", &m4ZZTightID2);
    //_ntuple1->SetBranchAddress("m4ZZTightIDNoVtx", &m4ZZTightIDNoVtx1);
    //_ntuple2->SetBranchAddress("m4ZZTightIDNoVtx", &m4ZZTightIDNoVtx2);
    //_ntuple1->SetBranchAddress("m4Charge", &m4Charge1);
    //_ntuple2->SetBranchAddress("m4Charge", &m4Charge2);
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m4GenCharge", &m4GenCharge1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m4GenCharge", &m4GenCharge2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m4GenPdgId", &m4GenPdgId1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m4GenPdgId", &m4GenPdgId2);
    }
    _ntuple1->SetBranchAddress("m4PdgId", &m4PdgId1);
    _ntuple2->SetBranchAddress("m4PdgId", &m4PdgId2);
    //_ntuple1->SetBranchAddress("m4BestTrackType", &m4BestTrackType1);
    //_ntuple2->SetBranchAddress("m4BestTrackType", &m4BestTrackType2);
    //_ntuple1->SetBranchAddress("m4CutBasedLoose", &m4CutBasedLoose1);
    //_ntuple2->SetBranchAddress("m4CutBasedLoose", &m4CutBasedLoose2);
    //_ntuple1->SetBranchAddress("m4CutBasedMedium", &m4CutBasedMedium1);
    //_ntuple2->SetBranchAddress("m4CutBasedMedium", &m4CutBasedMedium2);
    //_ntuple1->SetBranchAddress("m4CutBasedTight", &m4CutBasedTight1);
    //_ntuple2->SetBranchAddress("m4CutBasedTight", &m4CutBasedTight2);
    //_ntuple1->SetBranchAddress("m4IsTight", &m4IsTight1);
    //_ntuple2->SetBranchAddress("m4IsTight", &m4IsTight2);
    //_ntuple1->SetBranchAddress("m4MatchedStations", &m4MatchedStations1);
    //_ntuple2->SetBranchAddress("m4MatchedStations", &m4MatchedStations2);
    //_ntuple1->SetBranchAddress("m4NoOfMatches", &m4NoOfMatches1);
    //_ntuple2->SetBranchAddress("m4NoOfMatches", &m4NoOfMatches2);
    //_ntuple1->SetBranchAddress("m4PFIsoLoose", &m4PFIsoLoose1);
    //_ntuple2->SetBranchAddress("m4PFIsoLoose", &m4PFIsoLoose2);
    //_ntuple1->SetBranchAddress("m4PFIsoMedium", &m4PFIsoMedium1);
    //_ntuple2->SetBranchAddress("m4PFIsoMedium", &m4PFIsoMedium2);
    //_ntuple1->SetBranchAddress("m4PFIsoTight", &m4PFIsoTight1);
    //_ntuple2->SetBranchAddress("m4PFIsoTight", &m4PFIsoTight2);
    //_ntuple1->SetBranchAddress("m4PFIsoVTight", &m4PFIsoVTight1);
    //_ntuple2->SetBranchAddress("m4PFIsoVTight", &m4PFIsoVTight2);
  }

  if (_channel == "eeee" || _channel == "eemm"){
    //_ntuple1->SetBranchAddress("e1_e2_DR", &e1_e2_DR1);
    //_ntuple2->SetBranchAddress("e1_e2_DR", &e1_e2_DR2);
    //_ntuple1->SetBranchAddress("e1_e2_Energy", &e1_e2_Energy1);
    //_ntuple2->SetBranchAddress("e1_e2_Energy", &e1_e2_Energy2);
    //_ntuple1->SetBranchAddress("e1_e2_EnergyNoFSR", &e1_e2_EnergyNoFSR1);
    //_ntuple2->SetBranchAddress("e1_e2_EnergyNoFSR", &e1_e2_EnergyNoFSR2);
    //_ntuple1->SetBranchAddress("e1_e2_Eta", &e1_e2_Eta1);
    //_ntuple2->SetBranchAddress("e1_e2_Eta", &e1_e2_Eta2);
    //_ntuple1->SetBranchAddress("e1_e2_EtaNoFSR", &e1_e2_EtaNoFSR1);
    //_ntuple2->SetBranchAddress("e1_e2_EtaNoFSR", &e1_e2_EtaNoFSR2);
    _ntuple1->SetBranchAddress("e1_e2_Mass", &e1_e2_Mass1);
    _ntuple2->SetBranchAddress("e1_e2_Mass", &e1_e2_Mass2);
    //_ntuple1->SetBranchAddress("e1_e2_MassNoFSR", &e1_e2_MassNoFSR1);
    //_ntuple2->SetBranchAddress("e1_e2_MassNoFSR", &e1_e2_MassNoFSR2);
    //_ntuple1->SetBranchAddress("e1_e2_Mt", &e1_e2_Mt1);
    //_ntuple2->SetBranchAddress("e1_e2_Mt", &e1_e2_Mt2);
    //_ntuple1->SetBranchAddress("e1_e2_MtToMET", &e1_e2_MtToMET1);
    //_ntuple2->SetBranchAddress("e1_e2_MtToMET", &e1_e2_MtToMET2);
    //_ntuple1->SetBranchAddress("e1_e2_Phi", &e1_e2_Phi1);
    //_ntuple2->SetBranchAddress("e1_e2_Phi", &e1_e2_Phi2);
    //_ntuple1->SetBranchAddress("e1_e2_PhiNoFSR", &e1_e2_PhiNoFSR1);
    //_ntuple2->SetBranchAddress("e1_e2_PhiNoFSR", &e1_e2_PhiNoFSR2);
    //_ntuple1->SetBranchAddress("e1_e2_Pt", &e1_e2_Pt1);
    //_ntuple2->SetBranchAddress("e1_e2_Pt", &e1_e2_Pt2);
    //_ntuple1->SetBranchAddress("e1_e2_PtNoFSR", &e1_e2_PtNoFSR1);
    //_ntuple2->SetBranchAddress("e1_e2_PtNoFSR", &e1_e2_PtNoFSR2);
    //_ntuple1->SetBranchAddress("e1_e2_SS", &e1_e2_SS1);
    //_ntuple2->SetBranchAddress("e1_e2_SS", &e1_e2_SS2);
    //_ntuple1->SetBranchAddress("e1_e2_Charge", &e1_e2_Charge1);
    //_ntuple2->SetBranchAddress("e1_e2_Charge", &e1_e2_Charge2);
    //_ntuple1->SetBranchAddress("e1_e2_PdgId", &e1_e2_PdgId1);
    //_ntuple2->SetBranchAddress("e1_e2_PdgId", &e1_e2_PdgId2);
    //_ntuple1->SetBranchAddress("e1EffectiveArea", &e1EffectiveArea1);
    //_ntuple2->SetBranchAddress("e1EffectiveArea", &e1EffectiveArea2);
    _ntuple1->SetBranchAddress("e1Energy", &e1Energy1);
    _ntuple2->SetBranchAddress("e1Energy", &e1Energy2);
    _ntuple1->SetBranchAddress("e1Eta", &e1Eta1);
    _ntuple2->SetBranchAddress("e1Eta", &e1Eta2);
    //_ntuple1->SetBranchAddress("e1FSRDREt", &e1FSRDREt1);
    //_ntuple2->SetBranchAddress("e1FSRDREt", &e1FSRDREt2);
    //_ntuple1->SetBranchAddress("e1FSREt", &e1FSREt1);
    //_ntuple2->SetBranchAddress("e1FSREt", &e1FSREt2);
    //_ntuple1->SetBranchAddress("e1FSREta", &e1FSREta1);
    //_ntuple2->SetBranchAddress("e1FSREta", &e1FSREta2);
    //_ntuple1->SetBranchAddress("e1FSRPhi", &e1FSRPhi1);
    //_ntuple2->SetBranchAddress("e1FSRPhi", &e1FSRPhi2);
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e1GenEta", &e1GenEta1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e1GenEta", &e1GenEta2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e1GenPhi", &e1GenPhi1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e1GenPhi", &e1GenPhi2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e1GenPt", &e1GenPt1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e1GenPt", &e1GenPt2);
    }
    //_ntuple1->SetBranchAddress("e1IP2D", &e1IP2D1);
    //_ntuple2->SetBranchAddress("e1IP2D", &e1IP2D2);
    //_ntuple1->SetBranchAddress("e1IP2DUncertainty", &e1IP2DUncertainty1);
    //_ntuple2->SetBranchAddress("e1IP2DUncertainty", &e1IP2DUncertainty2);
    //_ntuple1->SetBranchAddress("e1IP3D", &e1IP3D1);
    //_ntuple2->SetBranchAddress("e1IP3D", &e1IP3D2);
    //_ntuple1->SetBranchAddress("e1IP3DUncertainty", &e1IP3DUncertainty1);
    //_ntuple2->SetBranchAddress("e1IP3DUncertainty", &e1IP3DUncertainty2);
    //_ntuple1->SetBranchAddress("e1Loose", &e1Loose1);
    //_ntuple2->SetBranchAddress("e1Loose", &e1Loose2);
    //_ntuple1->SetBranchAddress("e1Mass", &e1Mass1);
    //_ntuple2->SetBranchAddress("e1Mass", &e1Mass2);
    //_ntuple1->SetBranchAddress("e1Medium", &e1Medium1);
    //_ntuple2->SetBranchAddress("e1Medium", &e1Medium2);
    //_ntuple1->SetBranchAddress("e1MtToMET", &e1MtToMET1);
    //_ntuple2->SetBranchAddress("e1MtToMET", &e1MtToMET2);
    //_ntuple1->SetBranchAddress("e1PFChargedIso", &e1PFChargedIso1);
    //_ntuple2->SetBranchAddress("e1PFChargedIso", &e1PFChargedIso2);
    //_ntuple1->SetBranchAddress("e1PFNeutralIso", &e1PFNeutralIso1);
    //_ntuple2->SetBranchAddress("e1PFNeutralIso", &e1PFNeutralIso2);
    //_ntuple1->SetBranchAddress("e1PFPUIso", &e1PFPUIso1);
    //_ntuple2->SetBranchAddress("e1PFPUIso", &e1PFPUIso2);
    //_ntuple1->SetBranchAddress("e1PFPhotonIso", &e1PFPhotonIso1);
    //_ntuple2->SetBranchAddress("e1PFPhotonIso", &e1PFPhotonIso2);
    //_ntuple1->SetBranchAddress("e1PVDXY", &e1PVDXY1);
    //_ntuple2->SetBranchAddress("e1PVDXY", &e1PVDXY2);
    //_ntuple1->SetBranchAddress("e1PVDZ", &e1PVDZ1);
    //_ntuple2->SetBranchAddress("e1PVDZ", &e1PVDZ2);
    _ntuple1->SetBranchAddress("e1Phi", &e1Phi1);
    _ntuple2->SetBranchAddress("e1Phi", &e1Phi2);
    _ntuple1->SetBranchAddress("e1Pt", &e1Pt1);
    _ntuple2->SetBranchAddress("e1Pt", &e1Pt2);
    //_ntuple1->SetBranchAddress("e1RelPFIsoRho", &e1RelPFIsoRho1);
    //_ntuple2->SetBranchAddress("e1RelPFIsoRho", &e1RelPFIsoRho2);
    //_ntuple1->SetBranchAddress("e1Rho", &e1Rho1);
    //_ntuple2->SetBranchAddress("e1Rho", &e1Rho2);
    //_ntuple1->SetBranchAddress("e1SCEnergy", &e1SCEnergy1);
    //_ntuple2->SetBranchAddress("e1SCEnergy", &e1SCEnergy2);
    //_ntuple1->SetBranchAddress("e1SCEta", &e1SCEta1);
    //_ntuple2->SetBranchAddress("e1SCEta", &e1SCEta2);
    //_ntuple1->SetBranchAddress("e1SCPhi", &e1SCPhi1);
    //_ntuple2->SetBranchAddress("e1SCPhi", &e1SCPhi2);
    //_ntuple1->SetBranchAddress("e1SCRawEnergy", &e1SCRawEnergy1);
    //_ntuple2->SetBranchAddress("e1SCRawEnergy", &e1SCRawEnergy2);
    //_ntuple1->SetBranchAddress("e1SIP2D", &e1SIP2D1);
    //_ntuple2->SetBranchAddress("e1SIP2D", &e1SIP2D2);
    //_ntuple1->SetBranchAddress("e1SIP3D", &e1SIP3D1);
    //_ntuple2->SetBranchAddress("e1SIP3D", &e1SIP3D2);
    //_ntuple1->SetBranchAddress("e1ScaleGainDn", &e1ScaleGainDn1);
    //_ntuple2->SetBranchAddress("e1ScaleGainDn", &e1ScaleGainDn2);
    //_ntuple1->SetBranchAddress("e1ScaleGainUp", &e1ScaleGainUp1);
    //_ntuple2->SetBranchAddress("e1ScaleGainUp", &e1ScaleGainUp2);
    //_ntuple1->SetBranchAddress("e1ScaleStatDn", &e1ScaleStatDn1);
    //_ntuple2->SetBranchAddress("e1ScaleStatDn", &e1ScaleStatDn2);
    //_ntuple1->SetBranchAddress("e1ScaleStatUp", &e1ScaleStatUp1);
    //_ntuple2->SetBranchAddress("e1ScaleStatUp", &e1ScaleStatUp2);
    //_ntuple1->SetBranchAddress("e1ScaleSystDn", &e1ScaleSystDn1);
    //_ntuple2->SetBranchAddress("e1ScaleSystDn", &e1ScaleSystDn2);
    //_ntuple1->SetBranchAddress("e1ScaleSystUp", &e1ScaleSystUp1);
    //_ntuple2->SetBranchAddress("e1ScaleSystUp", &e1ScaleSystUp2);
    //_ntuple1->SetBranchAddress("e1ScaleTotDn", &e1ScaleTotDn1);
    //_ntuple2->SetBranchAddress("e1ScaleTotDn", &e1ScaleTotDn2);
    //_ntuple1->SetBranchAddress("e1ScaleTotUp", &e1ScaleTotUp1);
    //_ntuple2->SetBranchAddress("e1ScaleTotUp", &e1ScaleTotUp2);
    //_ntuple1->SetBranchAddress("e1SigmaPhiDn", &e1SigmaPhiDn1);
    //_ntuple2->SetBranchAddress("e1SigmaPhiDn", &e1SigmaPhiDn2);
    //_ntuple1->SetBranchAddress("e1SigmaPhiUp", &e1SigmaPhiUp1);
    //_ntuple2->SetBranchAddress("e1SigmaPhiUp", &e1SigmaPhiUp2);
    //_ntuple1->SetBranchAddress("e1SigmaRhoDn", &e1SigmaRhoDn1);
    //_ntuple2->SetBranchAddress("e1SigmaRhoDn", &e1SigmaRhoDn2);
    //_ntuple1->SetBranchAddress("e1SigmaRhoUp", &e1SigmaRhoUp1);
    //_ntuple2->SetBranchAddress("e1SigmaRhoUp", &e1SigmaRhoUp2);
    //_ntuple1->SetBranchAddress("e1SigmaTotDn", &e1SigmaTotDn1);
    //_ntuple2->SetBranchAddress("e1SigmaTotDn", &e1SigmaTotDn2);
    //_ntuple1->SetBranchAddress("e1SigmaTotUp", &e1SigmaTotUp1);
    //_ntuple2->SetBranchAddress("e1SigmaTotUp", &e1SigmaTotUp2);
    //_ntuple1->SetBranchAddress("e1Tight", &e1Tight1);
    //_ntuple2->SetBranchAddress("e1Tight", &e1Tight2);
    //_ntuple1->SetBranchAddress("e1UnCorrPt", &e1UnCorrPt1);
    //_ntuple2->SetBranchAddress("e1UnCorrPt", &e1UnCorrPt2);
    //_ntuple1->SetBranchAddress("e1Veto", &e1Veto1);
    //_ntuple2->SetBranchAddress("e1Veto", &e1Veto2);
    //_ntuple1->SetBranchAddress("e1ZZIso", &e1ZZIso1);
    //_ntuple2->SetBranchAddress("e1ZZIso", &e1ZZIso2);
    //_ntuple1->SetBranchAddress("e1IsEB", &e1IsEB1);
    //_ntuple2->SetBranchAddress("e1IsEB", &e1IsEB2);
    //_ntuple1->SetBranchAddress("e1IsGap", &e1IsGap1);
    //_ntuple2->SetBranchAddress("e1IsGap", &e1IsGap2);
    //_ntuple1->SetBranchAddress("e1ZZIsoPass", &e1ZZIsoPass1);
    //_ntuple2->SetBranchAddress("e1ZZIsoPass", &e1ZZIsoPass2);
    //_ntuple1->SetBranchAddress("e1ZZLooseID", &e1ZZLooseID1);
    //_ntuple2->SetBranchAddress("e1ZZLooseID", &e1ZZLooseID2);
    //_ntuple1->SetBranchAddress("e1ZZLooseIDNoVtx", &e1ZZLooseIDNoVtx1);
    //_ntuple2->SetBranchAddress("e1ZZLooseIDNoVtx", &e1ZZLooseIDNoVtx2);
    //_ntuple1->SetBranchAddress("e1ZZTightID", &e1ZZTightID1);
    //_ntuple2->SetBranchAddress("e1ZZTightID", &e1ZZTightID2);
    //_ntuple1->SetBranchAddress("e1ZZTightIDNoVtx", &e1ZZTightIDNoVtx1);
    //_ntuple2->SetBranchAddress("e1ZZTightIDNoVtx", &e1ZZTightIDNoVtx2);
    //_ntuple1->SetBranchAddress("e1Charge", &e1Charge1);
    //_ntuple2->SetBranchAddress("e1Charge", &e1Charge2);
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e1GenCharge", &e1GenCharge1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e1GenCharge", &e1GenCharge2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e1GenPdgId", &e1GenPdgId1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e1GenPdgId", &e1GenPdgId2);
    }
    _ntuple1->SetBranchAddress("e1PdgId", &e1PdgId1);
    _ntuple2->SetBranchAddress("e1PdgId", &e1PdgId2);
    //_ntuple1->SetBranchAddress("e1MissingHits", &e1MissingHits1);
    //_ntuple2->SetBranchAddress("e1MissingHits", &e1MissingHits2);
    //_ntuple1->SetBranchAddress("e2EffectiveArea", &e2EffectiveArea1);
    //_ntuple2->SetBranchAddress("e2EffectiveArea", &e2EffectiveArea2);
    _ntuple1->SetBranchAddress("e2Energy", &e2Energy1);
    _ntuple2->SetBranchAddress("e2Energy", &e2Energy2);
    _ntuple1->SetBranchAddress("e2Eta", &e2Eta1);
    _ntuple2->SetBranchAddress("e2Eta", &e2Eta2);
    //_ntuple1->SetBranchAddress("e2FSRDREt", &e2FSRDREt1);
    //_ntuple2->SetBranchAddress("e2FSRDREt", &e2FSRDREt2);
    //_ntuple1->SetBranchAddress("e2FSREt", &e2FSREt1);
    //_ntuple2->SetBranchAddress("e2FSREt", &e2FSREt2);
    //_ntuple1->SetBranchAddress("e2FSREta", &e2FSREta1);
    //_ntuple2->SetBranchAddress("e2FSREta", &e2FSREta2);
    //_ntuple1->SetBranchAddress("e2FSRPhi", &e2FSRPhi1);
    //_ntuple2->SetBranchAddress("e2FSRPhi", &e2FSRPhi2);
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e2GenEta", &e2GenEta1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e2GenEta", &e2GenEta2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e2GenPhi", &e2GenPhi1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e2GenPhi", &e2GenPhi2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e2GenPt", &e2GenPt1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e2GenPt", &e2GenPt2);
    }
    //_ntuple1->SetBranchAddress("e2IP2D", &e2IP2D1);
    //_ntuple2->SetBranchAddress("e2IP2D", &e2IP2D2);
    //_ntuple1->SetBranchAddress("e2IP2DUncertainty", &e2IP2DUncertainty1);
    //_ntuple2->SetBranchAddress("e2IP2DUncertainty", &e2IP2DUncertainty2);
    //_ntuple1->SetBranchAddress("e2IP3D", &e2IP3D1);
    //_ntuple2->SetBranchAddress("e2IP3D", &e2IP3D2);
    //_ntuple1->SetBranchAddress("e2IP3DUncertainty", &e2IP3DUncertainty1);
    //_ntuple2->SetBranchAddress("e2IP3DUncertainty", &e2IP3DUncertainty2);
    //_ntuple1->SetBranchAddress("e2Loose", &e2Loose1);
    //_ntuple2->SetBranchAddress("e2Loose", &e2Loose2);
    //_ntuple1->SetBranchAddress("e2Mass", &e2Mass1);
    //_ntuple2->SetBranchAddress("e2Mass", &e2Mass2);
    //_ntuple1->SetBranchAddress("e2Medium", &e2Medium1);
    //_ntuple2->SetBranchAddress("e2Medium", &e2Medium2);
    //_ntuple1->SetBranchAddress("e2MtToMET", &e2MtToMET1);
    //_ntuple2->SetBranchAddress("e2MtToMET", &e2MtToMET2);
    //_ntuple1->SetBranchAddress("e2PFChargedIso", &e2PFChargedIso1);
    //_ntuple2->SetBranchAddress("e2PFChargedIso", &e2PFChargedIso2);
    //_ntuple1->SetBranchAddress("e2PFNeutralIso", &e2PFNeutralIso1);
    //_ntuple2->SetBranchAddress("e2PFNeutralIso", &e2PFNeutralIso2);
    //_ntuple1->SetBranchAddress("e2PFPUIso", &e2PFPUIso1);
    //_ntuple2->SetBranchAddress("e2PFPUIso", &e2PFPUIso2);
    //_ntuple1->SetBranchAddress("e2PFPhotonIso", &e2PFPhotonIso1);
    //_ntuple2->SetBranchAddress("e2PFPhotonIso", &e2PFPhotonIso2);
    //_ntuple1->SetBranchAddress("e2PVDXY", &e2PVDXY1);
    //_ntuple2->SetBranchAddress("e2PVDXY", &e2PVDXY2);
    //_ntuple1->SetBranchAddress("e2PVDZ", &e2PVDZ1);
    //_ntuple2->SetBranchAddress("e2PVDZ", &e2PVDZ2);
    _ntuple1->SetBranchAddress("e2Phi", &e2Phi1);
    _ntuple2->SetBranchAddress("e2Phi", &e2Phi2);
    _ntuple1->SetBranchAddress("e2Pt", &e2Pt1);
    _ntuple2->SetBranchAddress("e2Pt", &e2Pt2);
    //_ntuple1->SetBranchAddress("e2RelPFIsoRho", &e2RelPFIsoRho1);
    //_ntuple2->SetBranchAddress("e2RelPFIsoRho", &e2RelPFIsoRho2);
    //_ntuple1->SetBranchAddress("e2Rho", &e2Rho1);
    //_ntuple2->SetBranchAddress("e2Rho", &e2Rho2);
    //_ntuple1->SetBranchAddress("e2SCEnergy", &e2SCEnergy1);
    //_ntuple2->SetBranchAddress("e2SCEnergy", &e2SCEnergy2);
    //_ntuple1->SetBranchAddress("e2SCEta", &e2SCEta1);
    //_ntuple2->SetBranchAddress("e2SCEta", &e2SCEta2);
    //_ntuple1->SetBranchAddress("e2SCPhi", &e2SCPhi1);
    //_ntuple2->SetBranchAddress("e2SCPhi", &e2SCPhi2);
    //_ntuple1->SetBranchAddress("e2SCRawEnergy", &e2SCRawEnergy1);
    //_ntuple2->SetBranchAddress("e2SCRawEnergy", &e2SCRawEnergy2);
    //_ntuple1->SetBranchAddress("e2SIP2D", &e2SIP2D1);
    //_ntuple2->SetBranchAddress("e2SIP2D", &e2SIP2D2);
    //_ntuple1->SetBranchAddress("e2SIP3D", &e2SIP3D1);
    //_ntuple2->SetBranchAddress("e2SIP3D", &e2SIP3D2);
    //_ntuple1->SetBranchAddress("e2ScaleGainDn", &e2ScaleGainDn1);
    //_ntuple2->SetBranchAddress("e2ScaleGainDn", &e2ScaleGainDn2);
    //_ntuple1->SetBranchAddress("e2ScaleGainUp", &e2ScaleGainUp1);
    //_ntuple2->SetBranchAddress("e2ScaleGainUp", &e2ScaleGainUp2);
    //_ntuple1->SetBranchAddress("e2ScaleStatDn", &e2ScaleStatDn1);
    //_ntuple2->SetBranchAddress("e2ScaleStatDn", &e2ScaleStatDn2);
    //_ntuple1->SetBranchAddress("e2ScaleStatUp", &e2ScaleStatUp1);
    //_ntuple2->SetBranchAddress("e2ScaleStatUp", &e2ScaleStatUp2);
    //_ntuple1->SetBranchAddress("e2ScaleSystDn", &e2ScaleSystDn1);
    //_ntuple2->SetBranchAddress("e2ScaleSystDn", &e2ScaleSystDn2);
    //_ntuple1->SetBranchAddress("e2ScaleSystUp", &e2ScaleSystUp1);
    //_ntuple2->SetBranchAddress("e2ScaleSystUp", &e2ScaleSystUp2);
    //_ntuple1->SetBranchAddress("e2ScaleTotDn", &e2ScaleTotDn1);
    //_ntuple2->SetBranchAddress("e2ScaleTotDn", &e2ScaleTotDn2);
    //_ntuple1->SetBranchAddress("e2ScaleTotUp", &e2ScaleTotUp1);
    //_ntuple2->SetBranchAddress("e2ScaleTotUp", &e2ScaleTotUp2);
    //_ntuple1->SetBranchAddress("e2SigmaPhiDn", &e2SigmaPhiDn1);
    //_ntuple2->SetBranchAddress("e2SigmaPhiDn", &e2SigmaPhiDn2);
    //_ntuple1->SetBranchAddress("e2SigmaPhiUp", &e2SigmaPhiUp1);
    //_ntuple2->SetBranchAddress("e2SigmaPhiUp", &e2SigmaPhiUp2);
    //_ntuple1->SetBranchAddress("e2SigmaRhoDn", &e2SigmaRhoDn1);
    //_ntuple2->SetBranchAddress("e2SigmaRhoDn", &e2SigmaRhoDn2);
    //_ntuple1->SetBranchAddress("e2SigmaRhoUp", &e2SigmaRhoUp1);
    //_ntuple2->SetBranchAddress("e2SigmaRhoUp", &e2SigmaRhoUp2);
    //_ntuple1->SetBranchAddress("e2SigmaTotDn", &e2SigmaTotDn1);
    //_ntuple2->SetBranchAddress("e2SigmaTotDn", &e2SigmaTotDn2);
    //_ntuple1->SetBranchAddress("e2SigmaTotUp", &e2SigmaTotUp1);
    //_ntuple2->SetBranchAddress("e2SigmaTotUp", &e2SigmaTotUp2);
    //_ntuple1->SetBranchAddress("e2Tight", &e2Tight1);
    //_ntuple2->SetBranchAddress("e2Tight", &e2Tight2);
    //_ntuple1->SetBranchAddress("e2UnCorrPt", &e2UnCorrPt1);
    //_ntuple2->SetBranchAddress("e2UnCorrPt", &e2UnCorrPt2);
    //_ntuple1->SetBranchAddress("e2Veto", &e2Veto1);
    //_ntuple2->SetBranchAddress("e2Veto", &e2Veto2);
    //_ntuple1->SetBranchAddress("e2ZZIso", &e2ZZIso1);
    //_ntuple2->SetBranchAddress("e2ZZIso", &e2ZZIso2);
    //_ntuple1->SetBranchAddress("e2IsEB", &e2IsEB1);
    //_ntuple2->SetBranchAddress("e2IsEB", &e2IsEB2);
    //_ntuple1->SetBranchAddress("e2IsGap", &e2IsGap1);
    //_ntuple2->SetBranchAddress("e2IsGap", &e2IsGap2);
    //_ntuple1->SetBranchAddress("e2ZZIsoPass", &e2ZZIsoPass1);
    //_ntuple2->SetBranchAddress("e2ZZIsoPass", &e2ZZIsoPass2);
    //_ntuple1->SetBranchAddress("e2ZZLooseID", &e2ZZLooseID1);
    //_ntuple2->SetBranchAddress("e2ZZLooseID", &e2ZZLooseID2);
    //_ntuple1->SetBranchAddress("e2ZZLooseIDNoVtx", &e2ZZLooseIDNoVtx1);
    //_ntuple2->SetBranchAddress("e2ZZLooseIDNoVtx", &e2ZZLooseIDNoVtx2);
    //_ntuple1->SetBranchAddress("e2ZZTightID", &e2ZZTightID1);
    //_ntuple2->SetBranchAddress("e2ZZTightID", &e2ZZTightID2);
    //_ntuple1->SetBranchAddress("e2ZZTightIDNoVtx", &e2ZZTightIDNoVtx1);
    //_ntuple2->SetBranchAddress("e2ZZTightIDNoVtx", &e2ZZTightIDNoVtx2);
    //_ntuple1->SetBranchAddress("e2Charge", &e2Charge1);
    //_ntuple2->SetBranchAddress("e2Charge", &e2Charge2);
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e2GenCharge", &e2GenCharge1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e2GenCharge", &e2GenCharge2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("e2GenPdgId", &e2GenPdgId1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("e2GenPdgId", &e2GenPdgId2);
    }
    _ntuple1->SetBranchAddress("e2PdgId", &e2PdgId1);
    _ntuple2->SetBranchAddress("e2PdgId", &e2PdgId2);
    //_ntuple1->SetBranchAddress("e2MissingHits", &e2MissingHits1);
    //_ntuple2->SetBranchAddress("e2MissingHits", &e2MissingHits2);
  }

  if (_channel == "mmmm" || _channel == "eemm"){
    //_ntuple1->SetBranchAddress("m1_m2_DR", &m1_m2_DR1);
    //_ntuple2->SetBranchAddress("m1_m2_DR", &m1_m2_DR2);
    //_ntuple1->SetBranchAddress("m1_m2_Energy", &m1_m2_Energy1);
    //_ntuple2->SetBranchAddress("m1_m2_Energy", &m1_m2_Energy2);
    //_ntuple1->SetBranchAddress("m1_m2_EnergyNoFSR", &m1_m2_EnergyNoFSR1);
    //_ntuple2->SetBranchAddress("m1_m2_EnergyNoFSR", &m1_m2_EnergyNoFSR2);
    //_ntuple1->SetBranchAddress("m1_m2_Eta", &m1_m2_Eta1);
    //_ntuple2->SetBranchAddress("m1_m2_Eta", &m1_m2_Eta2);
    //_ntuple1->SetBranchAddress("m1_m2_EtaNoFSR", &m1_m2_EtaNoFSR1);
    //_ntuple2->SetBranchAddress("m1_m2_EtaNoFSR", &m1_m2_EtaNoFSR2);
    _ntuple1->SetBranchAddress("m1_m2_Mass", &m1_m2_Mass1);
    _ntuple2->SetBranchAddress("m1_m2_Mass", &m1_m2_Mass2);
    //_ntuple1->SetBranchAddress("m1_m2_MassNoFSR", &m1_m2_MassNoFSR1);
    //_ntuple2->SetBranchAddress("m1_m2_MassNoFSR", &m1_m2_MassNoFSR2);
    //_ntuple1->SetBranchAddress("m1_m2_Mt", &m1_m2_Mt1);
    //_ntuple2->SetBranchAddress("m1_m2_Mt", &m1_m2_Mt2);
    //_ntuple1->SetBranchAddress("m1_m2_MtToMET", &m1_m2_MtToMET1);
    //_ntuple2->SetBranchAddress("m1_m2_MtToMET", &m1_m2_MtToMET2);
    //_ntuple1->SetBranchAddress("m1_m2_Phi", &m1_m2_Phi1);
    //_ntuple2->SetBranchAddress("m1_m2_Phi", &m1_m2_Phi2);
    //_ntuple1->SetBranchAddress("m1_m2_PhiNoFSR", &m1_m2_PhiNoFSR1);
    //_ntuple2->SetBranchAddress("m1_m2_PhiNoFSR", &m1_m2_PhiNoFSR2);
    //_ntuple1->SetBranchAddress("m1_m2_Pt", &m1_m2_Pt1);
    //_ntuple2->SetBranchAddress("m1_m2_Pt", &m1_m2_Pt2);
    //_ntuple1->SetBranchAddress("m1_m2_PtNoFSR", &m1_m2_PtNoFSR1);
    //_ntuple2->SetBranchAddress("m1_m2_PtNoFSR", &m1_m2_PtNoFSR2);
    //_ntuple1->SetBranchAddress("m1_m2_SS", &m1_m2_SS1);
    //_ntuple2->SetBranchAddress("m1_m2_SS", &m1_m2_SS2);
    //_ntuple1->SetBranchAddress("m1_m2_Charge", &m1_m2_Charge1);
    //_ntuple2->SetBranchAddress("m1_m2_Charge", &m1_m2_Charge2);
    //_ntuple1->SetBranchAddress("m1_m2_PdgId", &m1_m2_PdgId1);
    //_ntuple2->SetBranchAddress("m1_m2_PdgId", &m1_m2_PdgId2);
    //_ntuple1->SetBranchAddress("m1EffScaleFactor", &m1EffScaleFactor1);
    //_ntuple2->SetBranchAddress("m1EffScaleFactor", &m1EffScaleFactor2);
    //_ntuple1->SetBranchAddress("m1EffScaleFactorError", &m1EffScaleFactorError1);
    //_ntuple2->SetBranchAddress("m1EffScaleFactorError", &m1EffScaleFactorError2);
    _ntuple1->SetBranchAddress("m1Energy", &m1Energy1);
    _ntuple2->SetBranchAddress("m1Energy", &m1Energy2);
    _ntuple1->SetBranchAddress("m1Eta", &m1Eta1);
    _ntuple2->SetBranchAddress("m1Eta", &m1Eta2);
    //_ntuple1->SetBranchAddress("m1FSRDREt", &m1FSRDREt1);
    //_ntuple2->SetBranchAddress("m1FSRDREt", &m1FSRDREt2);
    //_ntuple1->SetBranchAddress("m1FSREt", &m1FSREt1);
    //_ntuple2->SetBranchAddress("m1FSREt", &m1FSREt2);
    //_ntuple1->SetBranchAddress("m1FSREta", &m1FSREta1);
    //_ntuple2->SetBranchAddress("m1FSREta", &m1FSREta2);
    //_ntuple1->SetBranchAddress("m1FSRPhi", &m1FSRPhi1);
    //_ntuple2->SetBranchAddress("m1FSRPhi", &m1FSRPhi2);
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m1GenEta", &m1GenEta1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m1GenEta", &m1GenEta2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m1GenPhi", &m1GenPhi1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m1GenPhi", &m1GenPhi2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m1GenPt", &m1GenPt1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m1GenPt", &m1GenPt2);
    }
    //_ntuple1->SetBranchAddress("m1IP2D", &m1IP2D1);
    //_ntuple2->SetBranchAddress("m1IP2D", &m1IP2D2);
    //_ntuple1->SetBranchAddress("m1IP2DUncertainty", &m1IP2DUncertainty1);
    //_ntuple2->SetBranchAddress("m1IP2DUncertainty", &m1IP2DUncertainty2);
    //_ntuple1->SetBranchAddress("m1IP3D", &m1IP3D1);
    //_ntuple2->SetBranchAddress("m1IP3D", &m1IP3D2);
    //_ntuple1->SetBranchAddress("m1IP3DUncertainty", &m1IP3DUncertainty1);
    //_ntuple2->SetBranchAddress("m1IP3DUncertainty", &m1IP3DUncertainty2);
    //_ntuple1->SetBranchAddress("m1Mass", &m1Mass1);
    //_ntuple2->SetBranchAddress("m1Mass", &m1Mass2);
    //_ntuple1->SetBranchAddress("m1MtToMET", &m1MtToMET1);
    //_ntuple2->SetBranchAddress("m1MtToMET", &m1MtToMET2);
    //_ntuple1->SetBranchAddress("m1PFChargedIso", &m1PFChargedIso1);
    //_ntuple2->SetBranchAddress("m1PFChargedIso", &m1PFChargedIso2);
    //_ntuple1->SetBranchAddress("m1PFNeutralIso", &m1PFNeutralIso1);
    //_ntuple2->SetBranchAddress("m1PFNeutralIso", &m1PFNeutralIso2);
    //_ntuple1->SetBranchAddress("m1PFPUIso", &m1PFPUIso1);
    //_ntuple2->SetBranchAddress("m1PFPUIso", &m1PFPUIso2);
    //_ntuple1->SetBranchAddress("m1PFPhotonIso", &m1PFPhotonIso1);
    //_ntuple2->SetBranchAddress("m1PFPhotonIso", &m1PFPhotonIso2);
    //_ntuple1->SetBranchAddress("m1PVDXY", &m1PVDXY1);
    //_ntuple2->SetBranchAddress("m1PVDXY", &m1PVDXY2);
    //_ntuple1->SetBranchAddress("m1PVDZ", &m1PVDZ1);
    //_ntuple2->SetBranchAddress("m1PVDZ", &m1PVDZ2);
    _ntuple1->SetBranchAddress("m1Phi", &m1Phi1);
    _ntuple2->SetBranchAddress("m1Phi", &m1Phi2);
    _ntuple1->SetBranchAddress("m1Pt", &m1Pt1);
    _ntuple2->SetBranchAddress("m1Pt", &m1Pt2);
    //_ntuple1->SetBranchAddress("m1PtErr", &m1PtErr1);
    //_ntuple2->SetBranchAddress("m1PtErr", &m1PtErr2);
    //_ntuple1->SetBranchAddress("m1PtErrUncorrected", &m1PtErrUncorrected1);
    //_ntuple2->SetBranchAddress("m1PtErrUncorrected", &m1PtErrUncorrected2);
    //_ntuple1->SetBranchAddress("m1PtUncorrected", &m1PtUncorrected1);
    //_ntuple2->SetBranchAddress("m1PtUncorrected", &m1PtUncorrected2);
    //_ntuple1->SetBranchAddress("m1RelPFIsoDBR03", &m1RelPFIsoDBR031);
    //_ntuple2->SetBranchAddress("m1RelPFIsoDBR03", &m1RelPFIsoDBR032);
    //_ntuple1->SetBranchAddress("m1RelPFIsoDBR04", &m1RelPFIsoDBR041);
    //_ntuple2->SetBranchAddress("m1RelPFIsoDBR04", &m1RelPFIsoDBR042);
    //_ntuple1->SetBranchAddress("m1RochesterScaleUncDn", &m1RochesterScaleUncDn1);
    //_ntuple2->SetBranchAddress("m1RochesterScaleUncDn", &m1RochesterScaleUncDn2);
    //_ntuple1->SetBranchAddress("m1RochesterScaleUncUp", &m1RochesterScaleUncUp1);
    //_ntuple2->SetBranchAddress("m1RochesterScaleUncUp", &m1RochesterScaleUncUp2);
    //_ntuple1->SetBranchAddress("m1RochesterSmearUncDn", &m1RochesterSmearUncDn1);
    //_ntuple2->SetBranchAddress("m1RochesterSmearUncDn", &m1RochesterSmearUncDn2);
    //_ntuple1->SetBranchAddress("m1RochesterSmearUncUp", &m1RochesterSmearUncUp1);
    //_ntuple2->SetBranchAddress("m1RochesterSmearUncUp", &m1RochesterSmearUncUp2);
    //_ntuple1->SetBranchAddress("m1SIP2D", &m1SIP2D1);
    //_ntuple2->SetBranchAddress("m1SIP2D", &m1SIP2D2);
    //_ntuple1->SetBranchAddress("m1SIP3D", &m1SIP3D1);
    //_ntuple2->SetBranchAddress("m1SIP3D", &m1SIP3D2);
    //_ntuple1->SetBranchAddress("m1TrackIso", &m1TrackIso1);
    //_ntuple2->SetBranchAddress("m1TrackIso", &m1TrackIso2);
    //_ntuple1->SetBranchAddress("m1ZZIso", &m1ZZIso1);
    //_ntuple2->SetBranchAddress("m1ZZIso", &m1ZZIso2);
    //_ntuple1->SetBranchAddress("m1HighPtID", &m1HighPtID1);
    //_ntuple2->SetBranchAddress("m1HighPtID", &m1HighPtID2);
    //_ntuple1->SetBranchAddress("m1HighPtIDNoVtx", &m1HighPtIDNoVtx1);
    //_ntuple2->SetBranchAddress("m1HighPtIDNoVtx", &m1HighPtIDNoVtx2);
    //_ntuple1->SetBranchAddress("m1IsGlobal", &m1IsGlobal1);
    //_ntuple2->SetBranchAddress("m1IsGlobal", &m1IsGlobal2);
    //_ntuple1->SetBranchAddress("m1IsLoose", &m1IsLoose1);
    //_ntuple2->SetBranchAddress("m1IsLoose", &m1IsLoose2);
    //_ntuple1->SetBranchAddress("m1IsMedium", &m1IsMedium1);
    //_ntuple2->SetBranchAddress("m1IsMedium", &m1IsMedium2);
    //_ntuple1->SetBranchAddress("m1IsPFMuon", &m1IsPFMuon1);
    //_ntuple2->SetBranchAddress("m1IsPFMuon", &m1IsPFMuon2);
    //_ntuple1->SetBranchAddress("m1IsTracker", &m1IsTracker1);
    //_ntuple2->SetBranchAddress("m1IsTracker", &m1IsTracker2);
    //_ntuple1->SetBranchAddress("m1PASTightID", &m1PASTightID1);
    //_ntuple2->SetBranchAddress("m1PASTightID", &m1PASTightID2);
    //_ntuple1->SetBranchAddress("m1PASTightIDNoVtx", &m1PASTightIDNoVtx1);
    //_ntuple2->SetBranchAddress("m1PASTightIDNoVtx", &m1PASTightIDNoVtx2);
    //_ntuple1->SetBranchAddress("m1PFID", &m1PFID1);
    //_ntuple2->SetBranchAddress("m1PFID", &m1PFID2);
    //_ntuple1->SetBranchAddress("m1PFIDNoVtx", &m1PFIDNoVtx1);
    //_ntuple2->SetBranchAddress("m1PFIDNoVtx", &m1PFIDNoVtx2);
    //_ntuple1->SetBranchAddress("m1ZZIsoPass", &m1ZZIsoPass1);
    //_ntuple2->SetBranchAddress("m1ZZIsoPass", &m1ZZIsoPass2);
    //_ntuple1->SetBranchAddress("m1ZZLooseID", &m1ZZLooseID1);
    //_ntuple2->SetBranchAddress("m1ZZLooseID", &m1ZZLooseID2);
    //_ntuple1->SetBranchAddress("m1ZZLooseIDNoVtx", &m1ZZLooseIDNoVtx1);
    //_ntuple2->SetBranchAddress("m1ZZLooseIDNoVtx", &m1ZZLooseIDNoVtx2);
    //_ntuple1->SetBranchAddress("m1ZZTightID", &m1ZZTightID1);
    //_ntuple2->SetBranchAddress("m1ZZTightID", &m1ZZTightID2);
    //_ntuple1->SetBranchAddress("m1ZZTightIDNoVtx", &m1ZZTightIDNoVtx1);
    //_ntuple2->SetBranchAddress("m1ZZTightIDNoVtx", &m1ZZTightIDNoVtx2);
    //_ntuple1->SetBranchAddress("m1Charge", &m1Charge1);
    //_ntuple2->SetBranchAddress("m1Charge", &m1Charge2);
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m1GenCharge", &m1GenCharge1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m1GenCharge", &m1GenCharge2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m1GenPdgId", &m1GenPdgId1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m1GenPdgId", &m1GenPdgId2);
    }
    _ntuple1->SetBranchAddress("m1PdgId", &m1PdgId1);
    _ntuple2->SetBranchAddress("m1PdgId", &m1PdgId2);
    //_ntuple1->SetBranchAddress("m1BestTrackType", &m1BestTrackType1);
    //_ntuple2->SetBranchAddress("m1BestTrackType", &m1BestTrackType2);
    //_ntuple1->SetBranchAddress("m1CutBasedLoose", &m1CutBasedLoose1);
    //_ntuple2->SetBranchAddress("m1CutBasedLoose", &m1CutBasedLoose2);
    //_ntuple1->SetBranchAddress("m1CutBasedMedium", &m1CutBasedMedium1);
    //_ntuple2->SetBranchAddress("m1CutBasedMedium", &m1CutBasedMedium2);
    //_ntuple1->SetBranchAddress("m1CutBasedTight", &m1CutBasedTight1);
    //_ntuple2->SetBranchAddress("m1CutBasedTight", &m1CutBasedTight2);
    //_ntuple1->SetBranchAddress("m1IsTight", &m1IsTight1);
    //_ntuple2->SetBranchAddress("m1IsTight", &m1IsTight2);
    //_ntuple1->SetBranchAddress("m1MatchedStations", &m1MatchedStations1);
    //_ntuple2->SetBranchAddress("m1MatchedStations", &m1MatchedStations2);
    //_ntuple1->SetBranchAddress("m1NoOfMatches", &m1NoOfMatches1);
    //_ntuple2->SetBranchAddress("m1NoOfMatches", &m1NoOfMatches2);
    //_ntuple1->SetBranchAddress("m1PFIsoLoose", &m1PFIsoLoose1);
    //_ntuple2->SetBranchAddress("m1PFIsoLoose", &m1PFIsoLoose2);
    //_ntuple1->SetBranchAddress("m1PFIsoMedium", &m1PFIsoMedium1);
    //_ntuple2->SetBranchAddress("m1PFIsoMedium", &m1PFIsoMedium2);
    //_ntuple1->SetBranchAddress("m1PFIsoTight", &m1PFIsoTight1);
    //_ntuple2->SetBranchAddress("m1PFIsoTight", &m1PFIsoTight2);
    //_ntuple1->SetBranchAddress("m1PFIsoVTight", &m1PFIsoVTight1);
    //_ntuple2->SetBranchAddress("m1PFIsoVTight", &m1PFIsoVTight2);
    //_ntuple1->SetBranchAddress("m2EffScaleFactor", &m2EffScaleFactor1);
    //_ntuple2->SetBranchAddress("m2EffScaleFactor", &m2EffScaleFactor2);
    //_ntuple1->SetBranchAddress("m2EffScaleFactorError", &m2EffScaleFactorError1);
    //_ntuple2->SetBranchAddress("m2EffScaleFactorError", &m2EffScaleFactorError2);
    _ntuple1->SetBranchAddress("m2Energy", &m2Energy1);
    _ntuple2->SetBranchAddress("m2Energy", &m2Energy2);
    _ntuple1->SetBranchAddress("m2Eta", &m2Eta1);
    _ntuple2->SetBranchAddress("m2Eta", &m2Eta2);
    //_ntuple1->SetBranchAddress("m2FSRDREt", &m2FSRDREt1);
    //_ntuple2->SetBranchAddress("m2FSRDREt", &m2FSRDREt2);
    //_ntuple1->SetBranchAddress("m2FSREt", &m2FSREt1);
    //_ntuple2->SetBranchAddress("m2FSREt", &m2FSREt2);
    //_ntuple1->SetBranchAddress("m2FSREta", &m2FSREta1);
    //_ntuple2->SetBranchAddress("m2FSREta", &m2FSREta2);
    //_ntuple1->SetBranchAddress("m2FSRPhi", &m2FSRPhi1);
    //_ntuple2->SetBranchAddress("m2FSRPhi", &m2FSRPhi2);
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m2GenEta", &m2GenEta1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m2GenEta", &m2GenEta2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m2GenPhi", &m2GenPhi1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m2GenPhi", &m2GenPhi2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m2GenPt", &m2GenPt1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m2GenPt", &m2GenPt2);
    }
    //_ntuple1->SetBranchAddress("m2IP2D", &m2IP2D1);
    //_ntuple2->SetBranchAddress("m2IP2D", &m2IP2D2);
    //_ntuple1->SetBranchAddress("m2IP2DUncertainty", &m2IP2DUncertainty1);
    //_ntuple2->SetBranchAddress("m2IP2DUncertainty", &m2IP2DUncertainty2);
    //_ntuple1->SetBranchAddress("m2IP3D", &m2IP3D1);
    //_ntuple2->SetBranchAddress("m2IP3D", &m2IP3D2);
    //_ntuple1->SetBranchAddress("m2IP3DUncertainty", &m2IP3DUncertainty1);
    //_ntuple2->SetBranchAddress("m2IP3DUncertainty", &m2IP3DUncertainty2);
    //_ntuple1->SetBranchAddress("m2Mass", &m2Mass1);
    //_ntuple2->SetBranchAddress("m2Mass", &m2Mass2);
    //_ntuple1->SetBranchAddress("m2MtToMET", &m2MtToMET1);
    //_ntuple2->SetBranchAddress("m2MtToMET", &m2MtToMET2);
    //_ntuple1->SetBranchAddress("m2PFChargedIso", &m2PFChargedIso1);
    //_ntuple2->SetBranchAddress("m2PFChargedIso", &m2PFChargedIso2);
    //_ntuple1->SetBranchAddress("m2PFNeutralIso", &m2PFNeutralIso1);
    //_ntuple2->SetBranchAddress("m2PFNeutralIso", &m2PFNeutralIso2);
    //_ntuple1->SetBranchAddress("m2PFPUIso", &m2PFPUIso1);
    //_ntuple2->SetBranchAddress("m2PFPUIso", &m2PFPUIso2);
    //_ntuple1->SetBranchAddress("m2PFPhotonIso", &m2PFPhotonIso1);
    //_ntuple2->SetBranchAddress("m2PFPhotonIso", &m2PFPhotonIso2);
    //_ntuple1->SetBranchAddress("m2PVDXY", &m2PVDXY1);
    //_ntuple2->SetBranchAddress("m2PVDXY", &m2PVDXY2);
    //_ntuple1->SetBranchAddress("m2PVDZ", &m2PVDZ1);
    //_ntuple2->SetBranchAddress("m2PVDZ", &m2PVDZ2);
    _ntuple1->SetBranchAddress("m2Phi", &m2Phi1);
    _ntuple2->SetBranchAddress("m2Phi", &m2Phi2);
    _ntuple1->SetBranchAddress("m2Pt", &m2Pt1);
    _ntuple2->SetBranchAddress("m2Pt", &m2Pt2);
    //_ntuple1->SetBranchAddress("m2PtErr", &m2PtErr1);
    //_ntuple2->SetBranchAddress("m2PtErr", &m2PtErr2);
    //_ntuple1->SetBranchAddress("m2PtErrUncorrected", &m2PtErrUncorrected1);
    //_ntuple2->SetBranchAddress("m2PtErrUncorrected", &m2PtErrUncorrected2);
    //_ntuple1->SetBranchAddress("m2PtUncorrected", &m2PtUncorrected1);
    //_ntuple2->SetBranchAddress("m2PtUncorrected", &m2PtUncorrected2);
    //_ntuple1->SetBranchAddress("m2RelPFIsoDBR03", &m2RelPFIsoDBR031);
    //_ntuple2->SetBranchAddress("m2RelPFIsoDBR03", &m2RelPFIsoDBR032);
    //_ntuple1->SetBranchAddress("m2RelPFIsoDBR04", &m2RelPFIsoDBR041);
    //_ntuple2->SetBranchAddress("m2RelPFIsoDBR04", &m2RelPFIsoDBR042);
    //_ntuple1->SetBranchAddress("m2RochesterScaleUncDn", &m2RochesterScaleUncDn1);
    //_ntuple2->SetBranchAddress("m2RochesterScaleUncDn", &m2RochesterScaleUncDn2);
    //_ntuple1->SetBranchAddress("m2RochesterScaleUncUp", &m2RochesterScaleUncUp1);
    //_ntuple2->SetBranchAddress("m2RochesterScaleUncUp", &m2RochesterScaleUncUp2);
    //_ntuple1->SetBranchAddress("m2RochesterSmearUncDn", &m2RochesterSmearUncDn1);
    //_ntuple2->SetBranchAddress("m2RochesterSmearUncDn", &m2RochesterSmearUncDn2);
    //_ntuple1->SetBranchAddress("m2RochesterSmearUncUp", &m2RochesterSmearUncUp1);
    //_ntuple2->SetBranchAddress("m2RochesterSmearUncUp", &m2RochesterSmearUncUp2);
    //_ntuple1->SetBranchAddress("m2SIP2D", &m2SIP2D1);
    //_ntuple2->SetBranchAddress("m2SIP2D", &m2SIP2D2);
    //_ntuple1->SetBranchAddress("m2SIP3D", &m2SIP3D1);
    //_ntuple2->SetBranchAddress("m2SIP3D", &m2SIP3D2);
    //_ntuple1->SetBranchAddress("m2TrackIso", &m2TrackIso1);
    //_ntuple2->SetBranchAddress("m2TrackIso", &m2TrackIso2);
    //_ntuple1->SetBranchAddress("m2ZZIso", &m2ZZIso1);
    //_ntuple2->SetBranchAddress("m2ZZIso", &m2ZZIso2);
    //_ntuple1->SetBranchAddress("m2HighPtID", &m2HighPtID1);
    //_ntuple2->SetBranchAddress("m2HighPtID", &m2HighPtID2);
    //_ntuple1->SetBranchAddress("m2HighPtIDNoVtx", &m2HighPtIDNoVtx1);
    //_ntuple2->SetBranchAddress("m2HighPtIDNoVtx", &m2HighPtIDNoVtx2);
    //_ntuple1->SetBranchAddress("m2IsGlobal", &m2IsGlobal1);
    //_ntuple2->SetBranchAddress("m2IsGlobal", &m2IsGlobal2);
    //_ntuple1->SetBranchAddress("m2IsLoose", &m2IsLoose1);
    //_ntuple2->SetBranchAddress("m2IsLoose", &m2IsLoose2);
    //_ntuple1->SetBranchAddress("m2IsMedium", &m2IsMedium1);
    //_ntuple2->SetBranchAddress("m2IsMedium", &m2IsMedium2);
    //_ntuple1->SetBranchAddress("m2IsPFMuon", &m2IsPFMuon1);
    //_ntuple2->SetBranchAddress("m2IsPFMuon", &m2IsPFMuon2);
    //_ntuple1->SetBranchAddress("m2IsTracker", &m2IsTracker1);
    //_ntuple2->SetBranchAddress("m2IsTracker", &m2IsTracker2);
    //_ntuple1->SetBranchAddress("m2PASTightID", &m2PASTightID1);
    //_ntuple2->SetBranchAddress("m2PASTightID", &m2PASTightID2);
    //_ntuple1->SetBranchAddress("m2PASTightIDNoVtx", &m2PASTightIDNoVtx1);
    //_ntuple2->SetBranchAddress("m2PASTightIDNoVtx", &m2PASTightIDNoVtx2);
    //_ntuple1->SetBranchAddress("m2PFID", &m2PFID1);
    //_ntuple2->SetBranchAddress("m2PFID", &m2PFID2);
    //_ntuple1->SetBranchAddress("m2PFIDNoVtx", &m2PFIDNoVtx1);
    //_ntuple2->SetBranchAddress("m2PFIDNoVtx", &m2PFIDNoVtx2);
    //_ntuple1->SetBranchAddress("m2ZZIsoPass", &m2ZZIsoPass1);
    //_ntuple2->SetBranchAddress("m2ZZIsoPass", &m2ZZIsoPass2);
    //_ntuple1->SetBranchAddress("m2ZZLooseID", &m2ZZLooseID1);
    //_ntuple2->SetBranchAddress("m2ZZLooseID", &m2ZZLooseID2);
    //_ntuple1->SetBranchAddress("m2ZZLooseIDNoVtx", &m2ZZLooseIDNoVtx1);
    //_ntuple2->SetBranchAddress("m2ZZLooseIDNoVtx", &m2ZZLooseIDNoVtx2);
    //_ntuple1->SetBranchAddress("m2ZZTightID", &m2ZZTightID1);
    //_ntuple2->SetBranchAddress("m2ZZTightID", &m2ZZTightID2);
    //_ntuple1->SetBranchAddress("m2ZZTightIDNoVtx", &m2ZZTightIDNoVtx1);
    //_ntuple2->SetBranchAddress("m2ZZTightIDNoVtx", &m2ZZTightIDNoVtx2);
    //_ntuple1->SetBranchAddress("m2Charge", &m2Charge1);
    //_ntuple2->SetBranchAddress("m2Charge", &m2Charge2);
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m2GenCharge", &m2GenCharge1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m2GenCharge", &m2GenCharge2);
    }
    if (_isT1MC){
      //_ntuple1->SetBranchAddress("m2GenPdgId", &m2GenPdgId1);
    }
    if (_isT2MC){
      //_ntuple2->SetBranchAddress("m2GenPdgId", &m2GenPdgId2);
    }
    _ntuple1->SetBranchAddress("m2PdgId", &m2PdgId1);
    _ntuple2->SetBranchAddress("m2PdgId", &m2PdgId2);
    //_ntuple1->SetBranchAddress("m2BestTrackType", &m2BestTrackType1);
    //_ntuple2->SetBranchAddress("m2BestTrackType", &m2BestTrackType2);
    //_ntuple1->SetBranchAddress("m2CutBasedLoose", &m2CutBasedLoose1);
    //_ntuple2->SetBranchAddress("m2CutBasedLoose", &m2CutBasedLoose2);
    //_ntuple1->SetBranchAddress("m2CutBasedMedium", &m2CutBasedMedium1);
    //_ntuple2->SetBranchAddress("m2CutBasedMedium", &m2CutBasedMedium2);
    //_ntuple1->SetBranchAddress("m2CutBasedTight", &m2CutBasedTight1);
    //_ntuple2->SetBranchAddress("m2CutBasedTight", &m2CutBasedTight2);
    //_ntuple1->SetBranchAddress("m2IsTight", &m2IsTight1);
    //_ntuple2->SetBranchAddress("m2IsTight", &m2IsTight2);
    //_ntuple1->SetBranchAddress("m2MatchedStations", &m2MatchedStations1);
    //_ntuple2->SetBranchAddress("m2MatchedStations", &m2MatchedStations2);
    //_ntuple1->SetBranchAddress("m2NoOfMatches", &m2NoOfMatches1);
    //_ntuple2->SetBranchAddress("m2NoOfMatches", &m2NoOfMatches2);
    //_ntuple1->SetBranchAddress("m2PFIsoLoose", &m2PFIsoLoose1);
    //_ntuple2->SetBranchAddress("m2PFIsoLoose", &m2PFIsoLoose2);
    //_ntuple1->SetBranchAddress("m2PFIsoMedium", &m2PFIsoMedium1);
    //_ntuple2->SetBranchAddress("m2PFIsoMedium", &m2PFIsoMedium2);
    //_ntuple1->SetBranchAddress("m2PFIsoTight", &m2PFIsoTight1);
    //_ntuple2->SetBranchAddress("m2PFIsoTight", &m2PFIsoTight2);
    //_ntuple1->SetBranchAddress("m2PFIsoVTight", &m2PFIsoVTight1);
    //_ntuple2->SetBranchAddress("m2PFIsoVTight", &m2PFIsoVTight2);
  }
}

void CompLooperBase::SetTitles(TH1 *hist, std::string x, std::string y) const{
  if (hist->InheritsFrom("TH3")) return;
  hist->GetXaxis()->SetTitle(x.c_str());
  hist->GetYaxis()->SetTitle(y.c_str());
  //hist->GetYaxis()->SetTitleOffset(1.4);
}

double CompLooperBase::GetDeltaPhi(double phi1, double phi2) const{
  double deltaphi = phi1 - phi2;
  while (deltaphi < 0) deltaphi += (2 * M_PI);
  while (deltaphi > M_PI) deltaphi -= (2*M_PI);
  return deltaphi;
}

double CompLooperBase::GetPolCosTheta(const ROOT::Math::PtEtaPhiEVector &lp, const ROOT::Math::PtEtaPhiEVector &ln) const{
  ROOT::Math::PtEtaPhiEVector z = lp+ln;
  ROOT::Math::Boost boost(z.BoostToCM());
  ROOT::Math::PtEtaPhiEVector lep_boost = boost(ln);
  return lep_boost.Vect().Dot(z.Vect()) / std::sqrt( lep_boost.Vect().Mag2() * z.Vect().Mag2() );
}

void CompLooperBase::MakeIndexFile(std::string filename) const{
  try{
    std::ifstream fin("inputs/index.php");
    std::ofstream fout((filename + "/index.php").c_str());
    fout << fin.rdbuf();
  }
  catch (...){
    std::cerr << "Error making index file " << filename << std::endl;
  }
}

#endif//COMPLOOPERBASE_H
