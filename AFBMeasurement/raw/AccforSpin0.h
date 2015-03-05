//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug 13 11:56:52 2013 by ROOT version 5.32/00
// from TTree tree/tree
// found on file: dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/lathomas/DYToEE_M-200_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/TOTAL.root
//////////////////////////////////////////////////////////

#ifndef AccforSpin0_h
#define AccforSpin0_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TString.h"
#include <TH1.h>
#include <TTree.h>
// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxnWasRun = 1;
const Int_t kMaxnAccept = 1;
const Int_t kMaxnErrors = 1;

class AccforSpin0 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   virtual void Initialize(TString samplefile, TString mysamplename,  float sampleweight);
   virtual void DrawPlot(vector <TH1F*> myhisto );
   virtual void ComputeAcc(TH1F* hnumS0 , TH1F* hdenS0,TH1F* hnumS1 , TH1F* hdenS1,TH1F* hnumS2 , TH1F* hdenS2);
   virtual double CalcInvariantMass(const int& n, const int& l);
   virtual void CalcCosthetaStar(float &costh, float &costhmeas, const int m);
   virtual void CalcCostheta(float &costh, float &costhmeas, const int m);
   virtual void ComputeAcc(TH1F* hnum , TH1F* hden, TString gphname);
   virtual void BuildTree(TTree * T1) ;
   virtual bool PassHEEP(const int &n);
   virtual float  CalcCosthetaCSData(const int &m, const int & n);
   virtual float  CalcRapidity(const int &m, const int & n);
   // Declaration of leaf types
   UInt_t          runnumber;
   UInt_t          eventnumber;
   UInt_t          luminosityBlock;
   UInt_t          eventcounter;
   Int_t           processid;
   Float_t         pthat;
   Float_t         alphaqcd;
   Float_t         alphaqed;
   Float_t         qscale;
   Float_t         weight;
   Int_t           hltCount;
   Int_t           PhysDecl_bool;
   Int_t           nWasRun_;
   Int_t           nAccept_;
   Int_t           nErrors_;
   Int_t           HLT_Mu15_eta2p1;
   Int_t           HLT_Mu24_eta2p1;
   Int_t           HLT_Mu30_eta2p1;
   Int_t           HLT_Mu40_eta2p1;
   Int_t           HLT_Mu50_eta2p1;
   Int_t           HLT_Mu22_TkMu22;
   Int_t           HLT_Mu22_Photon22_CaloIdL;
   Int_t           HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Int_t           HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Int_t           HLT_Ele8_CaloIdL_CaloIsoVL;
   Int_t           HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Int_t           HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Int_t           HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50;
   Int_t           HLT_DoubleEle33_CaloIdL;
   Int_t           HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;
   Int_t           HLT_DoubleEle33_CaloIdT;
   Int_t           HLT_Photon20_CaloIdVL_IsoL;
   Int_t           HLT_Photon30_CaloIdVL;
   Int_t           HLT_Photon50_CaloIdVL;
   Int_t           HLT_Photon50_CaloIdVL_IsoL;
   Int_t           HLT_Photon75_CaloIdVL;
   Int_t           HLT_Photon90_CaloIdVL;
   Int_t           HLT_Photon135;
   Int_t           HLT_Photon150;
   Int_t           HLT_Photon250_NoHE;
   Int_t           HLT_Photon300_NoHE;
   Int_t           HLT_Photon26_Photon18;
   Int_t           HLT_Photon36_Photon22;
   Int_t           HLT_DoublePhoton70;
   Int_t           HLT_DoublePhoton80;
   Int_t           prescale_HLT_Mu15_eta2p1;
   Int_t           prescale_HLT_Mu30_eta2p1;
   Int_t           prescale_HLT_Mu40_eta2p1;
   Int_t           prescale_HLT_Mu22_TkMu22;
   Int_t           prescale_HLT_Mu22_Photon22_CaloIdL;
   Int_t           prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Int_t           prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Int_t           prescale_HLT_Ele8_CaloIdL_CaloIsoVL;
   Int_t           prescale_HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Int_t           prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Int_t           prescale_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50;
   Int_t           prescale_HLT_DoubleEle33_CaloIdL;
   Int_t           prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;
   Int_t           prescale_HLT_DoubleEle33_CaloIdT;
   Int_t           prescale_HLT_Photon20_CaloIdVL_IsoL;
   Int_t           prescale_HLT_Photon30_CaloIdVL;
   Int_t           prescale_HLT_Photon50_CaloIdVL;
   Int_t           prescale_HLT_Photon50_CaloIdVL_IsoL;
   Int_t           prescale_HLT_Photon75_CaloIdVL;
   Int_t           prescale_HLT_Photon90_CaloIdVL;
   Int_t           prescale_HLT_Photon135;
   Int_t           prescale_HLT_Photon150;
   Int_t           prescale_HLT_Photon250_NoHE;
   Int_t           prescale_HLT_Photon300_NoHE;
   Int_t           prescale_HLT_Photon26_Photon18;
   Int_t           prescale_HLT_Photon36_Photon22;
   Int_t           prescale_HLT_DoublePhoton70;
   Int_t           prescale_HLT_DoublePhoton80;
   Float_t         rho;
   Float_t         calomet;
   Float_t         calomet_phi;
   Float_t         met;
   Float_t         pfmet;
   Float_t         pfmet_phi;
   Float_t         pfmetcor;
   Float_t         pfmetcor_phi;
   Float_t         sigmaZ;
   Float_t         sigmaZ0Error;
   Float_t         sq;
   Float_t         bsposx;
   Float_t         bsposy;
   Float_t         bsposz;
   Int_t           pvsize;
   Float_t         pvx[57];   //[pvsize]
   Float_t         pvy[57];   //[pvsize]
   Float_t         pvz[57];   //[pvsize]
   Bool_t          pv_isValid[57];   //[pvsize]
   Int_t           pv_ndof[57];   //[pvsize]
   Int_t           pv_nTracks[57];   //[pvsize]
   Float_t         pv_normChi2[57];   //[pvsize]
   Int_t           pv_totTrackSize[57];   //[pvsize]
   Int_t           pfJetColl_size;
   Float_t         pfJet_pt[49];   //[pfJetColl_size]
   Float_t         pfJet_eta[49];   //[pfJetColl_size]
   Float_t         pfJet_phi[49];   //[pfJetColl_size]
   Int_t           JetColl_size;
   Float_t         Jet_pt[29];   //[JetColl_size]
   Float_t         Jet_eta[29];   //[JetColl_size]
   Float_t         Jet_phi[29];   //[JetColl_size]
   Float_t         Jet_vx[29];   //[JetColl_size]
   Float_t         Jet_vy[29];   //[JetColl_size]
   Float_t         Jet_vz[29];   //[JetColl_size]
   Float_t         tCHighEffBTags[29];   //[JetColl_size]
   Float_t         tCHighPurBTags[29];   //[JetColl_size]
   Float_t         jetProbBTags[29];   //[JetColl_size]
   Float_t         jetBProbBTags[29];   //[JetColl_size]
   Float_t         sSecVertHighEffBTags[29];   //[JetColl_size]
   Float_t         sSecVertHighPurBTags[29];   //[JetColl_size]
   Float_t         cSecVertBTags[29];   //[JetColl_size]
   Float_t         cSecVertMVABTags[29];   //[JetColl_size]
   Float_t         ghostTrkBTags[29];   //[JetColl_size]
   Float_t         softEleIP3dBTags[29];   //[JetColl_size]
   Float_t         softElePtBTags[29];   //[JetColl_size]
   Float_t         softMuBTags[29];   //[JetColl_size]
   Float_t         softMuIP3dBTags[29];   //[JetColl_size]
   Float_t         softMuPtBTags[29];   //[JetColl_size]
   Int_t           muon_size;
   Float_t         muon_pt[3];   //[muon_size]
   Float_t         muon_ptError[3];   //[muon_size]
   Float_t         muon_gTrk_pt[3];   //[muon_size]
   Float_t         muon_gTrk_ptError[3];   //[muon_size]
   Float_t         muon_eta[3];   //[muon_size]
   Float_t         muon_etaError[3];   //[muon_size]
   Float_t         muon_theta[3];   //[muon_size]
   Float_t         muon_thetaError[3];   //[muon_size]
   Float_t         muon_phi[3];   //[muon_size]
   Float_t         muon_phiError[3];   //[muon_size]
   Float_t         muon_outerPt[3];   //[muon_size]
   Float_t         muon_outerEta[3];   //[muon_size]
   Float_t         muon_outerPhi[3];   //[muon_size]
   Float_t         muon_outerTheta[3];   //[muon_size]
   Float_t         muon_px[3];   //[muon_size]
   Float_t         muon_py[3];   //[muon_size]
   Float_t         muon_pz[3];   //[muon_size]
   Int_t           muon_charge[3];   //[muon_size]
   Int_t           muon_nhitspixel[3];   //[muon_size]
   Int_t           muon_nhitstrack[3];   //[muon_size]
   Int_t           muon_nhitsmuons[3];   //[muon_size]
   Int_t           muon_nhitstotal[3];   //[muon_size]
   Int_t           muon_nlayerswithhits[3];   //[muon_size]
   Int_t           muon_nlosthits[3];   //[muon_size]
   Int_t           muon_nSegmentMatch[3];   //[muon_size]
   Bool_t          muon_isTrackerMuon[3];   //[muon_size]
   Bool_t          muon_isPFMuon[3];   //[muon_size]
   Bool_t          muon_isPFIsolationValid[3];   //[muon_size]
   Float_t         muon_chi2[3];   //[muon_size]
   Int_t           muon_ndof[3];   //[muon_size]
   Float_t         muon_normChi2[3];   //[muon_size]
   Float_t         muon_d0[3];   //[muon_size]
   Float_t         muon_d0Error[3];   //[muon_size]
   Float_t         muon_dz_cmsCenter[3];   //[muon_size]
   Float_t         muon_dz_beamSpot[3];   //[muon_size]
   Float_t         muon_dz_firstPVtx[3];   //[muon_size]
   Float_t         muon_dz_firstPVtxwithBS[3];   //[muon_size]
   Float_t         muon_dzError[3];   //[muon_size]
   Float_t         muon_dxy_cmsCenter[3];   //[muon_size]
   Float_t         muon_dxy_beamSpot[3];   //[muon_size]
   Float_t         muon_dxy_firstPVtx[3];   //[muon_size]
   Float_t         muon_dxy_firstPVtxwithBS[3];   //[muon_size]
   Float_t         muon_dxyError[3];   //[muon_size]
   Float_t         muon_innerPosx[3];   //[muon_size]
   Float_t         muon_innerPosy[3];   //[muon_size]
   Float_t         muon_innerPosz[3];   //[muon_size]
   Float_t         muon_trackIso03[3];   //[muon_size]
   Float_t         muon_trackIso05[3];   //[muon_size]
   Float_t         muon_trackIso03_ptInVeto[3];   //[muon_size]
   Float_t         muon_trackIso05_ptInVeto[3];   //[muon_size]
   Float_t         muon_emIso03[3];   //[muon_size]
   Float_t         muon_emIso05[3];   //[muon_size]
   Float_t         muon_emIso03_ptInVeto[3];   //[muon_size]
   Float_t         muon_emIso05_ptInVeto[3];   //[muon_size]
   Float_t         muon_hadIso03[3];   //[muon_size]
   Float_t         muon_hadIso05[3];   //[muon_size]
   Float_t         muon_hadIso03_ptInVeto[3];   //[muon_size]
   Float_t         muon_hadIso05_ptInVeto[3];   //[muon_size]
   Int_t           scsize;
   Float_t         scenergy[7];   //[scsize]
   Float_t         sceta[7];   //[scsize]
   Float_t         scetacorr[7];   //[scsize]
   Float_t         scet[7];   //[scsize]
   Float_t         scphi[7];   //[scsize]
   Float_t         scpx[7];   //[scsize]
   Float_t         scpy[7];   //[scsize]
   Float_t         scpz[7];   //[scsize]
   Float_t         scx[7];   //[scsize]
   Float_t         scy[7];   //[scsize]
   Float_t         scz[7];   //[scsize]
   Int_t           pfele_size;
   Float_t         pfele_pt[4];   //[pfele_size]
   Float_t         pfele_eta[4];   //[pfele_size]
   Float_t         pfele_phi[4];   //[pfele_size]
   Int_t           pfele_charge[4];   //[pfele_size]
   Int_t           gsf_size;
   Bool_t          gsf_isEB[8];   //[gsf_size]
   Bool_t          gsf_isEE[8];   //[gsf_size]
   Float_t         gsf_px[8];   //[gsf_size]
   Float_t         gsf_py[8];   //[gsf_size]
   Float_t         gsf_pz[8];   //[gsf_size]
   Float_t         gsf_pt[8];   //[gsf_size]
   Float_t         gsf_eta[8];   //[gsf_size]
   Float_t         gsf_phi[8];   //[gsf_size]
   Float_t         gsf_theta[8];   //[gsf_size]
   Int_t           gsf_charge[8];   //[gsf_size]
   Float_t         gsf_deltaEtaATcalo[8];   //[gsf_size]
   Float_t         gsf_deltaPhiATcalo[8];   //[gsf_size]
   Float_t         gsf_sigmaetaeta[8];   //[gsf_size]
   Float_t         gsf_sigmaIetaIeta[8];   //[gsf_size]
   Float_t         gsf_ecalEnergy[8];   //[gsf_size]
   Float_t         gsf_eOVERp[8];   //[gsf_size]
   Float_t         gsf_dxy[8];   //[gsf_size]
   Float_t         gsf_dxy_beamSpot[8];   //[gsf_size]
   Float_t         gsf_dxy_firstPVtx[8];   //[gsf_size]
   Float_t         gsf_dxy_firstPVtxwithBS[8];   //[gsf_size]
   Float_t         gsf_dxyError[8];   //[gsf_size]
   Float_t         gsf_dz[8];   //[gsf_size]
   Float_t         gsf_dz_beamSpot[8];   //[gsf_size]
   Float_t         gsf_dz_firstPVtx[8];   //[gsf_size]
   Float_t         gsf_dz_firstPVtxwithBS[8];   //[gsf_size]
   Float_t         gsf_dzError[8];   //[gsf_size]
   Float_t         gsf_vz[8];   //[gsf_size]
   Int_t           gsf_nHits[8];   //[gsf_size]
   Int_t           gsf_nLostInnerHits[8];   //[gsf_size]
   Int_t           gsf_nLostOuterHits[8];   //[gsf_size]
   Int_t           gsf_convFlags[8];   //[gsf_size]
   Float_t         gsf_convDist[8];   //[gsf_size]
   Float_t         gsf_convDcot[8];   //[gsf_size]
   Float_t         gsf_convRadius[8];   //[gsf_size]
   Float_t         gsf_fBrem[8];   //[gsf_size]
   Float_t         gsf_e1x5[8];   //[gsf_size]
   Float_t         gsf_e2x5[8];   //[gsf_size]
   Float_t         gsf_e5x5[8];   //[gsf_size]
   Float_t         gsf_e1x3[8];   //[gsf_size]
   Float_t         gsf_p[8];   //[gsf_size]
   Float_t         gsf_e[8];   //[gsf_size]
   Float_t         gsf_deltaeta[8];   //[gsf_size]
   Float_t         gsf_deltaphi[8];   //[gsf_size]
   Float_t         gsf_hovere[8];   //[gsf_size]
   Float_t         gsf_hdepth1overe[8];   //[gsf_size]
   Float_t         gsf_hdepth2overe[8];   //[gsf_size]
   Float_t         gsf_hovere2012[8];   //[gsf_size]
   Float_t         gsf_hdepth1overe2012[8];   //[gsf_size]
   Float_t         gsf_hdepth2overe2012[8];   //[gsf_size]
   Float_t         gsf_trackiso[8];   //[gsf_size]
   Float_t         gsf_ecaliso[8];   //[gsf_size]
   Float_t         gsf_hcaliso1[8];   //[gsf_size]
   Float_t         gsf_hcaliso2[8];   //[gsf_size]
   Float_t         gsf_hcaliso12012[8];   //[gsf_size]
   Float_t         gsf_hcaliso22012[8];   //[gsf_size]
   Float_t         gsf_PFisocharged[8];   //[gsf_size]
   Float_t         gsf_PFisophoton[8];   //[gsf_size]
   Float_t         gsf_PFisoneutral[8];   //[gsf_size]
   Float_t         gsf_class[8];   //[gsf_size]
   Bool_t          gsf_isecaldriven[8];   //[gsf_size]
   Bool_t          gsf_istrackerdriven[8];   //[gsf_size]
   Float_t         gsfsc_e[8];   //[gsf_size]
   Float_t         gsfsc_pt[8];   //[gsf_size]
   Float_t         gsfsc_eta[8];   //[gsf_size]
   Float_t         gsfsc_phi[8];   //[gsf_size]
   Float_t         gsfsc_px[8];   //[gsf_size]
   Float_t         gsfsc_py[8];   //[gsf_size]
   Float_t         gsfsc_pz[8];   //[gsf_size]
   Float_t         gsf_e2x5overe5x5[8];   //[gsf_size]
   Float_t         gsf_e1x5overe5x5[8];   //[gsf_size]
   Float_t         gsf_gsfet[8];   //[gsf_size]
   UInt_t          gsf_hitsinfo[2][25];
   Int_t           scindexforgsf[8];   //[gsf_size]
   Int_t           gsftracksize;
   Float_t         gsftracketa[12];   //[gsftracksize]
   Float_t         gsftrackphi[12];   //[gsftracksize]
   Float_t         gsftrackp[12];   //[gsftracksize]
   Float_t         gsftrackpt[12];   //[gsftracksize]
   Float_t         gsftrackpx[12];   //[gsftracksize]
   Float_t         gsftrackpy[12];   //[gsftracksize]
   Float_t         gsftrackpz[12];   //[gsftracksize]
   Bool_t          gsfpass_ET[8];   //[gsf_size]
   Bool_t          gsfpass_PT[8];   //[gsf_size]
   Bool_t          gsfpass_DETETA[8];   //[gsf_size]
   Bool_t          gsfpass_CRACK[8];   //[gsf_size]
   Bool_t          gsfpass_DETAIN[8];   //[gsf_size]
   Bool_t          gsfpass_DPHIIN[8];   //[gsf_size]
   Bool_t          gsfpass_HADEM[8];   //[gsf_size]
   Bool_t          gsfpass_SIGMAIETAIETA[8];   //[gsf_size]
   Bool_t          gsfpass_E2X5OVER5X5[8];   //[gsf_size]
   Bool_t          gsfpass_ISOLEMHADDEPTH1[8];   //[gsf_size]
   Bool_t          gsfpass_ISOLHADDEPTH2[8];   //[gsf_size]
   Bool_t          gsfpass_ISOLPTTRKS[8];   //[gsf_size]
   Bool_t          gsfpass_ECALDRIVEN[8];   //[gsf_size]
   Bool_t          gsfpass_INVALID[8];   //[gsf_size]
   Bool_t          gsfpass_NOMISSINGHITS[8];   //[gsf_size]
   Bool_t          gsfpass_NOCONVERSION[8];   //[gsf_size]
   Bool_t          gsfpass_DXYFIRSTPV[8];   //[gsf_size]
   Bool_t          gsfpass_HEEP[8];   //[gsf_size]
   Bool_t          gsfpass_ID[8];   //[gsf_size]
   Bool_t          gsfpass_ISO[8];   //[gsf_size]
   Bool_t          gsfmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter[8];   //[gsf_size]
   Bool_t          gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter[8];   //[gsf_size]
   Bool_t          gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter[8];   //[gsf_size]
   Bool_t          gsfmatch_hltL1sL1SingleEG22[8];   //[gsf_size]
   Bool_t          gsfmatch_hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter[8];   //[gsf_size]
   Bool_t          gsfmatch_hltEle33CaloIdLPixelMatchFilter[8];   //[gsf_size]
   Bool_t          scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter[7];   //[scsize]
   Float_t         heepHeepMass;
   Int_t           scpixcharge[8];   //[gsf_size]
   Int_t           ctfcharge[8];   //[gsf_size]
   Int_t           gsfcharge[8];   //[gsf_size]
   Bool_t          gsfctfscpixconsistent[8];   //[gsf_size]
   Bool_t          gsfscpixconsistent[8];   //[gsf_size]
   Bool_t          gsfctfconsistent[8];   //[gsf_size]
   Int_t           genparticles_size;
   Int_t           genquarks_size;
   Int_t           gengluons_size;
   Int_t           conv_size;
   Float_t         conv_vtxProb[23];   //[conv_size]
   Float_t         conv_lxy[23];   //[conv_size]
   Int_t           conv_nHitsMax[23];   //[conv_size]
   Int_t           conv_eleind[23];   //[conv_size]
   Int_t           gsf0_crystal_size;
   Int_t           gsf0_crystal_ietaorix[252];   //[gsf0_crystal_size]
   Int_t           gsf0_crystal_iphioriy[252];   //[gsf0_crystal_size]
   Float_t         gsf0_crystal_energy[252];   //[gsf0_crystal_size]
   Float_t         gsf0_crystal_eta[252];   //[gsf0_crystal_size]
   Int_t           gsf1_crystal_size;
   Int_t           gsf1_crystal_ietaorix[288];   //[gsf1_crystal_size]
   Int_t           gsf1_crystal_iphioriy[288];   //[gsf1_crystal_size]
   Float_t         gsf1_crystal_energy[288];   //[gsf1_crystal_size]
   Float_t         gsf1_crystal_eta[288];   //[gsf1_crystal_size]
   Float_t         genquark_e[4];   //[genquarks_size]
   Float_t         genquark_eta[4];   //[genquarks_size]
   Float_t         genquark_phi[4];   //[genquarks_size]
   Float_t         genquark_pt[4];   //[genquarks_size]
   Float_t         genquark_px[4];   //[genquarks_size]
   Float_t         genquark_py[4];   //[genquarks_size]
   Float_t         genquark_pz[4];   //[genquarks_size]
   Int_t           genquark_status[4];   //[genquarks_size]
   Int_t           genquark_charge[4];   //[genquarks_size]
   Int_t           genquark_pdgid[4];   //[genquarks_size]
   Float_t         gengluon_e[2];   //[gengluons_size]
   Float_t         gengluon_eta[2];   //[gengluons_size]
   Float_t         gengluon_phi[2];   //[gengluons_size]
   Float_t         gengluon_pt[2];   //[gengluons_size]
   Float_t         gengluon_px[2];   //[gengluons_size]
   Float_t         gengluon_py[2];   //[gengluons_size]
   Float_t         gengluon_pz[2];   //[gengluons_size]
   Int_t           gengluon_status[2];   //[gengluons_size]
   Int_t           gengluon_charge[2];   //[gengluons_size]
   Int_t           gengluon_pdgid[2];   //[gengluons_size]
   Float_t         genele_e[2];   //[genparticles_size]
   Float_t         genele_eta[2];   //[genparticles_size]
   Float_t         genele_phi[2];   //[genparticles_size]
   Float_t         genele_pt[2];   //[genparticles_size]
   Float_t         genele_px[2];   //[genparticles_size]
   Float_t         genele_py[2];   //[genparticles_size]
   Float_t         genele_pz[2];   //[genparticles_size]
   Int_t           genele_charge[2];   //[genparticles_size]
   Float_t         unstableGenEle_e[2];   //[genparticles_size]
   Float_t         unstableGenEle_eta[2];   //[genparticles_size]
   Float_t         unstableGenEle_phi[2];   //[genparticles_size]
   Float_t         unstableGenEle_pt[2];   //[genparticles_size]
   Float_t         unstableGenEle_px[2];   //[genparticles_size]
   Float_t         unstableGenEle_py[2];   //[genparticles_size]
   Float_t         unstableGenEle_pz[2];   //[genparticles_size]
   Int_t           unstableGenEle_charge[2];   //[genparticles_size]
   Float_t         genelemom_e[2];   //[genparticles_size]
   Float_t         genelemom_eta[2];   //[genparticles_size]
   Float_t         genelemom_phi[2];   //[genparticles_size]
   Float_t         genelemom_pt[2];   //[genparticles_size]
   Float_t         genelemom_px[2];   //[genparticles_size]
   Float_t         genelemom_py[2];   //[genparticles_size]
   Float_t         genelemom_pz[2];   //[genparticles_size]
   Int_t           genelemom_charge[2];   //[genparticles_size]
   Int_t           genelemom_pdgid[2];   //[genparticles_size]
   Float_t         genelemom_mass[2];   //[genparticles_size]
   Float_t         x1quark[2];   //[genparticles_size]
   Float_t         x2quark[2];   //[genparticles_size]
   Float_t         genPair_mass;
   Int_t           trueNVtx;
   Int_t           nVtxBefore;
   Int_t           nVtxNow;
   Int_t           nVtxAfter;

   // List of branches
   TBranch        *b_runnumber;   //!
   TBranch        *b_eventnumber;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_eventcounter;   //!
   TBranch        *b_processid;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_alphaqcd;   //!
   TBranch        *b_alphaqed;   //!
   TBranch        *b_qscale;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_hltCount;   //!
   TBranch        *b_PhysDecl_bool;   //!
   TBranch        *b_nWasRun_;   //!
   TBranch        *b_nAccept_;   //!
   TBranch        *b_nErrors_;   //!
   TBranch        *b_HLT_Mu15_eta2p1;   //!
   TBranch        *b_HLT_Mu24_eta2p1;   //!
   TBranch        *b_HLT_Mu30_eta2p1;   //!
   TBranch        *b_HLT_Mu40_eta2p1;   //!
   TBranch        *b_HLT_Mu50_eta2p1;   //!
   TBranch        *b_HLT_Mu22_TkMu22;   //!
   TBranch        *b_HLT_Mu22_Photon22_CaloIdL;   //!
   TBranch        *b_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_HLT_Ele8_CaloIdL_CaloIsoVL;   //!
   TBranch        *b_HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50;   //!
   TBranch        *b_HLT_DoubleEle33_CaloIdL;   //!
   TBranch        *b_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;   //!
   TBranch        *b_HLT_DoubleEle33_CaloIdT;   //!
   TBranch        *b_HLT_Photon20_CaloIdVL_IsoL;   //!
   TBranch        *b_HLT_Photon30_CaloIdVL;   //!
   TBranch        *b_HLT_Photon50_CaloIdVL;   //!
   TBranch        *b_HLT_Photon50_CaloIdVL_IsoL;   //!
   TBranch        *b_HLT_Photon75_CaloIdVL;   //!
   TBranch        *b_HLT_Photon90_CaloIdVL;   //!
   TBranch        *b_HLT_Photon135;   //!
   TBranch        *b_HLT_Photon150;   //!
   TBranch        *b_HLT_Photon250_NoHE;   //!
   TBranch        *b_HLT_Photon300_NoHE;   //!
   TBranch        *b_HLT_Photon26_Photon18;   //!
   TBranch        *b_HLT_Photon36_Photon22;   //!
   TBranch        *b_HLT_DoublePhoton70;   //!
   TBranch        *b_HLT_DoublePhoton80;   //!
   TBranch        *b_prescale_HLT_Mu15_eta2p1;   //!
   TBranch        *b_prescale_HLT_Mu30_eta2p1;   //!
   TBranch        *b_prescale_HLT_Mu40_eta2p1;   //!
   TBranch        *b_prescale_HLT_Mu22_TkMu22;   //!
   TBranch        *b_prescale_HLT_Mu22_Photon22_CaloIdL;   //!
   TBranch        *b_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_prescale_HLT_Ele8_CaloIdL_CaloIsoVL;   //!
   TBranch        *b_prescale_HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_prescale_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50;   //!
   TBranch        *b_prescale_HLT_DoubleEle33_CaloIdL;   //!
   TBranch        *b_prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;   //!
   TBranch        *b_prescale_HLT_DoubleEle33_CaloIdT;   //!
   TBranch        *b_prescale_HLT_Photon20_CaloIdVL_IsoL;   //!
   TBranch        *b_prescale_HLT_Photon30_CaloIdVL;   //!
   TBranch        *b_prescale_HLT_Photon50_CaloIdVL;   //!
   TBranch        *b_prescale_HLT_Photon50_CaloIdVL_IsoL;   //!
   TBranch        *b_prescale_HLT_Photon75_CaloIdVL;   //!
   TBranch        *b_prescale_HLT_Photon90_CaloIdVL;   //!
   TBranch        *b_prescale_HLT_Photon135;   //!
   TBranch        *b_prescale_HLT_Photon150;   //!
   TBranch        *b_prescale_HLT_Photon250_NoHE;   //!
   TBranch        *b_prescale_HLT_Photon300_NoHE;   //!
   TBranch        *b_prescale_HLT_Photon26_Photon18;   //!
   TBranch        *b_prescale_HLT_Photon36_Photon22;   //!
   TBranch        *b_prescale_HLT_DoublePhoton70;   //!
   TBranch        *b_prescale_HLT_DoublePhoton80;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_calomet;   //!
   TBranch        *b_calomet_phi;   //!
   TBranch        *b_met;   //!
   TBranch        *b_pfmet;   //!
   TBranch        *b_pfmet_phi;   //!
   TBranch        *b_pfmetcor;   //!
   TBranch        *b_pfmetcor_phi;   //!
   TBranch        *b_sigmaZ;   //!
   TBranch        *b_sigmaZ0Error;   //!
   TBranch        *b_sq;   //!
   TBranch        *b_bsposx;   //!
   TBranch        *b_bsposy;   //!
   TBranch        *b_bsposz;   //!
   TBranch        *b_pvsize;   //!
   TBranch        *b_pvx;   //!
   TBranch        *b_pvy;   //!
   TBranch        *b_pvz;   //!
   TBranch        *b_pv_isValid;   //!
   TBranch        *b_pv_ndof;   //!
   TBranch        *b_pv_nTracks;   //!
   TBranch        *b_pv_normChi2;   //!
   TBranch        *b_pv_totTrackSize;   //!
   TBranch        *b_pfJetColl_size;   //!
   TBranch        *b_pfJet_pt;   //!
   TBranch        *b_pfJet_eta;   //!
   TBranch        *b_pfJet_phi;   //!
   TBranch        *b_JetColl_size;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_vx;   //!
   TBranch        *b_Jet_vy;   //!
   TBranch        *b_Jet_vz;   //!
   TBranch        *b_tCHighEffBTags;   //!
   TBranch        *b_tCHighPurBTags;   //!
   TBranch        *b_jetProbBTags;   //!
   TBranch        *b_jetBProbBTags;   //!
   TBranch        *b_sSecVertHighEffBTags;   //!
   TBranch        *b_sSecVertHighPurBTags;   //!
   TBranch        *b_cSecVertBTags;   //!
   TBranch        *b_cSecVertMVABTags;   //!
   TBranch        *b_ghostTrkBTags;   //!
   TBranch        *b_softEleIP3dBTags;   //!
   TBranch        *b_softElePtBTags;   //!
   TBranch        *b_softMuBTags;   //!
   TBranch        *b_softMuIP3dBTags;   //!
   TBranch        *b_softMuPtBTags;   //!
   TBranch        *b_muon_size;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_ptError;   //!
   TBranch        *b_muon_gTrk_pt;   //!
   TBranch        *b_muon_gTrk_ptError;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_etaError;   //!
   TBranch        *b_muon_theta;   //!
   TBranch        *b_muon_thetaError;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_phiError;   //!
   TBranch        *b_muon_outerPt;   //!
   TBranch        *b_muon_outerEta;   //!
   TBranch        *b_muon_outerPhi;   //!
   TBranch        *b_muon_outerTheta;   //!
   TBranch        *b_muon_px;   //!
   TBranch        *b_muon_py;   //!
   TBranch        *b_muon_pz;   //!
   TBranch        *b_muon_charge;   //!
   TBranch        *b_muon_nhitspixel;   //!
   TBranch        *b_muon_nhitstrack;   //!
   TBranch        *b_muon_nhitsmuons;   //!
   TBranch        *b_muon_nhitstotal;   //!
   TBranch        *b_muon_nlayerswithhits;   //!
   TBranch        *b_muon_nlosthits;   //!
   TBranch        *b_muon_nSegmentMatch;   //!
   TBranch        *b_muon_isTrackerMuon;   //!
   TBranch        *b_muon_isPFMuon;   //!
   TBranch        *b_muon_isPFIsolationValid;   //!
   TBranch        *b_muon_chi2;   //!
   TBranch        *b_muon_ndof;   //!
   TBranch        *b_muon_normChi2;   //!
   TBranch        *b_muon_d0;   //!
   TBranch        *b_muon_d0Error;   //!
   TBranch        *b_muon_dz_cmsCenter;   //!
   TBranch        *b_muon_dz_beamSpot;   //!
   TBranch        *b_muon_dz_firstPVtx;   //!
   TBranch        *b_muon_dz_firstPVtxwithBS;   //!
   TBranch        *b_muon_dzError;   //!
   TBranch        *b_muon_dxy_cmsCenter;   //!
   TBranch        *b_muon_dxy_beamSpot;   //!
   TBranch        *b_muon_dxy_firstPVtx;   //!
   TBranch        *b_muon_dxy_firstPVtxwithBS;   //!
   TBranch        *b_muon_dxyError;   //!
   TBranch        *b_muon_innerPosx;   //!
   TBranch        *b_muon_innerPosy;   //!
   TBranch        *b_muon_innerPosz;   //!
   TBranch        *b_muon_trackIso03;   //!
   TBranch        *b_muon_trackIso05;   //!
   TBranch        *b_muon_trackIso03_ptInVeto;   //!
   TBranch        *b_muon_trackIso05_ptInVeto;   //!
   TBranch        *b_muon_emIso03;   //!
   TBranch        *b_muon_emIso05;   //!
   TBranch        *b_muon_emIso03_ptInVeto;   //!
   TBranch        *b_muon_emIso05_ptInVeto;   //!
   TBranch        *b_muon_hadIso03;   //!
   TBranch        *b_muon_hadIso05;   //!
   TBranch        *b_muon_hadIso03_ptInVeto;   //!
   TBranch        *b_muon_hadIso05_ptInVeto;   //!
   TBranch        *b_scsize;   //!
   TBranch        *b_scenergy;   //!
   TBranch        *b_sceta;   //!
   TBranch        *b_scetacorr;   //!
   TBranch        *b_scet;   //!
   TBranch        *b_scphi;   //!
   TBranch        *b_scpx;   //!
   TBranch        *b_scpy;   //!
   TBranch        *b_scpz;   //!
   TBranch        *b_scx;   //!
   TBranch        *b_scy;   //!
   TBranch        *b_scz;   //!
   TBranch        *b_pfele_size;   //!
   TBranch        *b_pfele_pt;   //!
   TBranch        *b_pfele_eta;   //!
   TBranch        *b_pfele_phi;   //!
   TBranch        *b_pfele_charge;   //!
   TBranch        *b_gsf_size;   //!
   TBranch        *b_gsf_isEB;   //!
   TBranch        *b_gsf_isEE;   //!
   TBranch        *b_gsf_px;   //!
   TBranch        *b_gsf_py;   //!
   TBranch        *b_gsf_pz;   //!
   TBranch        *b_gsf_pt;   //!
   TBranch        *b_gsf_eta;   //!
   TBranch        *b_gsf_phi;   //!
   TBranch        *b_gsf_theta;   //!
   TBranch        *b_gsf_charge;   //!
   TBranch        *b_gsf_deltaEtaATcalo;   //!
   TBranch        *b_gsf_deltaPhiATcalo;   //!
   TBranch        *b_gsf_sigmaetaeta;   //!
   TBranch        *b_gsf_sigmaIetaIeta;   //!
   TBranch        *b_gsf_ecalEnergy;   //!
   TBranch        *b_gsf_eOVERp;   //!
   TBranch        *b_gsf_dxy;   //!
   TBranch        *b_gsf_dxy_beamSpot;   //!
   TBranch        *b_gsf_dxy_firstPVtx;   //!
   TBranch        *b_gsf_dxy_firstPVtxwithBS;   //!
   TBranch        *b_gsf_dxyError;   //!
   TBranch        *b_gsf_dz;   //!
   TBranch        *b_gsf_dz_beamSpot;   //!
   TBranch        *b_gsf_dz_firstPVtx;   //!
   TBranch        *b_gsf_dz_firstPVtxwithBS;   //!
   TBranch        *b_gsf_dzError;   //!
   TBranch        *b_gsf_vz;   //!
   TBranch        *b_gsf_nHits;   //!
   TBranch        *b_gsf_nLostInnerHits;   //!
   TBranch        *b_gsf_nLostOuterHits;   //!
   TBranch        *b_gsf_convFlags;   //!
   TBranch        *b_gsf_convDist;   //!
   TBranch        *b_gsf_convDcot;   //!
   TBranch        *b_gsf_convRadius;   //!
   TBranch        *b_gsf_fBrem;   //!
   TBranch        *b_gsf_e1x5;   //!
   TBranch        *b_gsf_e2x5;   //!
   TBranch        *b_gsf_e5x5;   //!
   TBranch        *b_gsf_e1x3;   //!
   TBranch        *b_gsf_p;   //!
   TBranch        *b_gsf_e;   //!
   TBranch        *b_gsf_deltaeta;   //!
   TBranch        *b_gsf_deltaphi;   //!
   TBranch        *b_gsf_hovere;   //!
   TBranch        *b_gsf_hdepth1overe;   //!
   TBranch        *b_gsf_hdepth2overe;   //!
   TBranch        *b_gsf_hovere2012;   //!
   TBranch        *b_gsf_hdepth1overe2012;   //!
   TBranch        *b_gsf_hdepth2overe2012;   //!
   TBranch        *b_gsf_trackiso;   //!
   TBranch        *b_gsf_ecaliso;   //!
   TBranch        *b_gsf_hcaliso1;   //!
   TBranch        *b_gsf_hcaliso2;   //!
   TBranch        *b_gsf_hcaliso12012;   //!
   TBranch        *b_gsf_hcaliso22012;   //!
   TBranch        *b_gsf_PFisocharged;   //!
   TBranch        *b_gsf_PFisophoton;   //!
   TBranch        *b_gsf_PFisoneutral;   //!
   TBranch        *b_gsf_class;   //!
   TBranch        *b_gsf_isecaldriven;   //!
   TBranch        *b_gsf_istrackerdriven;   //!
   TBranch        *b_gsfsc_e;   //!
   TBranch        *b_gsfsc_pt;   //!
   TBranch        *b_gsfsc_eta;   //!
   TBranch        *b_gsfsc_phi;   //!
   TBranch        *b_gsfsc_px;   //!
   TBranch        *b_gsfsc_py;   //!
   TBranch        *b_gsfsc_pz;   //!
   TBranch        *b_gsf_e2x5overe5x5;   //!
   TBranch        *b_gsf_e1x5overe5x5;   //!
   TBranch        *b_gsf_gsfet;   //!
   TBranch        *b_gsf_hitsinfo;   //!
   TBranch        *b_scindexforgsf;   //!
   TBranch        *b_gsftracksize;   //!
   TBranch        *b_gsftracketa;   //!
   TBranch        *b_gsftrackphi;   //!
   TBranch        *b_gsftrackp;   //!
   TBranch        *b_gsftrackpt;   //!
   TBranch        *b_gsftrackpx;   //!
   TBranch        *b_gsftrackpy;   //!
   TBranch        *b_gsftrackpz;   //!
   TBranch        *b_gsfpass_ET;   //!
   TBranch        *b_gsfpass_PT;   //!
   TBranch        *b_gsfpass_DETETA;   //!
   TBranch        *b_gsfpass_CRACK;   //!
   TBranch        *b_gsfpass_DETAIN;   //!
   TBranch        *b_gsfpass_DPHIIN;   //!
   TBranch        *b_gsfpass_HADEM;   //!
   TBranch        *b_gsfpass_SIGMAIETAIETA;   //!
   TBranch        *b_gsfpass_E2X5OVER5X5;   //!
   TBranch        *b_gsfpass_ISOLEMHADDEPTH1;   //!
   TBranch        *b_gsfpass_ISOLHADDEPTH2;   //!
   TBranch        *b_gsfpass_ISOLPTTRKS;   //!
   TBranch        *b_gsfpass_ECALDRIVEN;   //!
   TBranch        *b_gsfpass_INVALID;   //!
   TBranch        *b_gsfpass_NOMISSINGHITS;   //!
   TBranch        *b_gsfpass_NOCONVERSION;   //!
   TBranch        *b_gsfpass_DXYFIRSTPV;   //!
   TBranch        *b_gsfpass_HEEP;   //!
   TBranch        *b_gsfpass_ID;   //!
   TBranch        *b_gsfpass_ISO;   //!
   TBranch        *b_gsfmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter;   //!
   TBranch        *b_gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter;   //!
   TBranch        *b_gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter;   //!
   TBranch        *b_gsfmatch_hltL1sL1SingleEG22;   //!
   TBranch        *b_gsfmatch_hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter;   //!
   TBranch        *b_gsfmatch_hltEle33CaloIdLPixelMatchFilter;   //!
   TBranch        *b_scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter;   //!
   TBranch        *b_heepHeepMass;   //!
   TBranch        *b_scpixcharge;   //!
   TBranch        *b_ctfcharge;   //!
   TBranch        *b_gsfcharge;   //!
   TBranch        *b_gsfctfscpixconsistent;   //!
   TBranch        *b_gsfscpixconsistent;   //!
   TBranch        *b_gsfctfconsistent;   //!
   TBranch        *b_genparticles_size;   //!
   TBranch        *b_genquarks_size;   //!
   TBranch        *b_gengluons_size;   //!
   TBranch        *b_conv_size;   //!
   TBranch        *b_conv_vtxProb;   //!
   TBranch        *b_conv_lxy;   //!
   TBranch        *b_conv_nHitsMax;   //!
   TBranch        *b_conv_eleind;   //!
   TBranch        *b_gsf0_crystal_size;   //!
   TBranch        *b_gsf0_crystal_ietaorix;   //!
   TBranch        *b_gsf0_crystal_iphioriy;   //!
   TBranch        *b_gsf0_crystal_energy;   //!
   TBranch        *b_gsf0_crystal_eta;   //!
   TBranch        *b_gsf1_crystal_size;   //!
   TBranch        *b_gsf1_crystal_ietaorix;   //!
   TBranch        *b_gsf1_crystal_iphioriy;   //!
   TBranch        *b_gsf1_crystal_energy;   //!
   TBranch        *b_gsf1_crystal_eta;   //!
   TBranch        *b_genquark_e;   //!
   TBranch        *b_genquark_eta;   //!
   TBranch        *b_genquark_phi;   //!
   TBranch        *b_genquark_pt;   //!
   TBranch        *b_genquark_px;   //!
   TBranch        *b_genquark_py;   //!
   TBranch        *b_genquark_pz;   //!
   TBranch        *b_genquark_status;   //!
   TBranch        *b_genquark_charge;   //!
   TBranch        *b_genquark_pdgid;   //!
   TBranch        *b_gengluon_e;   //!
   TBranch        *b_gengluon_eta;   //!
   TBranch        *b_gengluon_phi;   //!
   TBranch        *b_gengluon_pt;   //!
   TBranch        *b_gengluon_px;   //!
   TBranch        *b_gengluon_py;   //!
   TBranch        *b_gengluon_pz;   //!
   TBranch        *b_gengluon_status;   //!
   TBranch        *b_gengluon_charge;   //!
   TBranch        *b_gengluon_pdgid;   //!
   TBranch        *b_genele_e;   //!
   TBranch        *b_genele_eta;   //!
   TBranch        *b_genele_phi;   //!
   TBranch        *b_genele_pt;   //!
   TBranch        *b_genele_px;   //!
   TBranch        *b_genele_py;   //!
   TBranch        *b_genele_pz;   //!
   TBranch        *b_genele_charge;   //!
   TBranch        *b_unstableGenEle_e;   //!
   TBranch        *b_unstableGenEle_eta;   //!
   TBranch        *b_unstableGenEle_phi;   //!
   TBranch        *b_unstableGenEle_pt;   //!
   TBranch        *b_unstableGenEle_px;   //!
   TBranch        *b_unstableGenEle_py;   //!
   TBranch        *b_unstableGenEle_pz;   //!
   TBranch        *b_unstableGenEle_charge;   //!
   TBranch        *b_genelemom_e;   //!
   TBranch        *b_genelemom_eta;   //!
   TBranch        *b_genelemom_phi;   //!
   TBranch        *b_genelemom_pt;   //!
   TBranch        *b_genelemom_px;   //!
   TBranch        *b_genelemom_py;   //!
   TBranch        *b_genelemom_pz;   //!
   TBranch        *b_genelemom_charge;   //!
   TBranch        *b_genelemom_pdgid;   //!
   TBranch        *b_genelemom_mass;   //!
   TBranch        *b_x1quark;   //!
   TBranch        *b_x2quark;   //!
   TBranch        *b_genPair_mass;   //!
   TBranch        *b_trueNVtx;   //!
   TBranch        *b_nVtxBefore;   //!
   TBranch        *b_nVtxNow;   //!
   TBranch        *b_nVtxAfter;   //!

   AccforSpin0(TTree *tree=0);
   virtual ~AccforSpin0();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AccforSpin0_cxx
AccforSpin0::AccforSpin0(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
/*    if (tree == 0) { */
/*       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/lathomas/DYToEE_M-200_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/TOTAL.root"); */
/*       if (!f || !f->IsOpen()) { */
/*          f = new TFile("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/lathomas/DYToEE_M-200_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/TOTAL.root"); */
/*       } */
/*       TDirectory * dir = (TDirectory*)f->Get("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/lathomas/DYToEE_M-200_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/TOTAL.root:/gsfcheckerjob"); */
/*       dir->GetObject("tree",tree); */

/*    } */
/*    Init(tree); */
}

AccforSpin0::~AccforSpin0()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AccforSpin0::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AccforSpin0::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AccforSpin0::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runnumber", &runnumber, &b_runnumber);
   fChain->SetBranchAddress("eventnumber", &eventnumber, &b_eventnumber);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("eventcounter", &eventcounter, &b_eventcounter);
   fChain->SetBranchAddress("processid", &processid, &b_processid);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("alphaqcd", &alphaqcd, &b_alphaqcd);
   fChain->SetBranchAddress("alphaqed", &alphaqed, &b_alphaqed);
   fChain->SetBranchAddress("qscale", &qscale, &b_qscale);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("hltCount", &hltCount, &b_hltCount);
   fChain->SetBranchAddress("PhysDecl_bool", &PhysDecl_bool, &b_PhysDecl_bool);
   fChain->SetBranchAddress("nWasRun_", &nWasRun_, &b_nWasRun_);
   fChain->SetBranchAddress("nAccept_", &nAccept_, &b_nAccept_);
   fChain->SetBranchAddress("nErrors_", &nErrors_, &b_nErrors_);
   fChain->SetBranchAddress("HLT_Mu15_eta2p1", &HLT_Mu15_eta2p1, &b_HLT_Mu15_eta2p1);
   fChain->SetBranchAddress("HLT_Mu24_eta2p1", &HLT_Mu24_eta2p1, &b_HLT_Mu24_eta2p1);
   fChain->SetBranchAddress("HLT_Mu30_eta2p1", &HLT_Mu30_eta2p1, &b_HLT_Mu30_eta2p1);
   fChain->SetBranchAddress("HLT_Mu40_eta2p1", &HLT_Mu40_eta2p1, &b_HLT_Mu40_eta2p1);
   fChain->SetBranchAddress("HLT_Mu50_eta2p1", &HLT_Mu50_eta2p1, &b_HLT_Mu50_eta2p1);
   fChain->SetBranchAddress("HLT_Mu22_TkMu22", &HLT_Mu22_TkMu22, &b_HLT_Mu22_TkMu22);
   fChain->SetBranchAddress("HLT_Mu22_Photon22_CaloIdL", &HLT_Mu22_Photon22_CaloIdL, &b_HLT_Mu22_Photon22_CaloIdL);
   fChain->SetBranchAddress("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdL_CaloIsoVL", &HLT_Ele8_CaloIdL_CaloIsoVL, &b_HLT_Ele8_CaloIdL_CaloIsoVL);
   fChain->SetBranchAddress("HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL", &HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50", &HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50, &b_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50);
   fChain->SetBranchAddress("HLT_DoubleEle33_CaloIdL", &HLT_DoubleEle33_CaloIdL, &b_HLT_DoubleEle33_CaloIdL);
   fChain->SetBranchAddress("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", &HLT_DoubleEle33_CaloIdL_GsfTrkIdVL, &b_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL);
   fChain->SetBranchAddress("HLT_DoubleEle33_CaloIdT", &HLT_DoubleEle33_CaloIdT, &b_HLT_DoubleEle33_CaloIdT);
   fChain->SetBranchAddress("HLT_Photon20_CaloIdVL_IsoL", &HLT_Photon20_CaloIdVL_IsoL, &b_HLT_Photon20_CaloIdVL_IsoL);
   fChain->SetBranchAddress("HLT_Photon30_CaloIdVL", &HLT_Photon30_CaloIdVL, &b_HLT_Photon30_CaloIdVL);
   fChain->SetBranchAddress("HLT_Photon50_CaloIdVL", &HLT_Photon50_CaloIdVL, &b_HLT_Photon50_CaloIdVL);
   fChain->SetBranchAddress("HLT_Photon50_CaloIdVL_IsoL", &HLT_Photon50_CaloIdVL_IsoL, &b_HLT_Photon50_CaloIdVL_IsoL);
   fChain->SetBranchAddress("HLT_Photon75_CaloIdVL", &HLT_Photon75_CaloIdVL, &b_HLT_Photon75_CaloIdVL);
   fChain->SetBranchAddress("HLT_Photon90_CaloIdVL", &HLT_Photon90_CaloIdVL, &b_HLT_Photon90_CaloIdVL);
   fChain->SetBranchAddress("HLT_Photon135", &HLT_Photon135, &b_HLT_Photon135);
   fChain->SetBranchAddress("HLT_Photon150", &HLT_Photon150, &b_HLT_Photon150);
   fChain->SetBranchAddress("HLT_Photon250_NoHE", &HLT_Photon250_NoHE, &b_HLT_Photon250_NoHE);
   fChain->SetBranchAddress("HLT_Photon300_NoHE", &HLT_Photon300_NoHE, &b_HLT_Photon300_NoHE);
   fChain->SetBranchAddress("HLT_Photon26_Photon18", &HLT_Photon26_Photon18, &b_HLT_Photon26_Photon18);
   fChain->SetBranchAddress("HLT_Photon36_Photon22", &HLT_Photon36_Photon22, &b_HLT_Photon36_Photon22);
   fChain->SetBranchAddress("HLT_DoublePhoton70", &HLT_DoublePhoton70, &b_HLT_DoublePhoton70);
   fChain->SetBranchAddress("HLT_DoublePhoton80", &HLT_DoublePhoton80, &b_HLT_DoublePhoton80);
   fChain->SetBranchAddress("prescale_HLT_Mu15_eta2p1", &prescale_HLT_Mu15_eta2p1, &b_prescale_HLT_Mu15_eta2p1);
   fChain->SetBranchAddress("prescale_HLT_Mu30_eta2p1", &prescale_HLT_Mu30_eta2p1, &b_prescale_HLT_Mu30_eta2p1);
   fChain->SetBranchAddress("prescale_HLT_Mu40_eta2p1", &prescale_HLT_Mu40_eta2p1, &b_prescale_HLT_Mu40_eta2p1);
   fChain->SetBranchAddress("prescale_HLT_Mu22_TkMu22", &prescale_HLT_Mu22_TkMu22, &b_prescale_HLT_Mu22_TkMu22);
   fChain->SetBranchAddress("prescale_HLT_Mu22_Photon22_CaloIdL", &prescale_HLT_Mu22_Photon22_CaloIdL, &b_prescale_HLT_Mu22_Photon22_CaloIdL);
   fChain->SetBranchAddress("prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("prescale_HLT_Ele8_CaloIdL_CaloIsoVL", &prescale_HLT_Ele8_CaloIdL_CaloIsoVL, &b_prescale_HLT_Ele8_CaloIdL_CaloIsoVL);
   fChain->SetBranchAddress("prescale_HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL", &prescale_HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_prescale_HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("prescale_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50", &prescale_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50, &b_prescale_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50);
   fChain->SetBranchAddress("prescale_HLT_DoubleEle33_CaloIdL", &prescale_HLT_DoubleEle33_CaloIdL, &b_prescale_HLT_DoubleEle33_CaloIdL);
   fChain->SetBranchAddress("prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", &prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL, &b_prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL);
   fChain->SetBranchAddress("prescale_HLT_DoubleEle33_CaloIdT", &prescale_HLT_DoubleEle33_CaloIdT, &b_prescale_HLT_DoubleEle33_CaloIdT);
   fChain->SetBranchAddress("prescale_HLT_Photon20_CaloIdVL_IsoL", &prescale_HLT_Photon20_CaloIdVL_IsoL, &b_prescale_HLT_Photon20_CaloIdVL_IsoL);
   fChain->SetBranchAddress("prescale_HLT_Photon30_CaloIdVL", &prescale_HLT_Photon30_CaloIdVL, &b_prescale_HLT_Photon30_CaloIdVL);
   fChain->SetBranchAddress("prescale_HLT_Photon50_CaloIdVL", &prescale_HLT_Photon50_CaloIdVL, &b_prescale_HLT_Photon50_CaloIdVL);
   fChain->SetBranchAddress("prescale_HLT_Photon50_CaloIdVL_IsoL", &prescale_HLT_Photon50_CaloIdVL_IsoL, &b_prescale_HLT_Photon50_CaloIdVL_IsoL);
   fChain->SetBranchAddress("prescale_HLT_Photon75_CaloIdVL", &prescale_HLT_Photon75_CaloIdVL, &b_prescale_HLT_Photon75_CaloIdVL);
   fChain->SetBranchAddress("prescale_HLT_Photon90_CaloIdVL", &prescale_HLT_Photon90_CaloIdVL, &b_prescale_HLT_Photon90_CaloIdVL);
   fChain->SetBranchAddress("prescale_HLT_Photon135", &prescale_HLT_Photon135, &b_prescale_HLT_Photon135);
   fChain->SetBranchAddress("prescale_HLT_Photon150", &prescale_HLT_Photon150, &b_prescale_HLT_Photon150);
   fChain->SetBranchAddress("prescale_HLT_Photon250_NoHE", &prescale_HLT_Photon250_NoHE, &b_prescale_HLT_Photon250_NoHE);
   fChain->SetBranchAddress("prescale_HLT_Photon300_NoHE", &prescale_HLT_Photon300_NoHE, &b_prescale_HLT_Photon300_NoHE);
   fChain->SetBranchAddress("prescale_HLT_Photon26_Photon18", &prescale_HLT_Photon26_Photon18, &b_prescale_HLT_Photon26_Photon18);
   fChain->SetBranchAddress("prescale_HLT_Photon36_Photon22", &prescale_HLT_Photon36_Photon22, &b_prescale_HLT_Photon36_Photon22);
   fChain->SetBranchAddress("prescale_HLT_DoublePhoton70", &prescale_HLT_DoublePhoton70, &b_prescale_HLT_DoublePhoton70);
   fChain->SetBranchAddress("prescale_HLT_DoublePhoton80", &prescale_HLT_DoublePhoton80, &b_prescale_HLT_DoublePhoton80);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("calomet", &calomet, &b_calomet);
   fChain->SetBranchAddress("calomet_phi", &calomet_phi, &b_calomet_phi);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("pfmet", &pfmet, &b_pfmet);
   fChain->SetBranchAddress("pfmet_phi", &pfmet_phi, &b_pfmet_phi);
   fChain->SetBranchAddress("pfmetcor", &pfmetcor, &b_pfmetcor);
   fChain->SetBranchAddress("pfmetcor_phi", &pfmetcor_phi, &b_pfmetcor_phi);
   fChain->SetBranchAddress("sigmaZ", &sigmaZ, &b_sigmaZ);
   fChain->SetBranchAddress("sigmaZ0Error", &sigmaZ0Error, &b_sigmaZ0Error);
   fChain->SetBranchAddress("sq", &sq, &b_sq);
   fChain->SetBranchAddress("bsposx", &bsposx, &b_bsposx);
   fChain->SetBranchAddress("bsposy", &bsposy, &b_bsposy);
   fChain->SetBranchAddress("bsposz", &bsposz, &b_bsposz);
   fChain->SetBranchAddress("pvsize", &pvsize, &b_pvsize);
   fChain->SetBranchAddress("pvx", pvx, &b_pvx);
   fChain->SetBranchAddress("pvy", pvy, &b_pvy);
   fChain->SetBranchAddress("pvz", pvz, &b_pvz);
   fChain->SetBranchAddress("pv_isValid", pv_isValid, &b_pv_isValid);
   fChain->SetBranchAddress("pv_ndof", pv_ndof, &b_pv_ndof);
   fChain->SetBranchAddress("pv_nTracks", pv_nTracks, &b_pv_nTracks);
   fChain->SetBranchAddress("pv_normChi2", pv_normChi2, &b_pv_normChi2);
   fChain->SetBranchAddress("pv_totTrackSize", pv_totTrackSize, &b_pv_totTrackSize);
   fChain->SetBranchAddress("pfJetColl_size", &pfJetColl_size, &b_pfJetColl_size);
   fChain->SetBranchAddress("pfJet_pt", pfJet_pt, &b_pfJet_pt);
   fChain->SetBranchAddress("pfJet_eta", pfJet_eta, &b_pfJet_eta);
   fChain->SetBranchAddress("pfJet_phi", pfJet_phi, &b_pfJet_phi);
   fChain->SetBranchAddress("JetColl_size", &JetColl_size, &b_JetColl_size);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_vx", Jet_vx, &b_Jet_vx);
   fChain->SetBranchAddress("Jet_vy", Jet_vy, &b_Jet_vy);
   fChain->SetBranchAddress("Jet_vz", Jet_vz, &b_Jet_vz);
   fChain->SetBranchAddress("tCHighEffBTags", tCHighEffBTags, &b_tCHighEffBTags);
   fChain->SetBranchAddress("tCHighPurBTags", tCHighPurBTags, &b_tCHighPurBTags);
   fChain->SetBranchAddress("jetProbBTags", jetProbBTags, &b_jetProbBTags);
   fChain->SetBranchAddress("jetBProbBTags", jetBProbBTags, &b_jetBProbBTags);
   fChain->SetBranchAddress("sSecVertHighEffBTags", sSecVertHighEffBTags, &b_sSecVertHighEffBTags);
   fChain->SetBranchAddress("sSecVertHighPurBTags", sSecVertHighPurBTags, &b_sSecVertHighPurBTags);
   fChain->SetBranchAddress("cSecVertBTags", cSecVertBTags, &b_cSecVertBTags);
   fChain->SetBranchAddress("cSecVertMVABTags", cSecVertMVABTags, &b_cSecVertMVABTags);
   fChain->SetBranchAddress("ghostTrkBTags", ghostTrkBTags, &b_ghostTrkBTags);
   fChain->SetBranchAddress("softEleIP3dBTags", softEleIP3dBTags, &b_softEleIP3dBTags);
   fChain->SetBranchAddress("softElePtBTags", softElePtBTags, &b_softElePtBTags);
   fChain->SetBranchAddress("softMuBTags", softMuBTags, &b_softMuBTags);
   fChain->SetBranchAddress("softMuIP3dBTags", softMuIP3dBTags, &b_softMuIP3dBTags);
   fChain->SetBranchAddress("softMuPtBTags", softMuPtBTags, &b_softMuPtBTags);
   fChain->SetBranchAddress("muon_size", &muon_size, &b_muon_size);
   fChain->SetBranchAddress("muon_pt", muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon_ptError", muon_ptError, &b_muon_ptError);
   fChain->SetBranchAddress("muon_gTrk_pt", muon_gTrk_pt, &b_muon_gTrk_pt);
   fChain->SetBranchAddress("muon_gTrk_ptError", muon_gTrk_ptError, &b_muon_gTrk_ptError);
   fChain->SetBranchAddress("muon_eta", muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon_etaError", muon_etaError, &b_muon_etaError);
   fChain->SetBranchAddress("muon_theta", muon_theta, &b_muon_theta);
   fChain->SetBranchAddress("muon_thetaError", muon_thetaError, &b_muon_thetaError);
   fChain->SetBranchAddress("muon_phi", muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_phiError", muon_phiError, &b_muon_phiError);
   fChain->SetBranchAddress("muon_outerPt", muon_outerPt, &b_muon_outerPt);
   fChain->SetBranchAddress("muon_outerEta", muon_outerEta, &b_muon_outerEta);
   fChain->SetBranchAddress("muon_outerPhi", muon_outerPhi, &b_muon_outerPhi);
   fChain->SetBranchAddress("muon_outerTheta", muon_outerTheta, &b_muon_outerTheta);
   fChain->SetBranchAddress("muon_px", muon_px, &b_muon_px);
   fChain->SetBranchAddress("muon_py", muon_py, &b_muon_py);
   fChain->SetBranchAddress("muon_pz", muon_pz, &b_muon_pz);
   fChain->SetBranchAddress("muon_charge", muon_charge, &b_muon_charge);
   fChain->SetBranchAddress("muon_nhitspixel", muon_nhitspixel, &b_muon_nhitspixel);
   fChain->SetBranchAddress("muon_nhitstrack", muon_nhitstrack, &b_muon_nhitstrack);
   fChain->SetBranchAddress("muon_nhitsmuons", muon_nhitsmuons, &b_muon_nhitsmuons);
   fChain->SetBranchAddress("muon_nhitstotal", muon_nhitstotal, &b_muon_nhitstotal);
   fChain->SetBranchAddress("muon_nlayerswithhits", muon_nlayerswithhits, &b_muon_nlayerswithhits);
   fChain->SetBranchAddress("muon_nlosthits", muon_nlosthits, &b_muon_nlosthits);
   fChain->SetBranchAddress("muon_nSegmentMatch", muon_nSegmentMatch, &b_muon_nSegmentMatch);
   fChain->SetBranchAddress("muon_isTrackerMuon", muon_isTrackerMuon, &b_muon_isTrackerMuon);
   fChain->SetBranchAddress("muon_isPFMuon", muon_isPFMuon, &b_muon_isPFMuon);
   fChain->SetBranchAddress("muon_isPFIsolationValid", muon_isPFIsolationValid, &b_muon_isPFIsolationValid);
   fChain->SetBranchAddress("muon_chi2", muon_chi2, &b_muon_chi2);
   fChain->SetBranchAddress("muon_ndof", muon_ndof, &b_muon_ndof);
   fChain->SetBranchAddress("muon_normChi2", muon_normChi2, &b_muon_normChi2);
   fChain->SetBranchAddress("muon_d0", muon_d0, &b_muon_d0);
   fChain->SetBranchAddress("muon_d0Error", muon_d0Error, &b_muon_d0Error);
   fChain->SetBranchAddress("muon_dz_cmsCenter", muon_dz_cmsCenter, &b_muon_dz_cmsCenter);
   fChain->SetBranchAddress("muon_dz_beamSpot", muon_dz_beamSpot, &b_muon_dz_beamSpot);
   fChain->SetBranchAddress("muon_dz_firstPVtx", muon_dz_firstPVtx, &b_muon_dz_firstPVtx);
   fChain->SetBranchAddress("muon_dz_firstPVtxwithBS", muon_dz_firstPVtxwithBS, &b_muon_dz_firstPVtxwithBS);
   fChain->SetBranchAddress("muon_dzError", muon_dzError, &b_muon_dzError);
   fChain->SetBranchAddress("muon_dxy_cmsCenter", muon_dxy_cmsCenter, &b_muon_dxy_cmsCenter);
   fChain->SetBranchAddress("muon_dxy_beamSpot", muon_dxy_beamSpot, &b_muon_dxy_beamSpot);
   fChain->SetBranchAddress("muon_dxy_firstPVtx", muon_dxy_firstPVtx, &b_muon_dxy_firstPVtx);
   fChain->SetBranchAddress("muon_dxy_firstPVtxwithBS", muon_dxy_firstPVtxwithBS, &b_muon_dxy_firstPVtxwithBS);
   fChain->SetBranchAddress("muon_dxyError", muon_dxyError, &b_muon_dxyError);
   fChain->SetBranchAddress("muon_innerPosx", muon_innerPosx, &b_muon_innerPosx);
   fChain->SetBranchAddress("muon_innerPosy", muon_innerPosy, &b_muon_innerPosy);
   fChain->SetBranchAddress("muon_innerPosz", muon_innerPosz, &b_muon_innerPosz);
   fChain->SetBranchAddress("muon_trackIso03", muon_trackIso03, &b_muon_trackIso03);
   fChain->SetBranchAddress("muon_trackIso05", muon_trackIso05, &b_muon_trackIso05);
   fChain->SetBranchAddress("muon_trackIso03_ptInVeto", muon_trackIso03_ptInVeto, &b_muon_trackIso03_ptInVeto);
   fChain->SetBranchAddress("muon_trackIso05_ptInVeto", muon_trackIso05_ptInVeto, &b_muon_trackIso05_ptInVeto);
   fChain->SetBranchAddress("muon_emIso03", muon_emIso03, &b_muon_emIso03);
   fChain->SetBranchAddress("muon_emIso05", muon_emIso05, &b_muon_emIso05);
   fChain->SetBranchAddress("muon_emIso03_ptInVeto", muon_emIso03_ptInVeto, &b_muon_emIso03_ptInVeto);
   fChain->SetBranchAddress("muon_emIso05_ptInVeto", muon_emIso05_ptInVeto, &b_muon_emIso05_ptInVeto);
   fChain->SetBranchAddress("muon_hadIso03", muon_hadIso03, &b_muon_hadIso03);
   fChain->SetBranchAddress("muon_hadIso05", muon_hadIso05, &b_muon_hadIso05);
   fChain->SetBranchAddress("muon_hadIso03_ptInVeto", muon_hadIso03_ptInVeto, &b_muon_hadIso03_ptInVeto);
   fChain->SetBranchAddress("muon_hadIso05_ptInVeto", muon_hadIso05_ptInVeto, &b_muon_hadIso05_ptInVeto);
   fChain->SetBranchAddress("scsize", &scsize, &b_scsize);
   fChain->SetBranchAddress("scenergy", scenergy, &b_scenergy);
   fChain->SetBranchAddress("sceta", sceta, &b_sceta);
   fChain->SetBranchAddress("scetacorr", scetacorr, &b_scetacorr);
   fChain->SetBranchAddress("scet", scet, &b_scet);
   fChain->SetBranchAddress("scphi", scphi, &b_scphi);
   fChain->SetBranchAddress("scpx", scpx, &b_scpx);
   fChain->SetBranchAddress("scpy", scpy, &b_scpy);
   fChain->SetBranchAddress("scpz", scpz, &b_scpz);
   fChain->SetBranchAddress("scx", scx, &b_scx);
   fChain->SetBranchAddress("scy", scy, &b_scy);
   fChain->SetBranchAddress("scz", scz, &b_scz);
   fChain->SetBranchAddress("pfele_size", &pfele_size, &b_pfele_size);
   fChain->SetBranchAddress("pfele_pt", pfele_pt, &b_pfele_pt);
   fChain->SetBranchAddress("pfele_eta", pfele_eta, &b_pfele_eta);
   fChain->SetBranchAddress("pfele_phi", pfele_phi, &b_pfele_phi);
   fChain->SetBranchAddress("pfele_charge", pfele_charge, &b_pfele_charge);
   fChain->SetBranchAddress("gsf_size", &gsf_size, &b_gsf_size);
   fChain->SetBranchAddress("gsf_isEB", gsf_isEB, &b_gsf_isEB);
   fChain->SetBranchAddress("gsf_isEE", gsf_isEE, &b_gsf_isEE);
   fChain->SetBranchAddress("gsf_px", gsf_px, &b_gsf_px);
   fChain->SetBranchAddress("gsf_py", gsf_py, &b_gsf_py);
   fChain->SetBranchAddress("gsf_pz", gsf_pz, &b_gsf_pz);
   fChain->SetBranchAddress("gsf_pt", gsf_pt, &b_gsf_pt);
   fChain->SetBranchAddress("gsf_eta", gsf_eta, &b_gsf_eta);
   fChain->SetBranchAddress("gsf_phi", gsf_phi, &b_gsf_phi);
   fChain->SetBranchAddress("gsf_theta", gsf_theta, &b_gsf_theta);
   fChain->SetBranchAddress("gsf_charge", gsf_charge, &b_gsf_charge);
   fChain->SetBranchAddress("gsf_deltaEtaATcalo", gsf_deltaEtaATcalo, &b_gsf_deltaEtaATcalo);
   fChain->SetBranchAddress("gsf_deltaPhiATcalo", gsf_deltaPhiATcalo, &b_gsf_deltaPhiATcalo);
   fChain->SetBranchAddress("gsf_sigmaetaeta", gsf_sigmaetaeta, &b_gsf_sigmaetaeta);
   fChain->SetBranchAddress("gsf_sigmaIetaIeta", gsf_sigmaIetaIeta, &b_gsf_sigmaIetaIeta);
   fChain->SetBranchAddress("gsf_ecalEnergy", gsf_ecalEnergy, &b_gsf_ecalEnergy);
   fChain->SetBranchAddress("gsf_eOVERp", gsf_eOVERp, &b_gsf_eOVERp);
   fChain->SetBranchAddress("gsf_dxy", gsf_dxy, &b_gsf_dxy);
   fChain->SetBranchAddress("gsf_dxy_beamSpot", gsf_dxy_beamSpot, &b_gsf_dxy_beamSpot);
   fChain->SetBranchAddress("gsf_dxy_firstPVtx", gsf_dxy_firstPVtx, &b_gsf_dxy_firstPVtx);
   fChain->SetBranchAddress("gsf_dxy_firstPVtxwithBS", gsf_dxy_firstPVtxwithBS, &b_gsf_dxy_firstPVtxwithBS);
   fChain->SetBranchAddress("gsf_dxyError", gsf_dxyError, &b_gsf_dxyError);
   fChain->SetBranchAddress("gsf_dz", gsf_dz, &b_gsf_dz);
   fChain->SetBranchAddress("gsf_dz_beamSpot", gsf_dz_beamSpot, &b_gsf_dz_beamSpot);
   fChain->SetBranchAddress("gsf_dz_firstPVtx", gsf_dz_firstPVtx, &b_gsf_dz_firstPVtx);
   fChain->SetBranchAddress("gsf_dz_firstPVtxwithBS", gsf_dz_firstPVtxwithBS, &b_gsf_dz_firstPVtxwithBS);
   fChain->SetBranchAddress("gsf_dzError", gsf_dzError, &b_gsf_dzError);
   fChain->SetBranchAddress("gsf_vz", gsf_vz, &b_gsf_vz);
   fChain->SetBranchAddress("gsf_nHits", gsf_nHits, &b_gsf_nHits);
   fChain->SetBranchAddress("gsf_nLostInnerHits", gsf_nLostInnerHits, &b_gsf_nLostInnerHits);
   fChain->SetBranchAddress("gsf_nLostOuterHits", gsf_nLostOuterHits, &b_gsf_nLostOuterHits);
   fChain->SetBranchAddress("gsf_convFlags", gsf_convFlags, &b_gsf_convFlags);
   fChain->SetBranchAddress("gsf_convDist", gsf_convDist, &b_gsf_convDist);
   fChain->SetBranchAddress("gsf_convDcot", gsf_convDcot, &b_gsf_convDcot);
   fChain->SetBranchAddress("gsf_convRadius", gsf_convRadius, &b_gsf_convRadius);
   fChain->SetBranchAddress("gsf_fBrem", gsf_fBrem, &b_gsf_fBrem);
   fChain->SetBranchAddress("gsf_e1x5", gsf_e1x5, &b_gsf_e1x5);
   fChain->SetBranchAddress("gsf_e2x5", gsf_e2x5, &b_gsf_e2x5);
   fChain->SetBranchAddress("gsf_e5x5", gsf_e5x5, &b_gsf_e5x5);
   fChain->SetBranchAddress("gsf_e1x3", gsf_e1x3, &b_gsf_e1x3);
   fChain->SetBranchAddress("gsf_p", gsf_p, &b_gsf_p);
   fChain->SetBranchAddress("gsf_e", gsf_e, &b_gsf_e);
   fChain->SetBranchAddress("gsf_deltaeta", gsf_deltaeta, &b_gsf_deltaeta);
   fChain->SetBranchAddress("gsf_deltaphi", gsf_deltaphi, &b_gsf_deltaphi);
   fChain->SetBranchAddress("gsf_hovere", gsf_hovere, &b_gsf_hovere);
   fChain->SetBranchAddress("gsf_hdepth1overe", gsf_hdepth1overe, &b_gsf_hdepth1overe);
   fChain->SetBranchAddress("gsf_hdepth2overe", gsf_hdepth2overe, &b_gsf_hdepth2overe);
   fChain->SetBranchAddress("gsf_hovere2012", gsf_hovere2012, &b_gsf_hovere2012);
   fChain->SetBranchAddress("gsf_hdepth1overe2012", gsf_hdepth1overe2012, &b_gsf_hdepth1overe2012);
   fChain->SetBranchAddress("gsf_hdepth2overe2012", gsf_hdepth2overe2012, &b_gsf_hdepth2overe2012);
   fChain->SetBranchAddress("gsf_trackiso", gsf_trackiso, &b_gsf_trackiso);
   fChain->SetBranchAddress("gsf_ecaliso", gsf_ecaliso, &b_gsf_ecaliso);
   fChain->SetBranchAddress("gsf_hcaliso1", gsf_hcaliso1, &b_gsf_hcaliso1);
   fChain->SetBranchAddress("gsf_hcaliso2", gsf_hcaliso2, &b_gsf_hcaliso2);
   fChain->SetBranchAddress("gsf_hcaliso12012", gsf_hcaliso12012, &b_gsf_hcaliso12012);
   fChain->SetBranchAddress("gsf_hcaliso22012", gsf_hcaliso22012, &b_gsf_hcaliso22012);
   fChain->SetBranchAddress("gsf_PFisocharged", gsf_PFisocharged, &b_gsf_PFisocharged);
   fChain->SetBranchAddress("gsf_PFisophoton", gsf_PFisophoton, &b_gsf_PFisophoton);
   fChain->SetBranchAddress("gsf_PFisoneutral", gsf_PFisoneutral, &b_gsf_PFisoneutral);
   fChain->SetBranchAddress("gsf_class", gsf_class, &b_gsf_class);
   fChain->SetBranchAddress("gsf_isecaldriven", gsf_isecaldriven, &b_gsf_isecaldriven);
   fChain->SetBranchAddress("gsf_istrackerdriven", gsf_istrackerdriven, &b_gsf_istrackerdriven);
   fChain->SetBranchAddress("gsfsc_e", gsfsc_e, &b_gsfsc_e);
   fChain->SetBranchAddress("gsfsc_pt", gsfsc_pt, &b_gsfsc_pt);
   fChain->SetBranchAddress("gsfsc_eta", gsfsc_eta, &b_gsfsc_eta);
   fChain->SetBranchAddress("gsfsc_phi", gsfsc_phi, &b_gsfsc_phi);
   fChain->SetBranchAddress("gsfsc_px", gsfsc_px, &b_gsfsc_px);
   fChain->SetBranchAddress("gsfsc_py", gsfsc_py, &b_gsfsc_py);
   fChain->SetBranchAddress("gsfsc_pz", gsfsc_pz, &b_gsfsc_pz);
   fChain->SetBranchAddress("gsf_e2x5overe5x5", gsf_e2x5overe5x5, &b_gsf_e2x5overe5x5);
   fChain->SetBranchAddress("gsf_e1x5overe5x5", gsf_e1x5overe5x5, &b_gsf_e1x5overe5x5);
   fChain->SetBranchAddress("gsf_gsfet", gsf_gsfet, &b_gsf_gsfet);
   fChain->SetBranchAddress("gsf_hitsinfo", gsf_hitsinfo, &b_gsf_hitsinfo);
   fChain->SetBranchAddress("scindexforgsf", scindexforgsf, &b_scindexforgsf);
   fChain->SetBranchAddress("gsftracksize", &gsftracksize, &b_gsftracksize);
   fChain->SetBranchAddress("gsftracketa", gsftracketa, &b_gsftracketa);
   fChain->SetBranchAddress("gsftrackphi", gsftrackphi, &b_gsftrackphi);
   fChain->SetBranchAddress("gsftrackp", gsftrackp, &b_gsftrackp);
   fChain->SetBranchAddress("gsftrackpt", gsftrackpt, &b_gsftrackpt);
   fChain->SetBranchAddress("gsftrackpx", gsftrackpx, &b_gsftrackpx);
   fChain->SetBranchAddress("gsftrackpy", gsftrackpy, &b_gsftrackpy);
   fChain->SetBranchAddress("gsftrackpz", gsftrackpz, &b_gsftrackpz);
   fChain->SetBranchAddress("gsfpass_ET", gsfpass_ET, &b_gsfpass_ET);
   fChain->SetBranchAddress("gsfpass_PT", gsfpass_PT, &b_gsfpass_PT);
   fChain->SetBranchAddress("gsfpass_DETETA", gsfpass_DETETA, &b_gsfpass_DETETA);
   fChain->SetBranchAddress("gsfpass_CRACK", gsfpass_CRACK, &b_gsfpass_CRACK);
   fChain->SetBranchAddress("gsfpass_DETAIN", gsfpass_DETAIN, &b_gsfpass_DETAIN);
   fChain->SetBranchAddress("gsfpass_DPHIIN", gsfpass_DPHIIN, &b_gsfpass_DPHIIN);
   fChain->SetBranchAddress("gsfpass_HADEM", gsfpass_HADEM, &b_gsfpass_HADEM);
   fChain->SetBranchAddress("gsfpass_SIGMAIETAIETA", gsfpass_SIGMAIETAIETA, &b_gsfpass_SIGMAIETAIETA);
   fChain->SetBranchAddress("gsfpass_E2X5OVER5X5", gsfpass_E2X5OVER5X5, &b_gsfpass_E2X5OVER5X5);
   fChain->SetBranchAddress("gsfpass_ISOLEMHADDEPTH1", gsfpass_ISOLEMHADDEPTH1, &b_gsfpass_ISOLEMHADDEPTH1);
   fChain->SetBranchAddress("gsfpass_ISOLHADDEPTH2", gsfpass_ISOLHADDEPTH2, &b_gsfpass_ISOLHADDEPTH2);
   fChain->SetBranchAddress("gsfpass_ISOLPTTRKS", gsfpass_ISOLPTTRKS, &b_gsfpass_ISOLPTTRKS);
   fChain->SetBranchAddress("gsfpass_ECALDRIVEN", gsfpass_ECALDRIVEN, &b_gsfpass_ECALDRIVEN);
   fChain->SetBranchAddress("gsfpass_INVALID", gsfpass_INVALID, &b_gsfpass_INVALID);
   fChain->SetBranchAddress("gsfpass_NOMISSINGHITS", gsfpass_NOMISSINGHITS, &b_gsfpass_NOMISSINGHITS);
   fChain->SetBranchAddress("gsfpass_NOCONVERSION", gsfpass_NOCONVERSION, &b_gsfpass_NOCONVERSION);
   fChain->SetBranchAddress("gsfpass_DXYFIRSTPV", gsfpass_DXYFIRSTPV, &b_gsfpass_DXYFIRSTPV);
   fChain->SetBranchAddress("gsfpass_HEEP", gsfpass_HEEP, &b_gsfpass_HEEP);
   fChain->SetBranchAddress("gsfpass_ID", gsfpass_ID, &b_gsfpass_ID);
   fChain->SetBranchAddress("gsfpass_ISO", gsfpass_ISO, &b_gsfpass_ISO);
   fChain->SetBranchAddress("gsfmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter", gsfmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter, &b_gsfmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter);
   fChain->SetBranchAddress("gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter", gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter, &b_gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter);
   fChain->SetBranchAddress("gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter", gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter, &b_gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter);
   fChain->SetBranchAddress("gsfmatch_hltL1sL1SingleEG22", gsfmatch_hltL1sL1SingleEG22, &b_gsfmatch_hltL1sL1SingleEG22);
   fChain->SetBranchAddress("gsfmatch_hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter", gsfmatch_hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter, &b_gsfmatch_hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter);
   fChain->SetBranchAddress("gsfmatch_hltEle33CaloIdLPixelMatchFilter", gsfmatch_hltEle33CaloIdLPixelMatchFilter, &b_gsfmatch_hltEle33CaloIdLPixelMatchFilter);
   fChain->SetBranchAddress("scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter", scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter, &b_scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter);
   fChain->SetBranchAddress("heepHeepMass", &heepHeepMass, &b_heepHeepMass);
   fChain->SetBranchAddress("scpixcharge", scpixcharge, &b_scpixcharge);
   fChain->SetBranchAddress("ctfcharge", ctfcharge, &b_ctfcharge);
   fChain->SetBranchAddress("gsfcharge", gsfcharge, &b_gsfcharge);
   fChain->SetBranchAddress("gsfctfscpixconsistent", gsfctfscpixconsistent, &b_gsfctfscpixconsistent);
   fChain->SetBranchAddress("gsfscpixconsistent", gsfscpixconsistent, &b_gsfscpixconsistent);
   fChain->SetBranchAddress("gsfctfconsistent", gsfctfconsistent, &b_gsfctfconsistent);
   fChain->SetBranchAddress("genparticles_size", &genparticles_size, &b_genparticles_size);
   fChain->SetBranchAddress("genquarks_size", &genquarks_size, &b_genquarks_size);
   fChain->SetBranchAddress("gengluons_size", &gengluons_size, &b_gengluons_size);
   fChain->SetBranchAddress("conv_size", &conv_size, &b_conv_size);
   fChain->SetBranchAddress("conv_vtxProb", conv_vtxProb, &b_conv_vtxProb);
   fChain->SetBranchAddress("conv_lxy", conv_lxy, &b_conv_lxy);
   fChain->SetBranchAddress("conv_nHitsMax", conv_nHitsMax, &b_conv_nHitsMax);
   fChain->SetBranchAddress("conv_eleind", conv_eleind, &b_conv_eleind);
   fChain->SetBranchAddress("gsf0_crystal_size", &gsf0_crystal_size, &b_gsf0_crystal_size);
   fChain->SetBranchAddress("gsf0_crystal_ietaorix", gsf0_crystal_ietaorix, &b_gsf0_crystal_ietaorix);
   fChain->SetBranchAddress("gsf0_crystal_iphioriy", gsf0_crystal_iphioriy, &b_gsf0_crystal_iphioriy);
   fChain->SetBranchAddress("gsf0_crystal_energy", gsf0_crystal_energy, &b_gsf0_crystal_energy);
   fChain->SetBranchAddress("gsf0_crystal_eta", gsf0_crystal_eta, &b_gsf0_crystal_eta);
   fChain->SetBranchAddress("gsf1_crystal_size", &gsf1_crystal_size, &b_gsf1_crystal_size);
   fChain->SetBranchAddress("gsf1_crystal_ietaorix", gsf1_crystal_ietaorix, &b_gsf1_crystal_ietaorix);
   fChain->SetBranchAddress("gsf1_crystal_iphioriy", gsf1_crystal_iphioriy, &b_gsf1_crystal_iphioriy);
   fChain->SetBranchAddress("gsf1_crystal_energy", gsf1_crystal_energy, &b_gsf1_crystal_energy);
   fChain->SetBranchAddress("gsf1_crystal_eta", gsf1_crystal_eta, &b_gsf1_crystal_eta);
   fChain->SetBranchAddress("genquark_e", genquark_e, &b_genquark_e);
   fChain->SetBranchAddress("genquark_eta", genquark_eta, &b_genquark_eta);
   fChain->SetBranchAddress("genquark_phi", genquark_phi, &b_genquark_phi);
   fChain->SetBranchAddress("genquark_pt", genquark_pt, &b_genquark_pt);
   fChain->SetBranchAddress("genquark_px", genquark_px, &b_genquark_px);
   fChain->SetBranchAddress("genquark_py", genquark_py, &b_genquark_py);
   fChain->SetBranchAddress("genquark_pz", genquark_pz, &b_genquark_pz);
   fChain->SetBranchAddress("genquark_status", genquark_status, &b_genquark_status);
   fChain->SetBranchAddress("genquark_charge", genquark_charge, &b_genquark_charge);
   fChain->SetBranchAddress("genquark_pdgid", genquark_pdgid, &b_genquark_pdgid);
   fChain->SetBranchAddress("gengluon_e", gengluon_e, &b_gengluon_e);
   fChain->SetBranchAddress("gengluon_eta", gengluon_eta, &b_gengluon_eta);
   fChain->SetBranchAddress("gengluon_phi", gengluon_phi, &b_gengluon_phi);
   fChain->SetBranchAddress("gengluon_pt", gengluon_pt, &b_gengluon_pt);
   fChain->SetBranchAddress("gengluon_px", gengluon_px, &b_gengluon_px);
   fChain->SetBranchAddress("gengluon_py", gengluon_py, &b_gengluon_py);
   fChain->SetBranchAddress("gengluon_pz", gengluon_pz, &b_gengluon_pz);
   fChain->SetBranchAddress("gengluon_status", gengluon_status, &b_gengluon_status);
   fChain->SetBranchAddress("gengluon_charge", gengluon_charge, &b_gengluon_charge);
   fChain->SetBranchAddress("gengluon_pdgid", gengluon_pdgid, &b_gengluon_pdgid);
   fChain->SetBranchAddress("genele_e", genele_e, &b_genele_e);
   fChain->SetBranchAddress("genele_eta", genele_eta, &b_genele_eta);
   fChain->SetBranchAddress("genele_phi", genele_phi, &b_genele_phi);
   fChain->SetBranchAddress("genele_pt", genele_pt, &b_genele_pt);
   fChain->SetBranchAddress("genele_px", genele_px, &b_genele_px);
   fChain->SetBranchAddress("genele_py", genele_py, &b_genele_py);
   fChain->SetBranchAddress("genele_pz", genele_pz, &b_genele_pz);
   fChain->SetBranchAddress("genele_charge", genele_charge, &b_genele_charge);
   fChain->SetBranchAddress("unstableGenEle_e", unstableGenEle_e, &b_unstableGenEle_e);
   fChain->SetBranchAddress("unstableGenEle_eta", unstableGenEle_eta, &b_unstableGenEle_eta);
   fChain->SetBranchAddress("unstableGenEle_phi", unstableGenEle_phi, &b_unstableGenEle_phi);
   fChain->SetBranchAddress("unstableGenEle_pt", unstableGenEle_pt, &b_unstableGenEle_pt);
   fChain->SetBranchAddress("unstableGenEle_px", unstableGenEle_px, &b_unstableGenEle_px);
   fChain->SetBranchAddress("unstableGenEle_py", unstableGenEle_py, &b_unstableGenEle_py);
   fChain->SetBranchAddress("unstableGenEle_pz", unstableGenEle_pz, &b_unstableGenEle_pz);
   fChain->SetBranchAddress("unstableGenEle_charge", unstableGenEle_charge, &b_unstableGenEle_charge);
   fChain->SetBranchAddress("genelemom_e", genelemom_e, &b_genelemom_e);
   fChain->SetBranchAddress("genelemom_eta", genelemom_eta, &b_genelemom_eta);
   fChain->SetBranchAddress("genelemom_phi", genelemom_phi, &b_genelemom_phi);
   fChain->SetBranchAddress("genelemom_pt", genelemom_pt, &b_genelemom_pt);
   fChain->SetBranchAddress("genelemom_px", genelemom_px, &b_genelemom_px);
   fChain->SetBranchAddress("genelemom_py", genelemom_py, &b_genelemom_py);
   fChain->SetBranchAddress("genelemom_pz", genelemom_pz, &b_genelemom_pz);
   fChain->SetBranchAddress("genelemom_charge", genelemom_charge, &b_genelemom_charge);
   fChain->SetBranchAddress("genelemom_pdgid", genelemom_pdgid, &b_genelemom_pdgid);
   fChain->SetBranchAddress("genelemom_mass", genelemom_mass, &b_genelemom_mass);
   fChain->SetBranchAddress("x1quark", x1quark, &b_x1quark);
   fChain->SetBranchAddress("x2quark", x2quark, &b_x2quark);
   fChain->SetBranchAddress("genPair_mass", &genPair_mass, &b_genPair_mass);
   fChain->SetBranchAddress("trueNVtx", &trueNVtx, &b_trueNVtx);
   fChain->SetBranchAddress("nVtxBefore", &nVtxBefore, &b_nVtxBefore);
   fChain->SetBranchAddress("nVtxNow", &nVtxNow, &b_nVtxNow);
   fChain->SetBranchAddress("nVtxAfter", &nVtxAfter, &b_nVtxAfter);







   tree->SetBranchStatus("*",0); //disable all branches
   tree->SetBranchStatus("gsf_gsfet",1); //enable used branches
   tree->SetBranchStatus("gsfsc_eta",1);
   tree->SetBranchStatus("gsfsc_phi",1);
   tree->SetBranchStatus("gsf_size",1);
   tree->SetBranchStatus("gsf_phi",1);
   tree->SetBranchStatus("gsf_eta",1);

   tree->SetBranchStatus("gsf_charge",1);
   tree->SetBranchStatus("gsf_isecaldriven",1);
   tree->SetBranchStatus("gsf_hovere",1);
   tree->SetBranchStatus("gsf_deltaeta",1);
   tree->SetBranchStatus("gsf_deltaphi",1);
   tree->SetBranchStatus("gsf_e2x5overe5x5",1);
   tree->SetBranchStatus("gsf_e1x5overe5x5",1);
   tree->SetBranchStatus("gsf_e5x5",1);
   tree->SetBranchStatus("gsf_ecaliso",1);
   tree->SetBranchStatus("gsf_hcaliso1",1);
   tree->SetBranchStatus("gsf_trackiso",1);
   tree->SetBranchStatus("gsf_nLostInnerHits",1);
   tree->SetBranchStatus("gsf_dxy_firstPVtx",1);
   tree->SetBranchStatus("gsf_sigmaIetaIeta",1);
  

   tree->SetBranchStatus("gsfpass_HEEP",1);
   tree->SetBranchStatus("genelemom_mass",1); 
 tree->SetBranchStatus("genquark_pdgid",1);
 tree->SetBranchStatus("genquark_pz",1);
  tree->SetBranchStatus("genquark_px",1);
   tree->SetBranchStatus("genquark_py",1);
   tree->SetBranchStatus("genquark_e",1);
    tree->SetBranchStatus("genquark_pt",1);
   tree->SetBranchStatus("genquark_phi",1);
   tree->SetBranchStatus("genquark_eta",1);
   
   tree->SetBranchStatus("gengluon_pz",1);
  tree->SetBranchStatus("gengluon_px",1);
   tree->SetBranchStatus("gengluon_py",1);
   tree->SetBranchStatus("gengluon_e",1);
 tree->SetBranchStatus("gengluon_pt",1);
   tree->SetBranchStatus("gengluon_eta",1);
tree->SetBranchStatus("gengluon_phi",1);

 tree->SetBranchStatus("genelemom_pz",1);
 tree->SetBranchStatus("genelemom_pt",1);
 tree->SetBranchStatus("unstableGenEle_charge",1);
 tree->SetBranchStatus("unstableGenEle_pz",1);
 tree->SetBranchStatus("unstableGenEle_pt",1);
 tree->SetBranchStatus("unstableGenEle_e",1);
 tree->SetBranchStatus("unstableGenEle_eta",1);
 tree->SetBranchStatus("unstableGenEle_phi",1);
 tree->SetBranchStatus("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL",1);
 tree->SetBranchStatus("gengluons_size",1);
 tree->SetBranchStatus("rho",1);
/*  tree->SetBranchStatus("",1); */
/*  tree->SetBranchStatus("",1); */
/*  tree->SetBranchStatus("",1); */
/*  tree->SetBranchStatus("",1); */

/*  tree->SetBranchStatus("",1); */
/*  tree->SetBranchStatus("",1); */
/*  tree->SetBranchStatus("",1); */
/*  tree->SetBranchStatus("",1); */
/*  tree->SetBranchStatus("",1); */
/*  tree->SetBranchStatus("",1); */
/*  tree->SetBranchStatus("",1); */
/*  tree->SetBranchStatus("",1); */






   Notify();
}

Bool_t AccforSpin0::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AccforSpin0::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AccforSpin0::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AccforSpin0_cxx
