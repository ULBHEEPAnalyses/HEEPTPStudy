#define AccforSpin0_cxx
#include "AccforSpin0.h"

#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include "TMultiGraph.h"
#include <iostream>
#include <stdio.h>
#include "TImage.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TLegend.h"
#include "TH1.h"
#include "TPaveLabel.h"
#include<sstream>
#include <TDCacheFile.h>
//#include "/beo5/charafbis/HomeSafe03062010/CMSSW_2_2_6/src/HEEPRoutine/pstyle.C"
#include "/user/lathomas/TandP/CMSSW_3_8_4/src/DataFormats/Math/interface/deltaR.h"

vector<TDCacheFile *> input;
vector <double> sweight;
vector <TString> samplename;
//TFile *f1 = new TFile("Test13TeV_ggvsqqbar.root","RECREATE");
//TFile *f1 = new TFile("Phase1noaging.root","RECREATE");
//TFile *f1 = new TFile("Phase11kaging.root","RECREATE");
TFile *f1 = new TFile("Phase2.root","RECREATE");

TPaveLabel *labelcmssimul = new TPaveLabel(0.6,0.7,0.8,0.8,"CMS simulation","brNDC");
TPaveLabel *labelCMSPrel = new TPaveLabel(0.4,0.7,0.6,0.8,"CMS Preliminary, #sqrt{s} = 14 TeV","brNDC");



vector<TH1F *> hcosthetastar;
vector<TH1F *> hfabscosthetastar;
vector<TH1F *> hcosthetastar_meas;
vector<TH1F *> hfabscosthetastar_meas;
vector<TH1F *> hcosthetastar0gluon;
vector<TH1F *> hfabscosthetastar0gluon;
vector<TH1F *> hcosthetastar1gluon;
vector<TH1F *> hfabscosthetastar1gluon;
vector<TH1F *> hcosthetastar2gluons;
vector<TH1F *> hfabscosthetastar2gluons;


float etabarrel = 1.442; 
float etaendcapslow = 1.56;
float etaendcapshigh = 3.0;
float ptcut = 35;

float rapidityheepheep = -1000; 
float costhetaCSheepheep = -1000;
float costhetastar = 0;
float costhetastar_meas = 0; 

//The following are the variables stored in the output trees
float genmass = 0;
float recomass = 0; 
int nbofgengluons = 0; 
float weightevt = 0; //Weight of each event, needed for the reweighted samples (e.g. Spin 0, AFB1,...). See below.
bool genelespt35 = false;
bool genelesBB = false;
bool genelesBE = false;
bool genelesEE = false;
bool genelesGap = false;
bool genelesBeamPipe = false;
bool gsfptcut = false;
bool gsfBB = false;
bool gsfBE = false; 
bool gsfEE = false;
bool hltfired = false; 
bool heepheep = false; 
bool mrecomgenmatching = false;



void AccforSpin0::Loop()
{
  //Some style attributes for the labels, not important
  labelcmssimul->SetFillColor(0);
  labelcmssimul->SetFillStyle(0);
  labelcmssimul->SetBorderSize(0);
  labelcmssimul->SetTextSize(0.35); 
  labelCMSPrel->SetFillColor(0);
  labelCMSPrel->SetFillStyle(0);
  labelCMSPrel->SetBorderSize(0);
  labelCMSPrel->SetTextSize(0.35); 

  /*Declaration of the trees corresponding to various processes studied
    TreeDY: Drell-Yan events, 
    TreeSpin0qqbar: A simulation of a spin 0 resonance, built by reweighting Z'Psi events to the right cos theta distribution
    TreeSpin0gg: A simulation of a spin 0 resonance, built by reweighting RS graviton (produced through gg) events to the right cos theta distribution
    TreeSpin0qqbarggmix: A simulation of a spin 0 resonance, built by reweighting RS graviton (produced through 50%gg-50%qqbar) events to the right cos theta distribution
    TreeSpin2qqbar, TreeSpin2gg, TreeSpin2qqbarggmix: same as spin 0, but for a spin 2. 
    TreeRSGrav: RS graviton sample
    TreeZprimePSI: Z' Psi sample
    TreeZprimeSSM: Z' SSM sample
    TreeADD: ADD graviton sample (used a long time ago).
    TreeAFB1: A simulation of a spin 1 with a maximally positive AFB (AFB=+0.75, the name is misleading...). Until recently I used Z'Psi reweighted events to generate this distribution. I recently had to move to Z'SSM events since only this sample is available for upgraded geomertries. 
    TreeAFBmin1: Same as TreeAFB1 for a maximally negative AFB (AFB=-0.75).

    N.B.: At the moment, only the following samples are filled: TreeZprimeSSM, TreeAFB1, TreeAFBmin1.
    The upgraded samples do not exist for the other samples...
  */ 
  TTree *TreeDY = new TTree("TreeDY","TreeDY");
  TTree *TreeSpin0qqbar = new TTree("TreeSpin0qqbar","TreeSpin0qqbar");
  TTree *TreeSpin0gg = new TTree("TreeSpin0gg","TreeSpin0gg");
  TTree *TreeSpin0qqbarggmix = new TTree("TreeSpin0qqbarggmix","TreeSpin0qqbarggmix");
  TTree *TreeSpin2qqbar = new TTree("TreeSpin2qqbar","TreeSpin2qqbar");
  TTree *TreeSpin2gg = new TTree("TreeSpin2gg","TreeSpin2gg");
  TTree *TreeSpin2qqbarggmix = new TTree("TreeSpin2qqbarggmix","TreeSpin2qqbarggmix");
  TTree *TreeRSGrav = new TTree("TreeRSGrav","TreeRSGrav");
  TTree *TreeZprimePSI = new TTree("TreeZprimePSI","TreePSI");
  TTree *TreeZprimeSSM = new TTree("TreeZprimeSSM","TreeZprimeSSM");
  TTree *TreeADD = new TTree("TreeADD","TreeADD");
  TTree *TreeAFB1 = new TTree("TreeAFB1","TreeAFB1");
  TTree *TreeAFBmin1 = new TTree("TreeAFBmin1","TreeAFBmin1");

 
  BuildTree(TreeDY); 
  BuildTree(TreeSpin0gg); 
  BuildTree(TreeSpin0qqbar); 
  BuildTree(TreeSpin0qqbarggmix); 
  BuildTree(TreeSpin2gg); 
  BuildTree(TreeSpin2qqbar); 
  BuildTree(TreeSpin2qqbarggmix); 
  BuildTree(TreeRSGrav); 
  BuildTree(TreeZprimePSI); 
  BuildTree(TreeZprimeSSM); 
  BuildTree(TreeADD);
  BuildTree(TreeAFB1); 
  BuildTree(TreeAFBmin1);
  
  //These are some histos used for acceptance studies vs M_gen
  TH1F *hdens0= new TH1F("hdens0", "M_{gen}", 60,0,6000);
  TH1F *hnums0= new TH1F("hnums0", "M_{gen}", 60,0,6000);
  TH1F *hdendy= new TH1F("hdendy", "M_{gen}", 60,0,6000);
  TH1F *hnumdy= new TH1F("hnumdy", "M_{gen}", 60,0,6000);
  TH1F *hdenrsgrav= new TH1F("hdenrsgrav", "M_{gen}", 60,0,6000);
  TH1F *hnumrsgrav= new TH1F("hnumrsgrav", "M_{gen}", 60,0,6000);
  TH1F *hdenzppsi= new TH1F("hdenzppsi", "M_{gen}", 60,0,6000);
  TH1F *hnumzppsi= new TH1F("hnumzppsi", "M_{gen}", 60,0,6000);
  TH1F *hdenzpssm= new TH1F("hdenzpssm", "M_{gen}", 60,0,6000);
  TH1F *hnumzpssm= new TH1F("hnumzpssm", "M_{gen}", 60,0,6000);
  hdens0->Sumw2();
  hnums0->Sumw2();
  hdendy->Sumw2();
  hnumdy->Sumw2();
  hdenrsgrav->Sumw2();
  hnumrsgrav->Sumw2();
  hdenzppsi->Sumw2();
  hnumzppsi->Sumw2();
  hdenzpssm->Sumw2();
  hnumzpssm->Sumw2();

  //Some time ago, I also used this macro for other stuff and also ran on ADD samples such as this one: 
// Initialize("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/lathomas/SignalMCnew/ADDdiLepton_LambdaT-1600_Tune4C_8TeV-pythia8_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/TOTAL.root","ADDLambdaT1600",1);

// I also initially used this macro with the following Phase 0 samples: 
/* Initialize("/user/lathomas/mcsamples/13TeV/DYToEE3000_Pythia8_CMSSW_6_2_0_patch1_PU_S10.root","DYM3000",1);
Initialize("/user/lathomas/mcsamples/13TeV/DYToEE4000_Pythia8_CMSSW_6_2_0_patch1_PU_S10.root","DYM4000",1);
Initialize("/user/lathomas/mcsamples/13TeV/ZprimeSSMToEE3000_Pythia8_CMSSW_6_2_0_patch1_PU_S10.root","ZprimeSSM3000",1);
Initialize("/user/lathomas/mcsamples/13TeV/ZprimeSSMToEE4000_Pythia8_CMSSW_6_2_0_patch1_PU_S10.root","ZprimeSSM4000",1);
Initialize("/user/lathomas/mcsamples/13TeV/ZprimeSSMToEE5000_Pythia8_CMSSW_6_2_0_patch1_PU_S10.root","ZprimeSSM5000",1);
Initialize("/user/lathomas/mcsamples/13TeV/ZprimePSIToEE3000_Pythia8_CMSSW_6_2_0_patch1_PU_S10.root","ZprimePSI3000",1);
Initialize("/user/lathomas/mcsamples/13TeV/ZprimePSIToEE4000_Pythia8_CMSSW_6_2_0_patch1_PU_S10.root","ZprimePSI4000",1);
Initialize("/user/lathomas/mcsamples/13TeV/ZprimePSIToEE5000_Pythia8_CMSSW_6_2_0_patch1_PU_S10.root","ZprimePSI5000",1);
Initialize("/user/lathomas/mcsamples/13TeV/RSGravToEE3000_Pythia8_CMSSW_6_2_0_patch1_PU_S10.root","RSGrav3000",1);
Initialize("/user/lathomas/mcsamples/13TeV/RSGravToEE4000_Pythia8_CMSSW_6_2_0_patch1_PU_S10.root","RSGrav4000",1);
Initialize("/user/lathomas/mcsamples/13TeV/RSGravToEE5000_Pythia8_CMSSW_6_2_0_patch1_PU_S10.root","RSGrav5000",1);*/


//The three samples to use now are those ones (pick one out of the three): 
//Initialize("/user/lathomas/mcsamples/ZprimeSSMToEE_M-2500_TuneZ2star_14TeV-pythia6_GEM2019Upg14DR-final_phase1_PU50bx25_DES19_62_V8-v1_AODSIM.root","ZprimeSSM2500",1);
//   Initialize("/user/lathomas/mcsamples/ZprimeSSMToEE_M-2500_TuneZ2star_14TeV-pythia6_GEM2019Upg14DR-age1k_PU140bx25_PH1_1K_FB_V1-v1_AODSIM_debug.root","ZprimeSSM2500",1);
Initialize("/user/lathomas/mcsamples/ZprimeSSMToEE_M-2500_TuneZ2star_14TeV-pythia6_Muon2023Upg14DR-PU140bx25_PH2_1K_FB_V2-v1_AODSIM.root","ZprimeSSM2500",1);

 cout<< "everything initialized"<< endl;
  const int nbFile = input.size();


  //First loop to get the cos theta at gen level, this is needed for the reweighting applied in some samples. 
  for (int p = 0; p<nbFile; p++) {
     cout<<"file "<<(input[p])->GetName()<<endl;
    (input[p])->cd();
    TTree *thetree = (TTree*)(input[p])->Get("gsfcheckerjob/tree");
    Init(thetree);
    Long64_t nentries = (*thetree).GetEntries();
    cout<<nentries<<" entries"<<endl;
    //Loop over the entries of the input tree
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      thetree->GetEntry(jentry);
      //  if( jentry>= 200 )break; //Uncomment here for quick tests
      if( (jentry%10000) == 0 ) cout<<"entry "<<jentry<<endl;
   
      CalcCosthetaStar(costhetastar,costhetastar_meas,p);// This function calculates costhetastar and costhetastar_meas at the gen level for each event. 
    
      //Filling some histos
      hcosthetastar[p]->Fill(costhetastar);
      hfabscosthetastar[p]->Fill(fabs(costhetastar));
      hcosthetastar_meas[p]->Fill(costhetastar_meas);
      hfabscosthetastar_meas[p]->Fill(fabs(costhetastar_meas)); 

      /*Split the case with 0, 1, 2 gluons in the initial state.
      For some reason, gluons are counted twice in the input tree at the moment.
      In principle, if you use Pythia, you should not get any event with 1 gluon in the initial state. 
       */
      if(gengluons_size==0){hcosthetastar0gluon[p]->Fill(costhetastar);hfabscosthetastar0gluon[p]->Fill(fabs(costhetastar));}
      if(gengluons_size==2){hcosthetastar1gluon[p]->Fill(costhetastar);hfabscosthetastar1gluon[p]->Fill(fabs(costhetastar));}
      if(gengluons_size==4){hcosthetastar2gluons[p]->Fill(costhetastar);hfabscosthetastar2gluons[p]->Fill(fabs(costhetastar));}

    }
  }

  //Drawing those histos. Can be commented.
  // DrawPlot(hcosthetastar);
  // DrawPlot(hfabscosthetastar);
  // DrawPlot(hcosthetastar_meas);
  // DrawPlot(hfabscosthetastar_meas);
  // DrawPlot(hcosthetastar0gluon);
  // DrawPlot(hfabscosthetastar0gluon); 
  // DrawPlot(hcosthetastar1gluon);
  // DrawPlot(hfabscosthetastar1gluon); 
  // DrawPlot(hcosthetastar2gluons);
  // DrawPlot(hfabscosthetastar2gluons);







 
  //Second loop 
  for (int p = 0; p<nbFile; p++) {
   
    cout<<"file "<<(input[p])->GetName()<<endl;
    (input[p])->cd();
    TTree *thetree = (TTree*)(input[p])->Get("gsfcheckerjob/tree");
    Init(thetree);
    Long64_t nentries = (*thetree).GetEntries();
    cout<<nentries<<" entries"<<endl;
    //Loop over the entries of the input tree
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      thetree->GetEntry(jentry);
      //if( jentry>= 200 )break; 
      if( (jentry%10000) == 0 ) cout<<"entry "<<jentry<<endl;

      //Calculate costheta at the gen level
      CalcCosthetaStar(costhetastar,costhetastar_meas,p);
  
      bool passevt = false; 
      genmass = genelemom_mass[0];

      // Old stuff: with the ADD sample, the variable genelemom_mass[0] is not filled.
      if(samplename[p].Index("ADD") !=kNPOS ) {
	TLorentzVector ele1;
	TLorentzVector ele2;
	ele1.SetPtEtaPhiE(unstableGenEle_pt[0],unstableGenEle_eta[0],unstableGenEle_phi[0], unstableGenEle_e[0]);
	ele2.SetPtEtaPhiE(unstableGenEle_pt[1],unstableGenEle_eta[1],unstableGenEle_phi[1], unstableGenEle_e[1]);
	genmass= (ele1+ele2).Mag();
	//	cout << "genmass for add " << genmass << endl;
      }

      //Now calculating the various variables to be stored in the tree. 
      nbofgengluons = gengluons_size/2; 
      weightevt = 1;
      if(samplename[p].Index("ADD") !=kNPOS)  nbofgengluons = gengluons_size;
      
      if(samplename[p].Index("ADD") !=kNPOS && gengluons_size+genquarks_size==0 )weightevt = 0;  


      costhetaCSheepheep =-1000;
      rapidityheepheep=-1000;
      recomass = 0;
      genelespt35 = false;
      genelesBB = false;
      genelesBE = false;
      genelesEE = false;
      genelesGap = false;
      genelesBeamPipe = false;
      gsfptcut = false;
      gsfBB = false;
      gsfBE = false; 
      gsfEE = false;
      hltfired = false;
      heepheep = false; 
      mrecomgenmatching = false;
      
      //Check where the gen electrons pass the Pt threshold and where they are emitted (Barrel-Barrel, Barrel-Endcaps, beampipe,...).
      if(fabs(unstableGenEle_pt[0]) >ptcut &&fabs(unstableGenEle_pt[1]) >ptcut ) genelespt35 =  true;
      if(fabs(unstableGenEle_eta[0]) <etabarrel &&fabs(unstableGenEle_eta[1]) <etabarrel ) genelesBB =  true;
      else if(fabs(unstableGenEle_eta[0])>etaendcapslow && fabs(unstableGenEle_eta[0]) <etaendcapshigh  &&fabs(unstableGenEle_eta[1]) <etabarrel  ) genelesBE = true; 
      else if(fabs(unstableGenEle_eta[1])>etaendcapslow && fabs(unstableGenEle_eta[1]) <etaendcapshigh  &&fabs(unstableGenEle_eta[0]) <etabarrel  ) genelesBE = true; 
      else if(fabs(unstableGenEle_eta[0])>etaendcapslow && fabs(unstableGenEle_eta[0]) <etaendcapshigh  && fabs(unstableGenEle_eta[1])>etaendcapslow && fabs(unstableGenEle_eta[1]) <etaendcapshigh ) genelesEE = true; 
      else if(fabs(unstableGenEle_eta[0]) <etaendcapshigh &&fabs(unstableGenEle_eta[1]) <etaendcapshigh ) genelesGap = true;
      else genelesBeamPipe = true; 

      if(HLT_DoubleEle33_CaloIdL_GsfTrkIdVL) hltfired = true;//N.B. No trigger menu in upgraded samples...

      /*Now, loop over the GSF
       Only select GSF matching a generate electron. 
       
       CAVEAT: The dR cut is not enough here. 
       There are so many jets from PU in some upgraded samples that it is not rare to have a GSF from a jet satisfying dR<0.3 with a gen ele by accident.
I therefore added the following cut 0.7 <E_T(reco)/E_T(gen) <1.3 
       */
      
      //First GSF 
      for(int i =0; i<gsf_size; i++){
	if(deltaR(gsf_eta[i],gsf_phi[i],unstableGenEle_eta[0],unstableGenEle_phi[0]) > 0.3 && deltaR(gsf_eta[i],gsf_phi[i],unstableGenEle_eta[1],unstableGenEle_phi[1]) >0.3) continue; 
	if(deltaR(gsf_eta[i],gsf_phi[i],unstableGenEle_eta[0],unstableGenEle_phi[0]) < 0.3 &&  (gsf_gsfet[i]/unstableGenEle_pt[0] <0.7 ||gsf_gsfet[i]/unstableGenEle_pt[0]   >1.3 ) ) continue; 
	if(deltaR(gsf_eta[i],gsf_phi[i],unstableGenEle_eta[1],unstableGenEle_phi[1]) < 0.3 &&  (gsf_gsfet[i]/unstableGenEle_pt[1] <0.7 ||gsf_gsfet[i]/unstableGenEle_pt[1]   >1.3 ) ) continue; 

	//Second GSF
	for(int j =  0; j< i; j++){
	  if(deltaR(gsf_eta[j],gsf_phi[j],unstableGenEle_eta[0],unstableGenEle_phi[0]) > 0.3 && deltaR(gsf_eta[j],gsf_phi[j],unstableGenEle_eta[1],unstableGenEle_phi[1]) >0.3) continue; 
	if(deltaR(gsf_eta[j],gsf_phi[j],unstableGenEle_eta[0],unstableGenEle_phi[0]) < 0.3 &&  (gsf_gsfet[j]/unstableGenEle_pt[0] <0.7 ||gsf_gsfet[j]/unstableGenEle_pt[0]   >1.3 ) ) continue; 
	if(deltaR(gsf_eta[j],gsf_phi[j],unstableGenEle_eta[1],unstableGenEle_phi[1]) < 0.3 &&  (gsf_gsfet[j]/unstableGenEle_pt[1] <0.7 ||gsf_gsfet[j]/unstableGenEle_pt[1]   >1.3 ) ) continue; 	

      
	if(gsf_gsfet[i]>ptcut && gsf_gsfet[j]>ptcut)gsfptcut = true; 
	if(fabs(gsfsc_eta[i])<etabarrel &&fabs(gsfsc_eta[j]) <etabarrel)gsfBB = true; 
	else if(fabs(gsfsc_eta[i])<etabarrel && fabs(gsfsc_eta[j])>etaendcapslow&&fabs(gsfsc_eta[j])<etaendcapshigh)gsfBE = true;
	else if(fabs(gsfsc_eta[j])<etabarrel && fabs(gsfsc_eta[i])>etaendcapslow&&fabs(gsfsc_eta[i])<etaendcapshigh)gsfBE = true;
	else if( fabs(gsfsc_eta[j])>etaendcapslow&&fabs(gsfsc_eta[j])<etaendcapshigh && fabs(gsfsc_eta[i])>etaendcapslow&&fabs(gsfsc_eta[i])<etaendcapshigh)gsfEE = true;
   

	if(PassHEEP(i) && PassHEEP(j))heepheep =true;
	recomass = CalcInvariantMass(i,j);
	if(recomass >0 && recomass/genmass<1.2 && recomass/genmass>0.8) mrecomgenmatching=true; 
	// Costheta can only be calculated if the two GSF have an opposite charge. 
	if(gsf_charge[i]*gsf_charge[j]==-1 &&PassHEEP(i) && PassHEEP(j)) costhetaCSheepheep = CalcCosthetaCSData(i,j); 
	if(gsf_charge[i]*gsf_charge[j]==-1 &&PassHEEP(i) && PassHEEP(j)) rapidityheepheep = CalcRapidity(i,j);
	
	}
      }
      
   
   
      for(int itgsf = 0; itgsf < gsf_size; itgsf++){
	//	if(!HLT_DoubleEle33_CaloIdL_GsfTrkIdVL) continue;  
	if(!gsfpass_HEEP[itgsf]) continue; 
    	for(int itgsf2 = 0; itgsf2 < gsf_size; itgsf2++){
	  if(itgsf2 == itgsf) continue;   
	  if(!gsfpass_HEEP[itgsf2]) continue; 
	  double heepheepmass = CalcInvariantMass(itgsf,itgsf2);
	  if (heepheepmass/genelemom_mass[0]<0.8 || heepheepmass/genelemom_mass[0]>1.2) continue; 
	  passevt =true; 
	  break;
	}
      }
      
      /*Reweighting based on the costhetastar distribution:
	This is where we use the histo   hcosthetastar filled earlier.
	For a spin 0, the weight is simply obtained such that the costheta distribution becomes flat. 
	For the maximally positive AFB, the distribution must be reweighted to follow (1+costheta)^2
	For the maximally negative AFB, the distribution must be reweighted to follow (1-cos^2theta)^2

	For the sample with a mixture 50% gg and 50% qqbar, events produced though each mode are reweighted to reach the desired fraction.
      */
      int mybin =hcosthetastar[p]->FindBin((costhetastar));
      
      double weights0=hcosthetastar[p]->Integral()/hcosthetastar[p]->GetBinContent(mybin)/200;
      double weights0gg= hcosthetastar2gluons[p]->Integral()/hcosthetastar2gluons[p]->GetBinContent(mybin)/200;
      
      double weightafb1 = (1+costhetastar)*(1+costhetastar)*hcosthetastar[p]->Integral()/hcosthetastar[p]->GetBinContent(mybin)/200;
      double weightafbmin1 = (1-costhetastar)*(1-costhetastar)*hcosthetastar[p]->Integral()/hcosthetastar[p]->GetBinContent(mybin)/200;

      double weights0qqbarggmix = 0; 
      double weights2qqbarggmix = 0; 
      if(nbofgengluons ==0){
	weights0qqbarggmix =hcosthetastar2gluons[p]->Integral()/hcosthetastar0gluon[p]->GetBinContent(mybin)/200;
	weights2qqbarggmix =hcosthetastar2gluons[p]->Integral()/hcosthetastar0gluon[p]->Integral();
      }
      if(nbofgengluons ==2){
	weights0qqbarggmix =hcosthetastar2gluons[p]->Integral()/hcosthetastar2gluons[p]->GetBinContent(mybin)/200;
	weights2qqbarggmix =1;
      }
      
      
      //Now filling the trees      
      if(samplename[p].Index("DY") !=kNPOS){
	hdendy->Fill(genelemom_mass[0]);
	if(passevt) hnumdy->Fill(genelemom_mass[0],1);
	TreeDY->Fill();
      }
      
      if(samplename[p].Index("RSGrav") !=kNPOS){
	hdenrsgrav->Fill(genelemom_mass[0]);
	if(passevt) hnumrsgrav->Fill(genelemom_mass[0],1); 
	TreeRSGrav->Fill();
	if(nbofgengluons ==0)TreeSpin2qqbar->Fill();
	if(nbofgengluons ==2)TreeSpin2gg->Fill();
	
	weightevt = weights0qqbarggmix ; 
	TreeSpin0qqbarggmix->Fill();
	
	weightevt = weights2qqbarggmix; 
	TreeSpin2qqbarggmix->Fill();
		
	if(nbofgengluons ==2)weightevt =weights0gg;
	if(nbofgengluons ==2)TreeSpin0gg->Fill();
      }
      if(samplename[p].Index("ZprimePSI") !=kNPOS){
	
	hdenzppsi->Fill(genelemom_mass[0]);
	if(passevt) hnumzppsi->Fill(genelemom_mass[0],1);
	TreeZprimePSI->Fill();

	 hdens0->Fill(genelemom_mass[0],weights0); 
	if(passevt)hnums0->Fill(genelemom_mass[0],weights0);
	weightevt = 1/(1+costhetastar*costhetastar);// weights0;
	TreeSpin0qqbar->Fill();


	
      }
      if(samplename[p].Index("ZprimeSSM") !=kNPOS){
	hdenzpssm->Fill(genelemom_mass[0]);
	if(passevt) hnumzpssm->Fill(genelemom_mass[0],1); 
	TreeZprimeSSM->Fill();

	weightevt =weightafb1;//test
	TreeAFB1->Fill();

	weightevt =weightafbmin1;//test
	TreeAFBmin1->Fill();
	
      }    
      if(samplename[p].Index("ADD") !=kNPOS){

	TreeADD->Fill();
      }  
      
    
     
    }
  }


  
  // TreeDY->Print();
  // TreeRSGrav->Print();
  // TreeZprimePSI->Print();
  // TreeZprimeSSM->Print();
  // TreeSpin0qqbar->Print();
  // TreeSpin0gg->Print();
  // TreeSpin0qqbarggmix->Print();
  // TreeSpin2qqbar->Print();
  // TreeSpin2gg->Print(); 
  // TreeSpin2qqbarggmix->Print();
  f1->cd(); 
  TreeDY->Write(); 
  TreeRSGrav->Write(); 
  TreeZprimePSI->Write(); 
  TreeZprimeSSM->Write(); 
  TreeSpin0qqbar->Write(); 
  TreeSpin0gg->Write(); 
  TreeSpin2qqbar->Write(); 
  TreeSpin2gg->Write(); 
  TreeSpin0qqbarggmix->Write(); 
  TreeSpin2qqbarggmix->Write(); 
  TreeAFB1->Write(); 
  TreeAFBmin1->Write(); 
  

}



void AccforSpin0::Initialize(TString samplefile, TString mysamplename,  float sampleweight)
{
    input.push_back(new TDCacheFile(samplefile,"open"));
    sweight.push_back(sampleweight);
    samplename.push_back(mysamplename);
    hcosthetastar.push_back(new TH1F(mysamplename+"hcosthetastar", "cos#theta^{*}",200,-1,1));
    hcosthetastar.back()->Sumw2();
    hfabscosthetastar.push_back(new TH1F(mysamplename+"hfabscosthetastar", "|cos#theta^{*}|",200,0,1));
    hfabscosthetastar.back()->Sumw2();
    hcosthetastar_meas.push_back(new TH1F(mysamplename+"hcosthetastar_meas", "cos#theta^{*} (meas.)",200,-1,1));
    hcosthetastar_meas.back()->Sumw2();
    hfabscosthetastar_meas.push_back(new TH1F(mysamplename+"hfabscosthetastar_meas", "|cos#theta^{*} (meas.)|",200,0,1));
    hfabscosthetastar_meas.back()->Sumw2();


    hcosthetastar0gluon.push_back(new TH1F(mysamplename+"hcosthetastar0gluon","cos#theta^{*} (meas.)",200,-1,1));
    hcosthetastar0gluon.back()->Sumw2();
    hfabscosthetastar0gluon.push_back(new TH1F(mysamplename+"hfabscosthetastar0gluon","cos#theta^{*} (meas.)",200,0,1));
    hfabscosthetastar0gluon.back()->Sumw2();
    hcosthetastar1gluon.push_back(new TH1F(mysamplename+"hcosthetastar1gluon","cos#theta^{*} (meas.)",200,-1,1));
    hcosthetastar1gluon.back()->Sumw2();
    hfabscosthetastar1gluon.push_back(new TH1F(mysamplename+"hfabscosthetastar1gluon","cos#theta^{*} (meas.)",200,0,1));
    hfabscosthetastar1gluon.back()->Sumw2();
    hcosthetastar2gluons.push_back(new TH1F(mysamplename+"hcosthetastar2gluons","cos#theta^{*} (meas.)",200,-1,1));
    hcosthetastar2gluons.back()->Sumw2();
    hfabscosthetastar2gluons.push_back(new TH1F(mysamplename+"hfabscosthetastar2gluons","cos#theta^{*} (meas.)",200,0,1));
    hfabscosthetastar2gluons.back()->Sumw2();



}

 void AccforSpin0::DrawPlot(vector <TH1F*> myhisto ){
 f1->cd();  
   for(int i =0;i<myhisto.size();i++){
     myhisto[i]->Write();
   }


}


void AccforSpin0::ComputeAcc(TH1F* hnumS0 , TH1F* hdenS0,TH1F* hnumS1 , TH1F* hdenS1,TH1F* hnumS2 , TH1F* hdenS2){
  
  //Remove bins with small stats (error must be smaller than 10% : 
  for(int i = 1; i<=hdenS0->GetNbinsX();i++)
    {
      if(hdenS0->GetBinError(i)/hdenS0->GetBinContent(i)>0.1){ hdenS0->SetBinError(i,0);hdenS0->SetBinContent(i,0);hnumS0->SetBinError(i,0);hnumS0->SetBinContent(i,0);}
    }

  for(int i = 1; i<=hdenS1->GetNbinsX();i++)
    {
      if(hdenS1->GetBinError(i)/hdenS1->GetBinContent(i)>0.1){ hdenS1->SetBinError(i,0);hdenS1->SetBinContent(i,0);hnumS1->SetBinError(i,0);hnumS1->SetBinContent(i,0);}
    }
  for(int i = 1; i<=hdenS2->GetNbinsX();i++)
    {
      if(hdenS2->GetBinError(i)/hdenS2->GetBinContent(i)>0.1){ hdenS2->SetBinError(i,0);hdenS2->SetBinContent(i,0);hnumS2->SetBinError(i,0);hnumS2->SetBinContent(i,0);}
    }


  TGraphAsymmErrors * mygphS0 = new TGraphAsymmErrors;
  TGraphAsymmErrors * mygphS1 = new TGraphAsymmErrors;
  TGraphAsymmErrors * mygphS2 = new TGraphAsymmErrors;
  mygphS0->BayesDivide(hnumS0,hdenS0);
  mygphS1->BayesDivide(hnumS1,hdenS1);
  mygphS2->BayesDivide(hnumS2,hdenS2);

  TString canvasname = hnumS0->GetName();
  TCanvas *c1 = new TCanvas(canvasname, "", 100, 100, 900, 650);
 
  c1->cd(); 
  gPad->SetGridx();
  gPad->SetGridy();
  mygphS0->SetLineColor(kBlack); 
  mygphS0->Draw("AP"); 
  mygphS1->SetLineColor(kRed); 
  mygphS2->SetLineColor(kBlue); 
  mygphS1->Draw("SAME,P"); 
   mygphS2->Draw("SAME,P"); 
  labelcmssimul->Draw("same");
  f1->cd();
  c1->Write();
}




void AccforSpin0::ComputeAcc(TH1F* hnum , TH1F* hden, TString gphname){
  
  //Remove bins with small stats (error must be smaller than 10% : 
  for(int i = 1; i<=hden->GetNbinsX();i++)
    {
      if(hden->GetBinError(i)/hden->GetBinContent(i)>0.1){ hden->SetBinError(i,0);hden->SetBinContent(i,0);hnum->SetBinError(i,0);hnum->SetBinContent(i,0);}
    }

  TGraphAsymmErrors * mygph = new TGraphAsymmErrors;

  mygph->BayesDivide(hnum,hden);
 
  mygph->SetName(gphname);

  f1->cd();
  mygph->Write(); 
}



double
AccforSpin0::CalcInvariantMass(const int& n, const int& l)
{
  bool etShift = false;
  float dataEtShiftFactorEB = 1.;
  float dataEtShiftFactorEE = 1.;

  if (etShift) {
    dataEtShiftFactorEB = 1.0036;
    dataEtShiftFactorEE = 1.0256;
  }

  TLorentzVector ele1;
  TLorentzVector ele2;

  Float_t et1 = gsf_gsfet[n];
  Float_t et2 = gsf_gsfet[l];
  // correct energy in data
  (fabs(gsfsc_eta[n]) > 1.566) ? et1 *= dataEtShiftFactorEE : et1 *= dataEtShiftFactorEB;
  (fabs(gsfsc_eta[l]) > 1.566) ? et2 *= dataEtShiftFactorEE : et2 *= dataEtShiftFactorEB;

  ele1.SetPtEtaPhiE(et1, gsf_eta[n], gsf_phi[n], (et1 * cosh(gsf_eta[n])));
  ele2.SetPtEtaPhiE(et2, gsf_eta[l], gsf_phi[l], (et2 * cosh(gsf_eta[l])));

  return (ele1+ele2).Mag();
}





void AccforSpin0::CalcCosthetaStar(float & costh, float & costhmeas, const int m){
  float quarkdir = (genquark_pdgid[0]>0) ? genquark_pz[0]/fabs(genquark_pz[0]) : - genquark_pz[0]/fabs(genquark_pz[0]) ;
  //  if(quarkdir ==0) quarkdir = 1;
      float Zpz = genelemom_pz[0];
      float Zpt = genelemom_pt[0];
      float mass = genelemom_mass[0];
      
      if(samplename[m].Index("ADD") !=kNPOS ) {

	TLorentzVector ele1;
	TLorentzVector ele2;
	ele1.SetPtEtaPhiE(unstableGenEle_pt[0],unstableGenEle_eta[0],unstableGenEle_phi[0], unstableGenEle_e[0]);
	ele2.SetPtEtaPhiE(unstableGenEle_pt[1],unstableGenEle_eta[1],unstableGenEle_phi[1], unstableGenEle_e[1]);
	mass= (ele1+ele2).Mag();
	Zpz=(ele1+ele2).Pz();
	Zpt=(ele1+ele2).Pt();
	
      }


      
      float elepz = (unstableGenEle_charge[0]<0)? unstableGenEle_pz[0] : unstableGenEle_pz[1];
      float posipz = (unstableGenEle_charge[0]<0)? unstableGenEle_pz[1] : unstableGenEle_pz[0];
      float eleenergy = (unstableGenEle_charge[0]<0)? unstableGenEle_e[0] : unstableGenEle_e[1];
      float posienergy = (unstableGenEle_charge[0]<0)? unstableGenEle_e[1] : unstableGenEle_e[0];
   
     
      costh = quarkdir*2/mass/sqrt(mass*mass+Zpt*Zpt)*(elepz*posienergy-posipz*eleenergy) ;
      costhmeas = Zpz/fabs(Zpz)*2/mass/sqrt(mass*mass+Zpt*Zpt)*(elepz*posienergy-posipz*eleenergy) ;

      

}

// Not to be used !
void AccforSpin0::CalcCostheta(float & costh, float & costhmeas, const int m){
  float quarkdir = (genquark_pdgid[0]>0) ? genquark_pz[0]/fabs(genquark_pz[0]) : - genquark_pz[0]/fabs(genquark_pz[0]) ;
  //  if(quarkdir ==0) quarkdir = 1;
      float Zpz = genelemom_pz[0];
      float Zpt = genelemom_pt[0];
      float mass = genelemom_mass[0];
      
      if(samplename[m].Index("ADD") !=kNPOS ) {

	TLorentzVector ele1;
	TLorentzVector ele2;
	ele1.SetPtEtaPhiE(unstableGenEle_pt[0],unstableGenEle_eta[0],unstableGenEle_phi[0], unstableGenEle_e[0]);
	ele2.SetPtEtaPhiE(unstableGenEle_pt[1],unstableGenEle_eta[1],unstableGenEle_phi[1], unstableGenEle_e[1]);
	mass= (ele1+ele2).Mag();
	Zpz=(ele1+ele2).Pz();
	Zpz=(ele1+ele2).Pt();
	
      }



      
      float elepz = (unstableGenEle_charge[0]<0)? unstableGenEle_pz[0] : unstableGenEle_pz[1];
      float posipz = (unstableGenEle_charge[0]<0)? unstableGenEle_pz[1] : unstableGenEle_pz[0];
      float eleenergy = (unstableGenEle_charge[0]<0)? unstableGenEle_e[0] : unstableGenEle_e[1];
      float posienergy = (unstableGenEle_charge[0]<0)? unstableGenEle_e[1] : unstableGenEle_e[0];
   
     
      costh = quarkdir*2/mass/sqrt(mass*mass+Zpt*Zpt)*(elepz*posienergy-posipz*eleenergy) ;
      costhmeas = Zpz/fabs(Zpz)*2/mass/sqrt(mass*mass+Zpt*Zpt)*(elepz*posienergy-posipz*eleenergy) ;

      

}


float  AccforSpin0::CalcCosthetaCSData(const int &m, const int & n){
  
  float cosinusthetaCollinsSoper = 0; 
  if(gsf_charge[m]*gsf_charge[n]!=-1 ){cout << "Same sign electrons " <<endl ; return -1000; }
   
    
  TLorentzVector ele1;
  TLorentzVector ele2;

  Float_t et1 = gsf_gsfet[m];
  Float_t et2 = gsf_gsfet[n];
 
  ele1.SetPtEtaPhiE(et1, gsf_eta[m], gsf_phi[m], (et1 * cosh(gsf_eta[m])));
  ele2.SetPtEtaPhiE(et2, gsf_eta[n], gsf_phi[n], (et2 * cosh(gsf_eta[n])));
  float elepz = (gsf_charge[m]<0)?   ele1.Pz(): ele2.Pz();
  float posipz = (gsf_charge[m]>0)?  ele1.Pz(): ele2.Pz();
  float eleenergy = (gsf_charge[m]<0)? ele1.E(): ele2.E(); 
  float posienergy = (gsf_charge[m]>0)? ele1.E(): ele2.E(); 
  
  Float_t pzpair = (ele1+ele2).Pz();
  Float_t mass = (ele1+ele2).Mag();
  Float_t Zpt = (ele1+ele2).Pt();
  cosinusthetaCollinsSoper = pzpair/fabs(pzpair)*2/mass/sqrt(mass*mass+Zpt*Zpt)*(elepz*posienergy-posipz*eleenergy) ;
  return cosinusthetaCollinsSoper;   

}




float  AccforSpin0::CalcRapidity(const int &m, const int & n){
  TLorentzVector ele1;
  TLorentzVector ele2;
  Float_t et1 = gsf_gsfet[m];
  Float_t et2 = gsf_gsfet[n];
 
  ele1.SetPtEtaPhiE(et1, gsf_eta[m], gsf_phi[m], (et1 * cosh(gsf_eta[m])));
  ele2.SetPtEtaPhiE(et2, gsf_eta[n], gsf_phi[n], (et2 * cosh(gsf_eta[n])));

  return (ele1+ele2).Rapidity();   

}






void AccforSpin0::BuildTree(TTree * T1) {
  T1->Branch("costhetastar", &costhetastar,"costhetastar/F");
  T1->Branch("costhetastar_meas", &costhetastar_meas,"costhetastar_meas/F"); 
  T1->Branch("genmass", &genmass, "genmass/F");
  T1->Branch("recomass", &recomass, "recomass/F");
  T1->Branch("nbofgengluons",&nbofgengluons,"nbofgengluons/I"); 
  T1->Branch("weightevt", &weightevt, "weightevt/F");
  T1->Branch("genelespt35", &genelespt35,"genelespt35/O");
  T1->Branch("genelesBB", &genelesBB,"genelesBB/O");
  T1->Branch("genelesBE", &genelesBE,"genelesBE/O");
  T1->Branch("genelesEE", &genelesEE,"genelesEE/O");
  T1->Branch("genelesGap", &genelesGap,"genelesGap/O");
  T1->Branch("genelesBeamPipe", &genelesBeamPipe,"genelesBeamPipe/O");
  T1->Branch("gsfptcut", &gsfptcut,"gsfptcut/O");
  T1->Branch("gsfBB", &gsfBB,"gsfBB/O");
  T1->Branch("gsfBE" , &gsfBE,"gsfBE/O");
  T1->Branch("gsfEE", &gsfEE,"gsfEE/O");
  T1->Branch("hltfired", &hltfired,"hltfired/O");
  T1->Branch("heepheep" , &heepheep,"heepheep/O");
  T1->Branch("mrecomgenmatching", &mrecomgenmatching,"mrecomgenmatching/O");
  T1->Branch("costhetaCSheepheep", &costhetaCSheepheep,"costhetaCSheepheep/F");
  T1->Branch("rapidityheepheep", &rapidityheepheep,"rapidityheepheep/F");
}



bool
AccforSpin0::PassHEEP(const int &n)
{

  // HEEP v4.1
  // barrel
  float bar_et = 35.;
  float bar_hoE = 0.05;
  float bar_DEta = 0.005;
  float bar_DPhi = 0.06;
  float bar_e2x5e5x5 = 0.94;
  float bar_e1x5e5x5 = 0.83;
  float bar_isoEcalHcal1_1 = 2.;
  float bar_isoEcalHcal1_2 = 0.03;
  float bar_isoEcalHcalRho = 0.28;
  float bar_isoTrack = 5.;
  int bar_missInnerHits = 1;
  float bar_gsfdxy = 0.02;

  // endcap
  float end_et = 35.;
  float end_hoE = 0.05;
  float end_DEta = 0.007 ;
  float end_DPhi = 0.06;
  float end_sigmaietaieta = 0.03;
  float end_isoEcalHcal1_1_1 = 2.5;
  float end_isoEcalHcal1_1_2 = 1.;
  float end_isoEcalHcal1_2 = 0.03;
  float end_isoEcalHcalRho = 0.28;
  float end_isoTrack = 5.;
  int end_missInnerHits = 1;
  float end_gsfdxy = 0.05;



  // HEEP v4.0
  // barrel
  //  if(gsf_eOVERp[n] >10) return false; 
  if (fabs(gsfsc_eta[n]) < 1.5
      && gsf_gsfet[n] > bar_et
      && gsf_isecaldriven[n]
      && gsf_hovere[n] < bar_hoE
      && fabs(gsf_deltaeta[n]) < bar_DEta
      && fabs(gsf_deltaphi[n]) < bar_DPhi
      && (gsf_e2x5overe5x5[n] > bar_e2x5e5x5 || gsf_e1x5overe5x5[n] > bar_e1x5e5x5)
      && (gsf_ecaliso[n] + gsf_hcaliso1[n]) < (bar_isoEcalHcal1_1 + bar_isoEcalHcal1_2 * gsf_gsfet[n] + bar_isoEcalHcalRho * rho)
      && gsf_trackiso[n] < bar_isoTrack
      && gsf_nLostInnerHits[n] <= bar_missInnerHits
      && fabs(gsf_dxy_firstPVtx[n])<= bar_gsfdxy
     
     ) return true;
  
  // endcap
  if ((fabs(gsfsc_eta[n]) > 1.5 && fabs(gsfsc_eta[n]) < 10000)
      && gsf_gsfet[n] > end_et
      && gsf_isecaldriven[n]
      && gsf_hovere[n] < end_hoE
      && fabs(gsf_deltaeta[n]) < end_DEta
      && fabs(gsf_deltaphi[n]) < end_DPhi
      && gsf_sigmaIetaIeta[n] < end_sigmaietaieta
      && ((gsf_gsfet[n] < 50. && (gsf_ecaliso[n] + gsf_hcaliso1[n]) < (end_isoEcalHcal1_1_1 + end_isoEcalHcalRho * rho))
         || (gsf_gsfet[n] >= 50. && (gsf_ecaliso[n] + gsf_hcaliso1[n]) < (end_isoEcalHcal1_1_2 + end_isoEcalHcal1_2 * gsf_gsfet[n] + end_isoEcalHcalRho * rho)))
      && gsf_trackiso[n] < end_isoTrack
      && gsf_nLostInnerHits[n] <= end_missInnerHits
      && fabs(gsf_dxy_firstPVtx[n])<= end_gsfdxy
     ) return true;
  return false;
}

