#include <TFile.h>
#include <TH2.h>
#include <TTree.h>
#include <TGraph.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TROOT.h>
#include "TStopwatch.h"
#include "TStyle.h"
#include "TLatex.h"
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TPaveLabel.h>
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TRandom3.h"
#include<fstream>
#include<sstream>
#include <vector>
#include "TPaveText.h"
#include <time.h>

#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/roofit/5.34.18-cms4/include/RooHist.h"
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/roofit/5.34.18-cms4/include/RooAbsDataStore.h"
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/roofit/5.34.18-cms4/include/RooAbsData.h"
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/roofit/5.34.18-cms4/include/RooGlobalFunc.h"
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/roofit/5.34.18-cms4/include/RooMsgService.h"
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/roofit/5.34.18-cms4/include/RooRealVar.h"
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/roofit/5.34.18-cms4/include/RooDataHist.h"
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/roofit/5.34.18-cms4/include/RooHistPdf.h"
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/roofit/5.34.18-cms4/include/RooDataSet.h"
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/roofit/5.34.18-cms4/include/RooRandom.h"

using namespace RooFit ;

class AFBMeasurer{
  private:
    // Members
    TFile* inputFile_  ;
    TFile* outputFile_ ;
    
    TH1F* hBase_cosThetaCS_ ;
    TH1F* h_cosThetaCS_likelihoodafbmin1_ ;
    TH1F* h_cosThetaCS_likelihoodafb1_ ;
    TH1F* h_cosThetaCS_genZprimeSSM_ ;
    TH1F* h_cosThetaCS_gen_ ;
    TH1F* h_between_ ;
    TH1D* h_results_ ;
    
    int nBins_quickscan_ ;
    int nBins_fullscan_ ;
    int NPE_ ;
    int events_per_PE_ ;
    float lumiScale_ ;
    int   nBins_cosThetaCS_1d_ ; // Binning used in templates
    float lower_cosThetaCS_ ;    // Range for CosThetaCS
    float upper_cosThetaCS_ ;    // Range for CosThetaCS
    TString theCuts_ ;
    TString thevars1d_ ;
    float massmin_ ;
    float massmax_ ;
    vector<double> PEpoints_ ;
    
    TTree * PEDataTree_ ;
    double costhtree_ ;
    double weighttree_ ;
    TBranch* b_costhtree_ ;
    TBranch* b_weighttree_ ;
    bool debug_ ;
    
    TCanvas* canvas_ ;
    
  public:
    void setNPseudoExperiments(int NPE){ NPE_ = NPE ; }
  
    AFBMeasurer(TString signaltype="ZprimePSI", TString massminString="2700", TString massmaxString="3300", TString lumi="3000fb", TString phasetype="Phase2"):
    inputFile_(0),
    outputFile_(0),
    hBase_cosThetaCS_(0),
    h_cosThetaCS_likelihoodafbmin1_(0),
    h_cosThetaCS_likelihoodafb1_(0),
    h_cosThetaCS_genZprimeSSM_(0),
    h_cosThetaCS_gen_(0),
    h_between_(0),
    h_results_(0),
    nBins_quickscan_(100),
    nBins_fullscan_(501),
    NPE_(1000),
    events_per_PE_(1),
    lumiScale_(1.0),
    nBins_cosThetaCS_1d_(100),
    lower_cosThetaCS_(-1),
    upper_cosThetaCS_( 1),
    theCuts_("1==1"),
    thevars1d_("costhetaCSheepheep"),
    PEDataTree_(0),
    costhtree_(0),
    weighttree_(0),
    b_costhtree_(0),
    b_weighttree_(0),
    debug_(false),
    canvas_(0){
      RooAbsData::setDefaultStorageType(RooAbsData::Tree);
      RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
      h_results_ = new TH1D("h_results","h_results",1001,-1.001,1.001); 
      
      /*For now, I give by hand the number of events passing the final selection in each case. 
      Yes it's ugly.*/
      if     (lumi=="100fb" ) lumiScale_ =  1.0/30.0 ;
      else if(lumi=="300fb" ) lumiScale_ =  0.1 ;
      else if(lumi=="3000fb") lumiScale_ =  1.0 ;
    
      if(signaltype=="ZprimePSI"){
        if(massminString=="2700") events_per_PE_ = round(1056*lumiScale_) ;
        if(massminString=="3600") events_per_PE_ = round( 128*lumiScale_) ;
      }
      else if(signaltype=="ZprimeSSM"){
        if(massminString=="2700") events_per_PE_ = round(2180*lumiScale_) ;
        if(massminString=="3600") events_per_PE_ = round( 369*lumiScale_) ;
        if     (lumi=="1000fb" &&massminString== "2000" && phasetype=="Phase1NoAging") events_per_PE_ = 3744 ; //events_per_PE_ = 3744 ; //Phase 1 no aging
        else if(lumi=="1000fb" &&massminString== "2000" && phasetype=="Phase11KAging") events_per_PE_ = 2598 ; //events_per_PE_ = 2598 ; //Phase 1 aging
        else if(lumi=="1000fb" &&massminString== "2000" && phasetype=="Phase2"       ) events_per_PE_ = 3541 ; //events_per_PE_ = 3541 ; //Phase 2
      }
      else if(signaltype=="ZprimeI"){
        if     (lumi=="100fb"  && massminString=="2700") events_per_PE_ =   49 ; //35 ;
        else if(lumi=="300fb"  && massminString=="2700") events_per_PE_ =  148 ;
        else if(lumi=="3000fb" && massminString=="2700") events_per_PE_ = 1485 ;
        else if(lumi=="100fb"  && massminString=="2700") events_per_PE_ =    5 ; //3 ;
        else if(lumi=="300fb"  && massminString=="2700") events_per_PE_ =   15 ;
        else if(lumi=="3000fb" && massminString=="3600") events_per_PE_ =  151 ;
      }
      else{
        // DY events
        events_per_PE_ = 11 ; //events_per_PE_=1056;
      }
      std::cout << "Events per psuedoexperiment = " << events_per_PE_ << std::endl ;
      
      for(int i=0 ; i<events_per_PE_ ; i++){ PEpoints_.push_back(0) ; }
      canvas_ = new TCanvas("canavs", "", 0, 0, 800, 600) ;
      
      // Call the input file, which is the output of the macro AccforSpin0
      TString inputFileName = Form("%s.root", phasetype.Data()) ;
      inputFile_ = new TFile(inputFileName, "OPEN") ;
      std::cout << "Input file name = " << inputFileName << std::endl ;
      
      inputFile_->cd() ;
      //TTree *thetreeZprimePSI = (TTree*)(inputFile_)->Get("TreeZprimePSI") ;
      //TTree *thetreeDY        = (TTree*)(inputFile_)->Get("TreeDY"       ) ;
      TTree *thetreeafb1      = (TTree*)(inputFile_)->Get("TreeAFB1"     ) ;
      TTree *thetreeafbmin1   = (TTree*)(inputFile_)->Get("TreeAFBmin1"  ) ;
      TTree *thetreeZprimeSSM = (TTree*)(inputFile_)->Get("TreeAFBmin1"  ) ;
      //TTree *thetreeZprimeSSM = (TTree*)(inputFile_)->Get("TreeZprimeSSM") ;
      
      TString outfilename = Form("output_Afb1D_%sm%sto%s%s_%s.root", signaltype.Data(), massminString.Data(), massmaxString.Data(), lumi.Data(), phasetype.Data()) ;
      outputFile_ = new TFile(outfilename,"RECREATE") ;
      outputFile_->cd() ;
  
      // This is the cuts you want to apply on the events
      // The last one removes the case where costhetaCSheepheep couldn't be calculated and has the initial value of -1000
      theCuts_ = Form("weightevt*(recomass>%s && recomass<%s && (gsfBB||gsfBE||gsfEE) && heepheep && costhetaCSheepheep>-10", massminString.Data(), massmaxString.Data()) ;
      thevars1d_ = "costhetaCSheepheep" ;
      
      // Declaration of some histo to generate the Pseudo Experiments and calculate the likelihood
      hBase_cosThetaCS_ = new TH1F("h_cosThetaCS", "", nBins_cosThetaCS_1d_, lower_cosThetaCS_, upper_cosThetaCS_) ;
      hBase_cosThetaCS_->Sumw2() ;
      h_cosThetaCS_likelihoodafbmin1_ = (TH1F*) hBase_cosThetaCS_->Clone("h_costhetaCS_likelihoodafbmin1") ;
      h_cosThetaCS_likelihoodafb1_    = (TH1F*) hBase_cosThetaCS_->Clone("h_costhetaCS_likelihoodafb1"   ) ;
      h_cosThetaCS_genZprimeSSM_      = (TH1F*) hBase_cosThetaCS_->Clone("h_costhetaCS_genZprimeSSM"     ) ;
      
      //Note the extra cut on recomass: events only enter one of the three histogram
      CreateHisto(thetreeafbmin1  , h_cosThetaCS_likelihoodafbmin1_, Form("%s && recomass%%3==1)", theCuts_.Data()), thevars1d_) ;
      CreateHisto(thetreeafb1     , h_cosThetaCS_likelihoodafb1_   , Form("%s && recomass%%3==2)", theCuts_.Data()), thevars1d_) ;
      CreateHisto(thetreeZprimeSSM, h_cosThetaCS_genZprimeSSM_     , Form("%s && recomass%%3==0)", theCuts_.Data()), thevars1d_) ;
      std::cout << h_cosThetaCS_likelihoodafbmin1_->GetEntries() << " " << Form("%s && recomass%%3==1)", theCuts_.Data()) << std::endl ;
      
      h_cosThetaCS_gen_ = (TH1F*) h_cosThetaCS_genZprimeSSM_->Clone("h_cosThetaCS_gen") ;
    }
    ~AFBMeasurer(){
      cout << "Destructing" << endl ;
      if(inputFile_ ) delete inputFile_  ;
      if(outputFile_) delete outputFile_ ;
    
      if(hBase_cosThetaCS_              ) delete hBase_cosThetaCS_ ;
      
      if(h_cosThetaCS_likelihoodafbmin1_) delete h_cosThetaCS_likelihoodafbmin1_ ;
      if(h_cosThetaCS_likelihoodafb1_   ) delete h_cosThetaCS_likelihoodafb1_ ;
      if(h_cosThetaCS_genZprimeSSM_     ) delete h_cosThetaCS_genZprimeSSM_ ;
      if(h_cosThetaCS_gen_              ) delete h_cosThetaCS_gen_ ;
      if(h_between_                     ) delete h_between_ ;
      if(h_results_                     ) delete h_results_ ;
    
      if(PEDataTree_  ) delete PEDataTree_ ;
      if(b_costhtree_ ) delete b_costhtree_ ;
      if(b_weighttree_) delete b_weighttree_ ;
      if(canvas_      ) delete canvas_ ;
    }
    
    void run(){
      //Repeat the following procedure 1000 times (this may be too small)
      for(int k=0 ; k<NPE_ ; ++k){
        cout << "k = " << k << endl ;
        // Step 1: generate 'events_per_PE_' events according to 'histcosthetaCSforgen' and fill 'PEpoints_'
        GeneratePseudoExperiment(lower_cosThetaCS_, upper_cosThetaCS_, h_cosThetaCS_gen_) ;
        
        // Step 2. First loop to localize roughly the maximum of the likelihood
        double logLmin = -1000000 ; //Minimal log likelihood value considered. If many events, it's maybe too high. 
        double minit = -10 ;
        std::vector<double> fracs ;
        std::vector<double> logLs ;
        for(int i=0 ; i<nBins_quickscan_+1 ; ++i){
          double frac = 1.0*i/(1.0*nBins_quickscan_) ;
          // Compute the likelihood for the linear combination (frac)*histoAfbMin1+(1-frac)histoAfb1
          // Build a pdf corresponding  to the current afb
          h_between_ = (TH1F*) h_cosThetaCS_likelihoodafbmin1_->Clone("h_between_") ;
          h_between_->Add(h_cosThetaCS_likelihoodafbmin1_,-(frac)) ;
          h_between_->Add(h_cosThetaCS_likelihoodafb1_   ,  frac ) ;
          if(k==50){
            TH1F* h_betweenClone = (TH1F*) h_between_->Clone(Form("h_between_%d", i)) ;
            h_betweenClone->Write() ;
          }
          
          double logL = computelikelihood(PEpoints_,h_between_) ;
          fracs.push_back(frac) ;
          logLs.push_back(logL) ;

          if(logL>logLmin){
            logLmin = logL ;
            minit = i ;
          }
          delete h_between_ ;
        }
        
        // Step 3. Find the interval such that (LogLmin - LogL)<MaxDev
        double MaxDev   = 0.5 ;
        double lowedge  = 0   ;
        double highedge = 1   ;
        for(int i=0 ; i<nBins_quickscan_+1 ; ++i){
          if(i<minit && (logLmin - logLs.at(i))<MaxDev){
            lowedge = fracs.at(i) ;
            i = minit ;
          }
          if(i>minit && (logLmin - logLs.at(i))>MaxDev){
            highedge = fracs.at(i) ;
            break ;
          }
        }
    
        // Step 4. Make accurate measurement
        const Int_t nBinsTmp = nBins_fullscan_ ;
        double  afbaccurate[nBinsTmp] ;
        double logLaccurate[nBinsTmp] ;
        double minLogFinal = -1000000000 ;
        double valueFinal = 0 ;
        for(int i=0 ; i<nBins_fullscan_ ; ++i){
          double frac= lowedge+(highedge-lowedge)*i/(nBins_fullscan_) ;
          h_between_ = (TH1F*) h_cosThetaCS_likelihoodafbmin1_->Clone("h_between_") ;
          h_between_->Add(h_cosThetaCS_likelihoodafbmin1_, -(frac)) ;
          h_between_->Add(h_cosThetaCS_likelihoodafb1_   ,   frac ) ;
      
          afbaccurate[i] = (2.0*frac-1.0)*3.0/4.0 ; //This is just the relation between frac and the afb (the demonstration is in the AN)
          logLaccurate[i] = computelikelihood(PEpoints_,h_between_) ;
          if(logLaccurate[i]>minLogFinal){
            minLogFinal = logLaccurate[i] ;
            valueFinal  = (2.0*frac-1.0)*3.0/4.0 ;
          }
          delete h_between_ ;
        }
  
        // Save the likelihood curve for the 4th PE generated.
        // Why the fourth? I didn't want the first one (in case something screws up at the end of the first iteration) and I wanted to have it even if running on
        // only 10 PE. 
        if(k==3){
          canvas_->cd() ;
          TGraph* mygraph = new TGraph(501, afbaccurate, logLaccurate) ;
          mygraph->Draw("AP") ;
          mygraph->Write() ;
        }
        if(k%25==0 || debug_){
          std::cout <<  "k = "<< k << std::endl ;
          std::cout << "best estimate " << valueFinal << " " << minLogFinal << std::endl ;
          std::cout << "range " << (2*lowedge-1)*3.0/4.0 << " , " << (2*highedge-1)*3.0/4.0 << std::endl ;
        }
        
        h_results_->Fill(valueFinal) ;
      }
      Float_t max = h_results_->GetMaximum() ;
      for(Int_t bin=1 ; bin<=h_results_->GetNbinsX() ; bin++){
        h_results_->SetBinContent(bin, h_results_->GetBinContent(bin)-max) ;
      }
      h_results_->Scale(-2.0) ;
      h_results_->Write() ;
    }
  
  private:
    
    // Methods
    void CreateHisto(TTree* mytree, TH1F* histo, TString cuts, TString vars){
      TString hNameTmp = histo->GetName() ;
      TString hName = Form("%s%s", hNameTmp.Data(),vars.Data()) ;
      std::vector<TString> patternsToRemove ;
      patternsToRemove.push_back(":") ;
      patternsToRemove.push_back("heep") ;
      patternsToRemove.push_back("(") ;
      patternsToRemove.push_back(")") ;
      for(unsigned int i=0 ; i<patternsToRemove.size() ; ++i){
        hName.ReplaceAll(patternsToRemove.at(i), "") ;
      }
      
      histo->SetName(hName) ;
      cout << "Creating histogram with name: " << hName << endl ;
      
      TString vartodraw = Form("%s>>%s", vars.Data(), hName.Data()) ;
      mytree->Draw(vartodraw,cuts,"LEP") ;
      histo->Scale(1/histo->Integral()) ;
      histo->Write() ;
    }
    
    void GeneratePseudoExperiment(float rrv_low, float rrv_high, TH1F* hgen){
      // This is just basic Roofit, one generates a PE according to a binned PDF obtained from a histogram
      RooRealVar varx("varx", "costhetaCS", rrv_low, rrv_high) ;
      
      // zzzz Here is where the segfaults are after compiling
      //RooDataHist roodatahist("roodatahist", "roodatahist", RooArgList(varx), Import((TH1&)*hgen)) ;
      RooDataHist roodatahist("roodatahist", "roodatahist", RooArgSet(varx), hgen) ;
      RooHistPdf    pdfforgen("pdfforgen"  , "pdf"        , RooArgSet(varx), roodatahist, 0) ;
      
      RooRandom::randomGenerator()->SetSeed() ;
      
      RooDataSet* data = pdfforgen.generate(RooArgList(varx), events_per_PE_) ;
      const RooAbsData* rad = (const RooAbsData*) data ;
      TTree* tree = (TTree*) rad->store()->tree() ;
      PEDataTree_ = tree->CloneTree() ;
      PEDataTree_->SetBranchAddress("varx"  , &costhtree_ , &b_costhtree_ ) ;
      PEDataTree_->SetBranchAddress("weight", &weighttree_, &b_weighttree_) ;
      
      //Loop over the TTree entries to count the nb of events in each bin
      PEpoints_.clear() ;
      for(int itree=0 ; itree<PEDataTree_->GetEntries() ; ++itree){
        PEDataTree_->GetEntry(itree) ;
        int eventsInTheBin = weighttree_ ;
        while(eventsInTheBin>0){
          PEpoints_.push_back(costhtree_) ;
          eventsInTheBin-- ;
        }
      }
      
      delete data;
      delete PEDataTree_;
    }
    
    double computelikelihood(vector<double> mypoints, TH1F* myh){
      double logLvalue = 0 ;
      for(unsigned int i=0 ; i<mypoints.size() ; ++i){
        double rrv_value = mypoints.at(i) ;
        int binnb = myh->FindBin(rrv_value) ;
        double proba = myh->GetBinContent(binnb)/myh->Integral() ;
        logLvalue += log(proba) ;
        //std::cout << i << " " << mypoints.size() << " " << proba << " " << logLvalue << std::endl ;
      }
      //std::cout << logLvalue << std::endl ;
      return logLvalue;
    }
} ;

void AfbMeasurement_1Dnew(){
  TString signaltype = "ZprimeSSM" ;
  TString massmin   =       "2000" ;
  TString massmax   =       "3000" ;
  TString lumi      =     "1000fb" ;
  TString phasetype =     "Phase2" ;
  
  AFBMeasurer* AFB = new AFBMeasurer(signaltype, massmin, massmax, lumi, phasetype) ;
  AFB->setNPseudoExperiments(1000) ;
  AFB->run() ;
  cout << "Tada!" << endl ;
}


