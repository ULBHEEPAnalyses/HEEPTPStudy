#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TLatex.h"
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include "TPaveLabel.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TRandom3.h"
using namespace RooFit;
#include<fstream>
#include<sstream>
#include <vector>
#include "TPaveText.h"
#include <time.h>
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

TFile *outputfile;  

double raptree ;
double costhtree;
double weighttree; 
TTree * mydatatree;
TBranch * b_raptree ;
TBranch * b_costhtree;
TBranch * b_weighttree;
bool debug = false;

void AfbMeasurement_1Dnew(TString signaltype = "ZprimePSI", TString massmin= "2700",TString massmax = "3300" , TString lumi ="3000fb", TString phasetype="Phase2"){
  RooAbsData::setDefaultStorageType(RooAbsData::Tree);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  TH1D * myh = new TH1D("myh","myh",1001,-1.001,1.001); 
  
  int evtsperPE =1;

  /*For now, I give by hand the number of events passing the final selection in each case. 
  Yes it's ugly.*/
  if(signaltype == "ZprimePSI" && lumi=="3000fb" &&massmin== "2700" ) evtsperPE =1056;
  else if(signaltype == "ZprimePSI" && lumi=="300fb" &&massmin== "2700") evtsperPE =106;
  else if(signaltype == "ZprimeSSM" && lumi=="3000fb" &&massmin== "2700") evtsperPE =3180;
  else if(signaltype == "ZprimeSSM" && lumi=="300fb" &&massmin== "2700") evtsperPE =318;
  else if(signaltype == "ZprimePSI" && lumi=="3000fb" &&massmin== "3600" ) evtsperPE =128;
  else if(signaltype == "ZprimePSI" && lumi=="300fb" &&massmin== "3600") evtsperPE =13;
  else if(signaltype == "ZprimeSSM" && lumi=="3000fb" &&massmin== "3600") evtsperPE =369;
  else if(signaltype == "ZprimeSSM" && lumi=="300fb" &&massmin== "3600") evtsperPE =37;
  else if(signaltype == "ZprimeI" && lumi=="100fb" &&massmin== "2700") evtsperPE =35;//13 TeV !
  else if(signaltype == "ZprimeI" && lumi=="300fb" &&massmin== "2700") evtsperPE =148;
  else if(signaltype == "ZprimeI" && lumi=="3000fb" &&massmin== "2700") evtsperPE =1485;
  else if(signaltype == "ZprimeI" && lumi=="100fb" &&massmin== "3600") evtsperPE =3;
  else if(signaltype == "ZprimeI" && lumi=="300fb" &&massmin== "3600") evtsperPE =15;
  else if(signaltype == "ZprimeI" && lumi=="3000fb" &&massmin== "3600") evtsperPE =151;//13 TeV!
  else if(signaltype != "ZprimePSI"&&signaltype != "ZprimeSSM" &&signaltype != "ZprimeI" ) evtsperPE=1056;
  else if(signaltype == "ZprimeSSM" && lumi=="1000fb" &&massmin== "2000" && phasetype=="Phase1NoAging") evtsperPE=3744; //Phase 1 no aging
  else if(signaltype == "ZprimeSSM" && lumi=="1000fb" &&massmin== "2000" && phasetype=="Phase11kAging" ) evtsperPE =2598; //Phase 1 aging
  else if(signaltype == "ZprimeSSM" && lumi=="1000fb" &&massmin== "2000"&& phasetype=="Phase2" ) evtsperPE =3541; //Phase 2 
  else {cout << "Problem in the name and nb of evts"<<endl; return;}
  
  TString outfilename ="output_Afb1D_"+signaltype+"m"+massmin+"to"+massmax+lumi+"_"+phasetype+".root";
  outputfile = new TFile(outfilename,"RECREATE"); 
 
 

  /*Set the binning used in the templates 
  Here, I chose 100 bins. This is probably too large.*/
  const int nbbinscosthetaCS_1d = 100; 
  const float lowbincosthetaCS = -1;  
  const float highbincosthetaCS = 1; 

  //Call the input file, which is the output of the macro AccforSpin0
  TFile * myfile;
  if(phasetype =="Phase1NoAging")myfile = new TFile("Phase1noaging.root","open"); 
  else if(phasetype =="Phase11kAging")myfile = new TFile("Phase11kaging.root","open"); 
  else  if(phasetype =="Phase2") myfile = new TFile("Phase2.root","open"); 
  else myfile=  new TFile("Test13TeV_ggvsqqbar.root","open"); 
  myfile->cd();
  TTree *thetreeZprimePSI= (TTree*)(myfile)->Get("TreeZprimePSI");
  TTree *thetreeDY = (TTree*)(myfile)->Get("TreeDY");
  TTree *thetreeafb1 = (TTree*)(myfile)->Get("TreeAFB1");
  TTree *thetreeafbmin1 = (TTree*)(myfile)->Get("TreeAFBmin1");
  TTree *thetreeZprimeSSM = (TTree*)(myfile)->Get("TreeZprimeSSM");
  
  outputfile->cd(); 

  /*This is the cuts you want to apply on the events
   The last one removes the case where costhetaCSheepheep couldn't be calculated and has the initial value of -1000;
   */
  TString thecuts = "weightevt*(recomass>"+massmin +"&&recomass<"+massmax+" &&(gsfBB||gsfBE||gsfEE)&&heepheep&&costhetaCSheepheep>-10";

  TString thevars1d = "costhetaCSheepheep";
  
  /*Declaration of some histo to generate the Pseudo Experiments and calculate the likelihood*/
 
  TH1F * histcosthetaCSforlikelihoodafbmin1= new TH1F("histcosthetaCSforlikelihoodafbmin1","",nbbinscosthetaCS_1d,lowbincosthetaCS,highbincosthetaCS); 
  TH1F * histcosthetaCSforlikelihoodafb1= new TH1F("histcosthetaCSforlikelihoodafb1","",nbbinscosthetaCS_1d,lowbincosthetaCS,highbincosthetaCS); 
  TH1F * histcosthetaCSforgenZprimeSSM= new TH1F("histcosthetaCSforgenZprimeSSM","",nbbinscosthetaCS_1d,lowbincosthetaCS,highbincosthetaCS); 


  //Note the extra cut on recomass: events only enter one of the three histogram
  CreateHisto(thetreeafbmin1,histcosthetaCSforlikelihoodafbmin1,thecuts+"&&recomass%3==1)", thevars1d); 
  CreateHisto(thetreeafb1,histcosthetaCSforlikelihoodafb1,thecuts+"&&recomass%3==2)", thevars1d); 
  CreateHisto(thetreeZprimeSSM,histcosthetaCSforgenZprimeSSM,thecuts+"&& recomass%3==0)", thevars1d); 
  

  vector<double> myPEpoints(evtsperPE);
 
  TH1F* histcosthetaCSforgen ; 
  
  histcosthetaCSforgen= (TH1F*)histcosthetaCSforgenZprimeSSM->Clone("histcosthetaCSforgen");
  

 
  //Repeat the following procedure 1000 times (this may be too small)
  for(int k = 0; k<1000; k ++){ 
    //Step 1: generate 'evtsperPE' events according to 'histcosthetaCSforgen' and fill 'myPEpoints'
    cout << "entering genpe" << endl;
    GenPE("costhetaCS",lowbincosthetaCS,highbincosthetaCS, histcosthetaCSforgen,myPEpoints);
    cout << "leavining genpe" << endl;
    //Step 2. First loop to localize roughly the maximum of the likelihood
    double logLmin = -1000000;//Minimal log likelihood value considered. If many events, it's maybe too high. 
    double minit = -10; 
    const int nbbinsquickscan = 100;
    double fracvslogLrough[nbbinsquickscan+1][2]; 
    TH1F * histbetween; 
    for(int i = 0; i<nbbinsquickscan+1 ; i ++){
      double frac= 1/(double)nbbinsquickscan*(double)i; 
      //Compute the likelihood for the linear combination (frac)*histoAfbMin1+(1-frac)histoAfb1
      //  Build a pdf corresponding  to the current afb
      histbetween=(TH1F*)histcosthetaCSforlikelihoodafbmin1->Clone("histbetween");
      histbetween->Add(histcosthetaCSforlikelihoodafbmin1,-(frac)); 
      histbetween->Add(histcosthetaCSforlikelihoodafb1,frac); 
      fracvslogLrough[i][0] = frac; 
      fracvslogLrough[i][1] = computelikelihood(myPEpoints,histbetween) ; 

      if(fracvslogLrough[i][1]>  logLmin){ logLmin  =  fracvslogLrough[i][1]   ; minit=i;}
      delete histbetween;
    }
    cout << "leavining genpe" << endl;
    //Step 3. Find the interval such that ( LogLmin - LogL)<MaxDev
    double MaxDev= 0.5; 
    double lowedge = 0;
    double highedge = 1;
    for(int i = 0; i<nbbinsquickscan+1 ; i ++){
      //double frac= -1+ 2/(double)nbbinsquickscan*(double)i; //orig
      double frac= 1/(double)nbbinsquickscan*(double)i; 
      if( i < minit&& (logLmin - fracvslogLrough[i][1] )<MaxDev){
	lowedge = fracvslogLrough[i][0]; 
	i = minit;
      }
      if( i>minit && (logLmin - fracvslogLrough[i][1] )>MaxDev){
	highedge = fracvslogLrough[i][0]; 
	break;
      }
    }
    
    //Step 4. Make accurate measurement  
    double afbaccurate[501];
    double logLaccurate[501];
    double minLogFinal = -1000000000;
    double valueFinal = 0;
    for(int i = 0; i<501;i++){
	double frac= lowedge+ (highedge-lowedge)/500*i; 
	
	histbetween=(TH1F*)histcosthetaCSforlikelihoodafbmin1->Clone("histbetween");
	histbetween->Add(histcosthetaCSforlikelihoodafbmin1,-(frac)); 
	histbetween->Add(histcosthetaCSforlikelihoodafb1,frac); 
	
	afbaccurate[i]=(double) 3/4*(2*frac-1); //This is just the relation between frac and the afb (the demonstration is in the AN)
	logLaccurate[i] = computelikelihood(myPEpoints,histbetween) ; 
	if(logLaccurate[i]>minLogFinal) {minLogFinal =logLaccurate[i]; valueFinal = (double) 3/4*(2*frac-1); }
	
	delete histbetween;
      } 
  
    //Save the likelihood curve for the 4th PE generated.
    //Why the fourth? I didn't want the first one (in case something screws up at the end of the first iteration) and I wanted to have it even if running on
    //only 10 PE. 
    if(k ==3)mygraph = new TGraph( 501, afbaccurate, logLaccurate); 
    if(k ==3)mygraph->Draw("AP");
    if(k ==3)mygraph->Write();
    if(k ==3)c1->Write(); 
    if(k %25==0||debug )  cout <<  "k = "<< k << endl;
    if(k %25==0||debug )cout << "best estimate " << valueFinal<< endl; 
    if(k %25==0||debug )cout << "range " << (double)3/4*(2*lowedge-1) << " , "<< (double)3/4*(2*highedge-1) << endl; 
    
    
    myh->Fill(valueFinal);

  }

  myh->Draw("");
  myh->Write();

  gApplication->Terminate();

}








void Initialize(TString samplefile, TString mysamplename, float  sampleweight)
{
  TFile *in =  TFile::Open(samplefile);
  input.push_back(in); 
  sweight.push_back(sampleweight);
  samplename.push_back(mysamplename);
}






void CreateHisto(TTree* mytree, TH2F* histo, TString cuts, TString vars){
  TString histoname = histo->GetName()+vars; 
  histoname.ReplaceAll(":",""); 
  histoname.ReplaceAll("heep",""); 
  histoname.ReplaceAll("(",""); 
  histoname.ReplaceAll(")",""); 
   

  histo->SetName(histoname); 
  cout << "hname " <<histoname << endl; 
  
  
  histo->Sumw2();
  TString vartodraw= vars+">>"+histoname;
  mytree->Draw(vartodraw,cuts,"ZCOL");

  histo->Scale(1/histo->Integral());  
  histo->Write();
}

void CreateHisto(TTree* mytree, TH1F* histo, TString cuts, TString vars){
  TString histoname = histo->GetName()+vars; 
  histoname.ReplaceAll(":",""); 
  histoname.ReplaceAll("heep",""); 
  histoname.ReplaceAll("(",""); 
  histoname.ReplaceAll(")",""); 
   

  histo->SetName(histoname); 
  cout << "hname " <<histoname << endl; 
  
  histo->Sumw2();
  TString vartodraw= vars+">>"+histoname;
  mytree->Draw(vartodraw,cuts,"LEP");
  histo->Scale(1/histo->Integral());  
  histo->Write();
  
}





void  GenPE(TString varxname, float varxlow, float varxhigh, TH1F* hgen, vector<double>   & costh_rap){
 
  //This is just basic Roofit, one generates a PE according to a binned PDF obtained from a histogram
  RooRealVar varx("varx",varxname,varxlow,varxhigh) ;
  RooDataHist roodatahist("roodatahist", "roodatahist", RooArgList(varx), hgen);
  RooHistPdf pdfforgen("pdfforgen","pdf",RooArgList(varx), roodatahist,0) ; 

  RooRandom::randomGenerator()->SetSeed(); 
  
  RooDataSet * data = pdfforgen.generate(RooArgList(varx),  costh_rap.size());
  //The pseudodata are stored in a TTree. 
  mydatatree = (TTree*)((RooTreeDataStore*)(data->store())->tree()->CloneTree());
  mydatatree->SetBranchAddress("varx", &costhtree, &b_costhtree);
  mydatatree->SetBranchAddress("weight", &weighttree, &b_weighttree);
  
  //Loop over the TTree entries to count the nb of events in each bin
  int j =0;

  for(int itree =0; itree< mydatatree->GetEntries();itree++){
    mydatatree->GetEntry(itree);
    int nbevtinthebin = weighttree;
    while(nbevtinthebin>0 ){
      costh_rap[j]=costhtree;
      j++; 
      nbevtinthebin--; 
    }
  }
  
  delete data;
  delete mydatatree;
 
}





double computelikelihood(vector<double> mypoints, TH1F * myh){

  double logLvalue = 0; 
  for(int i = 0; i< mypoints.size(); i++){
    double varxvalue = mypoints[i];
   
    int binnb = myh->FindBin(varxvalue);
    double proba = myh->GetBinContent(binnb)/myh->Integral();


    if(debug) cout << "i: " <<i << endl;
    if(debug)  cout << " costh , p " << mypoints[i]<<  " , " << myh->GetBinContent(binnb)/myh->Integral()<<endl; 
    logLvalue += log(proba); 
    if(debug) cout << "varxvalue, proba " << varxvalue << " , " << proba << endl; 
    if(debug) cout << "current logLvalue " << logLvalue << endl; 
  }
  
  return logLvalue;
}




