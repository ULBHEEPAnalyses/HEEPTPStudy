#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TGraph.h"
// #include "RooAbsPdf.h"
// #include "RooAddPdf.h"
// #include "RooDataHist.h"
// #include "RooArgList.h"
// #include "RooBreitWigner.h"
// #include "RooCBShape.h"
// #include "RooDataSet.h"
// #include "RooExponential.h"
// #include "RooFFTConvPdf.h"
// #include "RooGaussian.h"
// #include "RooGenericPdf.h"
// #include "RooPlot.h"
// #include "RooRealVar.h"
// #include "RooWorkspace.h"
// #include "RooFFTConvPdf.h"
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
//#include "RooFitResult.h"
using namespace RooFit;
//#include "tdrstyle.C"
#include<fstream>
#include<sstream>
#include <vector>
#include <iomanip>
#include <iostream>
//#include "AnaFuncs.hh"
#include "TPaveText.h"
#include <time.h>
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
TFile * outfile = new TFile("afbfor.root","RECREATE"); 
gStyle->SetOptStat(0);
TString labelname = "CMS Simulation, #sqrt{s} = 13 TeV" ;
TPaveLabel *labelmcordata = new TPaveLabel(0.1,0.86,0.8,1.0,labelname,"brNDC");
bool savepdf =false;

void ComputeStatUncertainty(){ 
  gStyle->SetOptStat(0);
  std::cout.precision(5); 

  //The function "calcuncertainty" calculates the stat uncertainty in the file "output_XXX"
  //Note that it also works for the qq/qqbar fraction
  calcuncertainty("output_Afb1D_ZprimeSSMm2000to30001000fb_Phase2.root"); 



  // return ;

  return;
}



void calcuncertainty(TString filename){
  TString labeldescst =""; 
  if(filename.Index("ZprimePSI")!=kNPOS) labeldescst += "$Z'_{\\Psi}$"; 
  if(filename.Index("ZprimeSSM")!=kNPOS) labeldescst += "$Z'_{SSM}$"; 
  if(filename.Index("DY")!=kNPOS) labeldescst += "Drell-Yan"; 
  if(filename.Index("RSGrav")!=kNPOS) labeldescst += "R-S graviton"; 

  if(filename.Index("m2700to3300")!=kNPOS) labeldescst += " (M = 3 TeV/c$^{2})$"; 
  if(filename.Index("m3600to4400")!=kNPOS) labeldescst += " (M = 4 TeV/c$^{2})$"; 
  if(filename.Index("m2000to2000")!=kNPOS) labeldescst += " (M = 2.5 TeV/c$^{2})$"; 

  if(filename.Index("1D")!=kNPOS)labeldescst += " 1D"; 
  if(filename.Index("2D")!=kNPOS)labeldescst += " 2D"; 
  
 


  TFile * myf = new TFile(filename,"open");
  TH1F* myh = (TH1F*)myf->Get("myh");
  
  double currentint = 0;
  double lowedge = -10; 
  double median =-10; 
  double highedge = -10; 
  for(int i =0;i<= myh->GetNbinsX()+1;i++){
    
    currentint += (double) myh->GetBinContent(i)/(double)myh->Integral(0,myh->GetNbinsX()+1); 
   
    if(currentint>=0.159 &&lowedge==-10 )lowedge = (double)myh->GetBinLowEdge(i);
    if(currentint>=0.5 && median ==-10)median  = (double)myh->GetBinCenter(i);
    if(currentint>=0.841&&highedge==-10 && i != myh->GetNbinsX()+1)highedge = myh->GetBinLowEdge(i+1);
  
  }
  lowedge -= median; 
  highedge -=median; 
  if(filename.Index("100fb") !=kNPOS ||filename.Index("DY")!=kNPOS ) cout << labeldescst ; 
  if(filename.Index("DY")!=kNPOS) cout << "& & " ; 
  cout <<"& $" <<fixed << median << "_{"<< fixed << setprecision(5)<< lowedge<< "}^{+"<< setprecision(5)<<highedge<<"}$ "; 
  if(filename.Index("3000fb") !=kNPOS  ||filename.Index("DY")!=kNPOS ) {  cout << "\\\\" <<endl; }
  
  
}



//Not used here
void calcsignif(TString filename, TString h0name, TString h1name, TString stringinfo)
{
  labelmcordata->SetFillColor(0);
  labelmcordata->SetFillStyle(0);
  labelmcordata->SetBorderSize(0);
  labelmcordata->SetTextSize(0.3); 
  labelmcordata->SetTextAlign(12);
  
  
  TString labeldescst =""; 
  TString h0label = "";
  TString h1label ="";
  TString labellumi = "";
  if(filename.Index("SignalZprimePSI") !=kNPOS) { labeldescst+="#splitline{Z'_{#Psi} "; h1label ="Z'_{#Psi}";}
  if(filename.Index("SignalRSGrav") !=kNPOS)  { labeldescst+="#splitline{G_{RS} "; h1label = "G_{RS}"; }
  if(filename.Index("SignalSpin0") !=kNPOS)  { labeldescst+="#splitline{S_{0} "; h1label ="S_{0}";}
  if(filename.Index("H0ZprimePSI") !=kNPOS) {  labeldescst+="vs S_{1} (q#bar{q})";h0label = "S_{1} (q#bar{q})";}  
  if(filename.Index("H0RSGrav") !=kNPOS)  { labeldescst+="vs G_{RS}";  h0label = "G_{RS}"; }
  if(filename.Index("H0Spin0qqbarggmix") !=kNPOS) {  labeldescst+="vs S_{0} (q#bar{q}-gg)";h0label ="S_{0} (q#bar{q}-gg)"; }
  else if(filename.Index("H0Spin0qqbar") !=kNPOS) {  labeldescst+="vs S_{0} (q#bar{q})";h0label ="S_{0} (q#bar{q})"; }
  else if(filename.Index("H0Spin0gg") !=kNPOS) {  labeldescst+="vs S_{0} (gg)";h0label ="S_{0} (gg)"; }
  if(filename.Index("H0Spin2qqbarggmix") !=kNPOS) {  labeldescst+="vs S_{2} (q#bar{q}-gg)";h0label ="S_{2} (q#bar{q}-gg)"; }
  else if(filename.Index("H0Spin2qqbar") !=kNPOS) {  labeldescst+="vs S_{2} (q#bar{q})";h0label ="S_{2} (q#bar{q})"; }
  else if(filename.Index("H0Spin2gg") !=kNPOS) {  labeldescst+="vs S_{2} (gg)";h0label ="S_{2} (gg)"; } 
  if(filename.Index("3300") !=kNPOS)  { labeldescst+=" (M = 3 TeV/c^{2})}";   }
  if(filename.Index("4400") !=kNPOS)  { labeldescst+=" (M = 4 TeV/c^{2})}";   }
  
  
  if(filename.Index("284evts") !=kNPOS ||filename.Index("30evts") !=kNPOS ||filename.Index("393evts") !=kNPOS ||filename.Index("35evts") !=kNPOS )  labellumi+=", 3000 fb^{-1}"; 
  else if(filename.Index("142evts") !=kNPOS ||filename.Index("15evts") !=kNPOS ||filename.Index("196evts") !=kNPOS ||filename.Index("18evts") !=kNPOS )  labellumi+=", 300 fb^{-1}";  
  else if(filename.Index("28evts") !=kNPOS ||filename.Index("3evts") !=kNPOS ||filename.Index("39evts") !=kNPOS ||filename.Index("4evts") !=kNPOS )  labellumi+=", 100 fb^{-1}";  

  labelmcordata->SetLabel(labelname+ labellumi );
  if(h0name.Index("2D") !=kNPOS)  labeldescst+="{2D likelihood}";  
  if(h0name.Index("1D") !=kNPOS) labeldescst +="{1D likelihood}";  
  TPaveLabel *labeldesc = new TPaveLabel(0.59,0.75,0.85,0.95,labeldescst,"brNDC");
  labeldesc->SetFillColor(0);
  labeldesc->SetFillStyle(0);
  labeldesc->SetBorderSize(0);
  labeldesc->SetTextColor(kGreen+3); 
  labeldesc->SetTextSize(0.2); 
  labeldesc->SetTextAlign(12);
  
  
  std::cout.precision(5); 
  if(labellumi.Index("100 fb") !=kNPOS&&labellumi.Index("3000 fb") ==kNPOS  ) cout << stringinfo<< " & " ;

  TFile * myfile = new TFile(filename,"open");
  TH1F * histoh0 = (TH1F*)myfile->Get(h0name); 
  TH1F * histoh1 = (TH1F*)myfile->Get(h1name);
 
  
  //cout << " nb of evts: " << histoh0->GetEntries()<< " &";
 // Compute expected significance under H0 
//   for(int i =0; i<= histoh0->GetNbinsX(); i++){
//     float signif_CLS = 0; 
//     if(histoh0->Integral(0,i)/histoh0->Integral(0,histoh0->GetNbinsX()+1)>=0.5  ){ signif_CLS  = histoh1->Integral(i+1,histoh0->GetNbinsX()+1)/ histoh0->Integral(i+1,histoh0->GetNbinsX()+1);
//       //cout << "exp. separation (H0 true) , sigmas : "<< signif_CLS << " , " << sqrt(2)*TMath::ErfInverse(1-signif_CLS) << endl;
//       cout <<  sqrt(2)*TMath::ErfInverse(1-signif_CLS) << " " ;
//       break; }
//   }
  
  // Same under H1 
  for(int i =histoh1->GetNbinsX()+1; i>=0; i--){
    float signif_CLS = 0; 
    if(histoh1->Integral(i,histoh1->GetNbinsX()+1)/histoh1->Integral(0,histoh0->GetNbinsX()+1)>=0.5  ) {
      
      signif_CLS  = histoh0->Integral(0,i)/ histoh0->Integral(i+1,histoh0->GetNbinsX()+1);
      if(signif_CLS  ==0)signif_CLS  = 1/ histoh0->Integral(i+1,histoh0->GetNbinsX()+1);
      //cout << "exp. separation (H1 true) , sigmas : "<< signif_CLS << " , " << sqrt(2)*TMath::ErfInverse(1-signif_CLS) << endl;
      cout  <<sqrt(2)*TMath::ErfInverse(1-signif_CLS) ; 
      if(labellumi.Index("3000 fb") !=kNPOS )cout << "\\\\"<< endl ;
      else cout << " & " ;
      //      if(i==0||i==1)cout << "end of histo reached" << endl;
      //if(histoh0->Integral(0,i)==0) cout <<"integral is zero" << endl;
   break; }
  }
  //cout << endl;
  TCanvas * c1 = new TCanvas(histoh0->GetName()+labellumi,""); 
  histoh0->Rebin(10);
  histoh1->Rebin(10);
  histoh1->SetLineColor(kRed);
  TLegend * mylegend = new TLegend(0.15,0.4,0.4,0.6);
  
  mylegend->SetTextFont(62);
  mylegend->SetTextSize(0.04);
  mylegend->SetFillStyle(0);
  mylegend->SetBorderSize(0);
  mylegend->AddEntry(histoh0,h0label,"lp"); 
  mylegend->AddEntry(histoh1,h1label,"lp"); 
  histoh0->GetXaxis()->SetTitle("ln L_{0}/L_{1}"); 
  histoh0->GetYaxis()->SetTitle("Normalized nb of entries"); 
  histoh0->GetYaxis()->SetTitleOffset(1.2);

//   if(filename.Index("20evts") !=kNPOS) histoh0->GetYaxis()->SetRangeUser(0,0.05*histoh0->GetEntries() );
//   if(filename.Index("10evts") !=kNPOS) histoh0->GetYaxis()->SetRangeUser(0,0.1*histoh0->GetEntries() );
  histoh0->DrawNormalized(); 
  histoh1->DrawNormalized("SAME"); 
 
  labelmcordata->Draw("same");
  labeldesc->Draw("same");mylegend->Draw("same");c1->Update();
  //  gPad->Redraw();
  outfile->cd(); 
  TString canvasname = labeldescst+labellumi;
  canvasname.ReplaceAll("^","");  canvasname.ReplaceAll("=",""); canvasname.ReplaceAll("/",""); 
   canvasname.ReplaceAll("(",""); canvasname.ReplaceAll(")","");
  canvasname.ReplaceAll(" ","");  canvasname.ReplaceAll("{",""); canvasname.ReplaceAll("}","");  canvasname.ReplaceAll("#",""); 
  canvasname.ReplaceAll(",","_");canvasname.ReplaceAll("splitline",""); canvasname.ReplaceAll("'","prime"); canvasname.ReplaceAll("-1","");canvasname.ReplaceAll("-","");
  c1->SetName(canvasname);
  if(savepdf) c1->SaveAs("plotsfinal/"+canvasname+".pdf"); 
  c1->Write(); 
  delete c1; 
  //cout << endl;
}















Drawsomehistos(TString filename, TString h2dname, TString stringinfo){

  labelmcordata->SetFillColor(0);
  labelmcordata->SetFillStyle(0);
  labelmcordata->SetBorderSize(0);
  labelmcordata->SetTextSize(0.3); 
  labelmcordata->SetTextAlign(12);

  TString labeldescst =""; 
  if(filename.Index("ZprimePSI") !=kNPOS&&h2dname.Index("h1") !=kNPOS)  labeldescst+="Z'_{#Psi} ";  
  if(filename.Index("RSGrav") !=kNPOS&&h2dname.Index("h1") !=kNPOS)  labeldescst+="G_{RS} ";  
  if(filename.Index("Spin0") !=kNPOS&&h2dname.Index("h1") !=kNPOS)  labeldescst+="S_{0} ";  
  if(h2dname.Index("h0") !=kNPOS)  labeldescst+="Drell-Yan ";  
  if(h2dname.Index("forgen") !=kNPOS)  labeldescst+="(generation pdf)";  
  if(h2dname.Index("forlikelihood") !=kNPOS)  labeldescst+="(likelihood calculation pdf)";  
 
  
  TPaveLabel *labeldesc = new TPaveLabel(0.3,0.1,0.55,0.3,labeldescst,"brNDC");
  labeldesc->SetFillColor(0);
  labeldesc->SetFillStyle(0);
  labeldesc->SetBorderSize(0);
  labeldesc->SetTextColor(kBlack); 
  labeldesc->SetTextSize(0.2); 
  labeldesc->SetTextAlign(12);
 
  TFile * myfile = new TFile(filename,"open");
  TH2F * histoh2d = (TH2F*)myfile->Get(h2dname); 
  TCanvas * c1 = new TCanvas("c1",""); c1->SetGridx(false);  c1->SetGridy(false); 
  histoh2d->GetXaxis()->SetTitle("|y_{e^{+}e^{-}}|"); 
  histoh2d->GetYaxis()->SetTitle("cos#theta_{CS_{meas}}"); 
  histoh2d->GetYaxis()->SetTitleOffset(1.2); 
  //histoh2d->GetZaxis()->SetTitle("Normalized nb of evts");
  //histoh2d->GetZaxis()->SetRangeUser(0,0.02);
  histoh2d->DrawNormalized("ZCOL"); 
  //histoh2d->SetMaximum(0.02);histoh2d->SetMinimum(0.000);
  //histoh2d->GetZaxis()->SetRangeUser(0,0.02);
  labelmcordata->Draw("same");
  labeldesc->Draw("same");  
  outfile->cd(); 
  TString canvasname = "plot2d_"+labeldescst;
  canvasname.ReplaceAll(" ","");  canvasname.ReplaceAll("{",""); canvasname.ReplaceAll("}","");  canvasname.ReplaceAll("#","");
  canvasname.ReplaceAll("(","");  canvasname.ReplaceAll(")","");
  canvasname.ReplaceAll(",","_");canvasname.ReplaceAll("splitline","");canvasname.ReplaceAll("'","prime"); canvasname.ReplaceAll("-","");
  c1->SetName(canvasname);c1->SaveAs("plotsfinal/"+canvasname+".pdf"); 
  c1->Write(); 
  delete c1; 


  TString hist1dname="";
  if(h2dname.Index("forgenh0") !=kNPOS)hist1dname = "histcosthetaCSforgenh0costhetaCS";
  if(h2dname.Index("forgenh1") !=kNPOS)hist1dname = "histcosthetaCSforgenh1costhetaCS";
  if(h2dname.Index("forlikelihoodh0") !=kNPOS) hist1dname = "histcosthetaCSforlikelihoodh0costhetaCS";
  if(h2dname.Index("forlikelihoodh1") !=kNPOS) hist1dname = "histcosthetaCSforlikelihoodh1costhetaCS"; 
  myfile->cd(); 
  TCanvas * c2 = new TCanvas("c2",""); c2->SetGridx(false);  c2->SetGridy(false); 
  TH1F * histoh1d = (TH1F*)myfile->Get(hist1dname); 
  histoh1d->GetXaxis()->SetTitle("cos#theta_{CS}_{meas}"); 
  histoh1d->GetYaxis()->SetTitleOffset(1.2); 
  histoh1d->GetYaxis()->SetTitle("Normalized nb of evts");
  histoh1d->DrawNormalized("ZCOL"); 
  
  labelmcordata->Draw("same");
  labeldesc->Draw("same");  
  outfile->cd(); 
  canvasname = "plot1d_"+labeldescst;
  canvasname.ReplaceAll(" ","");  canvasname.ReplaceAll("{",""); canvasname.ReplaceAll("}","");  canvasname.ReplaceAll("#",""); 
  canvasname.ReplaceAll(",","_");canvasname.ReplaceAll("splitline","");canvasname.ReplaceAll("'","prime"); canvasname.ReplaceAll("-","");
  canvasname.ReplaceAll("(","");  canvasname.ReplaceAll(")","");
  c2->SetName(canvasname); c2->SaveAs("plotsfinal/"+canvasname+".pdf"); 
  c2->Write(); 
  delete c2; 
}



DrawsomehistoswithPE(TString filename, TString h2dname, TString stringinfo)
{

  labelmcordata->SetFillColor(0);
  labelmcordata->SetFillStyle(0);
  labelmcordata->SetBorderSize(0);
  labelmcordata->SetTextSize(0.3); 
  labelmcordata->SetTextAlign(12);

  TString labeldescst =""; 
  if(filename.Index("ZprimePSI") !=kNPOS&&h2dname.Index("h1") !=kNPOS)  labeldescst+="Z'_{#Psi} ";  
  if(filename.Index("RSGrav") !=kNPOS&&h2dname.Index("h1") !=kNPOS)  labeldescst+="G_{RS} ";  
  if(filename.Index("Spin0") !=kNPOS&&h2dname.Index("h1") !=kNPOS)  labeldescst+="S_{0} ";  
  if(h2dname.Index("h0") !=kNPOS)  labeldescst+="Drell-Yan ";  
  if(h2dname.Index("forgen") !=kNPOS)  labeldescst+="(generation pdf)";  
  if(h2dname.Index("forlikelihood") !=kNPOS)  labeldescst+="(likelihood calculation pdf)";  
  
  TString labelPEstr =""; 
  if(filename.Index("20evts") !=kNPOS)labelPEstr+="PE: 20 evts"; 
  if(filename.Index("10evts") !=kNPOS)labelPEstr+="PE: 10 evts"; 
  TPaveLabel *labelPE = new TPaveLabel(0.2,0.17,0.6,0.37,"","brNDC");
  labelPE->SetFillColor(0);
  labelPE->SetFillStyle(0);
  labelPE->SetBorderSize(0);
  labelPE->SetTextColor(kBlack); 
  labelPE->SetTextSize(0.2); 
  labelPE->SetTextAlign(12);
  TPaveLabel *labeldesc = new TPaveLabel(0.3,0.1,0.55,0.3,labeldescst,"brNDC");
  labeldesc->SetFillColor(0);
  labeldesc->SetFillStyle(0);
  labeldesc->SetBorderSize(0);
  labeldesc->SetTextColor(kBlack); 
  labeldesc->SetTextSize(0.2); 
  labeldesc->SetTextAlign(12);
 
  TFile * myfile = new TFile(filename,"open");
  TH2F * histoh2d = (TH2F*)myfile->Get(h2dname); 
  TCanvas * c1 = new TCanvas("c1",""); c1->SetGridx(false);  c1->SetGridy(false); 
  histoh2d->GetXaxis()->SetTitle("|y_{e^{+}e^{-}}|"); 
  histoh2d->GetYaxis()->SetTitle("cos#theta_{CS}_{meas}"); 
  histoh2d->GetYaxis()->SetTitleOffset(1.2); 
  //histoh2d->GetZaxis()->SetTitle("Normalized nb of evts");
  //histoh2d->GetZaxis()->SetRangeUser(0,0.02);
  histoh2d->DrawNormalized("ZCOL"); 
  //histoh2d->SetMaximum(0.02);histoh2d->SetMinimum(0.000);
  //histoh2d->GetZaxis()->SetRangeUser(0,0.02);
  TGraph * mygph_h0 =  (TGraph*)myfile->Get("gph_h0_2d_2"); 
  labelPE->SetLabel(labelPEstr+" according to background");

  mygph_h0->SetMarkerStyle(20); 
  mygph_h0->SetMarkerSize(1.1); 
  mygph_h0->Draw("SAME,P");
  labelmcordata->Draw("same");labelPE->Draw("same");
  labeldesc->Draw("same");  
  outfile->cd(); 
  TString canvasname = "plot2dwithPEh0_"+labeldescst;
  canvasname.ReplaceAll(" ","");  canvasname.ReplaceAll("{",""); canvasname.ReplaceAll("}","");  canvasname.ReplaceAll("#",""); 
  canvasname.ReplaceAll(",","_");canvasname.ReplaceAll("splitline","");canvasname.ReplaceAll("'","prime"); canvasname.ReplaceAll("-","");
  canvasname.ReplaceAll("(","");  canvasname.ReplaceAll(")","");
  c1->SetName(canvasname);
  c1->Write();  c1->SaveAs("plotsfinal/"+canvasname+".pdf"); 
  delete c1; 
  delete mygph_h0; 


  TCanvas * c1 = new TCanvas("c1",""); c1->SetGridx(false);  c1->SetGridy(false); 
  histoh2d->DrawNormalized("ZCOL"); 
  TGraph * mygph_h1 =  (TGraph*)myfile->Get("gph_h1_2d_2"); 
  labelPE->SetLabel(labelPEstr+" according to signal");
  mygph_h1->SetMarkerStyle(20); 
  mygph_h1->SetMarkerSize(1.1); 
  
  mygph_h1->Draw("SAME,P");
  labelmcordata->Draw("same");labelPE->Draw("same");
  labeldesc->Draw("same");  
  outfile->cd(); 
  canvasname = "plot2dwithPEh1_"+labeldescst;
  canvasname.ReplaceAll(" ","");  canvasname.ReplaceAll("{",""); canvasname.ReplaceAll("}","");  canvasname.ReplaceAll("#",""); 
  canvasname.ReplaceAll(",","_");canvasname.ReplaceAll("splitline","");canvasname.ReplaceAll("'","prime"); canvasname.ReplaceAll("-","");
  canvasname.ReplaceAll("(","");  canvasname.ReplaceAll(")","");
  c1->SetName(canvasname);
  c1->Write();  c1->SaveAs("plotsfinal/"+canvasname+".pdf"); 
  delete c1; 
  delete mygph_h1; 

  TString hist1dname="";
  if(h2dname.Index("forgenh0") !=kNPOS)hist1dname = "histcosthetaCSforgenh0costhetaCS";
  if(h2dname.Index("forgenh1") !=kNPOS)hist1dname = "histcosthetaCSforgenh1costhetaCS";
  if(h2dname.Index("forlikelihoodh0") !=kNPOS) hist1dname = "histcosthetaCSforlikelihoodh0costhetaCS";
  if(h2dname.Index("forlikelihoodh1") !=kNPOS) hist1dname = "histcosthetaCSforlikelihoodh1costhetaCS"; 
  myfile->cd(); 
  TCanvas * c2 = new TCanvas("c2",""); c2->SetGridx(false);  c2->SetGridy(false); 
  TH1F * histoh1d = (TH1F*)myfile->Get(hist1dname); 
  histoh1d->GetXaxis()->SetTitle("cos#theta_{CS_{meas}}"); 
  histoh1d->GetYaxis()->SetTitleOffset(1.2); 
  histoh1d->GetYaxis()->SetTitle("Normalized nb of evts");
  histoh1d->DrawNormalized("ZCOL"); 
  TGraph * mygph_h0 =  (TGraph*)myfile->Get("gph_h0_2d_2");
  int nbpoints = mygph_h0->GetN();
  double * myvecx_h0 = mygph_h0->GetY();
  double * myvecy = mygph_h0->GetX();
  for(int j =0;j< nbpoints;j++) {myvecy[j]=0;} 
  labelPE->SetLabel(labelPEstr+" according to background");
  TGraph * mygph1d_h0 = new TGraph(nbpoints, myvecx_h0,myvecy); 
  mygph1d_h0->SetMarkerStyle(22); 
  mygph1d_h0->SetMarkerSize(1.1); 
  mygph1d_h0->Draw("SAME,P");
  labelmcordata->Draw("same");labelPE->Draw("same");
  labeldesc->Draw("same");  
  outfile->cd(); 
  canvasname = "plot1dwithPEh0_"+labeldescst;
  canvasname.ReplaceAll(" ","");  canvasname.ReplaceAll("{",""); canvasname.ReplaceAll("}","");  canvasname.ReplaceAll("#",""); 
  canvasname.ReplaceAll(",","_");canvasname.ReplaceAll("splitline","");canvasname.ReplaceAll("'","prime"); canvasname.ReplaceAll("-","");
  canvasname.ReplaceAll("(","");  canvasname.ReplaceAll(")","");
  c2->SetName(canvasname);
  c2->Write(); c2->SaveAs("plotsfinal/"+canvasname+".pdf");  
  delete c2; delete mygph_h0;  delete mygph1d_h0;  //delete  myvecx_h0 ; delete myvecy;


 TCanvas * c2 = new TCanvas("c2",""); c2->SetGridx(false);  c2->SetGridy(false); 
  TH1F * histoh1d = (TH1F*)myfile->Get(hist1dname); 
  histoh1d->GetXaxis()->SetTitle("cos#theta_{CS_{meas}}"); 
  histoh1d->GetYaxis()->SetTitleOffset(1.2); 
  histoh1d->GetYaxis()->SetTitle("Normalized nb of evts");
  histoh1d->DrawNormalized("ZCOL"); 
 
  TGraph * mygph_h1 =  (TGraph*)myfile->Get("gph_h1_2d_2");
  nbpoints = mygph_h1->GetN();
  double * myvecx_h1 = mygph_h1->GetY();
  myvecy = mygph_h1->GetX();
  for(int j =0;j< nbpoints;j++) {myvecy[j]=0;} 
  labelPE->SetLabel(labelPEstr+" according to signal");
  TGraph * mygph1d_h1 = new TGraph(nbpoints, myvecx_h1,myvecy); 
  mygph1d_h1->SetMarkerStyle(22); 
  mygph1d_h1->SetMarkerSize(1.1); 
  mygph1d_h1->Draw("SAME,P");

  labelmcordata->Draw("same");labelPE->Draw("same");
  labeldesc->Draw("same");  
  outfile->cd(); 
  canvasname = "plot1dwithPEh1_"+labeldescst;
  canvasname.ReplaceAll(" ","");  canvasname.ReplaceAll("{",""); canvasname.ReplaceAll("}","");  canvasname.ReplaceAll("#",""); 
  canvasname.ReplaceAll(",","_");canvasname.ReplaceAll("splitline","");canvasname.ReplaceAll("'","prime"); 
  canvasname.ReplaceAll("-",""); canvasname.ReplaceAll("(","");  canvasname.ReplaceAll(")","");
  c2->SetName(canvasname);
  c2->Write(); 
  c2->SaveAs("plotsfinal/"+canvasname+".pdf"); 
  delete c2; delete mygph_h1;  delete mygph1d_h1;  //delete  myvecx_h1 ; delete myvecy;



}
