
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TObjArray.h"
#include "TBranch.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TF1.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TPRegexp.h"

// RooFit includes
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"
#include <Math/MinimizerOptions.h>
#include "RooMinimizer.h"
#include "RooMinimizerFcn.h"

// standard includes
#include <iostream>




RooAbsReal *nll;
RooWorkspace *w;
RooStats::ModelConfig *mc_s;
RooDataSet *data;
//RooArgSet myArgs;

std::vector<TString> nuisListBestNLL;
std::vector<double> nuisValuesBestNLL;

// For LH Plots, n-sigma along x axis
// const int npoints = 1000;
const int npoints = 20;//10;
const int nsigma  = 4;
// Label size for Pull Summary
const double pullLabelSize = 0.028;
const int maxPullsPerPlot = 10;


TGraph *graphLH(TString nuisname, TString whichfit){

  //ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.000001);

  w->loadSnapshot(whichfit); // SetTo BestFit values as start

  // Get The parameter we want 
  RooRealVar *nuis =(RooRealVar*) w->var(nuisname);

  // Adjust the range on the signal strength
  RooRealVar *r = (RooRealVar *) w->var("r");
  r->setRange(-100,100);

  RooRealVar *jes = (RooRealVar *) w->var("CMS_scale_j");
  //jes->setRange(-2,-0.1);

  const double thres = 0.00001;
  double errLow = fabs(nuis->getAsymErrorLo());
  double errHigh = fabs(nuis->getAsymErrorHi());
  errLow = ( errLow<thres*nuis->getError() ) ? nuis->getError() : errLow;
  errHigh = ( errHigh<thres*nuis->getError() ) ? nuis->getError() : errHigh;

  nuis->setConstant(true);

  double bf = nuis->getVal();
  double nll_0=nll->getVal(); 

  //errLow = 1.; errHigh = 1.;
  //errLow = 0.75; errHigh = 0.75;
  //bf = 0;


  std::cout << " !! STARTING SCAN !! " << std::endl;

  printf("\t nuisance = %s,\t best-fit = %.3f, err = %.3f, errLow = %.3f, errHigh = %.3f,\t nll_0 = %.2f \n",
	 nuisname.Data(), bf, nuis->getError(), errLow, errHigh, nll_0 );


  //Configure the mimimizer for minimizing at each point along the scan.
  RooMinimizer m(*nll) ;
  m.setPrintEvalErrors(-1) ;
  m.setPrintLevel(-1) ;
  m.setVerbose(false) ;
  m.setStrategy(1) ;

  double minDeltaNLL = 99.;
  double nuis_minDeltaNLL = -99;

  //TGraph *gr = new TGraph(2*npoints+1);

  std::vector<double> xi, yi;
  xi.clear(); yi.clear();

  for (int i=-1*npoints;i<=npoints;i++){

    w->loadSnapshot("prefitparams"); // SetTo BestFit values as start
    nuis->setVal(bf+((i < 0) ? errLow : errHigh)*((float)i)*nsigma/npoints);
    nuis->setConstant(true);
    r->setRange(-100,100);
    //jes->setRange(-2,-0.1);

    //Save one running of the fit!
    m.minimize("Minuit2");
      
    double nll_v = nll->getVal();

    std::cout << "#### i="<<i
	      << ":x = " << nuis->getVal()
	      << ", y = " << nll_v-nll_0
	      << ", r = " << r->getVal()
	      << ", JES = " << jes->getVal()
	      << ", nll_v = " << nll_v
	      << ", nll_0 = " << nll_0
	      <<std::endl;

    if( (nll_v-nll_0)<minDeltaNLL ){
      minDeltaNLL = nll_v-nll_0;
      nuis_minDeltaNLL = nuis->getVal();
    }

    if( nll_v-nll_0 < 10000 ){
      xi.push_back(nuis->getVal());
      yi.push_back(nll_v-nll_0);
      //yi.push_back(jes->getVal());
    }
    else{
      std::cout << "   !!! WARNING !!! NLL = " << nll_v << " at " << nuis->getVal() << " for " <<  nuisname << std::endl;
    }
  }


  int NPoints = int( xi.size() );
  double x[NPoints], y[NPoints];
  for( int iPoint=0; iPoint<NPoints; iPoint++ ){
    x[iPoint] = xi[iPoint];
    y[iPoint] = yi[iPoint];
  }

  TGraph *gr = new TGraph(NPoints, x, y);

  std::cout << "\t value with minimal NLL at " << nuisname << " = " << nuis_minDeltaNLL << " for deltaNLL = " << minDeltaNLL << ",\t best-fit at " << bf << ",\t diff = " << nuis_minDeltaNLL - bf << std::endl;

  std::cout << " Minimize again for best fit " << std::endl;

  nuis->setVal( nuis_minDeltaNLL );
  nuis->setConstant( true );

  m.minimize("Minuit2");

  RooFitResult* res = m.save();

  RooArgSet useMyArgs = res->floatParsFinal();

  nuisListBestNLL.clear();
  nuisValuesBestNLL.clear();
  TIterator *nextArgPost(useMyArgs.createIterator());
  for (TObject *a = nextArgPost->Next(); a != 0; a = nextArgPost->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);      
    TString nuis_name = rrv->GetName();
    if( !nuis_name.Contains("ANNbin") ){
      nuisListBestNLL.push_back(nuis_name);
      nuisValuesBestNLL.push_back(rrv->getVal());
    }
  }

  nuisListBestNLL.push_back(nuisname);
  nuisValuesBestNLL.push_back(nuis_minDeltaNLL);


  gr->SetTitle("");
  gr->GetYaxis()->SetTitle("NLL - obs data");
  gr->GetYaxis()->SetTitleOffset(0.8);
  gr->GetXaxis()->SetTitleOffset(0.8);
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetTitle(nuisname);

  // 	gr->GetYaxis()->SetRangeUser(0.0,3);

  gr->SetLineColor(4);
  gr->SetLineWidth(2);
  gr->SetMarkerStyle(21);
  gr->SetMarkerSize(0.6);
	
  return gr;
	

}


void plotParametersNLL(TString nuisRE, std::string dataFits="mlfit.root", std::string workspace="wsTest.root", TString postfix="", TString histopostfix = "", int jobN=-1){

  // Some Global preferences
  gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so");
  gROOT->SetBatch(true);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  TString imgDir = "Images/Images_2014_01_21_parameterScan"+postfix;

  TString histoFileName = ( !nuisRE.EqualTo(".*") ) ? "HistoFiles/nll"+postfix+"_"+nuisRE+histopostfix+".root" : "HistoFiles/nll"+postfix+"_all"+histopostfix+".root";


  struct stat st;
  if( stat(imgDir.Data(),&st) != 0 )  mkdir(imgDir.Data(),0777);

  TFile *fd_=0;
  TFile *fw_=0;

  std::cout << "Getting fit to data from "<< dataFits <<std::endl;
  fd_ = TFile::Open(dataFits.c_str());
  
  // Toys are thrown from best fit to data (background only/mu=0) 
  RooFitResult *preFitFR = (RooFitResult*)fd_->Get("nuisances_prefit_res");  
  RooArgSet prefitargs = preFitFR->floatParsFinal();

  RooFitResult *bestfit=(RooFitResult*)fd_->Get("fit_b");
  RooArgSet fitargs = bestfit->floatParsFinal();
  
  RooFitResult *bestfit_s=(RooFitResult*)fd_->Get("fit_s");
  RooArgSet fitargs_s = bestfit_s->floatParsFinal();
//   // These are essentially the nuisances in the card (note, from toys file, they will be randomized)
//   // so need to use the data fit.
//   RooArgSet *prefitargs = (RooArgSet*)fd_->Get("nuisances_prefit");
  
  std::cout << "Getting the workspace from "<< workspace << std::endl;
  fw_ =  TFile::Open(workspace.c_str());
  w   = (RooWorkspace*) fw_->Get("w");
  data = (RooDataSet*) w->data("data_obs");
  mc_s = (RooStats::ModelConfig*)w->genobj("ModelConfig");
  std::cout << "make nll"<<std::endl;
  nll = mc_s->GetPdf()->createNLL(*data,
                                  RooFit::Constrain(*mc_s->GetNuisanceParameters()),
                                  RooFit::Extended(mc_s->GetPdf()->canBeExtended())
                                  );
  
  
//   // grab r (mu) from workspace to set to 0 for bonly fit since it wasnt floating 
//   RooRealVar *r = w->var("r"); r->setVal(0);fitargs.add(*r);
  
  w->saveSnapshot("prefitparams",prefitargs,true);
  w->saveSnapshot("bestfitparams",fitargs,true);	
  w->saveSnapshot("bestfitparams_sb",fitargs_s,true);	
  std::cout << "Workspace OK!"<<std::endl;
  


  TFile histofile(histoFileName,"recreate");
  histofile.cd();

  TCanvas *c = new TCanvas("c","",480,400);
  TGraph *gr = NULL;
  
  RooArgSet myArgs = bestfit_s->floatParsFinal();

  //This is a PERL regex.  The "^" and "$" require the whole string to match, not just part.
  //If you want to match anything that contains put ".*" before and after your RE in the function arguement.
  TPRegexp re("^"+nuisRE+"$");

  std::vector<TString> nuisListMaxFit;
  std::vector<double> nuisValuesMaxFit;
  std::vector<double> nuisValuesMaxFitErr;

  nuisListMaxFit.clear();
  nuisValuesMaxFit.clear();
  nuisValuesMaxFitErr.clear();

  int nNuis = 0;
  if( true ){
    TIterator *nextArgPre(myArgs.createIterator());
    for (TObject *a1 = nextArgPre->Next(); a1 != 0; a1 = nextArgPre->Next()) {
      RooRealVar *rrv = dynamic_cast<RooRealVar *>(a1);      
      TString nuis_name = rrv->GetName();
      if( !nuis_name.Contains("ANNbin") ){
	nNuis++;
	nuisListMaxFit.push_back(nuis_name);
	nuisValuesMaxFit.push_back(rrv->getVal());

	double thres = 0.00001;
	Double_t errValUp = fabs(rrv->getAsymErrorHi());
	Double_t errValDown = fabs(rrv->getAsymErrorLo());
	errValUp = ( errValUp<thres*rrv->getError() ) ? rrv->getError() : errValUp;
	errValDown = ( errValDown<thres*rrv->getError() ) ? rrv->getError() : errValDown;

	double err = 0.5 * ( errValUp + errValDown );
	nuisValuesMaxFitErr.push_back(err);
      }
    }
  }

  //nNuis = 0;
  TH1D* h_nuis_diff = new TH1D("h_nuis_diff",";nuisance;NLL Scan Best - Max Likelihood Fit Best", nNuis, 0, nNuis );
  TH1D* h_nuis_frac_diff = new TH1D("h_nuis_frac_diff",";nuisance;(NLL Scan Best - Max Likelihood Fit Best)/Max Fit Error", nNuis, 0, nNuis );

  for( int iBin=0; iBin<nNuis; iBin++ ){
    h_nuis_diff->GetXaxis()->SetBinLabel(iBin+1,nuisListMaxFit[iBin]);
    h_nuis_frac_diff->GetXaxis()->SetBinLabel(iBin+1,nuisListMaxFit[iBin]);
  }
  
  int bscounter=0;
  TIterator *nextArg(myArgs.createIterator());
  for (TObject *a = nextArg->Next(); a != 0; a = nextArg->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);      
    //RooRealVar *rrv = a;      
//     std::cout << " name = " << rrv->GetName() << ":\t value = " << rrv->getVal() << ",\t error = " << rrv->getError() << std::endl;
    TString nuis_name = rrv->GetName();
    
    if( nuis_name.Contains("ANNbin") ) continue;
    
    if( re.MatchB(nuis_name) ){
      //if( bscounter>2 ) break;
      bscounter++;
      if( jobN!=(bscounter-1) && jobN>=0 ) continue;

      gr = graphLH(nuis_name,"bestfitparams_sb");
      gr->Draw("ALP");
      
      c->Print(imgDir+"/nll"+postfix+"_"+nuis_name+".png");

      gr->Write("nll"+postfix+"_"+nuis_name);

      TCanvas *c2 = new TCanvas("c2","",960,500);
      //c2->SetBottomMargin(2);
      //c2->SetRightMargin(0.05);

      //std::cout << " Number of nuisances = " << nNuis << std::endl;
      for( int iNuis=0; iNuis<nNuis; iNuis++ ){
	if( nuisListMaxFit[iNuis].Contains("ANNbin") ) continue;

	//std::cout << "\t iNuis = " << iNuis << ",\t name = " << nuisListMaxFit[iNuis] << std::endl;
	for( int jNuis=0; jNuis<nNuis; jNuis++ ){
	  if( nuisListBestNLL[jNuis].Contains("ANNbin") ) continue;

	  //std::cout << "\t\t jNuis = " << jNuis << ",\t name = " << nuisListBestNLL[jNuis] << std::endl;
	
	  if( !nuisListMaxFit[iNuis].EqualTo(nuisListBestNLL[jNuis]) ) continue;

	  double preScan = nuisValuesMaxFit[iNuis];
	  double postScan = nuisValuesBestNLL[jNuis];
	  double preScanErr = nuisValuesMaxFitErr[iNuis];

	  //std::cout << "\t\t\t Found a match for " << nuisListMaxFit[iNuis] << "!  preScan = " << preScan << ",\t postScan = " << postScan << ",\t diff = " << postScan - preScan << std::endl;

	  h_nuis_diff->SetBinContent(iNuis+1,postScan-preScan);
	  h_nuis_frac_diff->SetBinContent(iNuis+1,(postScan-preScan)/preScanErr);

	  h_nuis_diff->Draw();
	  c2->Print(imgDir+"/nuisDiff"+postfix+"_"+nuis_name+".png");

	  h_nuis_frac_diff->Draw();
	  c2->Print(imgDir+"/nuisFracDiff"+postfix+"_"+nuis_name+".png");

	  //if( nuis_name.EqualTo(nuisListMaxFit[iNuis]) ) std::cout << " For nuisance " << nuis_name << " preFit = " << preScan << ",\t postFit = " << postScan << ",\t diff = " << postScan-preScan << std::endl;
	}
      }
      delete c2;
    }
  }

  histofile.Write();
  histofile.Close();

//   for( int iBin=0; iBin<nNuis; iBin++ ){
//     std::cout << " Content for bin " << iBin << "\t" << h_nuis_diff->GetXaxis()->GetBinLabel(iBin+1) << "\t is " << h_nuis_diff->GetBinContent(iBin+1) << std::endl;
//   }

}
