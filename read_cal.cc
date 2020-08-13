
#include<iostream>
#include<fstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<map>
#include<vector>
#include<set>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THStack.h"
#include "TF1.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"

#include "TRandom3.h"
#include "TGaxis.h"
#include "TStyle.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TROOT.h"

#include "TF1Convolution.h"


void func_canv_margin(TCanvas *canv, double left, double right, double top, double bot)
{
  canv->SetLeftMargin(left);
  canv->SetRightMargin(right);
  canv->SetTopMargin(top);
  canv->SetBottomMargin(bot);
}

/////////////////////////////////////////////////////////////////////////////////////////

double func_uniform(double *x, double *par)
{
  double result = 0;

  double xcur = x[0];

  if( xcur>=0 && xcur<=1 ) {
    result = 1;
  }
  
  return result;
}


double func_exp(double *x, double *par)
{
  double result = 0;

  double xcur = x[0];

  if( xcur>0 ) {
    result = exp(-xcur);
  }

  return result;
}
  
/////////////////////////////////////////////////////////////////////////////////////////

const double ana_low = 0;
const double ana_hgh = 1000;

double func_post_poisson(double *x, double *par)
{
  double val_meas = par[0];
  double val_weight = par[1];
  
  double val_true = x[0]/val_weight; 
  double numerator = exp( val_meas - val_true + val_meas*log(val_true/val_meas) );
  return numerator;
}


double pdf_func_post_poisson(double *x, double *par)
{
  double result = 0;
  
  double val_meas = par[0];
  double val_weight = par[1];
  
  double val_true = x[0];
  
  TF1 *ff = new TF1("ff", func_post_poisson, ana_low, ana_hgh, 2);
  ff->SetParameter(0, val_meas);
  ff->SetParameter(1, val_weight);
  double numerator = ff->Eval( val_true );  
  double denominator = ff->Integral(ana_low, ana_hgh);
  result = numerator/denominator;

  if( result<1e-5 ) result = 0;
  
  delete ff;

  return result;
}

/////////////////////////////////////////////////////////////////////////////////////////

double funcTypeA_post_poisson(double *x, double *par)
{
  double val_meas = 400;
  double val_weight = 0.25;
  
  double val_true = x[0]/val_weight; 
  double numerator = exp( val_meas - val_true + val_meas*log(val_true/val_meas) );
  return numerator;
}

double pdf_funcTypeA_post_poisson(double *x, double *par)
{
  double result = 0;
  double val_true = x[0];
  
  TF1 *ff = new TF1("ff", funcTypeA_post_poisson, ana_low, ana_hgh, 0);
  double numerator = ff->Eval( val_true );  
  double denominator = ff->Integral(ana_low, ana_hgh);
  result = numerator/denominator;
  if( result<1e-5 ) result = 0;  
  delete ff;
  return result;
}


double funcTypeB_post_poisson(double *x, double *par)
{
  double val_meas = 200;
  double val_weight = 1;
  
  double val_true = x[0]/val_weight; 
  double numerator = exp( val_meas - val_true + val_meas*log(val_true/val_meas) );
  return numerator;
}

double pdf_funcTypeB_post_poisson(double *x, double *par)
{
  double result = 0;
  double val_true = x[0];
  
  TF1 *ff = new TF1("ff", funcTypeB_post_poisson, ana_low, ana_hgh, 0);
  double numerator = ff->Eval( val_true );  
  double denominator = ff->Integral(ana_low, ana_hgh);
  result = numerator/denominator;
  if( result<1e-5 ) result = 0;  
  delete ff;
  return result;
}


/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// MAIN //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


void read_cal()
{
  TString roostr = "";
  
  //ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");
  
  ///////////////////////////////////////////
  /*
  TF1 *f_A = new TF1("f_A", func_uniform, 0,4,0);
  TF1 *f_B = new TF1("f_B", func_uniform, 0,4,0);
  
  TF1Convolution *f_mAB = new TF1Convolution(f_A, f_B, 0, 4, 1);
  TF1 *ff_mAB = new TF1("ff_mAB", f_mAB, 0, 2, 0);
  
  roostr = "canv_f_mAB";
  TCanvas *canv_f_mAB = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_f_mAB, 0.15, 0.2,0.1,0.15);
  ff_mAB->Draw();

  ///////////////////////////////////////////
  
  TF1 *f_C = new TF1("f_C", func_exp, 0,10,0);
  TF1 *f_D = new TF1("f_D", func_exp, 0,10,0);
  
  TF1Convolution *f_mCD = new TF1Convolution(f_C, f_D, 0, 10, 1);
  TF1 *ff_mCD = new TF1("ff_mCD", f_mCD, 0, 6, 0);
  
  roostr = "canv_f_mCD";
  TCanvas *canv_f_mCD = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_f_mCD, 0.15, 0.2,0.1,0.15);
  ff_mCD->Draw();
  */
  //////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  /*
  TF1 *func_AA = new TF1("func_AA", pdf_func_post_poisson, 0, 600, 2);
  func_AA->SetParameter(0, 400);
  func_AA->SetParameter(1, 0.25);
  func_AA->SetNpx(60000);
  
  TF1 *func_BB = new TF1("func_BB", pdf_func_post_poisson, 0, 600, 2);
  func_BB->SetParameter(0, 200);
  func_BB->SetParameter(1, 1);
  func_BB->SetNpx(60000);
  
  TF1 *func_CC = new TF1("func_CC", pdf_func_post_poisson, 0, 600, 2);
  func_CC->SetParameter(0, 300);
  func_CC->SetParameter(1, 0.5);
  func_CC->SetNpx(60000);

  cout<<endl<<" processing convolution"<<endl<<endl;
  
  TF1Convolution *f_uAB = new TF1Convolution("func_AA", "func_BB", 0.01, 600-0.01, 1);
  f_uAB->SetNofPointsFFT(10000);
  TF1 *func_uAB = new TF1("func_uAB", *f_uAB, 0, 600, f_uAB->GetNpar() );
  func_uAB->SetParameters(400, 0.25, 200, 1);
  
  cout<<func_uAB->Eval(300)<<endl;  
  // TF1Convolution *f_uABC = new TF1Convolution(func_uAB, func_CC, 0, 600, 1);
  // TF1 *func_uABC = new TF1("func_uABC", f_uABC, 0, 600, 0);
  // func_uABC->SetNpx(60000);  
  
  ////////////////////////////////////////////////////////////////
  
  roostr = "canv_func_AA";
  TCanvas *canv_func_AA = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_func_AA, 0.15, 0.2,0.1,0.15);
  
  func_AA->Draw();
  func_AA->SetLineColor(kBlue);

  func_BB->Draw("same");
  func_BB->SetLineColor(kGreen+1);

  func_CC->Draw("same");
  func_CC->SetLineColor(kRed);
  
  // func_uAB->Draw();
  // func_uAB->SetLineColor(kRed);
  
  */
  /////////////////////////////////////////////////////////////////////////////////////
  /*
  Double_t D_t_cm=237.*1e-4;//cm/sqrt(cm)
  Double_t L_drift = 30.;//cm
  Double_t v_drift=7.9*1e3;//cm/ns

  Double_t D_l_cm = 290.*1e-4;//cm/sqrt(cm)
	
  Double_t D_l=D_l_cm*D_l_cm*v_drift/2;//cm2/ns

  Double_t sigma_in=TMath::Sqrt(L_drift*L_drift/D_l);//;D_l*D_l*L_drift
	
  double tmax = 511*80;

	
  TF1* tf1_A = new TF1("tf1_A", "[2]*(x>[0])*TMath::Exp(-1*(x>[0])*(x-[0])*3/[1])*TMath::Power((x>[0])*(x-[0])/[1],3)*TMath::Sin((x>[0])*(x-[0])/[1])", 0, tmax);
  tf1_A->SetParameter(0, 9000.);
  tf1_A->SetParameter(1, 600.);
  tf1_A->SetParameter(2,1);

  TF1* tf1_L = new TF1("tf1_L", "gaus(0)", 0, tmax);
  tf1_L->SetParameter(0, 1);
  tf1_L->SetParameter(1, 9000);
  tf1_L->SetParameter(2, sigma_in);

  TF1Convolution* tf1conv_AL = new TF1Convolution("tf1_A", "tf1_L", 0.01, tmax, true);
  tf1conv_AL->SetNofPointsFFT(10000);
  cout<<" ---> "<< tf1conv_AL->GetNpar()<<endl;
  TF1   *f = new TF1("f",*tf1conv_AL, 0.01, tmax, tf1conv_AL->GetNpar());
  f->SetParameters(9000, 600, 1, 1, 0, sigma_in);

  TCanvas* c_tf = new TCanvas();
  //tf1_A->Draw();
  //tf1_L->Draw("same");
  f->Draw();  
 
  */

  
  
  TF1* tf1_A = new TF1("tf1_A", "exp( [0] - x + [0]*log(x/[0]) ) * (x>0) *(x<600)", 0.01, 600);
  tf1_A->SetParameter(0, 100);

  /*
  TF1* tf1_L = new TF1("tf1_L", "gaus(0)", 0, 600);
  tf1_L->SetParameter(0, 1);
  tf1_L->SetParameter(1, 300);
  tf1_L->SetParameter(2, 50);

  TF1Convolution* tf1conv_AL = new TF1Convolution("tf1_A", "tf1_L", 0.01, 600, true);
  tf1conv_AL->SetNofPointsFFT(10000);
  cout<<" ---> "<< tf1conv_AL->GetNpar()<<endl;
  TF1   *f = new TF1("f",*tf1conv_AL, 0.01, 600, tf1conv_AL->GetNpar());
  f->SetParameters(100, 1, 600, 100);
  */

  TF1* tf1_L = new TF1("tf1_L", "exp( [0] - x + [0]*log(x/[0]) ) * (x>0) *(x<600)", 0.01, 600);
  tf1_L->SetParameter(0, 200);

  TF1Convolution* tf1conv_AL = new TF1Convolution("tf1_A", "tf1_L", 0.01, 600, true);
  tf1conv_AL->SetNofPointsFFT(10000);
  cout<<" ---> "<< tf1conv_AL->GetNpar()<<endl;
  TF1   *f = new TF1("f",*tf1conv_AL, 0.01, 600, tf1conv_AL->GetNpar());
  f->SetParameters(100, 200);
  
  TCanvas* c_tf = new TCanvas();
  //tf1_A->Draw();
  //tf1_L->Draw("same");
  //f->Draw();

  f->Draw();

  //tf1_A->Draw("same");

  /*
  TF1* tf1_A = new TF1("tf1_A", "1",0, 1);
  //  tf1_A->SetParameter(0, 100);

  TF1* tf1_L = new TF1("tf1_L", "1", 0,1);
  //  tf1_L->SetParameter(0, 200);

  TF1Convolution* tf1conv_AL = new TF1Convolution("tf1_A", "tf1_L", 0,3, true);
  tf1conv_AL->SetNofPointsFFT(10000);
  cout<<" ---> "<< tf1conv_AL->GetNpar()<<endl;
  TF1   *f = new TF1("f",*tf1conv_AL, 0,3, tf1conv_AL->GetNpar());
  f->SetParameters(100,200);

  TCanvas* c_tf = new TCanvas();
  //tf1_A->Draw();
  //tf1_L->Draw("same");
  f->Draw();

  //tf1_A->Draw();
  */
  
}
