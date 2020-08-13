
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
const double eff_ana_low = 0;
const double eff_ana_hgh = 600;
double ratio_low = 0.6;
double ratio_hgh = 1.4;


/*
const double ana_low = 0;
const double ana_hgh = 20;
const double eff_ana_low = 0;
const double eff_ana_hgh = 20;
double ratio_low = 0;
double ratio_hgh = 4;
*/

double func_post_poisson(double *x, double *par)
{
  double result = 0;

  double xcur = x[0];
  
  double val_meas = par[0];
  double val_weight = par[1];
  
  double val_true = x[0]/val_weight; 
  double numerator = exp( val_meas - val_true + val_meas*log(val_true/val_meas) );
  result = numerator;

  if(xcur<ana_low || xcur>ana_hgh ) result = 0; // This line is important!!!
  
  return result;
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
  delete ff;

  return result;
}

/////////////////////////////////////////////////////////////////////////////////////////

TF1 *fX_func = NULL;// prediction
TF1 *fY_func = NULL;// data

double func_Y2X(double *x, double *par)
{
  double result = 0;
  
  double xcur = x[0];
  double par0 = par[0];
  
  if( xcur<ana_low || xcur>ana_hgh ) {
    result = 0;
  }
  else {
    double eval_fX = fX_func->Eval( xcur );
    double eval_fY = fY_func->Eval( par0 * xcur );
    double abs_X = fabs( xcur );
    result = abs_X * eval_fX * eval_fY;
  }
  
  return result;
}

double pdf_func_Y2X(double *x, double *par)
{
  double result = 0;

  double z = x[0];

  TF1 *roofunc_Y2X = new TF1("roofunc_Y2X", func_Y2X, eff_ana_low, eff_ana_hgh, 1);
  roofunc_Y2X->SetParameter(0, z);
  double val_integration = roofunc_Y2X->Integral( eff_ana_low, eff_ana_hgh );
  delete roofunc_Y2X;
  
  result = val_integration;
  
  return result;  
}

/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// MAIN //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


void read_err()
{
  TString roostr = "";
  
  //ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");
  
  ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-4);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////// validation

  /*
  TF1 *func_A = new TF1("func_A", func_uniform, 0,4,0);
  TF1 *func_B = new TF1("func_B", func_uniform, 0,4,0);
  
  TF1Convolution *fco_mAB = new TF1Convolution(func_A, func_B, 0, 4, 1);
  TF1 *func_mAB = new TF1("func_mAB", fco_mAB, 0, 2, 0);
  func_mAB->SetNpx(10000);
  
  roostr = "canv_func_mAB";
  TCanvas *canv_func_mAB = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_func_mAB, 0.15, 0.2,0.1,0.15);
  func_mAB->Draw();

  ///////////////////////////////////////////
  
  TF1 *func_C = new TF1("func_C", func_exp, 0,10,0);
  TF1 *func_D = new TF1("func_D", func_exp, 0,10,0);
  
  TF1Convolution *fco_mCD = new TF1Convolution(func_C, func_D, 0, 10, 1);
  TF1 *func_mCD = new TF1("func_mCD", fco_mCD, 0, 6, 0);
  func_mCD->SetNpx(10000);
 
  roostr = "canv_func_mCD";
  TCanvas *canv_func_mCD = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_func_mCD, 0.15, 0.2,0.1,0.15);
  func_mCD->Draw();
  */

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////// explore

  
  double E_obs = 400;
  double E_w   = 0.25;

  double F_obs = 200;
  double F_w   = 1;

  double G_obs = 300;
  double G_w   = 0.5;

  double OBS_obs = 460;
  double OBS_w   = 1;

  double central_ratio = 460./450;

  

  /*
  double E_obs = 40;
  double E_w   = 0.025;

  double F_obs = 10;
  double F_w   = 0.1;

  double G_obs = 1;
  double G_w   = 1;

  double OBS_obs = 3;
  double OBS_w   = 1;

  double central_ratio = 1.;
  */
  
  
  TF1 *func_E = new TF1("func_E", pdf_func_post_poisson, eff_ana_low, eff_ana_hgh, 2);
  func_E->SetParameters(E_obs, E_w);// Nobs, weight
  func_E->SetNpx(60000);
  
  // TF1 *func_gaus = new TF1("func_gaus", "gaus(0)", eff_ana_low, eff_ana_hgh);
  // func_gaus->SetParameters(0.08, 300, 20);
  // func_gaus->SetNpx(60000);  
  // TF1Convolution *fco_mEgaus = new TF1Convolution(func_E, func_gaus, eff_ana_low, eff_ana_hgh, true);
  // fco_mEgaus->SetNofPointsFFT(10000);  
  // TF1 *func_mEgaus = new TF1("func_mEgaus", fco_mEgaus, eff_ana_low, eff_ana_hgh, fco_mEgaus->GetNpar());
  // func_mEgaus->SetParameters(E_obs, E_w, 0.08, 300, 20);
  // func_mEgaus->SetNpx(60000);

  TF1 *func_F = new TF1("func_F", pdf_func_post_poisson, eff_ana_low, eff_ana_hgh, 2);
  func_F->SetParameters(F_obs, F_w);
  func_F->SetNpx(60000);  
  TF1Convolution *fco_mEF = new TF1Convolution(func_E, func_F, eff_ana_low, eff_ana_hgh, true);
  fco_mEF->SetNofPointsFFT(10000);  
  TF1 *func_mEF = new TF1("func_mEF", fco_mEF, eff_ana_low, eff_ana_hgh, fco_mEF->GetNpar());
  func_mEF->SetParameters(E_obs, E_w, F_obs, F_w);
  func_mEF->SetNpx(60000);
  
  TF1 *func_G = new TF1("func_G", pdf_func_post_poisson, eff_ana_low, eff_ana_hgh, 2);
  func_G->SetParameters(G_obs, G_w);
  func_G->SetNpx(60000);  
  TF1Convolution *fco_mEFG = new TF1Convolution(func_mEF, func_G, eff_ana_low, eff_ana_hgh, true);
  fco_mEFG->SetNofPointsFFT(10000);  
  TF1 *func_mEFG = new TF1("func_mEFG", fco_mEFG, eff_ana_low, eff_ana_hgh, fco_mEFG->GetNpar());
  func_mEFG->SetParameters(E_obs, E_w, F_obs, F_w, G_obs, G_w);
  func_mEFG->SetNpx(60000);
  
  TF1 *func_obs = new TF1("func_obs", pdf_func_post_poisson, eff_ana_low, eff_ana_hgh, 2);
  func_obs->SetParameters(OBS_obs, OBS_w);// Nobs, weight
  func_obs->SetNpx(60000);
  
  /////////////////////////// plotting
  
  roostr = "canv_obj";
  TCanvas *canv_obj = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_obj, 0.15, 0.2,0.1,0.15);

  func_E->Draw();
  func_E->SetLineColor(kBlue);
    
  // func_gaus->Draw("same");
  // func_gaus->SetLineColor(kBlack);  
  // func_mEgaus->Draw("same");
  // func_mEgaus->SetLineColor(kBlack);
  
  func_F->Draw("same");
  func_F->SetLineColor(kRed);
 
  func_G->Draw("same");
  func_G->SetLineColor(kGreen+1);
 
  func_mEFG->Draw("same");
  func_mEFG->SetLineColor(kGray+2);
 
  func_obs->Draw("same");
  func_obs->SetLineColor(kCyan);
 
  /////////////////////////////////////////////////////////////////////////////////////////
  
  // TF1 *f1_test = new TF1("f1_test","x", 0, 1);
  // cout<<endl<<" ---> f1 test: "<<f1_test->Integral(0,1)<<endl<<endl;

  fX_func = func_mEFG;
  fY_func = func_obs; 

  TF1 *roopdf_func_Y2X = new TF1("pdf_func_Y2X", pdf_func_Y2X, ratio_low, ratio_hgh, 0);
  roopdf_func_Y2X->SetNpx(400);
  
  // roostr = "canv_ratio";
  // TCanvas *canv_ratio = new TCanvas(roostr, roostr, 900, 650);
  // func_canv_margin(canv_ratio, 0.15, 0.2,0.1,0.15);  
  // roopdf_func_Y2X->Draw();
  
  roopdf_func_Y2X->SaveAs("file_roopdf_func_Y2X.root");
  //roopdf_func_Y2X->SaveAs("file_roopdf_func_Y2X.png");

  cout<<endl;
  cout<<" ------------------------------------------------------------------------------- "<<endl;
  cout<<" ------------------------------------------------------------------------------- "<<endl;
  cout<<" ------------------------------------------------------------------------------- "<<endl;
  cout<<" ------------------------------------------------------------------------------- "<<endl;
  cout<<endl;
  
  cout<<" 1sigma "<<1-TMath::Prob(1, 1)<<endl;
  double value_1sigma = 1-TMath::Prob(1, 1);
  
  TFile *file_obj = new TFile("file_roopdf_func_Y2X.root", "read");
  TF1 *fobj = (TF1*)file_obj->Get("pdf_func_Y2X");
  fobj->SetName("f_obj");
  
  
  int line = 0;

    
  double val_typeA = 0;
  double val_typeB = central_ratio;  

  ////////////////////////
  
  val_typeA = ratio_low;
  val_typeB = central_ratio;  
  while( fabs(val_typeA-val_typeB)>1e-4 ) {
    line++;   
    double val_mid = (val_typeA+val_typeB)/2;
    double y_mid = fobj->Integral(val_mid, central_ratio);
    
    if( y_mid > value_1sigma/2 ) {
      val_typeA = val_mid;
    }
    else {
      val_typeB = val_mid;
    }
    
    cout<<TString::Format(" -------> %3d, %10.6f %10.6f", line, val_typeA, val_typeB)<<endl;
    
    if( line>15 ) {
      break;
    }
    
  }
  double low_edge = (val_typeA+val_typeB)/2;
  line = 0;

  ////////////////////////
  
  val_typeA = central_ratio;
  val_typeB = ratio_hgh;  
  while( fabs(val_typeA-val_typeB)>1e-4 ) {
    line++;   
    double val_mid = (val_typeA+val_typeB)/2;
    double y_mid = fobj->Integral(central_ratio, val_mid);
    
    if( y_mid > value_1sigma/2 ) {
      val_typeB = val_mid;
    }
    else {
      val_typeA = val_mid;
    }
    
    cout<<TString::Format(" -------> %3d, %10.6f %10.6f", line, val_typeA, val_typeB)<<endl;
    
    if( line>15 ) {
      break;
    }
    
  }
  double hgh_edge = (val_typeA+val_typeB)/2;
  line = 0;

  cout<<endl<<TString::Format(" ---> %6.3f, 1sigma range %6.3f %6.3f", central_ratio, low_edge, hgh_edge)<<endl;

  
}
