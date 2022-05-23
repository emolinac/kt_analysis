#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <stdlib.h>
#include <gsl/gsl_integration.h>
#include "TF1.h"
#include "TF2.h"
#include "Math/WrappedTF1.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "TMath.h"
#include "Math/GaussLegendreIntegrator.h"
#include "Math/GSLIntegrator.h"
#include "Math/Integrator.h"
#include "Math/Functor.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TNtuple.h"
#include "TBenchmark.h"
using namespace std;

const double Pi = 3.14159265358979323846;
double kt2max  = 2.; // in GeV^2
double Pt2_MAX = 10.;

double meanPt2ddd(Double_t *x, Double_t *par) //x[0]:u2 , par[0] theta, par[1]:kt2, par[2]:z, par[3]:mkt2, par[4]:mpt2, par[5]=vm
{
  double result;
  
  result=exp(-par[1]/par[3])/(par[3]*(1-exp(-kt2max/par[3])));
  result=result*exp(-x[0])/(1.-exp(-par[5]));
  result=result/(2.*Pi);
  result=result*(par[2]*par[2]*par[1]+par[4]*x[0]+2.*par[2]*cos(par[0])*sqrt(par[4]*x[0]*par[1]));
  return result;
}

double meanPt2dd(Double_t *x, Double_t *par) // x[0]: theta, par[0]:kt2, par[1]:z, par[2]:mkt2, par[3]:mpt2, par[4]:Ptmax
{
  double vm,result;
  
  if(par[4]>2.*par[1]*cos(x[0])*sqrt(par[0])){
    vm=par[1]*par[1]*par[0]+par[4]*(par[4]-2.*par[1]*cos(x[0])*sqrt(par[0]));
  }
  else{
    vm=par[1]*par[1]*par[0];
  }
  vm=vm/par[3];
  TF1 f("func",meanPt2ddd,0,1,6);
  f.SetParameters(x[0],par[0],par[1],par[2],par[3],vm);
  ROOT::Math::WrappedTF1 wf1(f);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
  ig.SetFunction(wf1);
  if(vm<100.){
    result=ig.Integral(0.,vm);
  }
  else{
    result=ig.IntegralUp(0.);
  }
  if(ig.Status()!=0){cout<<"bad integration meanPt2dd"<<"  "<<ig.Status()<<endl;}

  return result;
}

double meanPt2d(Double_t *x, Double_t *par) // x[0]:kt2, par[0]:z, par[1]:mkt2, par[2]:mpt2, par[3]:Ptmax
{
  double result;

  TF1 f("func",meanPt2dd,0,1,5);
  f.SetParameters(x[0],par[0],par[1],par[2],par[3]);
  ROOT::Math::WrappedTF1 wf1(f);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
  ig.SetFunction(wf1);
  result=ig.Integral(0,2.*Pi);
  if(ig.Status()!=0){cout<<"bad integration meanPt2d"<<"  "<<ig.Status()<<endl;}
  return result;
}


double meanPt2(double z, double mkt2, double mpt2, double Ptmax)
{
  double result;
  TF1 f("func",meanPt2d,0,1,4);
  f.SetParameters(z,mkt2,mpt2,Ptmax);
  ROOT::Math::WrappedTF1 wf1(f);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
  ig.SetFunction(wf1);
  result=ig.Integral(0,kt2max);
  if(ig.Status()!=0){cout<<"bad integration meanPt2"<<"  "<<ig.Status()<<endl;}

  return result;
}



void MACRO2()
{
  TBenchmark bench;
  bench.Start("main");
  
  double z,mkt2,mpt2;
  z=0.99;
  mkt2=1;
  mpt2=0.01;
  
  std::cout<<"kt2_MAX="<<kt2max<<std::endl;
  std::cout<<"Pt2_MAX="<<Pt2_MAX<<std::endl;
  std::cout<<"meankt2="<<mkt2<<" meanpt2="<<mpt2<<" zh="<<z<<std::endl;
  
  double num_result = meanPt2(z,mkt2,mpt2,Pt2_MAX);
  double ana_result = z*z*mkt2+mpt2;
  double difference = (num_result - ana_result)*100/num_result;
  cout<<num_result<<"  "<<ana_result<<"   "<<difference<<"% difference"<<endl;

  bench.Show("main");
  return;
}

