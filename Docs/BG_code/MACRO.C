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
double kt2max=10.; // in GeV^2

double meanPt2ddd(Double_t *x, Double_t *par) //x[0]:u2 , par[0] theta, par[1]:kt2, par[2]:z, par[3]:mkt2, par[4]:mpt2
{
  double result;

  result=exp(-par[1]/par[3])/(par[3]*(1-exp(-kt2max/par[3])));                                     // Exponencial de kt con normalizacion para integracion hasta kt2max (sin Pi)
  result=result*exp(-x[0]);                                                                        // Exponencial de pt con cambio de variable
  result=result/(2.*Pi);                                                                           // Normalizacion
  result=result*(par[2]*par[2]*par[1]+par[4]*x[0]+2.*par[2]*cos(par[0])*sqrt(par[4]*x[0]*par[1])); // Multiplicacion por PT2
  return result;

  /*OBSERVACIONES                                                                                          */
  /*Normalizacion de f(x,kt) es 1/(Pi*mkt2*(1-exp(kt2MAX/mkt2)))                                           */
  /*Normalizacion de p(z,pt) es 1/(Pi*mpt2)                                                                */
  /*PT2 definido bajo relacion de PT = pt + zkt. Por lo tanto el angulo que se ocupa aca no es el correcto */
}

double meanPt2dd(Double_t *x, Double_t *par) // x[0]: theta, par[0]:kt2, par[1]:z, par[2]:mkt2, par[3]:mpt2, par[4]:Ptmax
{
  double vm,result;
  
  if(sqrt(kt2max)>par[4]){
    cout<<"sqrt(kt2max)>Ptmax: impossible"<<endl;
    exit(1);
  }
  if(par[4]>2.*par[1]*cos(x[0])*sqrt(par[0])){
    vm=par[1]*par[1]*par[0]+par[4]*(par[4]-2.*par[1]*cos(x[0])*sqrt(par[0]));
  }
  else{
    vm=par[1]*par[1]*par[0];
  }
  vm=vm/par[3];
  TF1 f("func",meanPt2ddd,0,1,5);
  f.SetParameters(x[0],par[0],par[1],par[2],par[3]);
  ROOT::Math::WrappedTF1 wf1(f);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
  ig.SetFunction(wf1);
  /*Eleccion debido a razones tecnicas de integracion*/
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



//int main(int argc,char *argv[])
void MACRO()
{
  TBenchmark bench;
  bench.Start("main");
  
  double z,mkt2,mpt2;
  z=0.8;
  mkt2=0.25;
  mpt2=0.2;
  cout<<meanPt2(z,mkt2,mpt2,10.)<<"  "<<z*z*mkt2+mpt2 <<endl;

  bench.Show("main");
  //  return 0;
  return;
}

