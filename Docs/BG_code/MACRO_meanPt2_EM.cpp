Double_t Pt2_MAX = 10.;
Double_t kt2_MAX = 2.;

Double_t Integrand(Double_t *x, Double_t *par)
{
  /* PARS */
  Double_t meankt2 = par[0]; Double_t meanpt2 = par[1] ; Double_t zh = par[2]; Double_t kt2 = par[3]; Double_t theta = par[4]; Double_t pt2_MAX = par[5];
  
  Double_t constants    = 1./(2.*TMath::Pi()*meankt2*meanpt2)/(1.-TMath::Exp(-kt2_MAX/meankt2))/(1.-TMath::Exp(-pt2_MAX/meanpt2));
  Double_t Pt2          = x[0] + zh*zh*kt2 + 2*zh*TMath::Sqrt(kt2*x[0])*TMath::Cos(theta);
  Double_t Exponentials = TMath::Exp(-x[0]/meanpt2)*TMath::Exp(-kt2/meankt2);

  return constants*Exponentials*Pt2;
}

Double_t pt2_integration(Double_t *x, Double_t *par)
{
  /* PARS */
  Double_t meankt2 = par[0]; Double_t meanpt2 = par[1] ; Double_t zh = par[2];
  /* VARS */
  Double_t kt2 = x[0] ; Double_t theta = x[1];

  Double_t pt2_MIN = 0.;
  Double_t pt2_MAX = Pt2_MAX + zh*zh*meankt2 - 2*zh*TMath::Sqrt(Pt2_MAX*kt2)*TMath::Cos(theta);
  
  //  TF1* IntFunc2 = new TF1("IntFunc2",Integrand,0,pt2_MAX,6);
  TF1* IntFunc2 = new TF1("IntFunc2",Integrand,0,1,6);
  IntFunc2->SetParameters(meankt2, meanpt2, zh, kt2, theta, pt2_MAX);

  /*INTEGRAL CONFIG*/
  ROOT::Math::WrappedTF1       wf2(*IntFunc2);
  ROOT::Math::GSLIntegrator ig2(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR,1.E-9,1.E-9,1000000);
  ig2.SetFunction(wf2);

  if(Pt2_MAX + zh*zh*meankt2 < 2*zh*TMath::Sqrt(Pt2_MAX*kt2)*TMath::Cos(theta))
    {std::cout<<"WARNING : Negative value of pt2_MAX. Re-stating value"<<std::endl;
      //std::cout<<ig2.IntegralUp(0.)<<std::endl;
      return ig2.Integral(pt2_MIN,zh*zh*meankt2);}
      //return 0;

  return ig2.Integral(pt2_MIN,pt2_MAX);
}

Double_t kt_theta_integration(Double_t meankt2, Double_t meanpt2, Double_t zh)
{
  Double_t min_lim[2] = {0       , 0            };
  Double_t max_lim[2] = {kt2_MAX , 2*TMath::Pi()};
  
  //  TF2* IntFunc = new TF2("IntFunc",pt2_integration,0,kt2_MAX,0,2*TMath::Pi(),3);
  TF2* IntFunc = new TF2("IntFunc",pt2_integration,0,1,0,1,3);
  IntFunc->SetParameters(meankt2,meanpt2,zh);

  /*INTEGRAL CONFIG*/
  ROOT::Math::Functor            wf(*IntFunc,2);
  ROOT::Math::AdaptiveIntegratorMultiDim ig(/*ROOT::Math::AdaptiveIntegratorMultiDim::kADAPTIVE,*/1.E-6,1.E-6,1000000);
  ig.SetFunction(wf);

  Double_t value   = ig.Integral(min_lim,max_lim);

  return value;
}

void MACRO_meanPt2_EM(Double_t meankt2, Double_t meanpt2, Double_t zh)
{
  TBenchmark bench;
  bench.Start("main");

  Double_t meanPt2 = kt_theta_integration(meankt2, meanpt2, zh);

  std::cout<<"Numerical  = "<<meanPt2<<std::endl;
  std::cout<<"Analytical = "<<zh*zh*meankt2+meanpt2<<std::endl;

  Double_t difference = (meanPt2-(zh*zh*meankt2+meanpt2))*100./meanPt2;
  
  std::cout<<"There is a "<<difference<<"% difference."<<std::endl;
  bench.Show("main");

}
