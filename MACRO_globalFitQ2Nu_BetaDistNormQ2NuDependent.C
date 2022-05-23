// Based on combinedFit.C from ROOT tutorials

// Description:
// MACRO that performs a global fit on the nine (Q2, Nu) bins in function of Zh

// Reqs:
// - File with TNtuple containing bin number and limits
// - File with Pt2 distributions (to calculate Pt2 MAX)
// - File with meanPt2 distributions


// Arrays that contains the integrals' limits
Double_t kt2_MAX[3][3], Pt2_MAX[3][3];

Double_t absTol=1.E-3;
Double_t relTol=1.E-3;

Double_t Q2_limits[4] = {1,1.3,1.8,4};
Double_t Nu_limits[4] = {2.2,3.2,3.7,4.26};

// Zh limits to define fit
Double_t zh_min = 0.4;
Double_t zh_max = 0.98;

// Parameters' index for every fitting function. Repeated indexes are shared.
// Index format : {meankt2, Norm, Beta Par 1, Beta Par 2}
int ipar00[4] = {0,1,2,3};   int ipar01[4] = {4,1,2,5};   int ipar02[4] = {6,1,2,7};
int ipar10[4] = {8,1,2,9};   int ipar11[4] = {10,1,2,11}; int ipar12[4] = {12,1,2,13};
int ipar20[4] = {14,1,2,15}; int ipar21[4] = {16,1,2,17}; int ipar22[4] = {18,1,2,19};

// FUNCTIONS
void Getkt2MAX();
void GetPt2MAX(TString material = "DC");

Double_t GetFirstEmptyPt2(TH1F* h);

double meanPt2ddd(Double_t *x, Double_t *par);
double meanPt2dd(Double_t *x, Double_t *par);
double meanPt2d(Double_t *x, Double_t *par);

//00
double meanPt200(double z, double mkt2, double mpt2, double Ptmax);
Double_t fitfunc00(Double_t* x, Double_t* par);
//01
double meanPt201(double z, double mkt2, double mpt2, double Ptmax);
Double_t fitfunc01(Double_t* x, Double_t* par);
//02
double meanPt202(double z, double mkt2, double mpt2, double Ptmax);
Double_t fitfunc02(Double_t* x, Double_t* par);
//10
double meanPt210(double z, double mkt2, double mpt2, double Ptmax);
Double_t fitfunc10(Double_t* x, Double_t* par);
//11
double meanPt211(double z, double mkt2, double mpt2, double Ptmax);
Double_t fitfunc11(Double_t* x, Double_t* par);
//12
double meanPt212(double z, double mkt2, double mpt2, double Ptmax);
Double_t fitfunc12(Double_t* x, Double_t* par);
//20
double meanPt220(double z, double mkt2, double mpt2, double Ptmax);
Double_t fitfunc20(Double_t* x, Double_t* par);
//21
double meanPt221(double z, double mkt2, double mpt2, double Ptmax);
Double_t fitfunc21(Double_t* x, Double_t* par);
//22
double meanPt222(double z, double mkt2, double mpt2, double Ptmax);
Double_t fitfunc22(Double_t* x, Double_t* par);

void SetQ2LimitsPads(TPad* p1, TPad* p2, TPad* p3, TLatex* t);
void SetNuLimitsPads(TPad* p1, TPad* p2, TPad* p3, TLatex* t);

void PrintTableLatex(TString material, Double_t* Q2_limits, Double_t* Nu_limits, ROOT::Fit::FitResult result, const Double_t* pars, const Double_t* pars_err);

// Create the GlobalCHi2 structure
struct GlobalChi2 {
  GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
	       ROOT::Math::IMultiGenFunction & f2,
	       ROOT::Math::IMultiGenFunction & f3,
	       ROOT::Math::IMultiGenFunction & f4,
	       ROOT::Math::IMultiGenFunction & f5,
	       ROOT::Math::IMultiGenFunction & f6,
	       ROOT::Math::IMultiGenFunction & f7,
	       ROOT::Math::IMultiGenFunction & f8,
	       ROOT::Math::IMultiGenFunction & f9) :
    fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3),fChi2_4(&f4), fChi2_5(&f5), fChi2_6(&f6),fChi2_7(&f7), fChi2_8(&f8), fChi2_9(&f9) {}

  double operator() (const double *par) const {

    double p00[4];    double p01[4];    double p02[4];
    double p10[4];    double p11[4];    double p12[4];
    double p20[4];    double p21[4];    double p22[4];
    for (int i = 0; i < 4; ++i){
      p00[i] = par[ ipar00[i] ];
      p01[i] = par[ ipar01[i] ];
      p02[i] = par[ ipar02[i] ];
      p10[i] = par[ ipar10[i] ];
      p11[i] = par[ ipar11[i] ];
      p12[i] = par[ ipar12[i] ];
      p20[i] = par[ ipar20[i] ];
      p21[i] = par[ ipar21[i] ];
      p22[i] = par[ ipar22[i] ];
    }
    
    return (*fChi2_1)(p00) + (*fChi2_2)(p01) + (*fChi2_3)(p02) + (*fChi2_4)(p10) + (*fChi2_5)(p11) + (*fChi2_6)(p12) + (*fChi2_7)(p20) + (*fChi2_8)(p21) + (*fChi2_9)(p22);
  }

  const  ROOT::Math::IMultiGenFunction * fChi2_1;  const  ROOT::Math::IMultiGenFunction * fChi2_2;  const  ROOT::Math::IMultiGenFunction * fChi2_3;
  const  ROOT::Math::IMultiGenFunction * fChi2_4;  const  ROOT::Math::IMultiGenFunction * fChi2_5;  const  ROOT::Math::IMultiGenFunction * fChi2_6;
  const  ROOT::Math::IMultiGenFunction * fChi2_7;  const  ROOT::Math::IMultiGenFunction * fChi2_8;  const  ROOT::Math::IMultiGenFunction * fChi2_9;
};

// Visual specs
Double_t lp_leftmargin  = 0.12;
Double_t rp_rightmargin = 0.1;

Double_t up_margin      = 0.1;
Double_t bottom_margin  = 0.1;

Double_t lp_length = 1./(2-lp_leftmargin+(1-lp_leftmargin)/(1-rp_rightmargin));
Double_t cp_length = (1-lp_leftmargin)*lp_length;
Double_t rp_length = (1-lp_leftmargin)*lp_length/(1-rp_rightmargin);

Double_t top_height    = 1./(2-up_margin+(1-up_margin)/(1-bottom_margin));
Double_t med_height    = (1-up_margin)*top_height;
Double_t bottom_height = (1-up_margin)*top_height/(1-bottom_margin);

Double_t pad_X1[3][3] = {{0 , lp_length , lp_length+cp_length},
			  {0 , lp_length , lp_length+cp_length},
			  {0 , lp_length , lp_length+cp_length}};

Double_t pad_X2[3][3] = {{lp_length , lp_length+cp_length , 1},
			  {lp_length , lp_length+cp_length , 1},
			  {lp_length , lp_length+cp_length , 1}};

Double_t pad_Y1[3][3] = {{bottom_height+med_height , bottom_height+med_height , bottom_height+med_height},
			  {bottom_height            , bottom_height            , bottom_height},
			  {0 , 0 , 0}};

Double_t pad_Y2[3][3] = {{1,1,1},
			  {bottom_height+med_height , bottom_height+med_height , bottom_height+med_height},
			  {bottom_height            , bottom_height            , bottom_height}};

Double_t b_margins[3][3] = {{0,0,0},
			    {0,0,0},
			    {bottom_margin,bottom_margin,bottom_margin}};

Double_t t_margins[3][3] = {{up_margin,up_margin,up_margin},
			    {0,0,0},
			    {0,0,0}};

Double_t l_margins[3][3] = {{lp_leftmargin,0,0},
			    {lp_leftmargin,0,0},
			    {lp_leftmargin,0,0}};

Double_t r_margins[3][3] = {{0,0,rp_rightmargin},
			    {0,0,rp_rightmargin},
			    {0,0,rp_rightmargin}};

void MACRO_globalFitQ2Nu_BetaDistNormQ2NuDependent(TString material = "DC")
{

  // Visualization 
  TCanvas* c = new TCanvas("c","",1200,1200);
  c->Draw();
  c->cd();
  TPad* p[3][3]; 
  TLatex* t = new TLatex();
  t->SetTextAlign(22);
  t->SetTextFont(62);
  t->SetTextSize(0.05);
  TLatex* t2 = new TLatex();
  t2->SetTextAlign(22);
  t2->SetTextFont(62);
  t2->SetTextSize(0.05);
  t2->SetTextAngle(-90);
  for(Int_t i = 0 ; i < 3 ; i++){
    for(Int_t j = 0 ; j < 3 ; j++){
      p[i][j] = new TPad(Form("p%i%i",i,j),Form("p%i%i",i,j),pad_X1[i][j],pad_Y1[i][j],pad_X2[i][j],pad_Y2[i][j]);
      p[i][j]->SetBottomMargin(b_margins[i][j]);
      p[i][j]->SetTopMargin(t_margins[i][j]);
      p[i][j]->SetLeftMargin(l_margins[i][j]);
      p[i][j]->SetRightMargin(r_margins[i][j]);

      p[i][j]->SetGridx(1);
      p[i][j]->SetGridy(1);
    }
  }

  // Open file with data
  TFile* f = new TFile("meanPt2_with_systerr.root","READ");

  // Determine integrals' upper limits
  Getkt2MAX();
  GetPt2MAX(material);
      
  // Get meanPt2 histos
  TH1F* h[3][3];
  for(Int_t Q2_bin = 0 ; Q2_bin < 3 ; Q2_bin++){
      for(Int_t Nu_bin = 0 ; Nu_bin < 3 ; Nu_bin++){
        h[Q2_bin][Nu_bin] = (TH1F*) f->Get(Form("meanPt2_"+material+"_%i%i_CLEAN_INTERPOLATED_SYST", Q2_bin, Nu_bin));
      }      
  }

  std::cout<<"Histos ready"<<std::endl;

  // Define fitting function
  TF1* func00 = new TF1("func00",fitfunc00,zh_min,zh_max,4,1);
  TF1* func01 = new TF1("func01",fitfunc01,zh_min,zh_max,4,1);
  TF1* func02 = new TF1("func02",fitfunc02,zh_min,zh_max,4,1);
  TF1* func10 = new TF1("func10",fitfunc10,zh_min,zh_max,4,1);
  TF1* func11 = new TF1("func11",fitfunc11,zh_min,zh_max,4,1);
  TF1* func12 = new TF1("func12",fitfunc12,zh_min,zh_max,4,1);
  TF1* func20 = new TF1("func20",fitfunc20,zh_min,zh_max,4,1);
  TF1* func21 = new TF1("func21",fitfunc21,zh_min,zh_max,4,1);
  TF1* func22 = new TF1("func22",fitfunc22,zh_min,zh_max,4,1);
  std::cout<<"Funcs ready"<<std::endl;

  // Wrap fitting functions

  ROOT::Math::WrappedMultiTF1 wfunc00(*func00,1);
  ROOT::Math::WrappedMultiTF1 wfunc01(*func01,1);
  ROOT::Math::WrappedMultiTF1 wfunc02(*func02,1);
  ROOT::Math::WrappedMultiTF1 wfunc10(*func10,1);
  ROOT::Math::WrappedMultiTF1 wfunc11(*func11,1);
  ROOT::Math::WrappedMultiTF1 wfunc12(*func12,1);
  ROOT::Math::WrappedMultiTF1 wfunc20(*func20,1);
  ROOT::Math::WrappedMultiTF1 wfunc21(*func21,1);
  ROOT::Math::WrappedMultiTF1 wfunc22(*func22,1);

  
  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange range;
  range.SetRange(zh_min,zh_max);

  ROOT::Fit::BinData data[3][3];
  for(Int_t Q2_bin = 0 ; Q2_bin < 3 ; Q2_bin++){
      for(Int_t Nu_bin = 0 ; Nu_bin < 3 ; Nu_bin++){
	data[Q2_bin][Nu_bin] = ROOT::Fit::BinData(opt,range);
	ROOT::Fit::FillData(data[Q2_bin][Nu_bin], h[Q2_bin][Nu_bin]);
      }
  }
  std::cout<<"BinData set"<<std::endl;
  
  ROOT::Fit::Chi2Function chi2_00(data[0][0], wfunc00);
  ROOT::Fit::Chi2Function chi2_01(data[0][1], wfunc01);
  ROOT::Fit::Chi2Function chi2_02(data[0][2], wfunc02);
  ROOT::Fit::Chi2Function chi2_10(data[1][0], wfunc10);
  ROOT::Fit::Chi2Function chi2_11(data[1][1], wfunc11);
  ROOT::Fit::Chi2Function chi2_12(data[1][2], wfunc12);
  ROOT::Fit::Chi2Function chi2_20(data[2][0], wfunc20);
  ROOT::Fit::Chi2Function chi2_21(data[2][1], wfunc21);
  ROOT::Fit::Chi2Function chi2_22(data[2][2], wfunc22);

  std::cout<<"Chi2Function set"<<std::endl;
  
  GlobalChi2 globalChi2(chi2_00, chi2_01, chi2_02,chi2_10, chi2_11, chi2_12,chi2_20, chi2_21, chi2_22);    

  std::cout<<"Global Chi2 set"<<std::endl;

  // Fit object declaration
  ROOT::Fit::Fitter fitter;

  // Parameters initialization
  const int Npar = 20;
  double par0[Npar] = {0.2,1.5,1.5,1,0.2,1,0.2,1,0.2,1,0.2,1,0.2,1,0.2,1,0.2,1,0.2,1};
  fitter.Config().SetParamsSettings(Npar,par0);

  // Parameters range setting
  fitter.Config().ParSettings(0).SetLimits(0.01,0.3);
  fitter.Config().ParSettings(1).SetLimits(0.5,8.5);
  fitter.Config().ParSettings(2).SetLimits(0.5,8.5);
  fitter.Config().ParSettings(3).SetLimits(0.001,5);
  fitter.Config().ParSettings(4).SetLimits(0.01,0.3);
  fitter.Config().ParSettings(5).SetLimits(0.001,5);
  fitter.Config().ParSettings(6).SetLimits(0.01,0.3);
  fitter.Config().ParSettings(7).SetLimits(0.001,5);
  fitter.Config().ParSettings(8).SetLimits(0.01,0.3);
  fitter.Config().ParSettings(9).SetLimits(0.001,5);
  fitter.Config().ParSettings(10).SetLimits(0.01,0.3);
  fitter.Config().ParSettings(11).SetLimits(0.001,5);
  fitter.Config().ParSettings(12).SetLimits(0.01,0.3);
  fitter.Config().ParSettings(13).SetLimits(0.001,5);
  fitter.Config().ParSettings(14).SetLimits(0.01,0.3);
  fitter.Config().ParSettings(15).SetLimits(0.001,5);
  fitter.Config().ParSettings(16).SetLimits(0.01,0.3);
  fitter.Config().ParSettings(17).SetLimits(0.001,5);
  fitter.Config().ParSettings(18).SetLimits(0.01,0.3);
  fitter.Config().ParSettings(19).SetLimits(0.001,5);

  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit","MigradImproved");

  std::cout<<"Fitter settings done"<<std::endl;
  
  // (specify optionally data size and flag to indicate that is a chi2 fit)
  fitter.FitFCN(Npar,globalChi2,0,data[0][0].Size()+data[0][1].Size()+data[0][2].Size()+data[1][0].Size()+data[1][1].Size()+data[1][2].Size()+data[2][0].Size()+data[2][1].Size()+data[2][2].Size(),true);

  std::cout<<"Fit done"<<std::endl;

  // Store results of the fit
  ROOT::Fit::FitResult result = fitter.Result();

  // Print results of fit
  result.Print(std::cout);

  // Parameter obtention
  const Double_t* pars     = result.GetParams();
  const Double_t* pars_err = result.GetErrors();
  
  std::cout<<"Chi2NDF = "<<result.Chi2()<<"/"<<result.Ndf()<<" = "<<result.Chi2()/result.Ndf()<<std::endl;
  std::cout<<"NCalls  = "<<result.NCalls()<<std::endl;  

  // Table output in LaTex format
  PrintTableLatex( material, Q2_limits, Nu_limits, result, pars, pars_err);

  // Plots customization
  c->cd();
  p[0][0]->Draw();
  p[0][0]->cd();

  func00->SetFitResult(result, ipar00);
  func00->SetRange(range().first, range().second);
  func00->SetLineColor(kBlue);
  h[0][0]->GetListOfFunctions()->Add(func00);

  h[0][0]->Draw();
  h[0][0]->GetXaxis()->SetRangeUser(zh_min,zh_max);
  h[0][0]->GetYaxis()->SetRangeUser(0.02,0.34);
  h[0][0]->SetLabelFont(62,"XY");
  h[0][0]->SetLabelSize(0.04,"XY");
  h[0][0]->SetTitleFont(62,"XY");
  h[0][0]->SetTitleSize(0.04,"XY");
  h[0][0]->SetTitle(";;<P^{2}_{T}>[GeV^{2}]");
  h[0][0]->GetYaxis()->CenterTitle();
  h[0][0]->Draw();
  
  c->cd();

  p[0][1]->Draw();
  p[0][1]->cd();
  
  func01->SetFitResult(result, ipar01);
  func01->SetRange(range().first, range().second);
  func01->SetLineColor(kBlue);
  h[0][1]->GetListOfFunctions()->Add(func01);
  h[0][1]->Draw();
  h[0][1]->GetXaxis()->SetRangeUser(zh_min,zh_max);
  h[0][1]->GetYaxis()->SetRangeUser(0.02,0.34);
  h[0][1]->SetLabelFont(62,"XY");
  h[0][1]->SetLabelSize(0.04,"XY");
  h[0][1]->SetTitleFont(62,"XY");
  h[0][1]->SetTitleSize(0.04,"XY");
  h[0][1]->Draw();
  c->cd();

  p[0][2]->Draw();
  p[0][2]->cd();
  func02->SetFitResult(result, ipar02);
  func02->SetRange(range().first, range().second);
  func02->SetLineColor(kBlue);
  h[0][2]->GetListOfFunctions()->Add(func02);
  h[0][2]->Draw();
  h[0][2]->GetXaxis()->SetRangeUser(zh_min,zh_max);
  h[0][2]->GetYaxis()->SetRangeUser(0.02,0.34);
  h[0][2]->SetLabelFont(62,"XY");
  h[0][2]->SetLabelSize(0.04,"XY");
  h[0][2]->SetTitleFont(62,"XY");
  h[0][2]->SetTitleSize(0.04,"XY");
  h[0][2]->Draw();
  c->cd();

  p[1][0]->Draw();
  p[1][0]->cd();

  gStyle->SetOptFit(0);
  func10->SetFitResult(result, ipar10);
  func10->SetRange(range().first, range().second);
  func10->SetLineColor(kBlue);
  h[1][0]->GetListOfFunctions()->Add(func10);
  h[1][0]->Draw();
  h[1][0]->GetXaxis()->SetRangeUser(zh_min,zh_max);
  h[1][0]->GetYaxis()->SetRangeUser(0.02,0.34);
  h[1][0]->SetLabelFont(62,"XY");
  h[1][0]->SetLabelSize(0.04,"XY");
  h[1][0]->SetTitleFont(62,"XY");
  h[1][0]->SetTitleSize(0.04,"XY");
  h[1][0]->SetTitle(";;<P^{2}_{T}>[GeV^{2}]");
  h[1][0]->GetYaxis()->CenterTitle();
  h[1][0]->Draw();
  c->cd();

  p[1][1]->Draw();
  p[1][1]->cd();
  func11->SetFitResult(result, ipar11);
  func11->SetRange(range().first, range().second);
  func11->SetLineColor(kBlue);
  h[1][1]->GetListOfFunctions()->Add(func11);
  h[1][1]->Draw();
  h[1][1]->GetXaxis()->SetRangeUser(zh_min,zh_max);
  h[1][1]->GetYaxis()->SetRangeUser(0.02,0.34);
  h[1][1]->SetLabelFont(62,"XY");
  h[1][1]->SetLabelSize(0.04,"XY");
  h[1][1]->SetTitleFont(62,"XY");
  h[1][1]->SetTitleSize(0.04,"XY");
  h[1][1]->Draw();
  c->cd();

  p[1][2]->Draw();
  p[1][2]->cd();
  func12->SetFitResult(result, ipar12);
  func12->SetRange(range().first, range().second);
  func12->SetLineColor(kBlue);
  h[1][2]->GetListOfFunctions()->Add(func12);
  h[1][2]->Draw();
  h[1][2]->GetXaxis()->SetRangeUser(zh_min,zh_max);
  h[1][2]->GetYaxis()->SetRangeUser(0.02,0.34);
  h[1][2]->SetLabelFont(62,"XY");
  h[1][2]->SetLabelSize(0.04,"XY");
  h[1][2]->SetTitleFont(62,"XY");
  h[1][2]->SetTitleSize(0.04,"XY");
  h[1][2]->Draw();
  c->cd();

  p[2][0]->Draw();
  p[2][0]->cd();
  gStyle->SetOptFit(0);
  func20->SetFitResult(result, ipar20);
  func20->SetRange(range().first, range().second);
  func20->SetLineColor(kBlue);
  h[2][0]->GetListOfFunctions()->Add(func20);
  h[2][0]->Draw();
  h[2][0]->GetXaxis()->SetRangeUser(zh_min,zh_max);
  h[2][0]->GetYaxis()->SetRangeUser(0.02,0.34);
  h[2][0]->SetLabelFont(62,"XY");
  h[2][0]->SetLabelSize(0.04,"XY");
  h[2][0]->SetTitleFont(62,"XY");
  h[2][0]->SetTitleSize(0.04,"XY");
  h[2][0]->SetTitle(";z_{h};<P^{2}_{T}>[GeV^{2}]");
  h[2][0]->GetYaxis()->CenterTitle();
  h[2][0]->GetXaxis()->CenterTitle();
  h[2][0]->Draw();
  c->cd();

  p[2][1]->Draw();
  p[2][1]->cd();
  func21->SetFitResult(result, ipar21);
  func21->SetRange(range().first, range().second);
  func21->SetLineColor(kBlue);
  h[2][1]->GetListOfFunctions()->Add(func21);
  h[2][1]->Draw();
  h[2][1]->GetXaxis()->SetRangeUser(zh_min,zh_max);
  h[2][1]->GetYaxis()->SetRangeUser(0.02,0.34);
  h[2][1]->SetLabelFont(62,"XY");
  h[2][1]->SetLabelSize(0.04,"XY");
  h[2][1]->SetTitleFont(62,"XY");
  h[2][1]->SetTitleSize(0.04,"XY");
  h[2][1]->SetTitle(";z_{h};");
  h[2][1]->GetXaxis()->CenterTitle();
  h[2][1]->Draw();
  c->cd();

  p[2][2]->Draw();
  p[2][2]->cd();
  func22->SetFitResult(result, ipar22);
  func22->SetRange(range().first, range().second);
  func22->SetLineColor(kBlue);
  h[2][2]->GetListOfFunctions()->Add(func22);
  h[2][2]->Draw();
  h[2][2]->GetXaxis()->SetRangeUser(zh_min,zh_max); 
  h[2][2]->GetYaxis()->SetRangeUser(0.02,0.34);
  h[2][2]->SetLabelFont(62,"XY");
  h[2][2]->SetLabelSize(0.04,"XY");
  h[2][2]->SetTitleFont(62,"XY");
  h[2][2]->SetTitleSize(0.04,"XY");
  h[2][2]->SetTitle(";z_{h};");
  h[2][2]->GetXaxis()->CenterTitle();
  h[2][2]->Draw();

  SetQ2LimitsPads(p[0][2],p[1][2],p[2][2],t2);
  SetNuLimitsPads(p[0][0],p[0][1],p[0][2],t);

  c->Print(Form("global_"+material+"_NormQ2NuDependent_%f_.pdf",result.Chi2()/result.Ndf()));
  c->Print(Form("global_"+material+"_NormQ2NuDependent_%f_.png",result.Chi2()/result.Ndf()));
}

void GetPt2MAX(TString material = "DC")
{
  TFile* f_Pt2 = new TFile("corr_data_Pt2_processed.root","READ");
  
  for(Int_t Q2_bin = 0 ; Q2_bin < 3 ; Q2_bin++){
    for(Int_t Nu_bin = 0 ; Nu_bin < 3 ; Nu_bin++){
      for(Int_t Zh_bin = 4 ; Zh_bin<5 ; Zh_bin++){

	TH1F* h = (TH1F*) f_Pt2->Get(Form("corr_data_Pt2_"+material+"_%i%i%i_CLEAN_INTERPOLATED",Q2_bin,Nu_bin,Zh_bin));  

	Pt2_MAX[Q2_bin][Nu_bin] = GetFirstEmptyPt2(h);
	
      }
    }
  }

  delete f_Pt2;
}

void Getkt2MAX()
{
  TFile* fbinning = new TFile("binning_3389012_Zh_fancy.root","READ");    

  TNtuple* limits_tuple = (TNtuple*) fbinning->Get("limits_tuple");
  Float_t Q2_min,Q2_max,Nu_min,Nu_max,Zh_min,Zh_max;
  Float_t Q2_min_local, Q2_max_local, Nu_min_local, Nu_max_local, Zh_min_local, Zh_max_local, Q2_bin_local, Nu_bin_local, Zh_bin_local;
  
  limits_tuple->SetBranchAddress("Q2_min",&Q2_min_local);
  limits_tuple->SetBranchAddress("Q2_max",&Q2_max_local);
  limits_tuple->SetBranchAddress("Nu_min",&Nu_min_local);
  limits_tuple->SetBranchAddress("Nu_max",&Nu_max_local);
  limits_tuple->SetBranchAddress("Zh_min",&Zh_min_local);
  limits_tuple->SetBranchAddress("Zh_max",&Zh_max_local);

  limits_tuple->SetBranchAddress("Q2_bin",&Q2_bin_local);
  limits_tuple->SetBranchAddress("Nu_bin",&Nu_bin_local);
  limits_tuple->SetBranchAddress("Zh_bin",&Zh_bin_local);

  for(Int_t Q2_bin = 0 ; Q2_bin < 3 ; Q2_bin++){
    for(Int_t Nu_bin = 0 ; Nu_bin < 3 ; Nu_bin++){

      for(Int_t tuple_entry = 0 ; tuple_entry < limits_tuple->GetEntries() ; tuple_entry++)
	{
	  limits_tuple->GetEntry(tuple_entry);
	  if(Q2_bin==Q2_bin_local&&Nu_bin==Nu_bin_local)
	    {
	      Q2_min = Q2_min_local;
	      Q2_max = Q2_max_local;
	      Nu_min = Nu_min_local;
	      Nu_max = Nu_max_local;
	      break;
	    }
	}
      Double_t Q2_center = (Q2_max+Q2_min)/2.;
      Double_t Nu_center = (Nu_max+Nu_min)/2.;
  
      Double_t xb = Q2_center/2./0.938/Nu_center;
  
      if(xb<0.5) {kt2_MAX[Q2_bin][Nu_bin] =  xb*(1-xb)*Q2_center/(1-2.*xb)/(1-2.*xb);}
      else       {kt2_MAX[Q2_bin][Nu_bin] = (2.-xb)*(1.-xb)*Q2_center;}
    }
  }  
}

double meanPt2ddd(Double_t *x, Double_t *par) //x[0]:u2 , par[0] theta, par[1]:kt2, par[2]:z, par[3]:mkt2, par[4]:mpt2, par[5]=vm, par[6]:kt2max
{
  double result;
  
  result=exp(-par[1]/par[3])/(par[3]*(1-exp(-par[6]/par[3])));
  result=result*exp(-x[0])/(1.-exp(-par[5]));
  result=result/(2.*TMath::Pi());
  result=result*(par[2]*par[2]*par[1]+par[4]*x[0]+2.*par[2]*cos(par[0])*sqrt(par[4]*x[0]*par[1]));
  return result;
}

double meanPt2dd(Double_t *x, Double_t *par) // x[0]: theta, par[0]:kt2, par[1]:z, par[2]:mkt2, par[3]:mpt2, par[4]:Ptmax, par[5]:kt2max
{
  double vm,result;
  
  if(par[4]>2.*par[1]*cos(x[0])*sqrt(par[0])){
    vm=par[1]*par[1]*par[0]+par[4]*(par[4]-2.*par[1]*cos(x[0])*sqrt(par[0]));
  }
  else{
    vm=par[1]*par[1]*par[0];
  }
  vm=vm/par[3];
  TF1 f("func",meanPt2ddd,0,1,7);
  f.SetParameters(x[0],par[0],par[1],par[2],par[3],vm, par[5]);
  ROOT::Math::WrappedTF1 wf1(f);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR,absTol,relTol);
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

double meanPt2d(Double_t *x, Double_t *par) // x[0]:kt2, par[0]:z, par[1]:mkt2, par[2]:mpt2, par[3]:Ptmax, par[4]:kt2max
{
  double result;

  TF1 f("func",meanPt2dd,0,1,6);
  f.SetParameters(x[0],par[0],par[1],par[2],par[3],par[4]);
  ROOT::Math::WrappedTF1 wf1(f);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR,absTol,relTol);
  ig.SetFunction(wf1);
  result=ig.Integral(0,2.*TMath::Pi());
  if(ig.Status()!=0){cout<<"bad integration meanPt2d"<<"  "<<ig.Status()<<endl;}
  return result;
}

//00
double meanPt200(double z, double mkt2, double mpt2, double Ptmax)
{
  double result;
  TF1 f("func",meanPt2d,0,1,5);
  f.SetParameters(z,mkt2,mpt2,Ptmax,kt2_MAX[0][0]);
  ROOT::Math::WrappedTF1 wf1(f);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR,absTol,relTol);
  ig.SetFunction(wf1);
  result=ig.Integral(0,kt2_MAX[0][0]);
  if(ig.Status()!=0){cout<<"bad integration meanPt2"<<"  "<<ig.Status()<<endl;}

  return result;
}
Double_t fitfunc00(Double_t* x, Double_t* par)
{
  Double_t zh      = x[0];
  Double_t meankt2 = par[0];  Double_t meanpt2_a = par[1];  Double_t meanpt2_b = par[2];  Double_t meanpt2_n = par[3];  
  Double_t meanpt2 = TMath::BetaDist(zh,meanpt2_a,meanpt2_b)*meanpt2_n;
  Double_t result  = meanPt200(zh,meankt2,meanpt2,Pt2_MAX[0][0]);
  return result;
}
//01
double meanPt201(double z, double mkt2, double mpt2, double Ptmax)
{
  double result;
  TF1 f("func",meanPt2d,0,1,5);
  f.SetParameters(z,mkt2,mpt2,Ptmax,kt2_MAX[0][1]);
  ROOT::Math::WrappedTF1 wf1(f);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR,absTol,relTol);
  ig.SetFunction(wf1);
  result=ig.Integral(0,kt2_MAX[0][1]);
  if(ig.Status()!=0){cout<<"bad integration meanPt2"<<"  "<<ig.Status()<<endl;}

  return result;
}
Double_t fitfunc01(Double_t* x, Double_t* par)
{
  Double_t zh      = x[0];
  Double_t meankt2 = par[0];  Double_t meanpt2_a = par[1];  Double_t meanpt2_b = par[2];  Double_t meanpt2_n = par[3];  
  Double_t meanpt2 = TMath::BetaDist(zh,meanpt2_a,meanpt2_b)*meanpt2_n;
  Double_t result  = meanPt201(zh,meankt2,meanpt2,Pt2_MAX[0][1]);
  return result;
}
//02
double meanPt202(double z, double mkt2, double mpt2, double Ptmax)
{
  double result;
  TF1 f("func",meanPt2d,0,1,5);
  f.SetParameters(z,mkt2,mpt2,Ptmax,kt2_MAX[0][2]);
  ROOT::Math::WrappedTF1 wf1(f);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR,absTol,relTol);
  ig.SetFunction(wf1);
  result=ig.Integral(0,kt2_MAX[0][2]);
  if(ig.Status()!=0){cout<<"bad integration meanPt2"<<"  "<<ig.Status()<<endl;}

  return result;
}
Double_t fitfunc02(Double_t* x, Double_t* par)
{
  Double_t zh      = x[0];
  Double_t meankt2 = par[0];  Double_t meanpt2_a = par[1];  Double_t meanpt2_b = par[2];  Double_t meanpt2_n = par[3];  
  Double_t meanpt2 = TMath::BetaDist(zh,meanpt2_a,meanpt2_b)*meanpt2_n;
  Double_t result  = meanPt202(zh,meankt2,meanpt2,Pt2_MAX[0][2]);
  return result;
}

//10
double meanPt210(double z, double mkt2, double mpt2, double Ptmax)
{
  double result;
  TF1 f("func",meanPt2d,0,1,5);
  f.SetParameters(z,mkt2,mpt2,Ptmax,kt2_MAX[1][0]);
  ROOT::Math::WrappedTF1 wf1(f);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR,absTol,relTol);
  ig.SetFunction(wf1);
  result=ig.Integral(0,kt2_MAX[1][0]);
  if(ig.Status()!=0){cout<<"bad integration meanPt2"<<"  "<<ig.Status()<<endl;}

  return result;
}
Double_t fitfunc10(Double_t* x, Double_t* par)
{
  Double_t zh      = x[0];
  Double_t meankt2 = par[0];  Double_t meanpt2_a = par[1];  Double_t meanpt2_b = par[2];  Double_t meanpt2_n = par[3];  
  Double_t meanpt2 = TMath::BetaDist(zh,meanpt2_a,meanpt2_b)*meanpt2_n;
  Double_t result  = meanPt210(zh,meankt2,meanpt2,Pt2_MAX[1][0]);
  return result;
}
//11
double meanPt211(double z, double mkt2, double mpt2, double Ptmax)
{
  double result;
  TF1 f("func",meanPt2d,0,1,5);
  f.SetParameters(z,mkt2,mpt2,Ptmax,kt2_MAX[1][1]);
  ROOT::Math::WrappedTF1 wf1(f);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR,absTol,relTol);
  ig.SetFunction(wf1);
  result=ig.Integral(0,kt2_MAX[1][1]);
  if(ig.Status()!=0){cout<<"bad integration meanPt2"<<"  "<<ig.Status()<<endl;}

  return result;
}
Double_t fitfunc11(Double_t* x, Double_t* par)
{
  Double_t zh      = x[0];
  Double_t meankt2 = par[0];  Double_t meanpt2_a = par[1];  Double_t meanpt2_b = par[2];  Double_t meanpt2_n = par[3];  
  Double_t meanpt2 = TMath::BetaDist(zh,meanpt2_a,meanpt2_b)*meanpt2_n;
  Double_t result  = meanPt211(zh,meankt2,meanpt2,Pt2_MAX[1][1]);
  return result;
}
//12
double meanPt212(double z, double mkt2, double mpt2, double Ptmax)
{
  double result;
  TF1 f("func",meanPt2d,0,1,5);
  f.SetParameters(z,mkt2,mpt2,Ptmax,kt2_MAX[1][2]);
  ROOT::Math::WrappedTF1 wf1(f);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR,absTol,relTol);
  ig.SetFunction(wf1);
  result=ig.Integral(0,kt2_MAX[1][2]);
  if(ig.Status()!=0){cout<<"bad integration meanPt2"<<"  "<<ig.Status()<<endl;}

  return result;
}
Double_t fitfunc12(Double_t* x, Double_t* par)
{
  Double_t zh      = x[0];
  Double_t meankt2 = par[0];  Double_t meanpt2_a = par[1];  Double_t meanpt2_b = par[2];  Double_t meanpt2_n = par[3];  
  Double_t meanpt2 = TMath::BetaDist(zh,meanpt2_a,meanpt2_b)*meanpt2_n;
  Double_t result  = meanPt212(zh,meankt2,meanpt2,Pt2_MAX[1][2]);
  return result;
}

//20
double meanPt220(double z, double mkt2, double mpt2, double Ptmax)
{
  double result;
  TF1 f("func",meanPt2d,0,1,5);
  f.SetParameters(z,mkt2,mpt2,Ptmax,kt2_MAX[2][0]);
  ROOT::Math::WrappedTF1 wf1(f);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR,absTol,relTol);
  ig.SetFunction(wf1);
  result=ig.Integral(0,kt2_MAX[2][0]);
  if(ig.Status()!=0){cout<<"bad integration meanPt2"<<"  "<<ig.Status()<<endl;}

  return result;
}
Double_t fitfunc20(Double_t* x, Double_t* par)
{
  Double_t zh      = x[0];
  Double_t meankt2 = par[0];  Double_t meanpt2_a = par[1];  Double_t meanpt2_b = par[2];  Double_t meanpt2_n = par[3];  
  Double_t meanpt2 = TMath::BetaDist(zh,meanpt2_a,meanpt2_b)*meanpt2_n;
  Double_t result  = meanPt220(zh,meankt2,meanpt2,Pt2_MAX[2][0]);
  return result;
}
//01
double meanPt221(double z, double mkt2, double mpt2, double Ptmax)
{
  double result;
  TF1 f("func",meanPt2d,0,1,5);
  f.SetParameters(z,mkt2,mpt2,Ptmax,kt2_MAX[2][1]);
  ROOT::Math::WrappedTF1 wf1(f);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR,absTol,relTol);
  ig.SetFunction(wf1);
  result=ig.Integral(0,kt2_MAX[2][1]);
  if(ig.Status()!=0){cout<<"bad integration meanPt2"<<"  "<<ig.Status()<<endl;}

  return result;
}
Double_t fitfunc21(Double_t* x, Double_t* par)
{
  Double_t zh      = x[0];
  Double_t meankt2 = par[0];  Double_t meanpt2_a = par[1];  Double_t meanpt2_b = par[2];  Double_t meanpt2_n = par[3];  
  Double_t meanpt2 = TMath::BetaDist(zh,meanpt2_a,meanpt2_b)*meanpt2_n;
  Double_t result  = meanPt221(zh,meankt2,meanpt2,Pt2_MAX[2][1]);
  return result;
}
//22
double meanPt222(double z, double mkt2, double mpt2, double Ptmax)
{
  double result;
  TF1 f("func",meanPt2d,0,1,5);
  f.SetParameters(z,mkt2,mpt2,Ptmax,kt2_MAX[2][2]);
  ROOT::Math::WrappedTF1 wf1(f);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR,absTol,relTol);
  ig.SetFunction(wf1);
  result=ig.Integral(0,kt2_MAX[2][2]);
  if(ig.Status()!=0){cout<<"bad integration meanPt2"<<"  "<<ig.Status()<<endl;}

  return result;
}
Double_t fitfunc22(Double_t* x, Double_t* par)
{
  Double_t zh      = x[0];
  Double_t meankt2 = par[0];  Double_t meanpt2_a = par[1];  Double_t meanpt2_b = par[2];  Double_t meanpt2_n = par[3];  
  Double_t meanpt2 = TMath::BetaDist(zh,meanpt2_a,meanpt2_b)*meanpt2_n;
  Double_t result  = meanPt222(zh,meankt2,meanpt2,Pt2_MAX[2][2]);
  return result;
}

Double_t GetFirstEmptyPt2(TH1F* h)
{
  for(Int_t i = 1 ; i <= h->GetNbinsX() ; i++){
    if(h->GetBinContent(i)==0){
      return h->GetBinCenter(i-1);
    }
  }
}

void SetQ2LimitsPads(TPad* p1, TPad* p2, TPad* p3, TLatex* t){
  p1->cd();
  t->DrawLatexNDC(1-p1->GetRightMargin()/2.,0.5,Form("%.1f<Q^{2}[GeV^{2}]<%.1f",Q2_limits[0],Q2_limits[1]));
  p2->cd();
  t->DrawLatexNDC(1-p1->GetRightMargin()/2.,0.5,Form("%.1f<Q^{2}[GeV^{2}]<%.1f",Q2_limits[1],Q2_limits[2]));
  p3->cd();
  t->DrawLatexNDC(1-p1->GetRightMargin()/2.,0.5,Form("%.1f<Q^{2}[GeV^{2}]<%.1f",Q2_limits[2],Q2_limits[3]));
}

void SetNuLimitsPads(TPad* p1, TPad* p2, TPad* p3, TLatex* t){
  p1->cd();
  t->DrawLatexNDC(0.5,1-p1->GetTopMargin()/2.,Form("%.1f<#nu[GeV]<%.1f",Nu_limits[0],Nu_limits[1]));
  p2->cd();
  t->DrawLatexNDC(0.5,1-p1->GetTopMargin()/2.,Form("%.1f<#nu[GeV]<%.1f",Nu_limits[1],Nu_limits[2]));
  p3->cd();
  t->DrawLatexNDC(0.5,1-p1->GetTopMargin()/2.,Form("%.1f<#nu[GeV]<%.1f",Nu_limits[2],Nu_limits[3]));
}

void PrintTableLatex(TString material, Double_t* Q2_limits, Double_t* Nu_limits, ROOT::Fit::FitResult result, const Double_t* pars, const Double_t* pars_err){
  std::cout<<"\\begin{table}[H]"<<std::endl;
  std::cout<<"\\centering"<<std::endl;
  std::cout<<"\\begin{tabular}{c|c|c|c|c|c}"<<std::endl;
  std::cout<<"\\hline"<<std::endl;
  std::cout<<"Target & $\\chi^2_{ndf}$ & $Q^2[GeV^2]$ & $\\nu[GeV]$ & $\\langle k^2_{\\perp}\\rangle[GeV^2]$  & $\\delta\\langle k^2_{\\perp}\\rangle[GeV^2]$\\\\ "<<std::endl;
  std::cout<<"\\hline"<<std::endl;
  std::cout<<"\\multirow{9}{*}{"+material+"} & "<<std::endl;
  std::cout<<"\\multirow{9}{*}{"<<std::setprecision(3)<<result.Chi2()/result.Ndf()<<"} & "<<std::endl;

  std::cout<<" "<<Q2_limits[0]<<"-"<<Q2_limits[1]<<" & "<<Nu_limits[0]<<"-"<<Nu_limits[1]<<" & "<<std::setprecision(3)<<pars[0]<<" & $\\pm$"<<std::setprecision(3)<<pars_err[0]<<" \\\\"<<std::endl;
  std::cout<<" "<<Q2_limits[0]<<"-"<<Q2_limits[1]<<" & "<<Nu_limits[1]<<"-"<<Nu_limits[2]<<" & "<<std::setprecision(3)<<pars[4]<<" & $\\pm$"<<std::setprecision(3)<<pars_err[4]<<" \\\\"<<std::endl;
  std::cout<<" "<<Q2_limits[0]<<"-"<<Q2_limits[1]<<" & "<<Nu_limits[2]<<"-"<<Nu_limits[3]<<" & "<<std::setprecision(3)<<pars[6]<<" & $\\pm$"<<std::setprecision(3)<<pars_err[6]<<" \\\\"<<std::endl;

  std::cout<<" "<<Q2_limits[1]<<"-"<<Q2_limits[2]<<" & "<<Nu_limits[0]<<"-"<<Nu_limits[1]<<" & "<<std::setprecision(3)<<pars[8]<<" & $\\pm$"<<std::setprecision(3)<<pars_err[8]<<" \\\\"<<std::endl;
  std::cout<<" "<<Q2_limits[1]<<"-"<<Q2_limits[2]<<" & "<<Nu_limits[1]<<"-"<<Nu_limits[2]<<" & "<<std::setprecision(3)<<pars[10]<<" & $\\pm$"<<std::setprecision(3)<<pars_err[10]<<" \\\\"<<std::endl;
  std::cout<<" "<<Q2_limits[1]<<"-"<<Q2_limits[2]<<" & "<<Nu_limits[2]<<"-"<<Nu_limits[3]<<" & "<<std::setprecision(3)<<pars[12]<<" & $\\pm$"<<std::setprecision(3)<<pars_err[12]<<" \\\\"<<std::endl;

  std::cout<<" "<<Q2_limits[2]<<"-"<<Q2_limits[3]<<" & "<<Nu_limits[0]<<"-"<<Nu_limits[1]<<" & "<<std::setprecision(3)<<pars[14]<<" & $\\pm$"<<std::setprecision(3)<<pars_err[14]<<" \\\\"<<std::endl;
  std::cout<<" "<<Q2_limits[2]<<"-"<<Q2_limits[3]<<" & "<<Nu_limits[1]<<"-"<<Nu_limits[2]<<" & "<<std::setprecision(3)<<pars[16]<<" & $\\pm$"<<std::setprecision(3)<<pars_err[16]<<" \\\\"<<std::endl;
  std::cout<<" "<<Q2_limits[2]<<"-"<<Q2_limits[3]<<" & "<<Nu_limits[2]<<"-"<<Nu_limits[3]<<" & "<<std::setprecision(3)<<pars[18]<<" & $\\pm$"<<std::setprecision(3)<<pars_err[18]<<" \\\\"<<std::endl;

  std::cout<<"\\hline"<<std::endl;
  std::cout<<"\\end{tabular}"<<std::endl;
  std::cout<<"\\caption{"<<material<<" global fit with Hypothesis 1. $x_f$ cut not applied.}"<<std::endl;
  std::cout<<"\\label{tab:"<<material<<"_global_integral_1}"<<std::endl;
  std::cout<<"\\end{table}"<<std::endl;
}
