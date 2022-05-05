#include "FitTools.h"


using namespace std;
using namespace blinding;
using namespace pocketfft;

const double Pi = 3.14159265358979323846264338327950288419716939937510582;

double time_scale = 1.0;
Blinders::fitType ftype = Blinders::kOmega_a;
Blinders * getBlinded = new Blinders(ftype, "Unexpected virtue of ignorance.");
TH1 * lost_muon;

void FillData(TH1* th1,std::vector<double> & data);
void FillHist(TH1* th1,const std::vector<double> & data,bool isAbs=true);
void FFT(TH1* hist, TH1* hist_fft, bool isAbs=true);

double func_5paras(double *x,double *p);
double func_9paras_cbo(double *x,double *p);
double func_10paras_cbo_lost(double *x, double *p);
double func_14paras_cbo_lost_vw(double *x, double *p);

Fitter::Fitter() {
}

void Fitter::SetTimeUnit(TimeUnit t_unit) {
    if(t_unit == Fitter::nano_second) {
        time_scale = 1e3;
    }
}

void Fitter::SetOutputDir(string output_dir) {
    m_output_dir = output_dir;
}

FitOutputInfo Fitter::doFit(const FitInput & fit_in) {
    TString wiggle_name, residual_name, func_name, res_name, fft_name, outf_name;
    wiggle_name.Form("wiggle_%s",fit_in.tag.Data());
    residual_name.Form("residual_%s",fit_in.tag.Data());
    func_name.Form("func_%s",fit_in.tag.Data());
    res_name.Form("result_%s",fit_in.tag.Data());
    fft_name.Form("fft_%s",fit_in.tag.Data());
    outf_name.Form("%s/result_%s.root",m_output_dir.c_str(),fit_in.tag.Data());


    lost_muon = fit_in.lost_muon;

    TH1 * fit_hist = (TH1*)fit_in.wiggle->Clone(wiggle_name.Data());
    TF1 * fit_func = new TF1(func_name,fit_in.func, fit_in.t_start, fit_in.t_end, fit_in.nvars);
    for(int n=0;n<fit_in.nvars;n++){
        fit_func->SetParName(n,fit_in.name_vars[n].c_str());
        fit_func->SetParameter(n,fit_in.init_values[n]);
    }
    fit_func->SetNpx(4500);
    TFitResultPtr fit_res = fit_hist->Fit(fit_func,"REMS");

    double chi2 = fit_func->GetChisquare();
    double ndf = fit_func->GetNDF();
    cout << "Chi2/NDF = " << chi2/ndf << "    chi2=" << chi2 << " ndf=" << ndf << endl;
    fit_res->SetName(res_name.Data());

    TH1 * residual = (TH1*)fit_in.wiggle->Clone(residual_name.Data());
    residual->Reset();
    int bin_t_start = residual->GetXaxis()->FindBin(fit_in.t_start) + 1;
    int bin_t_end = residual->GetXaxis()->FindBin(fit_in.t_end) - 1;
    for(int bin=bin_t_start;bin<=bin_t_end;bin++) {
        double time = residual->GetBinCenter(bin);
        double val = fit_hist->GetBinContent(bin);
        double err = fit_hist->GetBinError(bin);
        double fval = fit_func->Eval(time);
        residual->SetBinContent(bin,val - fval);
        residual->SetBinError(bin,err);
    }

    TH1D * fft = new TH1D(fft_name.Data(),fft_name.Data(),residual->GetNbinsX(),0,1./0.1492);
    FFT(residual,fft);


    TFile * fout = new TFile(outf_name.Data(),"recreate");
    fout->cd();
    fit_hist->Write();
    fit_func->Write();
    fit_res->Write();
    residual->Write();
    fft->Write();
    fout->Close();

    cout << "\nOutput file created : "<< outf_name << endl;
    cout << "   wiggle     : " << wiggle_name << endl;
    cout << "   residual   : " << residual_name << endl;
    cout << "   function   : " << func_name << endl;
    cout << "   FFT        : " << fft_name << endl;
    cout << "   FitResults : " << res_name << endl;

    FitOutputInfo info;
    info.file_name = outf_name;
    info.hist_name = wiggle_name;
    info.residual_name = residual_name;
    info.function_name = func_name;
    info.fft_name = fft_name;

    return info;
}

FitOutputInfo Fitter::Fit_5paras(string name, TH1* wiggle, double t_start, double t_end, vector<double> init_values) {
    TString tag;
    tag.Form("5paras_%s",name.c_str());

    FitInput fit_in;
    fit_in.tag = tag;
    fit_in.wiggle = wiggle;
    fit_in.t_start = t_start;
    fit_in.t_end = t_end;
    fit_in.init_values = init_values;
    fit_in.nvars = 5;
    string name_vars[] = {"N","#tau","A","R","#phi"};
    fit_in.name_vars = name_vars;
    std::function<double(double*,double*)> func = func_5paras;
    fit_in.func = func;

    return doFit(fit_in);
}

FitOutputInfo Fitter::Fit_9paras_cbo(string name, TH1* wiggle, double t_start, double t_end, vector<double> init_values) {
    TString tag;
    tag.Form("9paras_cbo_%s",name.c_str());

    FitInput fit_in;
    fit_in.tag = tag;
    fit_in.wiggle = wiggle;
    fit_in.t_start = t_start;
    fit_in.t_end = t_end;
    fit_in.init_values = init_values;
    fit_in.nvars = 9;
    string name_vars[] = {"N","#tau","A","R","#phi","#tau_{cbo}","A_{cbo}","#omega_{cbo}","#phi_{cbo}"};
    fit_in.name_vars = name_vars;
    std::function<double(double*,double*)> func = func_9paras_cbo;
    fit_in.func = func;

    return doFit(fit_in);    
}


FitOutputInfo Fitter::Fit_10paras_cbo_lost(string name, TH1* wiggle, double t_start, double t_end, vector<double> init_values,TH1* lm) {
    TString tag;
    tag.Form("10paras_cbo_lost_%s",name.c_str());

    FitInput fit_in;
    fit_in.tag = tag;
    fit_in.wiggle = wiggle;
    fit_in.t_start = t_start;
    fit_in.t_end = t_end;
    fit_in.init_values = init_values;
    fit_in.nvars = 10;
    string name_vars[] = {"N","#tau","A","R","#phi","#tau_{cbo}","A_{cbo}","#omega_{cbo}","#phi_{cbo}","k_{loss}"};
    fit_in.name_vars = name_vars;
    std::function<double(double*,double*)> func = func_10paras_cbo_lost;
    fit_in.func = func;
    fit_in.lost_muon = lm;

    return doFit(fit_in);
}

FitOutputInfo Fitter::Fit_14paras_cbo_lost_vw(string name, TH1* wiggle, double t_start, double t_end, vector<double> init_values,TH1* lm) {
    TString tag;
    tag.Form("14paras_cbo_lost_vw_%s",name.c_str());

    FitInput fit_in;
    fit_in.tag = tag;
    fit_in.wiggle = wiggle;
    fit_in.t_start = t_start;
    fit_in.t_end = t_end;
    fit_in.init_values = init_values;
    fit_in.nvars = 14;
    string name_vars[] = {
        "N","#tau","A","R","#phi",
        "#tau_{cbo}","A_{cbo}","#omega_{cbo}","#phi_{cbo}",
        "k_{loss}",
        "#tau_{vw}","A_{vw}","#omega_{vw}","#phi_{vw}"
    };
    fit_in.name_vars = name_vars;
    std::function<double(double*,double*)> func = func_14paras_cbo_lost_vw;

    fit_in.func = func;
    fit_in.lost_muon = lm;

    return doFit(fit_in);
}



// fft implement

void FillData(TH1* th1,std::vector<double> & data) {  
  int nbins = th1->GetNbinsX();
  for(int n=0;n<nbins;n++) {
    data[n] = th1->GetBinContent(n+1);
  }  
}

void FillHist(TH1* th1,const std::vector<double> & data,bool isAbs) {
  int nbins = th1->GetNbinsX();
  for(int n=0;n<nbins;n++){
    if(isAbs) {
      th1->SetBinContent(n+1,fabs(data[n]));  
    }
    else{
      th1->SetBinContent(n+1,data[n]);   
    }    
  }
}

void FFT(TH1* hist, TH1* hist_fft, bool isAbs) {
  size_t nbins = hist->GetNbinsX();
  std::vector<double> data_in(nbins);
  std::vector<double> data_out(nbins);
  FillData(hist,data_in);
  detail::shape_t shape{nbins};
  stride_t strided{sizeof(double)};
  shape_t axes{0}; 
  r2r_separable_hartley(shape, strided, strided, axes, data_in.data(), data_out.data(), 1.);
  FillHist(hist_fft,data_out,isAbs);
}  





//fit functions

double func_5paras(double *x,double *p) {
  double time = x[0] / time_scale;
  double norm = p[0];
  double life = p[1];
  double asym = p[2];  
  double omega = getBlinded->paramToFreq(p[3]);
  double phi = p[4];

  return norm * TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time + phi));
}

double func_9paras_cbo(double *x, double *p) {
    double time = x[0] / time_scale;

    // 5 paras
    double norm = p[0];
    double life = p[1];
    double asym = p[2];  
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // 4 paras: cbo
    double tau_cbo = p[5];
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];
    
    double cbo = 1-TMath::Exp(-time/tau_cbo)*asym_cbo*TMath::Cos(omega_cbo*time + phi_cbo);
    return  norm * TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time + phi)) *  cbo;
}


double func_10paras_cbo_lost(double *x, double *p) {
    double time = x[0] / time_scale;
    
    // 5 paras
    double norm = p[0];
    double life = p[1];
    double asym = p[2];  
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // 4 paras: cbo
    double tau_cbo = p[5];
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];

    // k_loss
    double k = p[9];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    double cbo = 1-TMath::Exp(-time/tau_cbo)*asym_cbo*TMath::Cos(omega_cbo*time + phi_cbo);

    return  norm *(1 - k*aloss)* TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time + phi)) *  cbo;
}


double func_14paras_cbo_lost_vw(double *x, double *p) {
    double time = x[0] / time_scale;
    
    // 5 paras
    double norm = p[0];
    double life = p[1];
    double asym = p[2];  
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // 4 paras: cbo
    double tau_cbo = p[5];
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];

    // k_loss
    double k = p[9];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // 4 paras: vw
    // double fcbo = 2.34
    // fc = 2*Pi/0.1492;
    // fvo = TMath::Sqrt(fcbo*(2*fc - fcbo));
    // fvw = fc - 2*fvo;
    double tau_vw = p[10];
    double asym_vw = p[11];
    double omega_vw = p[12];
    double phi_vw = p[13];

    double cbo = 1-TMath::Exp(-time/tau_cbo)*asym_cbo*TMath::Cos(omega_cbo*time + phi_cbo);
    double vw  = 1-TMath::Exp(-time/tau_vw)*asym_vw*TMath::Cos(omega_vw*time + phi_vw);
    return  norm *(1 - k*aloss)* TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time + phi)) *  cbo * vw;
}