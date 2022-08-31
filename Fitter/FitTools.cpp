#include "FitTools.h"
#include <sstream>

using namespace std;
using namespace blinding;
using namespace pocketfft;

const double Pi = 3.14159265358979323846264338327950288419716939937510582;

double time_scale = 1.0;
Blinders::fitType ftype = Blinders::kOmega_a;

#if data_version_major == 4
Blinders * getBlinded = new Blinders(ftype, "Unexpected virtue of ignorance.");
#elif data_version_major == 2
Blinders * getBlinded = new Blinders(ftype, "stay home, stay healthy!");
#elif data_version_major == 3
Blinders * getBlinded = new Blinders(ftype, "Bla Bla Bla!");
#endif

TH1 * lost_muon;


vector<double> GetInitialValuesD(string file_path, string index) {
  vector<double> init_values;
#ifdef USE_JSON
  Json::Value json_value;
  std::ifstream json_file(file_path.c_str());
  cout << "load parameters from json:" << file_path << endl;  
  json_file >> json_value;  
  for(auto val : json_value[index.c_str()]) {
      init_values.push_back(double(val.asFloat()));
  }  
#else 
  TFile * file = TFile::Open(file_path.c_str());
  TArrayD * arrayd = file->Get<TArrayD>(index.c_str());
  for(int n=0;n<arrayd->fN;n++) {
    init_values.push_back(arrayd->At(n));
  }
#endif
  return init_values;
}



Fitter::Fitter() : max_attempts(1) {
    cout << "Fitting on Run" << data_version_major << " data" << endl;
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

        
    Int_t fitStatus;    
    bool isValid = false;
    double chi2,ndf;

    TFitResultPtr fit_res;

    fit_res = fit_hist->Fit(fit_func,"REMS");


    fitStatus = fit_res;
    isValid = fit_res->IsValid() && (fitStatus%100==0);
    chi2 = fit_func->GetChisquare();
    ndf = fit_func->GetNDF();

    cout << "Chi2/NDF = " << chi2/ndf << "    Valid="<< isValid << "     Status="<< fitStatus << endl;

    int refit = 1;
    int max_attempts_ = max_attempts;
    while(max_attempts_>0 && !isValid) {
     cout << "Invalid fitting!!! Retry " << refit++ << endl;
     fit_res = fit_hist->Fit(fit_func,"QREMS");

     fitStatus = fit_res;
     isValid = fit_res->IsValid() && (fitStatus%100==0);
     chi2 = fit_func->GetChisquare();
     ndf = fit_func->GetNDF();

     cout << "Chi2/NDF = " << chi2/ndf << "    Valid="<< isValid << "     Status="<< fitStatus << endl;
     max_attempts_--;
    }

    cout << "\nResult-"<< fit_in.tag <<  "\tChi2/NDF=" << chi2/ndf << "\tValid="<< isValid << "\tStatus="<< fitStatus
    << "\tR=" << fit_func->GetParameter(3) << "+-" << fit_func->GetParError(3) << endl;

    vector<double> fit_values = fit_in.init_values;
    for(int n=0;n<fit_in.nvars;n++){
        fit_values[n] = fit_func->GetParameter(n);
    }

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
    // info.fitStatus = IsValid;
    info.fit_values = fit_values;

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
        "#tau_{vw}","A_{vw}","K_{vw}","#phi_{vw}"
    };
    fit_in.name_vars = name_vars;
    std::function<double(double*,double*)> func = func_14paras_cbo_lost_vw;

    fit_in.func = func;
    fit_in.lost_muon = lm;

    return doFit(fit_in);
}

FitOutputInfo Fitter::Fit_14paras_cbo_lost_vo(string name, TH1* wiggle, double t_start, double t_end, vector<double> init_values,TH1* lm) {
    TString tag;
    tag.Form("14paras_cbo_lost_vo_%s",name.c_str());

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
        "#tau_{y}","A_{y}","K_{y}","#phi_{y}",
    };
    fit_in.name_vars = name_vars;
    std::function<double(double*,double*)> func = func_14paras_cbo_lost_vo;

    fit_in.func = func;
    fit_in.lost_muon = lm;

    return doFit(fit_in);
}


FitOutputInfo Fitter::Fit_18paras_cbo_lost_vo_vw(string name, TH1* wiggle, double t_start, double t_end, vector<double> init_values,TH1* lm) {
    TString tag;
    tag.Form("18paras_cbo_lost_vo_vw_%s",name.c_str());

    FitInput fit_in;
    fit_in.tag = tag;
    fit_in.wiggle = wiggle;
    fit_in.t_start = t_start;
    fit_in.t_end = t_end;
    fit_in.init_values = init_values;
    fit_in.nvars = 18;
    string name_vars[] = {
        "N","#tau","A","R","#phi",
        "#tau_{cbo}","A_{cbo}","#omega_{cbo}","#phi_{cbo}",
        "k_{loss}",
        "#tau_{y}","A_{y}","K_{y}","#phi_{y}",
        "#tau_{vw}","A_{vw}","K_{vw}","#phi_{vw}",        
    };
    fit_in.name_vars = name_vars;
    std::function<double(double*,double*)> func = func_18paras_cbo_lost_vo_vw;

    fit_in.func = func;
    fit_in.lost_muon = lm;

    return doFit(fit_in);
}


FitOutputInfo Fitter::Fit_28paras_cbo_lost_vw_expansion(string name, TH1* wiggle, double t_start, double t_end, vector<double> init_values,TH1* lm) {
    TString tag;
    tag.Form("28paras_cbo_lost_vw_expansion_%s",name.c_str());

    FitInput fit_in;
    fit_in.tag = tag;
    fit_in.wiggle = wiggle;
    fit_in.t_start = t_start;
    fit_in.t_end = t_end;
    fit_in.init_values = init_values;
    fit_in.nvars = 28;
    string name_vars[] = {
        "N_{0}","#tau","A","R","#phi_{0}",
        "#tau_{cbo}","A_{cbo}","#omega_{cbo}","#phi_{cbo}",        
        "k_{loss}",
        "#tau_{vw}","A_{vw}","K_{vw}","#phi_{vw}",
        "A_{2cbo}","#phi_{2cbo}","A_{cbo,A}","#phi_{cbo,A}","A_{cbo,#phi}","#phi_{cbo,#phi}",
        "#tau_{y}","A_{y}","K_{y}","#phi_{y}",
        "A_{VW-cbo}","#phi_{VW-cbo}","A_{VW+cbo}","#phi_{VW+cbo}"
    };

    fit_in.name_vars = name_vars;
    std::function<double(double*,double*)> func = func_28paras_cbo_lost_vw_expansion;

    fit_in.func = func;
    fit_in.lost_muon = lm;

    return doFit(fit_in);
}

FitOutputInfo Fitter::Fit_22paras_cbo_lost_vw_expansion_lite(string name, TH1* wiggle, double t_start, double t_end, vector<double> init_values,TH1* lm) {
    TString tag;
    tag.Form("22paras_cbo_lost_vw_expansion_lite_%s",name.c_str());

    FitInput fit_in;
    fit_in.tag = tag;
    fit_in.wiggle = wiggle;
    fit_in.t_start = t_start;
    fit_in.t_end = t_end;
    fit_in.init_values = init_values;
    fit_in.nvars = 22;
    string name_vars[] = {
        "N_{0}","#tau","A","R","#phi_{0}",
        "#tau_{cbo}","A_{cbo}","#omega_{cbo}","#phi_{cbo}",        
        "k_{loss}",
        "#tau_{vw}","A_{vw}","K_{vw}","#phi_{vw}",
        "A_{VW+cbo}","#phi_{VW+cbo}",
        "#tau_{y}","A_{y}","K_{y}","#phi_{y}",
        "A_{cbo,A}","#phi_{cbo,A}",        
    };

    fit_in.name_vars = name_vars;
    std::function<double(double*,double*)> func = func_22paras_cbo_lost_vw_expansion_lite;

    fit_in.func = func;
    fit_in.lost_muon = lm;

    return doFit(fit_in);
}

FitOutputInfo Fitter::Fit_12paras_changing_cbo(string name, TH1* wiggle, double t_start, double t_end, vector<double> init_values,TH1* lm) {
    TString tag;
    tag.Form("12paras_changing_cbo_%s",name.c_str());

    FitInput fit_in;
    fit_in.tag = tag;
    fit_in.wiggle = wiggle;
    fit_in.t_start = t_start;
    fit_in.t_end = t_end;
    fit_in.init_values = init_values;
    fit_in.nvars = 12;
    string name_vars[] = {
        "N_{0}","#tau","A","R","#phi_{0}",
        "#tau_{cbo}","A_{cbo}","#omega_{cbo}","#phi_{cbo}",
        "k_{loss}"
        "A_{1}", "#tau_{1}",
        
    };

    fit_in.name_vars = name_vars;
    std::function<double(double*,double*)> func = func_12paras_changing_cbo;

    fit_in.func = func;
    fit_in.lost_muon = lm;

    return doFit(fit_in);
}

FitOutputInfo Fitter::Fit_11paras_changing_cbo(string name, TH1* wiggle, double t_start, double t_end, vector<double> init_values) {
    TString tag;
    tag.Form("11paras_changing_cbo_%s",name.c_str());

    FitInput fit_in;
    fit_in.tag = tag;
    fit_in.wiggle = wiggle;
    fit_in.t_start = t_start;
    fit_in.t_end = t_end;
    fit_in.init_values = init_values;
    fit_in.nvars = 11;
    
    string name_vars[] = {
        "N_{0}","#tau","A","R","#phi_{0}",
        "#tau_{cbo}","A_{cbo}","#omega_{cbo}","#phi_{cbo}",        
        "A_{1}", "#tau_{1}",        
    };

    fit_in.name_vars = name_vars;
    std::function<double(double*,double*)> func = func_11paras_changing_cbo;

    fit_in.func = func;    

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
    double k = p[9]*1e-9;
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
    double k = p[9]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // 4 paras: vw
    double fcbo = 2.32657;
    double fc = 2*M_PI/0.1492;
    double fvo = sqrt(fcbo*(2*fc - fcbo));
    double fvw = fc - 2*fvo;

    double tau_vw = p[10];
    double asym_vw = p[11];
    double omega_vw = p[12]*fvw;
    double phi_vw = p[13];




    double cbo = 1-TMath::Exp(-time/tau_cbo)*asym_cbo*TMath::Cos(omega_cbo*time + phi_cbo);
    double vw  = 1-TMath::Exp(-time/tau_vw)*asym_vw*TMath::Cos(omega_vw*time + phi_vw);
    return  norm *(1 - k*aloss)* TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time + phi)) *  cbo * vw;
}


double func_14paras_cbo_lost_vo(double *x, double *p) {
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
    double k = p[9]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // 4 paras: vw
    double fcbo = 2.34;
    double fc = 2*M_PI/0.1492;
    double fvo = sqrt(fcbo*(2*fc - fcbo));
    double fvw = fc - 2*fvo;

    double tau_vo = p[10];
    double asym_vo = p[11];
    double omega_vo = p[12]*fvo;
    double phi_vo = p[13];


    

    double cbo = 1-TMath::Exp(-time/tau_cbo)*asym_cbo*TMath::Cos(omega_cbo*time + phi_cbo);
    double vo  = vo = 1 - exp(-time/tau_vo)*asym_vo*cos(omega_vo*time + phi_vo);
    return  norm *(1 - k*aloss)* TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time + phi)) *  cbo * vo;
}



double func_18paras_cbo_lost_vo_vw(double *x, double *p) {
    double time = x[0] / time_scale;
    // 5-par
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double R = p[3];
    double phi = p[4];
    double omega = getBlinded->paramToFreq(R);
    
    // cbo-par
    double tau_cbo = p[5];
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];
    double cbo = 1-TMath::Exp(-time/tau_cbo)*asym_cbo*TMath::Cos(omega_cbo*time + phi_cbo);
    
    double fcbo = 2.32657;
    double fc = 2*M_PI/0.1492;
    double fvo = sqrt(fcbo*(2*fc - fcbo));
    double fvw = fc - 2*fvo;
    
    // k_loss
    double k = p[9]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    
    // vo-par
    double tau_vo = p[10];
    double asym_vo = p[11];
    double omega_vo = p[12]*fvo;
    double phi_vo = p[13];
    double vo = 1 - exp(-time/tau_vo)*asym_vo*cos(omega_vo*time + phi_vo);

    // vw-par
    double tau_vw = p[14];
    double asym_vw = p[15];
    double omega_vw = p[16]*fvw;
    double phi_vw = p[17];
    double vw = 1-exp(-time/tau_vw)*asym_vw*cos(omega_vw*time + phi_vw);

    return (1 - k*aloss) * norm * exp(-time/life) * (1 - asym*cos(omega*time + phi)) * cbo * vw * vo;
}


//1.9MHz func, 28 paras
double func_28paras_cbo_lost_vw_expansion(double *x, double *p) {
    double time = x[0] / time_scale;
    // 5-par
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double R = p[3];
    double phi = p[4];
    double omega = getBlinded->paramToFreq(R);
    
    // cbo-par
    double tau_cbo = p[5];
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];
    
    double fcbo = 2.32657;
    double fc = 2*M_PI/0.1492;
    double fvo = sqrt(fcbo*(2*fc - fcbo));
    double fvw = fc - 2*fvo;
    
    // k_loss
    double k = p[9]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // vw-par
    double tau_vw = p[10];
    double asym_vw = p[11];
    double omega_vw = p[12]*fvw;
    double phi_vw = p[13];

    // expansion-par
    double asym_vwcbo = p[24];
    double phi_vwcbo = p[25];
    double asym_vw_cbo = p[26];
    double phi_vw_cbo = p[27];
    double expan = (1 - exp(-time/tau_cbo)*asym_cbo*cos(omega_cbo*time + phi_cbo) - exp(-time/tau_vw)*asym_vw*cos(omega_vw*time + phi_vw) + exp(-time/tau_cbo - time/tau_vw)*(asym_vwcbo*cos((omega_vw + omega_cbo)*time + phi_vwcbo) + asym_vw_cbo*cos((omega_vw - omega_cbo)*time + phi_vw_cbo)));
    
    // dcbo-par
    double asym_dcbo = p[14];
    double phi_dcbo = p[15];
    double dcbo = 1 - exp(-2*time/tau_cbo)*asym_dcbo*cos(2*omega_cbo*time + phi_dcbo);
    
    // vo-par
    double tau_vo = p[20];
    double asym_vo = p[21];
    double omega_vo = p[22]*fvo;
    double phi_vo = p[23];
    double vo = 1 - exp(-time/tau_vo)*asym_vo*cos(omega_vo*time + phi_vo);
    
    // modification of A and phi
    double A1_cbo = p[16];
    double phi1_cbo = p[17];
    double A2_cbo = p[18];
    double phi2_cbo = p[19];
    double At = 1 - A1_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi1_cbo);
    double phit = 1 - A2_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi2_cbo);
    return (1 - k*aloss) * norm * exp(-time/life) * (1 - asym*At*cos(omega*time + phi*phit)) * expan * dcbo * vo;
}

//1.9MHz func, 28 paras
double func_22paras_cbo_lost_vw_expansion_lite(double *x, double *p) {
    
    double time = x[0] / time_scale;

    int nvar = 0;
    // 5-par
    double norm = p[nvar++];
    double life = p[nvar++];
    double asym = p[nvar++];
    double R = p[nvar++];
    double phi = p[nvar++];
    double omega = getBlinded->paramToFreq(R);
    
    // cbo-par
    double tau_cbo = p[nvar++];
    double asym_cbo = p[nvar++];
    double omega_cbo = p[nvar++];
    double phi_cbo = p[nvar++];
    
    double fcbo = 2.32657;
    double fc = 2*M_PI/0.1492;
    double fvo = sqrt(fcbo*(2*fc - fcbo));
    double fvw = fc - 2*fvo;
    
    // k_loss
    double k = p[nvar++]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // vw-par
    double tau_vw = p[nvar++];
    double asym_vw = p[nvar++];
    double omega_vw = p[nvar++]*fvw;
    double phi_vw = p[nvar++];

    // expansion-par
    // double asym_vwcbo = p[24];
    // double phi_vwcbo = p[25];
    double asym_vw_cbo = p[nvar++];
    double phi_vw_cbo = p[nvar++];
    double expan = (1 - exp(-time/tau_cbo)*asym_cbo*cos(omega_cbo*time + phi_cbo) - exp(-time/tau_vw)*asym_vw*cos(omega_vw*time + phi_vw) + exp(-time/tau_cbo - time/tau_vw)*(asym_vw_cbo*cos((omega_vw - omega_cbo)*time + phi_vw_cbo)));
    
    // // dcbo-par
    // // double asym_dcbo = p[14];
    // // double phi_dcbo = p[15];
    double dcbo = 1; // - exp(-2*time/tau_cbo)*asym_dcbo*cos(2*omega_cbo*time + phi_dcbo);
    
    // vo-par
    double tau_vo = p[nvar++];
    double asym_vo = p[nvar++];
    double omega_vo = p[nvar++]*fvo;
    double phi_vo = p[nvar++];
    double vo = 1 - exp(-time/tau_vo)*asym_vo*cos(omega_vo*time + phi_vo);
    
    // modification of A and phi
    double A1_cbo = p[nvar++];
    double phi1_cbo = p[nvar++];
    // double A2_cbo = p[18]; 
    // double phi2_cbo = p[19];
    double At = 1 - A1_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi1_cbo);
    double phit = 1; //- A2_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi2_cbo);

    return (1 - k*aloss) * norm * exp(-time/life) * (1 - asym*At*cos(omega*time + phi*phit)) * expan * dcbo * vo;
}


//frequency changing cbo per calo fit
// 5 + 4 cbo + 2 exp + 1 kloss
double func_12paras_changing_cbo(double *x, double *p) {
    double time = x[0] / time_scale;

    int nvar = 0;
    // 5-par
    double norm = p[nvar++];
    double life = p[nvar++];
    double asym = p[nvar++];
    double R = p[nvar++];
    double phi = p[nvar++];
    double omega = getBlinded->paramToFreq(R);
    
    // cbo-par
    double tau_cbo = p[nvar++];
    double asym_cbo = p[nvar++];
    double omega_cbo = p[nvar++];
    double phi_cbo = p[nvar++];
    

    // k_loss
    double k = p[nvar++]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // time chaging cbo terms
    double a1 = p[nvar++];
    double tau1 = p[nvar++];

    // double fcbo = 2.32657;
    // double fc = 2*M_PI/0.1492;
    // double fvo = sqrt(fcbo*(2*fc - fcbo));
    // double fvw = fc - 2*fvo;
    double cbo = 1-TMath::Exp(-time/tau_cbo)*asym_cbo*TMath::Cos(omega_cbo*time + a1*exp(-time/tau1)+ phi_cbo);
    return (1 - k*aloss) * norm * exp(-time/life) * (1 - asym*cos(omega*time + phi)) *cbo;
}

double func_11paras_changing_cbo(double *x, double *p) {
    double time = x[0] / time_scale;

    int nvar = 0;
    // 5-par
    double norm = p[nvar++];
    double life = p[nvar++];
    double asym = p[nvar++];
    double R = p[nvar++];
    double phi = p[nvar++];
    double omega = getBlinded->paramToFreq(R);
    
    // cbo-par
    double tau_cbo = p[nvar++];
    double asym_cbo = p[nvar++];
    double omega_cbo = p[nvar++];
    double phi_cbo = p[nvar++];
    

    // k_loss
    // double k = p[nvar++]*1e-9;
    // double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // time chaging cbo terms
    double a1 = p[nvar++];
    double tau1 = p[nvar++];

    // double fcbo = 2.32657;
    // double fc = 2*M_PI/0.1492;
    // double fvo = sqrt(fcbo*(2*fc - fcbo));
    // double fvw = fc - 2*fvo;
    double cbo = 1-TMath::Exp(-time/tau_cbo)*asym_cbo*TMath::Cos(omega_cbo*time + a1*exp(-time/tau1)+ phi_cbo);
    return norm * exp(-time/life) * (1 - asym*cos(omega*time + phi)) *cbo;
}