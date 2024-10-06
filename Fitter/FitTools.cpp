#include "FitTools.h"
#include <sstream>
// #include "new_functions.h"
#include "functions_run456_new.h"
#include "TMath.h"
#include <Math/MinimizerOptions.h>
#include "TVirtualFitter.h"

using namespace pocketfft;

vector<double> GetInitialValuesD(string file_path, string index){
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
  TFile* file = TFile::Open(file_path.c_str());  
  TObject* obj = file->Get(index.c_str());
  if(!obj->InheritsFrom(TF1::Class())){
    cout << "load parameters from Root::TArrayD    " << file_path << endl;
    TArrayD* arrayd = (TArrayD*)obj;
    for(int n=0; n<arrayd->fN; n++){
        // cout << "var" << n << " = "<< arrayd->At(n) << endl;
        init_values.push_back(arrayd->At(n));
    }
  }
  else{
    cout << "load parameters from Root::TF1    " << file_path << endl;
    TF1* func = (TF1*)obj;
    for(int n=0; n<func->GetNpar(); n++){
        init_values.push_back(func->GetParameter(n));
    }    
  }
#endif  
  return init_values;
}


void FillData(TH1* th1, vector<double> & data){  
  int nbins = th1->GetNbinsX();
  for(int n=0;n<nbins;n++){
    data[n] = th1->GetBinContent(n+1);
  }  
}

void FillHist(TH1* th1,const vector<double> & data,bool isAbs){
  int nbins = th1->GetNbinsX();
  for(int n=0; n<nbins; n++){
    if(isAbs) {
      th1->SetBinContent(n+1, fabs(data[n]));  
    }
    else{
      th1->SetBinContent(n+1, data[n]);   
    }    
  }
}

void FFT(TH1* hist, TH1* hist_fft, bool isAbs){
  size_t nbins = hist->GetNbinsX();
  vector<double> data_in(nbins);
  vector<double> data_out(nbins);
  FillData(hist, data_in);
  detail::shape_t shape{nbins};
  stride_t strided{sizeof(double)};
  shape_t axes{0}; 
  r2r_separable_hartley(shape, strided, strided, axes, data_in.data(), data_out.data(), 1.);
  FillHist(hist_fft, data_out, isAbs);
}  


void Fitter::SetTimeUnit(TimeUnit t_unit){
    if(t_unit == Fitter::nano_second){
        time_scale = 1e3;
    } else if(t_unit == Fitter::micro_second){
        time_scale = 1.0;
    } else{
        cout << "Unknown time unit" << endl;
        exit(1);
    }
    cout << "Time scale = " << time_scale << endl;
}

void Fitter::SetBlindedString(string _which_run){
    if(_which_run == "run4"){
        blindedString = "Unexpected virtue of ignorance.";
    } else if(_which_run == "run5"){
        blindedString = "No bold guesses, no great discoveries.";
    } else if(_which_run == "run6"){
        blindedString = "Bad times make a good man.";
    } else if(_which_run == "run456"){
        blindedString = "Sow nothing, reap nothing.";
    } else if(_which_run == "run2" || _which_run == "run3"){
        blindedString = "Random blind!";
    } else{
        cout << "Unknown run, please use run2, run3, run4, run5, run6 or run456!" << endl;
        exit(1);
    }
    getBlinded = new Blinders(ftype, blindedString.c_str());
    cout << "############# " << _which_run << " is now blinded! #############" << endl;
}

void Fitter::SetOutputDir(string _output_dir){
    output_dir = _output_dir;
}

FitOutputInfo Fitter::doFit(const FitInput & fit_in){
    TVirtualFitter::SetMaxIterations(100000);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
    TString wiggle_name, residual_name, func_name, res_name, fft_name, outf_name;
    wiggle_name.Form("wiggle_%s", fit_in.tag.Data());
    residual_name.Form("residual_%s", fit_in.tag.Data());
    func_name.Form("func_%s", fit_in.tag.Data());
    res_name.Form("result_%s", fit_in.tag.Data());
    fft_name.Form("fft_%s", fit_in.tag.Data());
    outf_name.Form("%s/result_%s.root", output_dir.c_str(), fit_in.tag.Data());


    lost_muon = fit_in.lost_muon;
    splines = fit_in.splines;
    // itp = fit_in.itp;

    TH1 * fit_hist = (TH1*)fit_in.wiggle->Clone(wiggle_name.Data());
    TF1 * fit_func = new TF1(func_name, fit_in.func, fit_in.t_start, fit_in.t_end, fit_in.nvars);
    

    for(int n=0; n<fit_in.nvars; n++){        
        fit_func->SetParName(n, fit_in.name_vars[n].c_str());
        if(n<fit_in.init_values.size()){
            fit_func->SetParameter(n, fit_in.init_values[n]);
        }
    }
    
    for(auto fix : fix_parameters){
        if(fix.first<fit_in.nvars){
            fit_func->FixParameter(fix.first,fix.second);
        }
    }

    for(auto range : range_parameters){
        if(range.first<fit_in.nvars){
            fit_func->SetParLimits(range.first, range.second.first, range.second.second);
        }
    }

    fit_func->SetNpx(45000);

        
    Int_t fitStatus;    
    bool isValid = false;
    double chi2, ndf;

    TFitResultPtr fit_res;

    cout << "Performing fit histogram with option '" << fit_option << "'" << endl;
    // fit_res = fit_hist->Fit(fit_func,"REMS");
    fit_res = fit_hist->Fit(fit_func, this->fit_option);


    fitStatus = fit_res;
    isValid = fit_res->IsValid();// && (fitStatus%100==0);
    chi2 = fit_func->GetChisquare();
    ndf = fit_func->GetNDF();

    cout << "Chi2/NDF = " << chi2/ndf << "    Valid = " << isValid << "     Status = " << fitStatus << endl;

    int refit = 1;
    int max_attempts_ = max_attempts;
    while(max_attempts_>0 && !isValid){
     cout << "Invalid fitting!!! Retry " << refit++ << endl;
     fit_res = fit_hist->Fit(fit_func, this->fit_option);

     fitStatus = fit_res;     
     isValid = fit_res->IsValid();// && (fitStatus%100==0);
     chi2 = fit_func->GetChisquare();
     ndf = fit_func->GetNDF();

     cout << "Chi2/NDF = " << chi2/ndf << "    Valid = " << isValid << "     Status = "<< fitStatus << endl;
     max_attempts_--;
    }

    cout << "\nResult-" << fit_in.tag <<  "\tChi2/NDF = " << chi2/ndf << "\tValid = "<< isValid << "\tStatus = "<< fitStatus
    << "\tR = " << fit_func->GetParameter(3) << "+-" << fit_func->GetParError(3) << endl;

    vector<double> fit_values = fit_in.init_values;
    for(int n=0; n<fit_in.nvars; n++){
        if(n<fit_values.size()){
            fit_values[n] = fit_func->GetParameter(n);
        }
        else{
            fit_values.push_back(fit_func->GetParameter(n));
        }
    }

    fit_res->SetName(res_name.Data());
    TH1* residual = (TH1*)fit_in.wiggle->Clone(residual_name.Data());
    residual->Reset();
    int bin_t_start = residual->GetXaxis()->FindBin(fit_in.t_start) + 1;
    int bin_t_end = residual->GetXaxis()->FindBin(fit_in.t_end) - 1;
    for(int bin=bin_t_start; bin<=bin_t_end; bin++){
        double time = residual->GetBinCenter(bin);
        double val = fit_hist->GetBinContent(bin);
        double err = fit_hist->GetBinError(bin);
        double fval = fit_func->Eval(time);
        residual->SetBinContent(bin,val - fval);
        residual->SetBinError(bin,err);
    }

    TH1D * fft = new TH1D(fft_name.Data(), fft_name.Data(), residual->GetNbinsX(), 0,1./0.1492);
    FFT(residual, fft);


    TFile * fout = new TFile(outf_name.Data(), "recreate");
    fout->cd();
    fit_hist->Write();
    fit_func->Write();
    fit_res->Write();
    residual->Write();
    fft->Write();
    fout->Close();



    cout << "\nOutput file created : "<< outf_name << endl;
    cout << "   wiggle    : " << wiggle_name << endl;
    cout << "   residual  : " << residual_name << endl;
    cout << "   function  : " << func_name << endl;
    cout << "   FFT       : " << fft_name << endl;
    cout << "   FitResults: " << res_name << endl;

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

Fitter::Fitter() : max_attempts(1) {
    cout << "Fitting on Run" << data_version_major << " data" << endl;

    name_vars["5paras"] = {"N", "#tau", "A", "R", "#phi"};
    name_vars["9paras"] = {
        "N", "#tau", "A", "R", "#phi", 
        "#tau_{cbo}", "#alpha_{cbo}", "#beta_{cbo}", "#omega_{cbo}"
        };
    name_vars["10paras"] = {
        "N", "#tau", "A", "R", "#phi", 
        "#tau_{cbo}", "#alpha_{cbo}", "#beta_{cbo}", "#omega_{cbo}", 
        "k_{loss}"
        };
    name_vars["28paras"] = {
        "N", "#tau", "A", "R", "#phi", 
        "#tau_{cbo}", "#alpha_{cbo}", "#beta_{cbo}", "#omega_{cbo}", 
        "k_{loss}", 
        "#tau_{vw}", "#alpha_{vw}", "#beta_{vw}", "#omega_{vw}", 
        "#tau_{vo}", "#alpha_{vo}", "#beta_{vo}", "#omega_{vo}", 
        "#alpha_{2cbo}", "#beta_{2cbo}",
        "#alpha_{cbo,A}", "#beta_{cbo,A}", "#alpha_{cbo,#phi}", "#beta_{cbo,#phi}",
        "#alpha_{vw+cbo}", "#beta_{vw+cbo}", "#alpha_{vw-cbo}", "#beta_{vw-cbo}"
        };
    name_vars["29paras"] = {
        "N", "#tau", "A", "R", "#phi", 
        "#tau_{cbo}", "#alpha_{cbo}", "#beta_{cbo}", "#omega_{cbo}", 
        "k_{loss}", 
        "#tau_{vw}", "#alpha_{vw}", "#beta_{vw}", "#omega_{vw}", 
        "#tau_{vo}", "#alpha_{vo}", "#beta_{vo}", "#omega_{vo}", 
        "#alpha_{2cbo}", "#beta_{2cbo}",
        "#alpha_{cbo,A}", "#beta_{cbo,A}", "#alpha_{cbo,#phi}", "#beta_{cbo,#phi}",
        "#alpha_{vw+cbo}", "#beta_{vw+cbo}", "#alpha_{vw-cbo}", "#beta_{vw-cbo}",
        "c_{cbo}"
        };

    name_vars["40paras"] = {
        "N", "#tau", "A", "R", "#phi",
        "#tau_{cbo}", "#alpha_{cbo}", "#beta_{cbo}", "#omega_{cbo}",
        "k_{loss}",
        "#tau_{vw}", "#alpha_{vw}", "#beta_{vw}", "#omega_{vw}",
        "#tau_{vo}", "#alpha_{vo}", "#beta_{vo}", "#omega_{vo}",
        "#alpha_{2cbo}", "#beta_{2cbo}",
        "#alpha_{cbo+a}", "#beta_{cbo+a}", "#alpha_{cbo-a}", "#beta_{cbo-a}",
        "#alpha_{vw+a}", "#beta_{vw+a}", "#alpha_{vw-a}", "#beta_{vw-a}",
        "#alpha_{vo+a}", "#beta_{vo+a}", "#alpha_{vo-a}", "#beta_{vo-a}",
        "#alpha_{vw+cbo}", "#beta_{vw+cbo}", "#alpha_{vw-cbo}", "#beta_{vw-cbo}",
        "#alpha_{vo+cbo}", "#beta_{vo+cbo}", "#alpha_{vo-cbo}", "#beta_{vo-cbo}"
        };
    name_vars["41paras"] = {
        "N", "#tau", "A", "R", "#phi",
        "#tau_{cbo}", "#alpha_{cbo}", "#beta_{cbo}", "#omega_{cbo}",
        "k_{loss}",
        "#tau_{vw}", "#alpha_{vw}", "#beta_{vw}", "#omega_{vw}",
        "#tau_{vo}", "#alpha_{vo}", "#beta_{vo}", "#omega_{vo}",
        "#alpha_{2cbo}", "#beta_{2cbo}",
        "#alpha_{cbo+a}", "#beta_{cbo+a}", "#alpha_{cbo-a}", "#beta_{cbo-a}",
        "#alpha_{vw+a}", "#beta_{vw+a}", "#alpha_{vw-a}", "#beta_{vw-a}",
        "#alpha_{vo+a}", "#beta_{vo+a}", "#alpha_{vo-a}", "#beta_{vo-a}",
        "#alpha_{vw+cbo}", "#beta_{vw+cbo}", "#alpha_{vw-cbo}", "#beta_{vw-cbo}",
        "#alpha_{vo+cbo}", "#beta_{vo+cbo}", "#alpha_{vo-cbo}", "#beta_{vo-cbo}",
        "c_{cbo}"
        };

    name_vars["42paras"] = {
        "N", "#tau", "A", "R", "#phi",
        "#tau_{cbo}", "#alpha_{cbo}", "#beta_{cbo}", "#omega_{cbo}",
        "k_{loss}",
        "#tau_{vw}", "#alpha_{vw}", "#beta_{vw}", "#omega_{vw}",
        "#tau_{vo}", "#alpha_{vo}", "#beta_{vo}", "#omega_{vo}",
        "#alpha_{2cbo}", "#beta_{2cbo}",
        "#tau_{vw-cbo}", "#alpha_{vw-cbo}", "#beta_{vw-cbo}", "#omega_{vw-cbo}",
        "#alpha_{cbo+a}", "#beta_{cbo+a}", "#alpha_{cbo-a}", "#beta_{cbo-a}",
        "#alpha_{vw+a}", "#beta_{vw+a}", "#alpha_{vw-a}", "#beta_{vw-a}",
        "#alpha_{vo+a}", "#beta_{vo+a}", "#alpha_{vo-a}", "#beta_{vo-a}",
        "#alpha_{vw-cbo-a}", "#beta_{vw-cbo-a}",
        "#alpha_{vo+cbo}", "#beta_{vo+cbo}", "#alpha_{vo-cbo}", "#beta_{vo-cbo}"
        };
    name_vars["43paras"] = {
        "N", "#tau", "A", "R", "#phi",
        "#tau_{cbo}", "#alpha_{cbo}", "#beta_{cbo}", "#omega_{cbo}",
        "k_{loss}",
        "#tau_{vw}", "#alpha_{vw}", "#beta_{vw}", "#omega_{vw}",
        "#tau_{vo}", "#alpha_{vo}", "#beta_{vo}", "#omega_{vo}",
        "#alpha_{2cbo}", "#beta_{2cbo}",
        "#tau_{vw-cbo}", "#alpha_{vw-cbo}", "#beta_{vw-cbo}", "#omega_{vw-cbo}",
        "#alpha_{cbo+a}", "#beta_{cbo+a}", "#alpha_{cbo-a}", "#beta_{cbo-a}",
        "#alpha_{vw+a}", "#beta_{vw+a}", "#alpha_{vw-a}", "#beta_{vw-a}",
        "#alpha_{vo+a}", "#beta_{vo+a}", "#alpha_{vo-a}", "#beta_{vo-a}",
        "#alpha_{vw-cbo-a}", "#beta_{vw-cbo-a}",
        "#alpha_{vo+cbo}", "#beta_{vo+cbo}", "#alpha_{vo-cbo}", "#beta_{vo-cbo}",
        "c_{cbo}"
    };

    name_vars["44paras"] = {
        "N", "#tau", "A", "R", "#phi",
        "#tau_{cbo}", "#alpha_{cbo}", "#beta_{cbo}", "#omega_{cbo}",
        "k_{loss}",
        "#tau_{vw}", "#alpha_{vw}", "#beta_{vw}", "#omega_{vw}",
        "#tau_{vo}", "#alpha_{vo}", "#beta_{vo}", "#omega_{vo}",
        "#alpha_{2cbo}", "#beta_{2cbo}",
        "#tau_{vw-cbo}", "#alpha_{vw-cbo}", "#beta_{vw-cbo}", "#omega_{vw-cbo}",
        "#alpha_{cbo+a}", "#beta_{cbo+a}", "#alpha_{cbo-a}", "#beta_{cbo-a}",
        "#alpha_{vw+a}", "#beta_{vw+a}", "#alpha_{vw-a}", "#beta_{vw-a}",
        "#alpha_{vo+a}", "#beta_{vo+a}", "#alpha_{vo-a}", "#beta_{vo-a}",
        "#alpha_{vw-cbo+a}", "#beta_{vw-cbo+a}", "#alpha_{vw-cbo-a}", "#beta_{vw-cbo-a}",
        "#alpha_{vo+cbo}", "#beta_{vo+cbo}", "#alpha_{vo-cbo}", "#beta_{vo-cbo}"
        };
    name_vars["45paras"] = {
        "N", "#tau", "A", "R", "#phi",
        "#tau_{cbo}", "#alpha_{cbo}", "#beta_{cbo}", "#omega_{cbo}",
        "k_{loss}",
        "#tau_{vw}", "#alpha_{vw}", "#beta_{vw}", "#omega_{vw}",
        "#tau_{vo}", "#alpha_{vo}", "#beta_{vo}", "#omega_{vo}",
        "#alpha_{2cbo}", "#beta_{2cbo}",
        "#tau_{vw-cbo}", "#alpha_{vw-cbo}", "#beta_{vw-cbo}", "#omega_{vw-cbo}",
        "#alpha_{cbo+a}", "#beta_{cbo+a}", "#alpha_{cbo-a}", "#beta_{cbo-a}",
        "#alpha_{vw+a}", "#beta_{vw+a}", "#alpha_{vw-a}", "#beta_{vw-a}",
        "#alpha_{vo+a}", "#beta_{vo+a}", "#alpha_{vo-a}", "#beta_{vo-a}",
        "#alpha_{vw-cbo+a}", "#beta_{vw-cbo+a}", "#alpha_{vw-cbo-a}", "#beta_{vw-cbo-a}",
        "#alpha_{vo+cbo}", "#beta_{vo+cbo}", "#alpha_{vo-cbo}", "#beta_{vo-cbo}",
        "c_{cbo}"
        };

    name_vars["48paras"] = {
        "N", "#tau", "A", "R", "#phi",
        "#tau_{cbo}", "#alpha_{cbo}", "#beta_{cbo}", "#omega_{cbo}",
        "k_{loss}",
        "#tau_{vw}", "#alpha_{vw}", "#beta_{vw}", "#omega_{vw}",
        "#tau_{vo}", "#alpha_{vo}", "#beta_{vo}", "#omega_{vo}",
        "#alpha_{2cbo}", "#beta_{2cbo}",
        "#tau_{vw-cbo}", "#alpha_{vw-cbo}", "#beta_{vw-cbo}", "#omega_{vw-cbo}",
        "#tau_{vw+cbo}", "#alpha_{vw+cbo}", "#beta_{vw+cbo}", "#omega_{vw+cbo}",
        "#alpha_{cbo+a}", "#beta_{cbo+a}", "#alpha_{cbo-a}", "#beta_{cbo-a}",
        "#alpha_{vw+a}", "#beta_{vw+a}", "#alpha_{vw-a}", "#beta_{vw-a}",
        "#alpha_{vo+a}", "#beta_{vo+a}", "#alpha_{vo-a}", "#beta_{vo-a}",
        "#alpha_{vw-cbo+a}", "#beta_{vw-cbo+a}", "#alpha_{vw-cbo-a}", "#beta_{vw-cbo-a}",
        "#alpha_{vo+cbo}", "#beta_{vo+cbo}", "#alpha_{vo-cbo}", "#beta_{vo-cbo}"
        };
    name_vars["49paras"] = {
        "N", "#tau", "A", "R", "#phi",
        "#tau_{cbo}", "#alpha_{cbo}", "#beta_{cbo}", "#omega_{cbo}",
        "k_{loss}",
        "#tau_{vw}", "#alpha_{vw}", "#beta_{vw}", "#omega_{vw}",
        "#tau_{vo}", "#alpha_{vo}", "#beta_{vo}", "#omega_{vo}",
        "#alpha_{2cbo}", "#beta_{2cbo}",
        "#tau_{vw-cbo}", "#alpha_{vw-cbo}", "#beta_{vw-cbo}", "#omega_{vw-cbo}",
        "#tau_{vw+cbo}", "#alpha_{vw+cbo}", "#beta_{vw+cbo}", "#omega_{vw+cbo}",
        "#alpha_{cbo+a}", "#beta_{cbo+a}", "#alpha_{cbo-a}", "#beta_{cbo-a}",
        "#alpha_{vw+a}", "#beta_{vw+a}", "#alpha_{vw-a}", "#beta_{vw-a}",
        "#alpha_{vo+a}", "#beta_{vo+a}", "#alpha_{vo-a}", "#beta_{vo-a}",
        "#alpha_{vw-cbo+a}", "#beta_{vw-cbo+a}", "#alpha_{vw-cbo-a}", "#beta_{vw-cbo-a}",
        "#alpha_{vo+cbo}", "#beta_{vo+cbo}", "#alpha_{vo-cbo}", "#beta_{vo-cbo}",
        "c_{cbo}"
        };

    name_vars["42parasGPR"] = {
        "N", "#tau", "A", "R", "#phi",
        "#tau_{cbo}", "#omega_{cbo}",
        "k_{loss}",
        "#tau_{vw}", "#alpha_{vw}", "#beta_{vw}", "#omega_{vw}",
        "#tau_{vo}", "#alpha_{vo}", "#beta_{vo}", "#omega_{vo}",
        "#alpha_{2cbo}", "#beta_{2cbo}",
        "#tau_{vw-cbo}", "#alpha_{vw-cbo}", "#beta_{vw-cbo}", "#omega_{vw-cbo}",
        "#alpha_{cbo+a}", "#beta_{cbo+a}", "#alpha_{cbo-a}", "#beta_{cbo-a}",
        "#alpha_{vw+a}", "#beta_{vw+a}", "#alpha_{vw-a}", "#beta_{vw-a}",
        "#alpha_{vo+a}", "#beta_{vo+a}", "#alpha_{vo-a}", "#beta_{vo-a}",
        "#alpha_{vw-cbo+a}", "#beta_{vw-cbo+a}", "#alpha_{vw-cbo-a}", "#beta_{vw-cbo-a}",
        "#alpha_{vo+cbo}", "#beta_{vo+cbo}", "#alpha_{vo-cbo}", "#beta_{vo-cbo}"
        };
    name_vars["40parasGPR"] = {
        "N", "#tau", "A", "R", "#phi",
        "#tau_{2cbo}", "#omega_{cbo}",
        "k_{loss}",
        "#tau_{vw}", "#alpha_{vw}", "#beta_{vw}", "#omega_{vw}",
        "#tau_{vo}", "#alpha_{vo}", "#beta_{vo}", "#omega_{vo}",
        "#alpha_{2cbo}", "#beta_{2cbo}",
        "#tau_{vw-cbo}", "#alpha_{vw-cbo}", "#beta_{vw-cbo}", "#omega_{vw-cbo}",
        "#alpha_{vw+a}", "#beta_{vw+a}", "#alpha_{vw-a}", "#beta_{vw-a}",
        "#alpha_{vo+a}", "#beta_{vo+a}", "#alpha_{vo-a}", "#beta_{vo-a}",
        "#alpha_{vw-cbo+a}", "#beta_{vw-cbo+a}", "#alpha_{vw-cbo-a}", "#beta_{vw-cbo-a}",
        "#tau_{vo+cbo}", "#alpha_{vo+cbo}", "#beta_{vo+cbo}", "#tau_{vo-cbo}", "#alpha_{vo-cbo}", "#beta_{vo-cbo}"
        };
    name_vars["44parasGPR"] = {
        "N", "#tau", "A", "R", "#phi",
        "#tau_{cbo}", "#omega_{cbo}",
        "k_{loss}",
        "#tau_{vw}", "#alpha_{vw}", "#beta_{vw}", "#omega_{vw}",
        "#tau_{vo}", "#alpha_{vo}", "#beta_{vo}", "#omega_{vo}",
        "#alpha_{2cbo}", "#beta_{2cbo}",
        "#tau_{vw-cbo}", "#alpha_{vw-cbo}", "#beta_{vw-cbo}", "#omega_{vw-cbo}",
        "#alpha_{cbo+a}", "#beta_{cbo+a}", "#alpha_{cbo-a}", "#beta_{cbo-a}",
        "#alpha_{vw+a}", "#beta_{vw+a}", "#alpha_{vw-a}", "#beta_{vw-a}",
        "#alpha_{vo+a}", "#beta_{vo+a}", "#alpha_{vo-a}", "#beta_{vo-a}",
        "#alpha_{vw-cbo+a}", "#beta_{vw-cbo+a}", "#alpha_{vw-cbo-a}", "#beta_{vw-cbo-a}",
        "#alpha_{vo+cbo}", "#beta_{vo+cbo}", "#alpha_{vo-cbo}", "#beta_{vo-cbo}",
        "a_{cbo}^{1}", "b_{cbo}^{1}"
        };
    name_vars["46parasGPR"] = {
        "N", "#tau", "A", "R", "#phi",
        "#tau_{2cbo}", "#omega_{cbo}",
        "k_{loss}",
        "#tau_{vw}", "#alpha_{vw}", "#beta_{vw}", "#omega_{vw}",
        "#tau_{vo}", "#alpha_{vo}", "#beta_{vo}", "#omega_{vo}",
        "#alpha_{2cbo}", "#beta_{2cbo}",
        "#tau_{vw-cbo}", "#alpha_{vw-cbo}", "#beta_{vw-cbo}", "#omega_{vw-cbo}",
        "#alpha_{vw+a}", "#beta_{vw+a}", "#alpha_{vw-a}", "#beta_{vw-a}",
        "#alpha_{vo+a}", "#beta_{vo+a}", "#alpha_{vo-a}", "#beta_{vo-a}",
        "#alpha_{vw-cbo+a}", "#beta_{vw-cbo+a}", "#alpha_{vw-cbo-a}", "#beta_{vw-cbo-a}",
        "#tau_{vo+cbo}", "#alpha_{vo+cbo}", "#beta_{vo+cbo}", "#tau_{vo-cbo}", "#alpha_{vo-cbo}", "#beta_{vo-cbo}",
        "a_{cbo}", "b_{cbo}", "a_{cbopa}", "b_{cbopa}", "a_{cboma}", "b_{cboma}"
        };
    name_vars["50parasGPR"] = {
        "N", "#tau", "A", "R", "#phi",
        "#tau_{cbo}", "#omega_{cbo}",
        "k_{loss}",
        "#tau_{vw}", "#alpha_{vw}", "#beta_{vw}", "#omega_{vw}",
        "#tau_{vo}", "#alpha_{vo}", "#beta_{vo}", "#omega_{vo}",
        "#alpha_{2cbo}", "#beta_{2cbo}",
        "#tau_{vw-cbo}", "#alpha_{vw-cbo}", "#beta_{vw-cbo}", "#omega_{vw-cbo}",
        "#tau_{vw+cbo}", "#alpha_{vw+cbo}", "#beta_{vw+cbo}", "#omega_{vw+cbo}",
        "#alpha_{cbo+a}", "#beta_{cbo+a}", "#alpha_{cbo-a}", "#beta_{cbo-a}",
        "#alpha_{vw+a}", "#beta_{vw+a}", "#alpha_{vw-a}", "#beta_{vw-a}",
        "#alpha_{vo+a}", "#beta_{vo+a}", "#alpha_{vo-a}", "#beta_{vo-a}",
        "#alpha_{vw-cbo+a}", "#beta_{vw-cbo+a}", "#alpha_{vw-cbo-a}", "#beta_{vw-cbo-a}",
        "#alpha_{vo+cbo}", "#beta_{vo+cbo}", "#alpha_{vo-cbo}", "#beta_{vo-cbo}",
        "a_{cbo}^{1}", "a_{cbo}^{2}", "b_{cbo}^{1}", "b_{cbo}^{2}"
        };
    name_vars["2paras"] = {
        "N", "#tau"
        };
    name_vars["3paras_kloss"] = {
        "N", "#tau",
        "k_{loss}"
        };
    name_vars["6paras_kloss"] = {
        "N", "#tau", "A", "R", "#phi", 
        "k_{loss}"
        };
    name_vars["11paras_cbo"] = {
        "N", "#tau", "A", "R", "#phi", 
        "#tau_{cbo}", "c_{cbo}", "A_{cbo}^{c}", "A_{cbo}^{s}", "#omega_{cbo}", 
        "k_{loss}"
        };
    name_vars["13paras_res"] = {
        "N", "#tau", "A", "R", "#phi", 
        "#tau_{cbo}", "c_{cbo}", "A_{cbo}^{c}", "A_{cbo}^{s}", "#omega_{cbo}", 
        "k_{loss}", "N_{res}", "#tau_{res}"
        };
    
}