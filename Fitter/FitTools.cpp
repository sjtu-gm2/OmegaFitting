#include "FitTools.h"
#include <sstream>
#include "functions.h"
#include "TMath.h"
#include <Math/MinimizerOptions.h>


using namespace pocketfft;


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


void Fitter::SetTimeUnit(TimeUnit t_unit) {
    if(t_unit == Fitter::nano_second) {
        time_scale = 1e3;
    } else if(t_unit == Fitter::micro_second){
        time_scale = 1.0;
    } else {
        cout <<"unknown time unit" << endl;
        exit(1);
    }
    cout << "time scale = " << time_scale << endl;
}

void Fitter::SetOutputDir(string output_dir) {
    m_output_dir = output_dir;
}

FitOutputInfo Fitter::doFit(const FitInput & fit_in) {
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);
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
    
    for(auto fix : fix_parameters) {
        if(fix.first<fit_in.nvars) {
            fit_func->FixParameter(fix.first,fix.second);
        }
    }

    for(auto range : range_parameters) {
        if(range.first<fit_in.nvars) {
            fit_func->SetParLimits(range.first, range.second.first, range.second.second);
        }
    }

    fit_func->SetNpx(4500);

        
    Int_t fitStatus;    
    bool isValid = false;
    double chi2,ndf;

    TFitResultPtr fit_res;

    fit_res = fit_hist->Fit(fit_func,"REMS");



    fitStatus = fit_res;
    isValid = fit_res->IsValid();// && (fitStatus%100==0);
    chi2 = fit_func->GetChisquare();
    ndf = fit_func->GetNDF();

    cout << "Chi2/NDF = " << chi2/ndf << "    Valid="<< isValid << "     Status="<< fitStatus << endl;

    int refit = 1;
    int max_attempts_ = max_attempts;
    while(max_attempts_>0 && !isValid) {
     cout << "Invalid fitting!!! Retry " << refit++ << endl;
     fit_res = fit_hist->Fit(fit_func,"QREMS");

     fitStatus = fit_res;     
     isValid = fit_res->IsValid();// && (fitStatus%100==0);
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

Fitter::Fitter() : max_attempts(1) {
    cout << "Fitting on Run" << data_version_major << " data" << endl;
    name_vars["5paras"] = {"N","#tau","A","R","#phi"};
    name_vars["9paras_cbo"] = {"N","#tau","A","R","#phi","#tau_{cbo}","A_{cbo}","#omega_{cbo}","#phi_{cbo}"};
    name_vars["10paras_cbo_lost"] = {
        "N","#tau","A","R","#phi",
        "#tau_{cbo}","A_{cbo}","#omega_{cbo}","#phi_{cbo}",
        "k_{loss}",
        "#tau_{vw}","A_{vw}","K_{vw}","#phi_{vw}"
    };
    name_vars["14paras_cbo_lost_vo"] = {
        "N","#tau","A","R","#phi",
        "#tau_{cbo}","A_{cbo}","#omega_{cbo}","#phi_{cbo}",
        "k_{loss}",
        "#tau_{y}","A_{y}","K_{y}","#phi_{y}",
    };

    name_vars["18paras_cbo_lost_vo_vw"] = {
        "N","#tau","A","R","#phi",
        "#tau_{cbo}","A^{N}_{cbo}","#omega_{cbo}","#phi_{cbo}",
        "k_{loss}",
        "#tau_{y}","A_{y}","K_{y}","#phi_{y}",
        "#tau_{vw}","A_{vw}","K_{vw}","#phi_{vw}",        
    };

    name_vars["28paras_cbo_lost_vw_expansion"] = {
        "N_{0}","#tau","A","R","#phi_{0}",
        "#tau_{cbo}","A^{N}_{cbo}","#omega_{cbo}","#phi_{cbo}",        
        "k_{loss}",
        "#tau_{vw}","A_{vw}","K_{vw}","#phi_{vw}",
        "A_{2cbo}","#phi_{2cbo}","A_{cbo,A}","#phi_{cbo,A}","A_{cbo,#phi}","#phi_{cbo,#phi}",
        "#tau_{y}","A_{y}","K_{y}","#phi_{y}",
        "A_{VW+cbo}","#phi_{VW+cbo}","A_{VW-cbo}","#phi_{VW-cbo}",
    };

    name_vars["22paras_cbo_lost_vw_expansion_lite"] = {
        "N_{0}","#tau","A","R","#phi_{0}",
        "#tau_{cbo}","A^{N}_{cbo}","#omega_{cbo}","#phi_{cbo}",        
        "k_{loss}",
        "#tau_{vw}","A_{vw}","K_{vw}","#phi_{vw}",
        "A_{VW-cbo}","#phi_{VW-cbo}",
        "#tau_{y}","A_{y}","K_{y}","#phi_{y}",
        "A_{cbo,A}","#phi_{cbo,A}",
    };

    name_vars["12paras_changing_cbo"] = {
        "N_{0}","#tau","A","R","#phi_{0}",
        "#tau_{cbo}","A_{cbo}","#omega_{cbo}","#phi_{cbo}",
        "k_{loss}"
        "A_{1}", "#tau_{1}",
    };

    name_vars["11paras_changing_cbo"] = {
        "N_{0}","#tau","A","R","#phi_{0}",
        "#tau_{cbo}","A_{cbo}","#omega_{cbo}","#phi_{cbo}",
        "A_{1}", "#tau_{1}",        
    };

    name_vars["13paras_cbo_vo"] = {
        "N","#tau","A","R","#phi",
        "#tau_{cbo}","A_{cbo}","#omega_{cbo}","#phi_{cbo}",        
        "#tau_{y}","A_{y}","K_{y}","#phi_{y}",
    };

    name_vars["15paras_changing_cbo_vo"] = {
        "N","#tau","A","R","#phi",
        "#tau_{cbo}","A_{cbo}","#omega_{cbo}","#phi_{cbo}",        
        "#tau_{y}","A_{y}","K_{y}","#phi_{y}",
        "A_{1}", "#tau_{1}",
    };
    name_vars["29paras_cbo_envelope_C"] = {
        "N_{0}","#tau","A","R","#phi_{0}",
        "#tau_{cbo}","A^{N}_{cbo}","#omega_{cbo}","#phi_{cbo}",        
        "k_{loss}",
        "#tau_{vw}","A_{vw}","K_{vw}","#phi_{vw}",
        "A_{2cbo}","#phi_{2cbo}","A_{cbo,A}","#phi_{cbo,A}","A_{cbo,#phi}","#phi_{cbo,#phi}",
        "#tau_{y}","A_{y}","K_{y}","#phi_{y}",
        "A_{VW+cbo}","#phi_{VW+cbo}","A_{VW-cbo}","#phi_{VW-cbo}",
        "C_{cbo}",
    };
    name_vars["30paras_cbo_freq"] = {
        "N_{0}","#tau","A","R","#phi_{0}",
        "#tau_{cbo}","A^{N}_{cbo}","#omega_{cbo}","#phi_{cbo}",        
        "k_{loss}",
        "#tau_{vw}","A_{vw}","K_{vw}","#phi_{vw}",
        "A_{2cbo}","#phi_{2cbo}","A_{cbo,A}","#phi_{cbo,A}","A_{cbo,#phi}","#phi_{cbo,#phi}",
        "#tau_{y}","A_{y}","K_{y}","#phi_{y}",
        "A_{VW+cbo}","#phi_{VW+cbo}","A_{VW-cbo}","#phi_{VW-cbo}",
        "A_{cbo}^{residual}","#tau_{cbo}^{residual}"
    };
    name_vars["31paras_cbo_time"] = {
        "N_{0}","#tau","A","R","#phi_{0}",
        "#tau_{cbo}","A^{N}_{cbo}","#omega_{cbo}","#phi_{cbo}",        
        "k_{loss}",
        "#tau_{vw}","A_{vw}","K_{vw}","#phi_{vw}",
        "A_{2cbo}","#phi_{2cbo}","A_{cbo,A}","#phi_{cbo,A}","A_{cbo,#phi}","#phi_{cbo,#phi}",
        "#tau_{y}","A_{y}","K_{y}","#phi_{y}",
        "A_{VW+cbo}","#phi_{VW+cbo}","A_{VW-cbo}","#phi_{VW-cbo}",
        "#tau_{cbo}^{a}","#tau_{cbo}^{phi}", "#tau_{cbo}^{2a}",
    };

    name_vars["simplified_22paras_calos"] = {
        "N_{0}","#tau","A","R","#phi_{0}",
        "#tau_{cbo}","A^{N}_{cbo}","#omega_{cbo}","#phi_{cbo}",        
        "k_{loss}",
        "#tau_{vw}","A_{vw}","K_{vw}","#phi_{vw}",
        "A_{2cbo}","#phi_{2cbo}",
        "#tau_{y}","A_{y}","K_{y}","#phi_{y}",
        "A_{VW-cbo}","#phi_{VW-cbo}",
    };
    
    name_vars["run1_24paras"] = {
        "N_{0}","#tau","A","R","#phi_{0}",
        "#tau_{cbo}","A^{N}_{cbo}","#omega_{cbo}","#phi_{cbo}",        
        "k_{loss}",
        "#tau_{vw}","A_{vw}","K_{vw}","#phi_{vw}",
        "A_{2cbo}","#phi_{2cbo}","A_{cbo,A}","#phi_{cbo,A}","A_{cbo,#phi}","#phi_{cbo,#phi}",
        "#tau_{y}","A_{y}","K_{y}","#phi_{y}",        
    };
    name_vars["calos_cbo"] = {
        "N_{0}","#tau","A","R","#phi_{0}",
        "#tau_{cbo}","A^{N}_{cbo}","#omega_{cbo}","#phi_{cbo}",        
        "k_{loss}",
        "A_{cbo,A}","#phi_{cbo,A}",
        "C_{cbo}","A_{cbo}^{ch}","#tau_{cbo}^{ch}"

    };
    name_vars["31paras_calos_1"] = {
        "N_{0}","#tau","A","R","#phi_{0}",
        "#tau_{cbo}","A^{N}_{cbo}","#omega_{cbo}","#phi_{cbo}",        
        "k_{loss}",
        "#tau_{vw}","A_{vw}","K_{vw}","#phi_{vw}",
        "A_{2cbo}","#phi_{2cbo}","A_{cbo,A}","#phi_{cbo,A}","A_{cbo,#phi}","#phi_{cbo,#phi}",
        "#tau_{y}","A_{y}","K_{y}","#phi_{y}",
        "A_{VW+cbo}","#phi_{VW+cbo}","A_{VW-cbo}","#phi_{VW-cbo}",
        "C_{cbo}","A_{cbo}^{residual}","#tau_{cbo}^{residual}"
    };    
    name_vars["31paras_calos_2"] = {
        "N_{0}","#tau","A","R","#phi_{0}",
        "#tau_{cbo}","A^{N}_{cbo}","#omega_{cbo}","#phi_{cbo}",        
        "k_{loss}",
        "#tau_{vw}","A_{vw}","K_{vw}","#phi_{vw}",
        "A_{2cbo}","#phi_{2cbo}","A_{cbo,A}","#phi_{cbo,A}","A_{cbo,#phi}","#phi_{cbo,#phi}",
        "#tau_{y}","A_{y}","K_{y}","#phi_{y}",
        "A_{VW+cbo}","#phi_{VW+cbo}","A_{VW-cbo}","#phi_{VW-cbo}",
        "C_{cbo}","A_{cbo}^{residual}","#tau_{cbo}^{residual}"
    };    
}