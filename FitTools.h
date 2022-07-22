#ifndef FITTOOLS_H
#define FITTOOLS_H

#include <string>
#include <iostream>
#include <functional>
#include <iomanip>

#include "TH1D.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TMath.h"
#include "TString.h"

#include "pocketfft_hdronly.h"
#include "Blinders.hh"


using namespace std;

struct FitInput {
    TString tag;    
    TH1 * wiggle;
    double t_start;
    double t_end;
    std::function<double (double *,double *)> func;
    int nvars;
    string * name_vars;
    vector<double> init_values;
    TH1 * lost_muon;    
};

struct FitOutputInfo {
    string file_name;
    string hist_name;
    string function_name;
    string residual_name;
    string fft_name;
};

class Fitter {
  public:
    enum TimeUnit {nano_second, micro_second};
    
    Fitter();
    void SetOutputDir(string output_dir);
    void SetTimeUnit(TimeUnit t_unit);

    FitOutputInfo Fit_5paras(string name, TH1* wiggle, double t_start, double t_end, vector<double> init_values);
    FitOutputInfo Fit_9paras_cbo(string name, TH1* wiggle, double t_start, double t_end, vector<double> init_values);
    FitOutputInfo Fit_10paras_cbo_lost(string name, TH1* wiggle, double t_start, double t_end, vector<double> init_values,TH1* lm);
    FitOutputInfo Fit_14paras_cbo_lost_vw(string name, TH1* wiggle, double t_start, double t_end, vector<double> init_values,TH1* lm);
    FitOutputInfo Fit_28paras_run23_official(string name, TH1* wiggle, double t_start, double t_end, vector<double> init_values,TH1* lm);


  private:    
    FitOutputInfo doFit(const FitInput & fit_in);
    string m_output_dir;
};
#endif