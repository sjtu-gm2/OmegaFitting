#ifndef FITTOOLS_H
#define FITTOOLS_H

#include "configure.h"

#include <string>
#include <iostream>
#include <functional>
#include <iomanip>
#include <fstream>
#include <vector>

#ifdef USE_JSON
    #include <json/json.h>
#endif

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


void FillData(TH1* th1, vector<double> & data);
void FillHist(TH1* th1, const vector<double> & data, bool isAbs=true);
void FFT(TH1* hist, TH1* hist_fft, bool isAbs=true);
vector<double> GetInitialValuesD(string file_path, string index);

struct FitInput {
    TString tag;    
    TH1* wiggle;
    double t_start;
    double t_end;
    std::function<double (double*, double*)> func;
    int nvars;
    vector<string> name_vars;
    vector<double> init_values;
    TH1* lost_muon;    
};
struct FitOutputInfo {
    string file_name;
    string hist_name;
    string function_name;
    string residual_name;
    string fft_name;
    int fitStatus;
    int nvars;
    vector<double> fit_values;    
};

DECLARE_FUNC(_5paras)
DECLARE_FUNC(_9paras)
DECLARE_FUNC(_10paras)
DECLARE_FUNC(_14paras_vo)
DECLARE_FUNC(_14paras_vw)
DECLARE_FUNC(_18paras)
DECLARE_FUNC(_20paras)
DECLARE_FUNC(_24paras)
DECLARE_FUNC(_28paras)

class Fitter {
  public:
    enum TimeUnit {nano_second, micro_second};
    int max_attempts;
    
    Fitter();

    void SetOutputDir(string _output_dir);
    void SetTimeUnit(TimeUnit t_unit);
    void SetBlindedString(string _which_run);
    void SetMaxAttempts(int _max_attempts){max_attempts = _max_attempts;};
    void SetFixParameters(map<int, double> _fix_parameters){fix_parameters=_fix_parameters;};
    void SetRangeParameters(map<int, pair<double,double> > _range_parameters){range_parameters = _range_parameters;};

    REGISTER_FUNC(_5paras, 5)
    REGISTER_FUNC(_9paras, 9)
    REGISTER_FUNC(_10paras, 10)
    REGISTER_FUNC(_14paras_vo, 14)
    REGISTER_FUNC(_14paras_vw, 14)
    REGISTER_FUNC(_18paras, 18)
    REGISTER_FUNC(_20paras, 20)
    REGISTER_FUNC(_24paras, 24)
    REGISTER_FUNC(_28paras, 28)

  private:    
    map<string,vector<string> > name_vars;
    FitOutputInfo doFit(const FitInput & fit_in);
    string output_dir;
    map<int, double> fix_parameters;
    map<int, pair<double,double> > range_parameters;
};

#endif
