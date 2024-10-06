#ifndef FITTOOLS_H
#define FITTOOLS_H

#include "configure.h"

#include <string>
#include <map>
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
#include "TSpline.h"

#include "pocketfft_hdronly.h"
#include "Blinders.hh"

// #include <Math/Interpolator.h>

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
    map<string, TSpline3*>* splines;
    // ROOT::Math::Interpolator * itp;  
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
DECLARE_FUNC(_28paras)
DECLARE_FUNC(_29paras)

DECLARE_FUNC(_40paras)
DECLARE_FUNC(_42paras)
DECLARE_FUNC(_44paras)
DECLARE_FUNC(_48paras)

DECLARE_FUNC(_41paras)
DECLARE_FUNC(_43paras)
DECLARE_FUNC(_45paras)
DECLARE_FUNC(_49paras)

DECLARE_FUNC(_42parasGPR)
DECLARE_FUNC(_40parasGPR)
DECLARE_FUNC(_44parasGPR)
DECLARE_FUNC(_46parasGPR)
DECLARE_FUNC(_50parasGPR)

DECLARE_FUNC(_2paras)
DECLARE_FUNC(_3paras_kloss)
DECLARE_FUNC(_6paras_kloss)
DECLARE_FUNC(_11paras_cbo)
DECLARE_FUNC(_13paras_res)

class Fitter {
  public:
    enum TimeUnit {nano_second, micro_second};
    int max_attempts;
    
    Fitter();

    void SetOutputDir(string _output_dir);
    void SetTimeUnit(TimeUnit t_unit);
    void SetBlindedString(string _which_run);
    void SetMaxAttempts(int _max_attempts){max_attempts = _max_attempts;};
    void SetFitOption(TString _fit_option){fit_option = _fit_option;};
    void SetFixParameters(map<int, double> _fix_parameters){fix_parameters=_fix_parameters;};
    void SetRangeParameters(map<int, pair<double,double> > _range_parameters){range_parameters = _range_parameters;};

    REGISTER_FUNC(_5paras, 5)
    REGISTER_FUNC(_28paras, 28)
    REGISTER_FUNC(_29paras, 29)

    REGISTER_FUNC(_40paras, 40)
    REGISTER_FUNC(_42paras, 42)
    REGISTER_FUNC(_44paras, 44)
    REGISTER_FUNC(_48paras, 48)

    REGISTER_FUNC(_41paras, 41)
    REGISTER_FUNC(_43paras, 43)
    REGISTER_FUNC(_45paras, 45)
    REGISTER_FUNC(_49paras, 49)

    REGISTER_FUNC(_42parasGPR, 42)
    REGISTER_FUNC(_40parasGPR, 40)
    REGISTER_FUNC(_44parasGPR, 44)
    REGISTER_FUNC(_46parasGPR, 46)
    REGISTER_FUNC(_50parasGPR, 50)

    REGISTER_FUNC(_2paras, 2)
    REGISTER_FUNC(_3paras_kloss, 3)
    REGISTER_FUNC(_6paras_kloss, 6)
    REGISTER_FUNC(_11paras_cbo, 11)
    REGISTER_FUNC(_13paras_res, 13)

  private:    
    map<string,vector<string> > name_vars;
    FitOutputInfo doFit(const FitInput & fit_in);
    string output_dir;
    TString fit_option = "RSML";
    map<int, double> fix_parameters;
    map<int, pair<double,double> > range_parameters;
};

#endif
