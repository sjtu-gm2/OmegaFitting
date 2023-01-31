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


void FillData(TH1* th1,std::vector<double> & data);
void FillHist(TH1* th1,const std::vector<double> & data,bool isAbs=true);
void FFT(TH1* hist, TH1* hist_fft, bool isAbs=true);
vector<double> GetInitialValuesD(string file_path, string index);

struct FitInput {
    TString tag;    
    TH1 * wiggle;
    double t_start;
    double t_end;
    std::function<double (double *,double *)> func;
    int nvars;
    vector<string> name_vars;
    vector<double> init_values;
    TH1 * lost_muon;    
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
DECLARE_FUNC(_9paras_cbo)
DECLARE_FUNC(_10paras_cbo_lost)
DECLARE_FUNC(_13paras_cbo_vo)
DECLARE_FUNC(_14paras_cbo_lost_vw)
DECLARE_FUNC(_14paras_cbo_lost_vo)
DECLARE_FUNC(_18paras_cbo_lost_vo_vw)
DECLARE_FUNC(_22paras_cbo_lost_vw_expansion_lite)
DECLARE_FUNC(_28paras_cbo_lost_vw_expansion)
DECLARE_FUNC(_11paras_changing_cbo)
DECLARE_FUNC(_12paras_changing_cbo)
DECLARE_FUNC(_15paras_changing_cbo_vo)

DECLARE_FUNC(_29paras_cbo_envelope_C)
DECLARE_FUNC(_30paras_cbo_freq)
DECLARE_FUNC(_31paras_cbo_time)
DECLARE_FUNC(_simplified_22paras_calos)


//calorimeter fit
DECLARE_FUNC(_calos_cbo)
DECLARE_FUNC(_31paras_calos_1)
DECLARE_FUNC(_31paras_calos_2)




class Fitter {
  public:
    enum TimeUnit {nano_second, micro_second};
    int max_attempts;
    
    Fitter();

    void SetOutputDir(string output_dir);
    void SetTimeUnit(TimeUnit t_unit);
    void SetMaxAttempts(int attempts) {max_attempts = attempts;};
    void SetFixParameters(map<int,double> _fix_parameters) {fix_parameters=_fix_parameters;};


    REGISTER_FUNC(_5paras,5)
    REGISTER_FUNC(_9paras_cbo,9)
    REGISTER_FUNC(_10paras_cbo_lost,10)

    REGISTER_FUNC(_14paras_cbo_lost_vw,14)
    REGISTER_FUNC(_14paras_cbo_lost_vo,14)    
    
    REGISTER_FUNC(_18paras_cbo_lost_vo_vw,18)

    REGISTER_FUNC(_22paras_cbo_lost_vw_expansion_lite,22)
    REGISTER_FUNC(_28paras_cbo_lost_vw_expansion,28)

    REGISTER_FUNC(_11paras_changing_cbo,11)
    REGISTER_FUNC(_12paras_changing_cbo,12)
    REGISTER_FUNC(_13paras_cbo_vo,13)

    REGISTER_FUNC(_15paras_changing_cbo_vo,15)

    
    REGISTER_FUNC(_simplified_22paras_calos,22)
    REGISTER_FUNC(_29paras_cbo_envelope_C,29)
    REGISTER_FUNC(_30paras_cbo_freq,30)
    REGISTER_FUNC(_31paras_cbo_time,31)

    // REGISTER_FUNC(_run1_24paras,24)

    REGISTER_FUNC(_calos_cbo,15)
    REGISTER_FUNC(_31paras_calos_1,31)
    REGISTER_FUNC(_31paras_calos_2,31)
    
    

  private:    
    map<string,vector<string>> name_vars;
    FitOutputInfo doFit(const FitInput & fit_in);
    string m_output_dir;
    map<int,double> fix_parameters;
};

#endif
