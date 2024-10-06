#include <string>
#include <sys/stat.h>
#include <ctype.h>
#include "TFile.h"
#include "TH1D.h"
#include "FitTools.h"


using namespace std;

void FullFit(string which_run, TH1* wiggle, int start_bin, int end_bin, TH1* lm, string outputDir, string method, 
            vector<double> init_values, vector<int> fit_chain, int attempts, TString fit_option, map<int,double> fix_parameters, 
            map<int, pair<double, double> > range_parameters, map<string, TSpline3*>* splines){
    int status = mkdir(outputDir.c_str(),0777);

    Fitter fitter;
    fitter.SetMaxAttempts(attempts);
    fitter.SetOutputDir(outputDir);
    fitter.SetFitOption(fit_option);
    fitter.SetFixParameters(fix_parameters);
    fitter.SetRangeParameters(range_parameters);
    fitter.SetBlindedString(which_run);

    double binW = wiggle->GetBinWidth(1);
    double start_time = wiggle->GetBinLowEdge(start_bin);
    double end_time   = wiggle->GetBinLowEdge(end_bin) + binW;

    cout << " Time range from " << start_time << " to " << end_time << endl;

    if(fabs(binW-0.1492)<0.0001){           
        fitter.SetTimeUnit(Fitter::micro_second);        
        cout << "wiggle time unit: micro second"<<endl;
    }
    else if(fabs(binW-149.2)<0.0001){
        fitter.SetTimeUnit(Fitter::nano_second);
        cout << "wiggle time unit: nano second"<<endl;
    }
    else{
        cout <<"unknown bin width = " << binW << endl;
        exit(1);
    }

    FitOutputInfo info;
    info.fit_values = init_values;
    for (int fit_mode : fit_chain) {
        if(fit_mode == 1002) info = fitter.Fit_2paras(method, wiggle, start_time, end_time, info.fit_values);
        if(fit_mode == 1003) info = fitter.Fit_3paras_kloss(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 1006) info = fitter.Fit_6paras_kloss(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 1011) info = fitter.Fit_11paras_cbo(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 1013) info = fitter.Fit_13paras_res(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 1005) info = fitter.Fit_5paras(method, wiggle, start_time, end_time, info.fit_values);
        if(fit_mode == 1028) info = fitter.Fit_28paras(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 1029) info = fitter.Fit_29paras(method, wiggle, start_time, end_time, info.fit_values, lm);

        if(fit_mode == 1040) info = fitter.Fit_40paras(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 1042) info = fitter.Fit_42paras(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 1044) info = fitter.Fit_44paras(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 1048) info = fitter.Fit_48paras(method, wiggle, start_time, end_time, info.fit_values, lm);

        if(fit_mode == 1041) info = fitter.Fit_41paras(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 1043) info = fitter.Fit_43paras(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 1045) info = fitter.Fit_45paras(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 1049) info = fitter.Fit_49paras(method, wiggle, start_time, end_time, info.fit_values, lm);

        if(fit_mode == 3042) info = fitter.Fit_42parasGPR(method, wiggle, start_time, end_time, info.fit_values, lm, splines);
        if(fit_mode == 3040) info = fitter.Fit_40parasGPR(method, wiggle, start_time, end_time, info.fit_values, lm, splines);
        if(fit_mode == 3044) info = fitter.Fit_44parasGPR(method, wiggle, start_time, end_time, info.fit_values, lm, splines);
        if(fit_mode == 3046) info = fitter.Fit_46parasGPR(method, wiggle, start_time, end_time, info.fit_values, lm, splines);
        if(fit_mode == 3050) info = fitter.Fit_50parasGPR(method, wiggle, start_time, end_time, info.fit_values, lm, splines);
    }
}

// argv[1]: which run
// argv[2]: wiggle file
// argv[3]: wiggle name
// argv[4]: lm file
// argv[5]: lm name

// argv[6]: initial values json/root file
// argv[7]: initial values name

// argv[8]: output directory
// argv[9]: version
// argv[10]: fit mode

// argv[11]: maximum attempts

// argv[12]: time bin start
// argv[13]: time bin end

// argv[14]: fit option string

// argv[i] start with '--fix': fix parameters
// argv[i] start with '--range': range parameters
int main(int argc, char **argv) {
    string which_run = argv[1];

    cout << "file " << argv[2] << endl;
    TFile* file = TFile::Open(argv[2]);

    char* hname = argv[3];
    cout << "Using wiggle: " << hname << endl;
    TH1D* wiggle = (TH1D*)file->Get(hname)->Clone();

    TFile* file_lm = TFile::Open(argv[4]);
    TH1* lm = (TH1*)file_lm->Get(argv[5]);

    cout << "Get initial values from " << argv[6] << "  " << argv[7] << endl;
    vector<double> init_values = GetInitialValuesD(argv[6], argv[7]);    


    string outputDir(argv[8]);
    string name = Form("%s", argv[9]);

    int mode = atoi(argv[10]);
    cout << "mode " << mode << endl;
    vector<int> fit_chain;
    if(mode==0) fit_chain = {1005};
    if(mode==1) fit_chain = {1028};
    if(mode==2) fit_chain = {1005, 1028};
    if(mode==3) fit_chain = {1005, 1009, 1010};
    if(mode==4) fit_chain = {1010};
    if(mode==5) fit_chain = {1005, 1010, 1015};
    if(mode==6) fit_chain = {1005, 1010, 1018};
    if(mode==7) fit_chain = {1005, 1010, 1024};
    if(mode==8) fit_chain = {1005, 1010, 1028};
    if(mode==41) fit_chain = {1005, 1010, 1041};
    if(mode==42) fit_chain = {1005, 1010, 1042};
    if(mode==43) fit_chain = {1005, 1010, 1043};
    if(mode==44) fit_chain = {1005, 1010, 1044};
    if(mode==45) fit_chain = {1005, 1010, 1045};
    if(mode==46) fit_chain = {1005, 1010, 1046};
    if(mode==48) fit_chain = {1005, 1010, 1048};
    if(mode==49) fit_chain = {1005, 1010, 1049};

    if(mode==1002) fit_chain = {1002};
    if(mode==1003) fit_chain = {1002, 1003};
    if(mode==1006) fit_chain = {5, 1006};
    if(mode==1013) fit_chain = {5, 1011, 1013};

    if(mode==30) fit_chain = {3042};
    if(mode==31) fit_chain = {3044};
    if(mode==32) fit_chain = {3046};
    if(mode==33) fit_chain = {3050};
    if(mode==34) fit_chain = {3040};

    int attempts;
    attempts = atoi(argv[11]);

    int start_bin, end_bin;

    start_bin = atoi(argv[12]);
    end_bin   = atoi(argv[13]);

    TString fit_option = argv[14];

    TString filename_GPRAmp = argv[15];
    TString folder = argv[16];
    map<string, TSpline3*>* splines = new map<string, TSpline3*>();

    if (filename_GPRAmp != "None" && folder != "None") {
        TFile* file_GPRAmp = TFile::Open(filename_GPRAmp);
        (*splines)["alpha_cbo"] = (TSpline3*)file_GPRAmp->Get(folder + "/alpha_cbo");
        (*splines)["beta_cbo"] = (TSpline3*)file_GPRAmp->Get(folder + "/beta_cbo");
        // (*splines)["alpha_vw"] = (TSpline3*)file_GPRAmp->Get(folder + "/alpha_vw");
        // (*splines)["beta_vw"] = (TSpline3*)file_GPRAmp->Get(folder + "/beta_vw");
        // (*splines)["alpha_vo"] = (TSpline3*)file_GPRAmp->Get(folder + "/alpha_vo");
        // (*splines)["beta_vo"] = (TSpline3*)file_GPRAmp->Get(folder + "/beta_vo");
        // (*splines)["alpha_2cbo"] = (TSpline3*)file_GPRAmp->Get(folder + "/alpha_2cbo");
        // (*splines)["beta_2cbo"] = (TSpline3*)file_GPRAmp->Get(folder + "/beta_2cbo");
        // (*splines)["alpha_vwmcbo"] = (TSpline3*)file_GPRAmp->Get(folder + "/alpha_vwmcbo");
        // (*splines)["beta_vwmcbo"] = (TSpline3*)file_GPRAmp->Get(folder + "/beta_vwmcbo");
        // (*splines)["alpha_vwpcbo"] = (TSpline3*)file_GPRAmp->Get(folder + "/alpha_vwpcbo");
        // (*splines)["beta_vwpcbo"] = (TSpline3*)file_GPRAmp->Get(folder + "/beta_vwpcbo");
        (*splines)["alpha_cbopa"] = (TSpline3*)file_GPRAmp->Get(folder + "/alpha_cbopa");
        (*splines)["beta_cbopa"] = (TSpline3*)file_GPRAmp->Get(folder + "/beta_cbopa");
        (*splines)["alpha_cboma"] = (TSpline3*)file_GPRAmp->Get(folder + "/alpha_cboma");
        (*splines)["beta_cboma"] = (TSpline3*)file_GPRAmp->Get(folder + "/beta_cboma");
        // (*splines)["alpha_vwpa"] = (TSpline3*)file_GPRAmp->Get(folder + "/alpha_vwpa");
        // (*splines)["beta_vwpa"] = (TSpline3*)file_GPRAmp->Get(folder + "/beta_vwpa");
        // (*splines)["alpha_vwma"] = (TSpline3*)file_GPRAmp->Get(folder + "/alpha_vwma");
        // (*splines)["beta_vwma"] = (TSpline3*)file_GPRAmp->Get(folder + "/beta_vwma");
        // (*splines)["alpha_vopa"] = (TSpline3*)file_GPRAmp->Get(folder + "/alpha_vopa");
        // (*splines)["beta_vopa"] = (TSpline3*)file_GPRAmp->Get(folder + "/beta_vopa");
        // (*splines)["alpha_voma"] = (TSpline3*)file_GPRAmp->Get(folder + "/alpha_voma");
        // (*splines)["beta_voma"] = (TSpline3*)file_GPRAmp->Get(folder + "/beta_voma");
        // (*splines)["alpha_vwmcboma"] = (TSpline3*)file_GPRAmp->Get(folder + "/alpha_vwmcboma");
        // (*splines)["beta_vwmcboma"] = (TSpline3*)file_GPRAmp->Get(folder + "/beta_vwmcboma");
        // (*splines)["alpha_vopcbo"] = (TSpline3*)file_GPRAmp->Get(folder + "/alpha_vopcbo");
        // (*splines)["beta_vopcbo"] = (TSpline3*)file_GPRAmp->Get(folder + "/beta_vopcbo");
        // (*splines)["alpha_vomcbo"] = (TSpline3*)file_GPRAmp->Get(folder + "/alpha_vomcbo");
        // (*splines)["beta_vomcbo"] = (TSpline3*)file_GPRAmp->Get(folder + "/beta_vomcbo");
    }

    map<int, double> fix_parameters;
    map<int, pair<double, double> > range_parameters;
    int i = 16;
    while(++i < argc){        
        if (strcmp(argv[i], "--fix") == 0 && strcmp(argv[i+1], "None") != 0){
            ++i;
            while(i < argc && isdigit(argv[i][0])){                
                int npar = atoi(argv[i++]);
                char* s_fix_value = argv[i++];
                double fix_value;
                if(strcmp(s_fix_value, "nan") == 0){
                    fix_value = init_values[npar];
                } else{                    
                    fix_value = atof(s_fix_value);
                }
                fix_parameters[npar] = fix_value;
                cout << "Fix parameter " << npar << " to " << fix_value << endl;
            }
        }        
        if (strcmp(argv[i], "--range") == 0 && strcmp(argv[i+1], "None") != 0){
            ++i;
            while (i < argc && isdigit(argv[i][0])) {
                int npar = atoi(argv[i++]);
                auto range = make_pair<double, double>(0, 0);
                range.first = atof(argv[i++]);
                range.second = atof(argv[i++]);
                range_parameters[npar] = range;
                cout << "Set parameter " << npar << " range from " << range.first << " to " << range.second << endl;
            }            
        }        
    }

    FullFit(which_run, wiggle, start_bin, end_bin, lm, outputDir, name, init_values, fit_chain, attempts, fit_option, fix_parameters, range_parameters, splines);
    
    return 0;
}