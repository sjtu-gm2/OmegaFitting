#include <string>
#include <sys/stat.h>
#include <ctype.h>
#include "TFile.h"
#include "TH1D.h"
#include "FitTools.h"


using namespace std;

void FullFit(string which_run, TH1* wiggle, int start_bin, int end_bin, TH1* lm, string outputDir, string method, 
            vector<double> init_values, vector<int> fit_chain, int attempts, TString fit_option, map<int,double> fix_parameters, 
            map<int, pair<double, double> > range_parameters){
    int status = mkdir(outputDir.c_str(),0777);

    Fitter fitter;
    fitter.SetMaxAttempts(attempts);
    fitter.SetOutputDir(outputDir);
    fitter.SetFitOption(fit_option);
    fitter.SetFixParameters(fix_parameters);
    fitter.SetRangeParameters(range_parameters);
    fitter.SetBlindedString(which_run);

    double binW = wiggle->GetBinWidth(1);
    double start_time = start_bin * binW;
    double end_time   = end_bin * binW;

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
    for(int fit_mode : fit_chain){
        if(fit_mode == 5) info = fitter.Fit_5paras(method, wiggle, start_time, end_time, info.fit_values);
        if(fit_mode == 9) info = fitter.Fit_9paras(method, wiggle, start_time, end_time, info.fit_values);
        if(fit_mode == 10) info = fitter.Fit_10paras(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 14) info = fitter.Fit_14paras_vo(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 15) info = fitter.Fit_14paras_vw(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 18) info = fitter.Fit_18paras(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 20) info = fitter.Fit_20paras(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 24) info = fitter.Fit_24paras(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 28) info = fitter.Fit_28paras(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 42) info = fitter.Fit_42paras(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 46) info = fitter.Fit_46paras(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 1002) info = fitter.Fit_2paras(method, wiggle, start_time, end_time, info.fit_values);
        if(fit_mode == 1003) info = fitter.Fit_3paras_kloss(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 1006) info = fitter.Fit_6paras_kloss(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 1011) info = fitter.Fit_11paras_cbo(method, wiggle, start_time, end_time, info.fit_values, lm);
        if(fit_mode == 1013) info = fitter.Fit_13paras_res(method, wiggle, start_time, end_time, info.fit_values, lm);
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
    if(mode==0) fit_chain = {5};
    if(mode==1) fit_chain = {28};
    if(mode==2) fit_chain = {5, 28};
    if(mode==3) fit_chain = {5, 9, 10};
    if(mode==4) fit_chain = {10};
    if(mode==5) fit_chain = {5, 9, 10, 14, 15};
    if(mode==6) fit_chain = {5, 9, 10, 14, 15, 18};
    if(mode==7) fit_chain = {5, 9, 10, 14, 15, 18, 20, 24};
    if(mode==8) fit_chain = {5, 9, 10, 28};
    if(mode==1002) fit_chain = {1002};
    if(mode==1003) fit_chain = {1002, 1003};
    if(mode==1006) fit_chain = {5, 1006};
    if(mode==1013) fit_chain = {5, 1011, 1013};
    if(mode==1042) fit_chain = {5, 42};
    if(mode==1046) fit_chain = {5, 46};

    int attempts;
    attempts = atoi(argv[11]);

    int start_bin, end_bin;

    start_bin = atoi(argv[12]);
    end_bin   = atoi(argv[13]);

    TString fit_option = argv[14];

    map<int, double> fix_parameters;
    map<int, pair<double, double> > range_parameters;
    int i = 14;
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
    FullFit(which_run, wiggle, start_bin, end_bin, lm, outputDir, name, init_values, fit_chain, attempts, fit_option, fix_parameters, range_parameters);
    return 0;
}