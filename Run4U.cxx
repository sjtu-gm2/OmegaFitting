#include <string>
#include <sys/stat.h>
#include "TFile.h"
#include "TH1D.h"
#include "FitTools.h"

float start_time = 30.1384e3;
float end_time = 650e3;

using namespace std;

vector<double> read_parameters(string file_name, string func_name, int nvars) {
    TFile * f = TFile::Open(file_name.c_str());
    TF1 * func = (TF1*)f->Get(func_name.c_str());
    vector<double> vals;
    for(int n=0;n<nvars;n++) {
        vals.push_back(func->GetParameter(n));
    }    
    return vals;
}

void fill_parameters_chain(const FitOutputInfo &info, vector<double> & init_values_run, int nvars) {
    string file_name = info.file_name;
    string func_name = info.function_name;
    auto vals = read_parameters(file_name,func_name,nvars);
    for(int n=0;n<nvars;n++) {
        init_values_run[n] = vals[n];
    }
}


void FullFit(TH1* wiggle, TH1 *lm, string outputDir, string method,vector<double> init_values,vector<int> fit_chain,int attempts) {
    int status = mkdir(outputDir.c_str(),0777);

    Fitter fitter;
    fitter.SetMaxAttempts(attempts);
    fitter.SetOutputDir(outputDir);
    fitter.SetTimeUnit(Fitter::nano_second);

    FitOutputInfo info;
    info.fit_values = init_values;
    for(int fit_mode : fit_chain){
        if(fit_mode == 5) info = fitter.Fit_5paras(method,wiggle,start_time,end_time,info.fit_values);
        if(fit_mode == 9) info = fitter.Fit_9paras_cbo(method,wiggle,start_time,end_time,info.fit_values);
        if(fit_mode == 10) info = fitter.Fit_10paras_cbo_lost(method,wiggle,start_time,end_time,info.fit_values,lm);
        
        //changing cbo mode
        if(fit_mode == 11) info = fitter.Fit_11paras_changing_cbo(method,wiggle,start_time,end_time,info.fit_values);
        if(fit_mode == 12) info = fitter.Fit_12paras_changing_cbo(method,wiggle,start_time,end_time,info.fit_values,lm);
        if(fit_mode == 13) info = fitter.Fit_13paras_cbo_vo(method,wiggle,start_time,end_time,info.fit_values,lm);
        if(fit_mode == 15) info = fitter.Fit_15paras_changing_cbo_vo(method,wiggle,start_time,end_time,info.fit_values,lm);



        if(fit_mode == 14) info = fitter.Fit_14paras_cbo_lost_vo(method,wiggle,start_time,end_time,info.fit_values,lm);
        if(fit_mode == 18) info = fitter.Fit_18paras_cbo_lost_vo_vw(method,wiggle,start_time,end_time,info.fit_values,lm);
        if(fit_mode == 22) info = fitter.Fit_22paras_cbo_lost_vw_expansion_lite(method,wiggle,start_time,end_time,info.fit_values,lm);
        if(fit_mode == 28) info = fitter.Fit_28paras_cbo_lost_vw_expansion(method,wiggle,start_time,end_time,info.fit_values,lm);
    }
}

// argv[1]: wiggle file
// argv[2]: wiggle name
// argv[3]: lm file
// argv[4]: lm name
// argv[5]: initial values json/root file
// argv[6]: initial values name
// argv[7]: output directory
// argv[8]: version
// argv[9]: fit mode

// argv[10]: maximum attempts (optional)
// argv[11]: time bin shift (optional)
int main(int argc,char **argv) {
    cout << "file " << argv[1] << endl;
    TFile * file = TFile::Open(argv[1]);

    char *hname = argv[2];
    cout << "Using wiggle: " << hname << endl;
    TH1D * wiggle = (TH1D*) file->Get(hname)->Clone();

    TFile * file_lm = TFile::Open(argv[3]);
    TH1 * lm = (TH1*)file_lm->Get(argv[4]);

    cout << "Get initialValues from " << argv[5] << "  " << argv[6] << endl;
    vector<double> init_values = GetInitialValuesD(argv[5],argv[6]);


    string outputDir(argv[7]);
    string name = Form("%s",argv[8]);

    int mode = atoi(argv[9]);
    vector<int> fit_chain;
    if(mode==0) fit_chain = {28};
    if(mode==1) fit_chain = {5,28};
    if(mode==2) fit_chain = {5};
    if(mode==3) fit_chain = {5,9,10};
    if(mode==4) fit_chain = {10};
    if(mode==5) fit_chain = {5,9,10,14,28};
    if(mode==6) fit_chain = {5,9,10,14,18};
    if(mode==7) fit_chain = {5,9,10,22};

    //frequency changing cbo
    if(mode==8) fit_chain = {5,9,10,12};
    if(mode==9) fit_chain = {5,9,11};
    if(mode==10) fit_chain = {5,9,13};
    if(mode==11) fit_chain = {5,9,13,15};

    int attempts = 1;
    int time_shift = 0;
    if(argc>10) {
        attempts = atoi(argv[10]);
    }

    if(argc>11) {
        time_shift = atoi(argv[11]);
        start_time += 0.1492e3 * time_shift;
    }


    cout << " Time range from " << start_time << " to " << end_time << endl;
    FullFit(wiggle,lm,outputDir,name,init_values,fit_chain,attempts);

    return 0;
}
