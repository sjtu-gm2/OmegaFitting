#include <string>
#include <sys/stat.h>
#include <json/json.h>
#include <fstream>
#include <tuple>

#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"


#include "FitTools.h"

float start_time = 30.1384;
float end_time = 650;

using namespace std;

struct EnerySlice {
    string name;
    int bin_start;
    int bin_end;
    int e_start; //MeV
    int e_end;
};

vector<double> read_parameters(string file_name, string func_name, int nvars) {
    TFile * f = TFile::Open(file_name.c_str());
    TF1 * func = (TF1*)f->Get(func_name.c_str());

    vector<double> vals;
    for(int n=0;n<nvars;n++) {
        vals.push_back(func->GetParameter(n));
    }    
    return vals;
}

void FullFit(TH1* wiggle, TH1 *lm, string outputDir, string method,vector<double> init_values) {
    int status = mkdir(outputDir.c_str(),0777);


    Fitter fitter;

    fitter.SetOutputDir(outputDir);
    fitter.SetTimeUnit(Fitter::micro_second);

    // 5 paras fit
    vector<double> init_values_5paras;
    init_values_5paras.insert(init_values_5paras.end(),init_values.begin(),init_values.begin()+5);
    auto info_5pars = fitter.Fit_5paras(method,wiggle,start_time,end_time,init_values_5paras);

    vector<double> init_values_28paras = read_parameters(info_5pars.file_name,info_5pars.function_name,5);
    init_values_28paras.insert(init_values_28paras.end(),init_values.begin()+5,init_values.end());

    fitter.Fit_28paras_run23_official(method,wiggle,start_time,end_time,init_values_28paras,lm);
}

EnerySlice ESliceParse(char * arg) {
    EnerySlice slice;
    int nslice = atoi(arg);

    slice.bin_start = 2*nslice+21;
    slice.bin_end = 2*nslice+22;
    slice.e_start = (slice.bin_start -1)*50;
    slice.e_end = slice.bin_end * 50;     

    TString name;
    name.Form("Slice_%d_%d_MeV",slice.e_start,slice.e_end);
    slice.name = name.Data();
    return slice;
}

EnerySlice ESliceParseTMethod(char * arg) {
    EnerySlice slice;
    int nslice = atoi(arg);

    slice.bin_start = nslice;
    slice.bin_end = 62;
    slice.e_start = (slice.bin_start -1)*50;
    slice.e_end = slice.bin_end * 50;     
    cout <<"nslice = "<<nslice<<endl;
    cout <<"e_start = "<<slice.e_start<<endl;
    cout <<"e_end = "<<slice.e_end<<endl;
    TString name;
    name.Form("Slice_%d_%d_MeV",slice.e_start,slice.e_end);
    slice.name = name.Data();
    return slice;
}

// argv[1]: input wiggle file
// argv[2]: input lm file
// argv[3]: output directory
// argv[4]: method
// argv[5]: syst scan gain_A,gain_T,stdp_A,stdp_T
// argv[6]: initial values json file
// argv[7]: dataset
// argv[8]: scan point
int main(int argc,char **argv) {
    TFile * file = TFile::Open(argv[1]);

    int scan_number = atoi(argv[8]);
    float scan = float(scan_number) * 0.1;
    char *hname = Form("wiggle1700/%ssys_%s_%.1f",argv[4],argv[5],scan);
    cout << "Using wiggle: " << hname << endl;
    TH1D * wiggle = (TH1D*) file->Get(hname)->Clone();

    TFile * file_lm = TFile::Open(argv[2]);
    TH1 * lm = (TH1*)file_lm->Get("topDir/Iter0/LostMuons/Cuts/Triples/Losses/triple_losses_spectra_integral");

    string outputDir(argv[3]);
    // EnerySlice slice = ESliceParseTMethod(argv[4]);

    // cout << "Processing " << slice.name << endl;


    Json::Value json_value;
    std::ifstream json_file(argv[6]);
    cout << "load parameters from " << argv[6] << endl;
    json_file >> json_value;
    vector<double> init_values;

    for(auto val : json_value["0"]) {
        init_values.push_back(double(val.asFloat()));
    }

    cout << "Using end time up to " << end_time << "micro second" << endl;
    string name = Form("%s_%smethod_%s_%s",argv[7],argv[4],argv[5],argv[8]);
    FullFit(wiggle,lm,outputDir,name,init_values);
    return 0;
}