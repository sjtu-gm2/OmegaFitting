#include <string>
#include <sys/stat.h>
#include <json/json.h>
#include <fstream>
#include <tuple>

#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"


#include "FitTools.h"

float end_time;

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

void FullFit(TH1* wiggle, TH1 *lm, string outputDir, string method,vector<float> init_values) {
    int status = mkdir(outputDir.c_str(),0777);

    vector<double> init_values_5paras;
    init_values_5paras.insert(init_values_5paras.end(),init_values.begin(),init_values.begin()+5);
    
    Fitter fitter;

    fitter.SetOutputDir(outputDir);
    fitter.SetTimeUnit(Fitter::nano_second);    

    float start_time = 32e3;
    // float end_time = 300e3;
    //5 paras fit
    auto info_5pars = fitter.Fit_5paras(method,wiggle,start_time,end_time,init_values_5paras);

    //9 paras fit
    vector<double> init_values_9paras = read_parameters(info_5pars.file_name,info_5pars.function_name,5);
    init_values_9paras.insert(init_values_9paras.end(),init_values.begin()+5,init_values.begin()+9);
    auto info_9paras = fitter.Fit_9paras_cbo(method,wiggle,start_time,end_time,init_values_9paras);

    //10 paras fit
    vector<double> init_values_10paras = read_parameters(info_9paras.file_name,info_9paras.function_name,9);
    init_values_10paras.insert(init_values_10paras.end(),init_values.begin()+9,init_values.begin()+10);
    fitter.Fit_10paras_cbo_lost(method,wiggle,start_time,end_time,init_values_9paras,lm);
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
// argv[4]: slice index
// argv[5]: initial values json file
// argv[6]: end time
int main(int argc,char **argv) {
    TFile * file = TFile::Open(argv[1]);
    TH2D * hist_ET = (TH2D*) file->Get("rhoc3_d");
    TFile * file_lm = TFile::Open(argv[2]);    
    TH1 * lm = (TH1*)file_lm->Get("topDir/Iter0/LostMuons/Cuts/Triples/Losses/triple_losses_spectra_integral");

    string outputDir(argv[3]);    
    EnerySlice slice = ESliceParseTMethod(argv[4]);

    cout << "Processing " << slice.name << endl;
    TH1D * wiggle = hist_ET->ProjectionX(slice.name.c_str(),slice.bin_start,slice.bin_end);

    Json::Value json_value;
    std::ifstream json_file(argv[5]);
    json_file >> json_value;
    vector<float> init_values;

    for(auto val : json_value["0"]) {
        init_values.push_back(val.asFloat());
    }

    int i_end_time = atoi(argv[6]);
    end_time = float(i_end_time) * 1.e3;
    cout << "Using end time up to " << end_time << endl;
    FullFit(wiggle,lm,outputDir,slice.name,init_values);
    return 0;
}