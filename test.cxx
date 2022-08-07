// #include <string>
// #include <sys/stat.h>
// #include <json/json.h>
// #include <fstream>
// #include <tuple>

// #include "TFile.h"
// #include "TH2D.h"
// #include "TH1D.h"


// #include "FitTools.h"

// using namespace std;

// int main() {

//     TFile * file = new TFile("test.root","RECREATE");
//     file->Write();
//     file->Close();

//     Json::Value json_value;
//     std::ifstream json_file("/Users/cheng/WorkRun4/run2_nominal_T.json");
//     json_file >> json_value;
//     vector<float> init_values;

//     for(auto val : json_value["0"]) {
//         init_values.push_back(val.asFloat());
//         cout << val.asFloat() << endl;
//     }
//     Fitter fitter;
//     return 0;
// }

#include <string>
#include <sys/stat.h>
#include <json/json.h>
#include <fstream>
#include <tuple>

#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"


#include "FitTools.h"

#define var 4

float start_time = 30.1384e3;
float end_time = 650e3;

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

void FullFit(TH1* wiggle, TH1 *lm, string outputDir, string method,vector<double> init_values,int chain_fit) {
    int status = mkdir(outputDir.c_str(),0777);


    Fitter fitter;

    fitter.SetOutputDir(outputDir);
    fitter.SetTimeUnit(Fitter::nano_second);
    if(chain_fit==0) {
        fitter.Fit_28paras_run23_official(method,wiggle,start_time,end_time,init_values,lm);
    }
    else if(chain_fit==1) {
        // 5 paras fit 
        vector<double> init_values_5paras;
        init_values_5paras.insert(init_values_5paras.end(),init_values.begin(),init_values.begin()+5);
        auto info_5pars = fitter.Fit_5paras(method,wiggle,start_time,end_time,init_values_5paras);

        // 28 paras fit
        vector<double> init_values_28paras = read_parameters(info_5pars.file_name,info_5pars.function_name,5);
        init_values_28paras.insert(init_values_28paras.end(),init_values.begin()+5,init_values.end());

        fitter.Fit_28paras_run23_official(method,wiggle,start_time,end_time,init_values_28paras,lm);
          
    }

    
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
// argv[2]: input wiggle name
// argv[3]: input lm file
// argv[4]: output directory
// argv[5]: initial values json file
// argv[6]: start time shift
// argv[7]: version
// argv[8]: chain fit 5 par -> 28 par

int main(int argc,char **argv) {
    cout << "file " << argv[1] << endl;
    TFile * file = TFile::Open(argv[1]);

    
    char *hname = argv[2];

    cout << "Using wiggle: " << hname << endl;

    int triple = atoi(argv[2]);

    TH2D * raw = (TH2D*)file->Get("topDir/Iter0/RawHists/RawHist/rawTimesAndEnergies");
    TH2D * dou = (TH2D*)file->Get("topDir/Iter0/EmpiricalPileupHists/Added/addedPileupTimesAndEnergies_1stOrder");
    TH2D * tri = (TH2D*)file->Get("topDir/Iter0/EmpiricalPileupHists/Added/addedPileupTimesAndEnergies_2ndOrder");
    TH2D * ET_corr = (TH2D*) raw->Clone();
    ET_corr->Add(dou,-1);
    if(triple==1) {
        ET_corr->Add(tri,-1);
    }

    TH1D * wiggle = ET_corr->ProjectionX("wiggle",ET_corr->GetYaxis()->FindBin(1700),10000);
    cout << wiggle->Integral() << endl;


    // TH1D * wiggle = (TH1D*) file->Get(hname)->Clone();

    TFile * file_lm = TFile::Open(argv[3]);

    TH1 * lm = (TH1*)file_lm->Get("integral_hist");

    string outputDir(argv[4]);
    
    Json::Value json_value;
    std::ifstream json_file(argv[5]);
    cout << "load parameters from " << argv[5] << endl;
    
    json_file >> json_value;
    vector<double> init_values;

    for(auto val : json_value["0"]) {
        init_values.push_back(double(val.asFloat()));
    }

    int start_time_shift = atoi(argv[6]);
    start_time =  29.84e3 + 0.1492e3 * start_time_shift;
    cout << " Time range from " << start_time << " to " << end_time << endl;

    string name = Form("%s",argv[7]);

    int chain_fit = atoi(argv[8]);
    // FullFit(wiggle,lm,outputDir,name,init_values,chain_fit);

    return 0;
}