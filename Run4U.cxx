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

void FullFit(TH1* wiggle, TH1 *lm, string outputDir, string method,vector<double> init_values,int chain_fit,int attempts) {
    int status = mkdir(outputDir.c_str(),0777);


    Fitter fitter;
    fitter.SetMaxAttempts(attempts);
    fitter.SetOutputDir(outputDir);
    fitter.SetTimeUnit(Fitter::nano_second);
    if(chain_fit==0) {
        fitter.Fit_28paras_cbo_lost_vw_expansion(method,wiggle,start_time,end_time,init_values,lm);
    }
    else if(chain_fit==1) {
        // 5 paras fit 
        vector<double> init_values_5paras;
        init_values_5paras.insert(init_values_5paras.end(),init_values.begin(),init_values.begin()+5);
        auto info_5pars = fitter.Fit_5paras(method,wiggle,start_time,end_time,init_values_5paras);

        // 28 paras fit
        vector<double> init_values_28paras = read_parameters(info_5pars.file_name,info_5pars.function_name,5);
        init_values_28paras.insert(init_values_28paras.end(),init_values.begin()+5,init_values.end());

        fitter.Fit_28paras_cbo_lost_vw_expansion(method,wiggle,start_time,end_time,init_values_28paras,lm);        
    }
    else if(chain_fit==2) {
        // 5 paras fit 
        vector<double> init_values_5paras;
        init_values_5paras.insert(init_values_5paras.end(),init_values.begin(),init_values.begin()+5);
        auto info_5pars = fitter.Fit_5paras(method,wiggle,start_time,end_time,init_values_5paras);
    }
    else if(chain_fit==3) {
        // 5 paras fit
        vector<double> init_values_5paras;
        init_values_5paras.insert(init_values_5paras.end(),init_values.begin(),init_values.begin()+5);
        auto info_5pars = fitter.Fit_5paras(method,wiggle,start_time,end_time,init_values_5paras);

        // 9 paras fit
        vector<double> init_values_9paras = read_parameters(info_5pars.file_name,info_5pars.function_name,5);
        init_values_9paras.insert(init_values_9paras.end(),init_values.begin()+5,init_values.begin()+9);
        auto info_9pars = fitter.Fit_9paras_cbo(method,wiggle,start_time,end_time,init_values_9paras);

        // 10 paras fit
        vector<double> init_values_10paras = read_parameters(info_9pars.file_name,info_9pars.function_name,9);
        if(init_values.size()<10) {
            init_values_10paras.push_back(0.);    
        }
        else {
            init_values_10paras.insert(init_values_10paras.end(),init_values.begin()+9,init_values.begin()+10);    
        }
        
        fitter.Fit_10paras_cbo_lost(method,wiggle,start_time,end_time,init_values_10paras,lm);
    }
    else if(chain_fit==4) {
        fitter.Fit_10paras_cbo_lost(method,wiggle,start_time,end_time,init_values,lm);
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
int main(int argc,char **argv) {
    cout << "file " << argv[1] << endl;
    TFile * file = TFile::Open(argv[1]);

    char *hname = argv[2];
    cout << "Using wiggle: " << hname << endl;
    TH1D * wiggle = (TH1D*) file->Get(hname)->Clone();
    

    TFile * file_lm = TFile::Open(argv[3]);
    TH1 * lm = (TH1*)file_lm->Get(argv[4]);

    vector<double> init_values = GetInitialValuesD(argv[5],argv[6]);

    cout << " Time range from " << start_time << " to " << end_time << endl;

    string outputDir(argv[7]);
    string name = Form("%s",argv[8]);

    int chain_fit = atoi(argv[9]);

    int attempts = 1;
    if(argc>10) {
        attempts = atoi(argv[10]);
    }
    FullFit(wiggle,lm,outputDir,name,init_values,chain_fit,attempts);

    return 0;
}