#include <string>

#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"

#include "FitTools.h"

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

int main() {
    TFile * file = TFile::Open("/Users/cheng/workspace/Data/Run4U/oldfitter/hists_PU_All.root");
    TH2D * hist_ET = (TH2D*) file->Get("rhoc3_d");

    TFile * file_lm = TFile::Open("/Users/cheng/workspace/Data/Run4U/LM_Run4U.root");
    
    TH1 * lm = (TH1*)file_lm->Get("topDir/Iter5/LostMuons/Cuts/Triples/Losses/triple_losses_spectra_integral");

    vector<double> init_values_5paras = {400e3, 6.44234e+01,-3.76183e-01, 1.50075e+03,2.16642e+00};

    int e_start = hist_ET->GetYaxis()->FindBin(1700.);
    int e_end = hist_ET->GetYaxis()->GetNbins();

    TH1D * wiggle = hist_ET->ProjectionX("wiggle_1p7GeV",e_start,e_end);

    string outputDir("./output");
    Fitter fitter;

    fitter.SetOutputDir(outputDir);
    fitter.SetTimeUnit(Fitter::nano_second);    

    string method = "test";
    //5 paras fit
    fitter.Fit_5paras(method,wiggle,30e3,600e3,init_values_5paras);

    //9 paras fit
    vector<double> init_values_9paras = read_parameters("./output/result_5paras_test.root","func_5paras_test",5);
    vector<double> cbos = {3.18528e+02,3.45506e-03,2.34045e+00,3.40741e+00};
    init_values_9paras.insert(init_values_9paras.end(),cbos.begin(),cbos.end());
    fitter.Fit_9paras_cbo(method,wiggle,30e3,600e3,init_values_9paras);

    //10 paras fit
    vector<double> init_values_10paras = read_parameters("./output/result_9paras_cbo_test.root","func_9paras_cbo_test",9);    
    init_values_10paras.push_back(0.);
    fitter.Fit_10paras_cbo_lost(method,wiggle,30e3,600e3,init_values_9paras,lm);

    //14 paras fit
    vector<double> init_values_14paras = read_parameters("./output/result_10paras_cbo_lost_test.root","func_10paras_cbo_lost_test",10);
    vector<double> vws = {88.44, 1.08532e-04, 2.73870, -4.95202e+00};
    init_values_14paras.insert(init_values_14paras.end(),vws.begin(),vws.end());
    fitter.Fit_14paras_cbo_lost_vw(method,wiggle,30e3,600e3,init_values_14paras,lm);

};