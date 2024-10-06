#include "TH1.h"
#include "TMath.h"
#include "TSpline.h"
#include <map>

using namespace std;
using namespace blinding;

const double Pi = 3.14159265358979323846264338327950288419716939937510582;
const double pipy = 3.141592653589793;

double time_scale = 1.0;
Blinders::fitType ftype = Blinders::kOmega_a;

string blindedString = "Unexpected virtue of ignorance.";
Blinders *getBlinded;

TH1 *lost_muon;
map<string, TSpline3*>* splines;

// ROOT::Math::Interpolator * itp;

// fit function

double func_5paras(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];
    double wiggle = norm * TMath::Exp(-time/life) * (1 + asym*TMath::Cos(omega*time-phi));

    return wiggle;
}

double func_28paras(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_cbo = p[5];
    double alpha_cbo = p[6];
    double beta_cbo = p[7];
    double omega_cbo = p[8];

    // kloss
    double kloss = p[9];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vw
    double tau_vw = p[10];
    double alpha_vw = p[11];
    double beta_vw = p[12];
    double omega_vw = p[13]*fvw;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo;
    double vo = 1 - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time));

    // dcbo
    double alpha_2cbo = p[18];
    double beta_2cbo = p[19];
    double dcbo = 1 - TMath::Exp(-2*time/tau_cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time));

    // At and phit
    double alpha_At = p[20];
    double beta_At = p[21];
    double At = 1 - TMath::Exp(-time/tau_cbo) * (alpha_At*TMath::Cos(omega_cbo*time)+beta_At*TMath::Sin(omega_cbo*time));
    double alpha_phit = p[22];
    double beta_phit = p[23];
    double phit = 1 - TMath::Exp(-time/tau_cbo) * (alpha_phit*TMath::Cos(omega_cbo*time)+beta_phit*TMath::Sin(omega_cbo*time));

    double wiggle = norm * TMath::Exp(-time/life) * (1 + asym*At*TMath::Cos(omega*time-phi*phit));

    // expansion
    double alpha_vwpcbo = p[24];
    double beta_vwpcbo = p[25];
    double alpha_vwmcbo = p[26];
    double beta_vwmcbo = p[27];

    double expansion = 1 - TMath::Exp(-time/tau_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_cbo-time/tau_vw) * (alpha_vwpcbo*TMath::Cos((omega_vw+omega_cbo)*time)+beta_vwpcbo*TMath::Sin((omega_vw+omega_cbo)*time)
                         + alpha_vwmcbo*TMath::Cos((omega_vw-omega_cbo)*time)+beta_vwmcbo*TMath::Sin((omega_vw-omega_cbo)*time));

    return wiggle * (1-kloss*aloss) * expansion * vo * dcbo;
}

double func_29paras(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_cbo = p[5];
    double alpha_cbo = p[6];
    double beta_cbo = p[7];
    double omega_cbo = p[8];

    // kloss
    double kloss = p[9];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vw
    double tau_vw = p[10];
    double alpha_vw = p[11];
    double beta_vw = p[12];
    double omega_vw = p[13]*fvw;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo;
    double vo = 1 - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time));

    // dcbo
    double alpha_2cbo = p[18];
    double beta_2cbo = p[19];
    double dcbo = 1 - TMath::Exp(-2*time/tau_cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time));

    // At and phit
    double alpha_At = p[20];
    double beta_At = p[21];
    double At = 1 - TMath::Exp(-time/tau_cbo) * (alpha_At*TMath::Cos(omega_cbo*time)+beta_At*TMath::Sin(omega_cbo*time));
    double alpha_phit = p[22];
    double beta_phit = p[23];
    double phit = 1 - TMath::Exp(-time/tau_cbo) * (alpha_phit*TMath::Cos(omega_cbo*time)+beta_phit*TMath::Sin(omega_cbo*time));

    double wiggle = norm * TMath::Exp(-time/life) * (1 + asym*At*TMath::Cos(omega*time-phi*phit));

    // expansion
    double alpha_vwpcbo = p[24];
    double beta_vwpcbo = p[25];
    double alpha_vwmcbo = p[26];
    double beta_vwmcbo = p[27];

    // CBO envelope
    double c_cbo = p[28];

    double expansion = 1 - (TMath::Exp(-time/tau_cbo)+c_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_cbo-time/tau_vw) * (alpha_vwpcbo*TMath::Cos((omega_vw+omega_cbo)*time)+beta_vwpcbo*TMath::Sin((omega_vw+omega_cbo)*time)
                         + alpha_vwmcbo*TMath::Cos((omega_vw-omega_cbo)*time)+beta_vwmcbo*TMath::Sin((omega_vw-omega_cbo)*time));

    return wiggle * (1-kloss*aloss) * expansion * vo * dcbo;
}

double func_40paras(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_cbo = p[5];
    double alpha_cbo = p[6];
    double beta_cbo = p[7];
    double omega_cbo = p[8];

    // kloss
    double kloss = p[9];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vw
    double tau_vw = p[10];
    double alpha_vw = p[11];
    double beta_vw = p[12];
    double omega_vw = p[13]*fvw;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo;

    // dcbo
    double alpha_2cbo = p[18];
    double beta_2cbo = p[19];

    // cbo +- a
    double alpha_cbopa = p[20];
    double beta_cbopa = p[21];
    double alpha_cboma = p[22];
    double beta_cboma = p[23];

    // vw +- a
    double alpha_vwpa = p[24];
    double beta_vwpa = p[25];
    double alpha_vwma = p[26];
    double beta_vwma = p[27];

    // vo +- a
    double alpha_vopa = p[28];
    double beta_vopa = p[29];
    double alpha_voma = p[30];
    double beta_voma = p[31];

    // vw +- cbo
    double alpha_vwpcbo = p[32];
    double beta_vwpcbo = p[33];
    double alpha_vwmcbo = p[34];
    double beta_vwmcbo = p[35];

    // vo +- cbo
    double alpha_vopcbo = p[36];
    double beta_vopcbo = p[37];
    double alpha_vomcbo = p[38];
    double beta_vomcbo = p[39];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - TMath::Exp(-time/tau_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - TMath::Exp(-2*time/tau_cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - TMath::Exp(-time/tau_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         - TMath::Exp(-time/tau_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         - TMath::Exp(-time/tau_vw-time/tau_cbo) * (alpha_vwpcbo*TMath::Cos((omega_vw+omega_cbo)*time)+beta_vwpcbo*TMath::Sin((omega_vw+omega_cbo)*time))
                         - TMath::Exp(-time/tau_vw-time/tau_cbo) * (alpha_vwmcbo*TMath::Cos((omega_vw-omega_cbo)*time)+beta_vwmcbo*TMath::Sin((omega_vw-omega_cbo)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_41paras(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_cbo = p[5];
    double alpha_cbo = p[6];
    double beta_cbo = p[7];
    double omega_cbo = p[8];

    // kloss
    double kloss = p[9];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vw
    double tau_vw = p[10];
    double alpha_vw = p[11];
    double beta_vw = p[12];
    double omega_vw = p[13]*fvw;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo;

    // dcbo
    double alpha_2cbo = p[18];
    double beta_2cbo = p[19];

    // cbo +- a
    double alpha_cbopa = p[20];
    double beta_cbopa = p[21];
    double alpha_cboma = p[22];
    double beta_cboma = p[23];

    // vw +- a
    double alpha_vwpa = p[24];
    double beta_vwpa = p[25];
    double alpha_vwma = p[26];
    double beta_vwma = p[27];

    // vo +- a
    double alpha_vopa = p[28];
    double beta_vopa = p[29];
    double alpha_voma = p[30];
    double beta_voma = p[31];

    // vw +- cbo
    double alpha_vwpcbo = p[32];
    double beta_vwpcbo = p[33];
    double alpha_vwmcbo = p[34];
    double beta_vwmcbo = p[35];

    // vo +- cbo
    double alpha_vopcbo = p[36];
    double beta_vopcbo = p[37];
    double alpha_vomcbo = p[38];
    double beta_vomcbo = p[39];

    // envelope
    double c_cbo = p[40];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - (TMath::Exp(-time/tau_cbo)+c_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - TMath::Exp(-2*time/tau_cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - (TMath::Exp(-time/tau_cbo)+c_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         - (TMath::Exp(-time/tau_cbo)+c_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         - TMath::Exp(-time/tau_vw-time/tau_cbo) * (alpha_vwpcbo*TMath::Cos((omega_vw+omega_cbo)*time)+beta_vwpcbo*TMath::Sin((omega_vw+omega_cbo)*time))
                         - TMath::Exp(-time/tau_vw-time/tau_cbo) * (alpha_vwmcbo*TMath::Cos((omega_vw-omega_cbo)*time)+beta_vwmcbo*TMath::Sin((omega_vw-omega_cbo)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_42paras(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_cbo = p[5];
    double alpha_cbo = p[6];
    double beta_cbo = p[7];
    double omega_cbo = p[8];

    // kloss
    double kloss = p[9];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vw
    double tau_vw = p[10];
    double alpha_vw = p[11];
    double beta_vw = p[12];
    double omega_vw = p[13]*fvw;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo;

    // dcbo
    double alpha_2cbo = p[18];
    double beta_2cbo = p[19];

    // vw-cbo
    double tau_vwmcbo = p[20];
    double alpha_vwmcbo = p[21];
    double beta_vwmcbo = p[22];
    double omega_vwmcbo = p[23];

    // cbo +- a
    double alpha_cbopa = p[24];
    double beta_cbopa = p[25];
    double alpha_cboma = p[26];
    double beta_cboma = p[27];

    // vw +- a
    double alpha_vwpa = p[28];
    double beta_vwpa = p[29];
    double alpha_vwma = p[30];
    double beta_vwma = p[31];

    // vo +- a
    double alpha_vopa = p[32];
    double beta_vopa = p[33];
    double alpha_voma = p[34];
    double beta_voma = p[35];

    // vw-cbo-a
    double alpha_vwmcboma = p[36];
    double beta_vwmcboma = p[37];

    // vo +- cbo
    double alpha_vopcbo = p[38];
    double beta_vopcbo = p[39];
    double alpha_vomcbo = p[40];
    double beta_vomcbo = p[41];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - TMath::Exp(-time/tau_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - TMath::Exp(-2*time/tau_cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbo*TMath::Cos(omega_vwmcbo*time)+beta_vwmcbo*TMath::Sin(omega_vwmcbo*time))
                         - TMath::Exp(-time/tau_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         - TMath::Exp(-time/tau_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_43paras(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_cbo = p[5];
    double alpha_cbo = p[6];
    double beta_cbo = p[7];
    double omega_cbo = p[8];

    // kloss
    double kloss = p[9];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vw
    double tau_vw = p[10];
    double alpha_vw = p[11];
    double beta_vw = p[12];
    double omega_vw = p[13]*fvw;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo;

    // dcbo
    double alpha_2cbo = p[18];
    double beta_2cbo = p[19];

    // vw-cbo
    double tau_vwmcbo = p[20];
    double alpha_vwmcbo = p[21];
    double beta_vwmcbo = p[22];
    double omega_vwmcbo = p[23];

    // cbo +- a
    double alpha_cbopa = p[24];
    double beta_cbopa = p[25];
    double alpha_cboma = p[26];
    double beta_cboma = p[27];

    // vw +- a
    double alpha_vwpa = p[28];
    double beta_vwpa = p[29];
    double alpha_vwma = p[30];
    double beta_vwma = p[31];

    // vo +- a
    double alpha_vopa = p[32];
    double beta_vopa = p[33];
    double alpha_voma = p[34];
    double beta_voma = p[35];

    // vw-cbo-a
    double alpha_vwmcboma = p[36];
    double beta_vwmcboma = p[37];

    // vo +- cbo
    double alpha_vopcbo = p[38];
    double beta_vopcbo = p[39];
    double alpha_vomcbo = p[40];
    double beta_vomcbo = p[41];

    // envelope
    double c_cbo = p[42];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - (TMath::Exp(-time/tau_cbo)+c_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - TMath::Exp(-2*time/tau_cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbo*TMath::Cos(omega_vwmcbo*time)+beta_vwmcbo*TMath::Sin(omega_vwmcbo*time))
                         - (TMath::Exp(-time/tau_cbo)+c_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         - (TMath::Exp(-time/tau_cbo)+c_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_44paras(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_cbo = p[5];
    double alpha_cbo = p[6];
    double beta_cbo = p[7];
    double omega_cbo = p[8];

    // kloss
    double kloss = p[9];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vw
    double tau_vw = p[10];
    double alpha_vw = p[11];
    double beta_vw = p[12];
    double omega_vw = p[13]*fvw;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo;

    // dcbo
    double alpha_2cbo = p[18];
    double beta_2cbo = p[19];

    // vw-cbo
    double tau_vwmcbo = p[20];
    double alpha_vwmcbo = p[21];
    double beta_vwmcbo = p[22];
    double omega_vwmcbo = p[23];

    // cbo +- a
    double alpha_cbopa = p[24];
    double beta_cbopa = p[25];
    double alpha_cboma = p[26];
    double beta_cboma = p[27];

    // vw +- a
    double alpha_vwpa = p[28];
    double beta_vwpa = p[29];
    double alpha_vwma = p[30];
    double beta_vwma = p[31];

    // vo +- a
    double alpha_vopa = p[32];
    double beta_vopa = p[33];
    double alpha_voma = p[34];
    double beta_voma = p[35];

    // vw-cbo +- a
    double alpha_vwmcbopa = p[36];
    double beta_vwmcbopa = p[37];
    double alpha_vwmcboma = p[38];
    double beta_vwmcboma = p[39];

    // vo +- cbo
    double alpha_vopcbo = p[40];
    double beta_vopcbo = p[41];
    double alpha_vomcbo = p[42];
    double beta_vomcbo = p[43];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - TMath::Exp(-time/tau_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - TMath::Exp(-2*time/tau_cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbo*TMath::Cos(omega_vwmcbo*time)+beta_vwmcbo*TMath::Sin(omega_vwmcbo*time))
                         - TMath::Exp(-time/tau_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         - TMath::Exp(-time/tau_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbopa*TMath::Cos((omega_vwmcbo+omega)*time)+beta_vwmcbopa*TMath::Sin((omega_vwmcbo+omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_45paras(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_cbo = p[5];
    double alpha_cbo = p[6];
    double beta_cbo = p[7];
    double omega_cbo = p[8];

    // kloss
    double kloss = p[9];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vw
    double tau_vw = p[10];
    double alpha_vw = p[11];
    double beta_vw = p[12];
    double omega_vw = p[13]*fvw;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo;

    // dcbo
    double alpha_2cbo = p[18];
    double beta_2cbo = p[19];

    // vw-cbo
    double tau_vwmcbo = p[20];
    double alpha_vwmcbo = p[21];
    double beta_vwmcbo = p[22];
    double omega_vwmcbo = p[23];

    // cbo +- a
    double alpha_cbopa = p[24];
    double beta_cbopa = p[25];
    double alpha_cboma = p[26];
    double beta_cboma = p[27];

    // vw +- a
    double alpha_vwpa = p[28];
    double beta_vwpa = p[29];
    double alpha_vwma = p[30];
    double beta_vwma = p[31];

    // vo +- a
    double alpha_vopa = p[32];
    double beta_vopa = p[33];
    double alpha_voma = p[34];
    double beta_voma = p[35];

    // vw-cbo +- a
    double alpha_vwmcbopa = p[36];
    double beta_vwmcbopa = p[37];
    double alpha_vwmcboma = p[38];
    double beta_vwmcboma = p[39];

    // vo +- cbo
    double alpha_vopcbo = p[40];
    double beta_vopcbo = p[41];
    double alpha_vomcbo = p[42];
    double beta_vomcbo = p[43];

    // envelope
    double c_cbo = p[44];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - (TMath::Exp(-time/tau_cbo)+c_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - TMath::Exp(-2*time/tau_cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbo*TMath::Cos(omega_vwmcbo*time)+beta_vwmcbo*TMath::Sin(omega_vwmcbo*time))
                         - (TMath::Exp(-time/tau_cbo)+c_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         - (TMath::Exp(-time/tau_cbo)+c_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbopa*TMath::Cos((omega_vwmcbo+omega)*time)+beta_vwmcbopa*TMath::Sin((omega_vwmcbo+omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_48paras(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_cbo = p[5];
    double alpha_cbo = p[6];
    double beta_cbo = p[7];
    double omega_cbo = p[8];

    // kloss
    double kloss = p[9];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vw
    double tau_vw = p[10];
    double alpha_vw = p[11];
    double beta_vw = p[12];
    double omega_vw = p[13]*fvw;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo;

    // dcbo
    double alpha_2cbo = p[18];
    double beta_2cbo = p[19];

    // vw-cbo
    double tau_vwmcbo = p[20];
    double alpha_vwmcbo = p[21];
    double beta_vwmcbo = p[22];
    double omega_vwmcbo = p[23];

    // vw+cbo
    double tau_vwpcbo = p[24];
    double alpha_vwpcbo = p[25];
    double beta_vwpcbo = p[26];
    double omega_vwpcbo = p[27];

    // cbo +- a
    double alpha_cbopa = p[28];
    double beta_cbopa = p[29];
    double alpha_cboma = p[30];
    double beta_cboma = p[31];

    // vw +- a
    double alpha_vwpa = p[32];
    double beta_vwpa = p[33];
    double alpha_vwma = p[34];
    double beta_vwma = p[35];

    // vo +- a
    double alpha_vopa = p[36];
    double beta_vopa = p[37];
    double alpha_voma = p[38];
    double beta_voma = p[39];

    // vw-cbo +- a
    double alpha_vwmcbopa = p[40];
    double beta_vwmcbopa = p[41];
    double alpha_vwmcboma = p[42];
    double beta_vwmcboma = p[43];

    // vo +- cbo
    double alpha_vopcbo = p[44];
    double beta_vopcbo = p[45];
    double alpha_vomcbo = p[46];
    double beta_vomcbo = p[47];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - TMath::Exp(-time/tau_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - TMath::Exp(-2*time/tau_cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbo*TMath::Cos(omega_vwmcbo*time)+beta_vwmcbo*TMath::Sin(omega_vwmcbo*time))
                         - TMath::Exp(-time/tau_vwpcbo) * (alpha_vwpcbo*TMath::Cos(omega_vwpcbo*time)+beta_vwpcbo*TMath::Sin(omega_vwpcbo*time))
                         - TMath::Exp(-time/tau_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         - TMath::Exp(-time/tau_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbopa*TMath::Cos((omega_vwmcbo+omega)*time)+beta_vwmcbopa*TMath::Sin((omega_vwmcbo+omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_49paras(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_cbo = p[5];
    double alpha_cbo = p[6];
    double beta_cbo = p[7];
    double omega_cbo = p[8];

    // kloss
    double kloss = p[9];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vw
    double tau_vw = p[10];
    double alpha_vw = p[11];
    double beta_vw = p[12];
    double omega_vw = p[13]*fvw;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo;

    // dcbo
    double alpha_2cbo = p[18];
    double beta_2cbo = p[19];

    // vw-cbo
    double tau_vwmcbo = p[20];
    double alpha_vwmcbo = p[21];
    double beta_vwmcbo = p[22];
    double omega_vwmcbo = p[23];

    // vw+cbo
    double tau_vwpcbo = p[24];
    double alpha_vwpcbo = p[25];
    double beta_vwpcbo = p[26];
    double omega_vwpcbo = p[27];

    // cbo +- a
    double alpha_cbopa = p[28];
    double beta_cbopa = p[29];
    double alpha_cboma = p[30];
    double beta_cboma = p[31];

    // vw +- a
    double alpha_vwpa = p[32];
    double beta_vwpa = p[33];
    double alpha_vwma = p[34];
    double beta_vwma = p[35];

    // vo +- a
    double alpha_vopa = p[36];
    double beta_vopa = p[37];
    double alpha_voma = p[38];
    double beta_voma = p[39];

    // vw-cbo +- a
    double alpha_vwmcbopa = p[40];
    double beta_vwmcbopa = p[41];
    double alpha_vwmcboma = p[42];
    double beta_vwmcboma = p[43];

    // vo +- cbo
    double alpha_vopcbo = p[44];
    double beta_vopcbo = p[45];
    double alpha_vomcbo = p[46];
    double beta_vomcbo = p[47];

    // envelope
    double c_cbo = p[48];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - (TMath::Exp(-time/tau_cbo)+c_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - TMath::Exp(-2*time/tau_cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbo*TMath::Cos(omega_vwmcbo*time)+beta_vwmcbo*TMath::Sin(omega_vwmcbo*time))
                         - TMath::Exp(-time/tau_vwpcbo) * (alpha_vwpcbo*TMath::Cos(omega_vwpcbo*time)+beta_vwpcbo*TMath::Sin(omega_vwpcbo*time))
                         - (TMath::Exp(-time/tau_cbo)+c_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         - (TMath::Exp(-time/tau_cbo)+c_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbopa*TMath::Cos((omega_vwmcbo+omega)*time)+beta_vwmcbopa*TMath::Sin((omega_vwmcbo+omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_42parasGPR(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_cbo = p[5];
    double alpha_cbo = (*splines)["alpha_cbo"]->Eval(time);
    double beta_cbo = (*splines)["beta_cbo"]->Eval(time);
    double omega_cbo = p[6];

    // kloss
    double kloss = p[7];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vw
    double tau_vw = p[8];
    double alpha_vw = p[9];
    double beta_vw = p[10];
    double omega_vw = p[11]*fvw;

    // vo
    double tau_vo = p[12];
    double alpha_vo = p[13];
    double beta_vo = p[14];
    double omega_vo = p[15]*fvo;

    // dcbo
    double alpha_2cbo = p[16];
    double beta_2cbo = p[17];

    // vw-cbo
    double tau_vwmcbo = p[18];
    double alpha_vwmcbo = p[19];
    double beta_vwmcbo = p[20];
    double omega_vwmcbo = p[21];

    // cbo +- a
    double alpha_cbopa = p[22];
    double beta_cbopa = p[23];
    double alpha_cboma = p[24];
    double beta_cboma = p[25];

    // vw +- a
    double alpha_vwpa = p[26];
    double beta_vwpa = p[27];
    double alpha_vwma = p[28];
    double beta_vwma = p[29];

    // vo +- a
    double alpha_vopa = p[30];
    double beta_vopa = p[31];
    double alpha_voma = p[32];
    double beta_voma = p[33];

    // vw-cbo +- a
    double alpha_vwmcbopa = p[34];
    double beta_vwmcbopa = p[35];
    double alpha_vwmcboma = p[36];
    double beta_vwmcboma = p[37];

    // vo +- cbo
    double alpha_vopcbo = p[38];
    double beta_vopcbo = p[39];
    double alpha_vomcbo = p[40];
    double beta_vomcbo = p[41];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - TMath::Exp(-2*time/tau_cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbo*TMath::Cos(omega_vwmcbo*time)+beta_vwmcbo*TMath::Sin(omega_vwmcbo*time))
                         - TMath::Exp(-time/tau_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         - TMath::Exp(-time/tau_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbopa*TMath::Cos((omega_vwmcbo+omega)*time)+beta_vwmcbopa*TMath::Sin((omega_vwmcbo+omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_40parasGPR(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_2cbo = p[5];
    double alpha_cbo = (*splines)["alpha_cbo"]->Eval(time);
    double beta_cbo = (*splines)["beta_cbo"]->Eval(time);
    double omega_cbo = p[6];

    // kloss
    double kloss = p[7];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vw
    double tau_vw = p[8];
    double alpha_vw = p[9];
    double beta_vw = p[10];
    double omega_vw = p[11]*fvw;

    // vo
    double tau_vo = p[12];
    double alpha_vo = p[13];
    double beta_vo = p[14];
    double omega_vo = p[15]*fvo;

    // dcbo
    double alpha_2cbo = p[16];
    double beta_2cbo = p[17];

    // vw-cbo
    double tau_vwmcbo = p[18];
    double alpha_vwmcbo = p[19];
    double beta_vwmcbo = p[20];
    double omega_vwmcbo = p[21];

    // cbo +- a
    double alpha_cbopa = (*splines)["alpha_cbopa"]->Eval(time);;
    double beta_cbopa = (*splines)["beta_cbopa"]->Eval(time);
    double alpha_cboma = (*splines)["alpha_cboma"]->Eval(time);
    double beta_cboma = (*splines)["beta_cboma"]->Eval(time);

    // vw +- a
    double alpha_vwpa = p[22];
    double beta_vwpa = p[23];
    double alpha_vwma = p[24];
    double beta_vwma = p[25];

    // vo +- a
    double alpha_vopa = p[26];
    double beta_vopa = p[27];
    double alpha_voma = p[28];
    double beta_voma = p[29];

    // vw-cbo +- a
    double alpha_vwmcbopa = p[30];
    double beta_vwmcbopa = p[31];
    double alpha_vwmcboma = p[32];
    double beta_vwmcboma = p[33];

    // vo +- cbo
    double tau_vopcbo = p[34];
    double alpha_vopcbo = p[35];
    double beta_vopcbo = p[36];
    double tau_vomcbo = p[37];
    double alpha_vomcbo = p[38];
    double beta_vomcbo = p[39];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - TMath::Exp(-time/tau_2cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbo*TMath::Cos(omega_vwmcbo*time)+beta_vwmcbo*TMath::Sin(omega_vwmcbo*time))
                         - (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         - (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbopa*TMath::Cos((omega_vwmcbo+omega)*time)+beta_vwmcbopa*TMath::Sin((omega_vwmcbo+omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         - TMath::Exp(-time/tau_vopcbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         - TMath::Exp(-time/tau_vomcbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_46parasGPR(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_2cbo = p[5];
    double alpha_cbo = p[40]*(*splines)["alpha_cbo"]->Eval(time);
    double beta_cbo = p[41]*(*splines)["beta_cbo"]->Eval(time);
    double omega_cbo = p[6];

    // kloss
    double kloss = p[7];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vw
    double tau_vw = p[8];
    double alpha_vw = p[9];
    double beta_vw = p[10];
    double omega_vw = p[11]*fvw;

    // vo
    double tau_vo = p[12];
    double alpha_vo = p[13];
    double beta_vo = p[14];
    double omega_vo = p[15]*fvo;

    // dcbo
    double alpha_2cbo = p[16];
    double beta_2cbo = p[17];

    // vw-cbo
    double tau_vwmcbo = p[18];
    double alpha_vwmcbo = p[19];
    double beta_vwmcbo = p[20];
    double omega_vwmcbo = p[21];

    // cbo +- a
    double alpha_cbopa = p[42]*(*splines)["alpha_cbopa"]->Eval(time);;
    double beta_cbopa = p[43]*(*splines)["beta_cbopa"]->Eval(time);
    double alpha_cboma = p[44]*(*splines)["alpha_cboma"]->Eval(time);
    double beta_cboma = p[45]*(*splines)["beta_cboma"]->Eval(time);

    // vw +- a
    double alpha_vwpa = p[22];
    double beta_vwpa = p[23];
    double alpha_vwma = p[24];
    double beta_vwma = p[25];

    // vo +- a
    double alpha_vopa = p[26];
    double beta_vopa = p[27];
    double alpha_voma = p[28];
    double beta_voma = p[29];

    // vw-cbo +- a
    double alpha_vwmcbopa = p[30];
    double beta_vwmcbopa = p[31];
    double alpha_vwmcboma = p[32];
    double beta_vwmcboma = p[33];

    // vo +- cbo
    double tau_vopcbo = p[34];
    double alpha_vopcbo = p[35];
    double beta_vopcbo = p[36];
    double tau_vomcbo = p[37];
    double alpha_vomcbo = p[38];
    double beta_vomcbo = p[39];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - TMath::Exp(-time/tau_2cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbo*TMath::Cos(omega_vwmcbo*time)+beta_vwmcbo*TMath::Sin(omega_vwmcbo*time))
                         - (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         - (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbopa*TMath::Cos((omega_vwmcbo+omega)*time)+beta_vwmcbopa*TMath::Sin((omega_vwmcbo+omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         - TMath::Exp(-time/tau_vopcbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         - TMath::Exp(-time/tau_vomcbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_44parasGPR(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_cbo = p[5];
    double alpha_cbo = p[42]*(*splines)["alpha_cbo"]->Eval(time);
    double beta_cbo = p[43]*(*splines)["beta_cbo"]->Eval(time);
    double omega_cbo = p[6];

    // kloss
    double kloss = p[7];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vw
    double tau_vw = p[8];
    double alpha_vw = p[9];
    double beta_vw = p[10];
    double omega_vw = p[11]*fvw;

    // vo
    double tau_vo = p[12];
    double alpha_vo = p[13];
    double beta_vo = p[14];
    double omega_vo = p[15]*fvo;

    // dcbo
    double alpha_2cbo = p[16];
    double beta_2cbo = p[17];

    // vw-cbo
    double tau_vwmcbo = p[18];
    double alpha_vwmcbo = p[19];
    double beta_vwmcbo = p[20];
    double omega_vwmcbo = p[21];

    // cbo +- a
    double alpha_cbopa = p[22];
    double beta_cbopa = p[23];
    double alpha_cboma = p[24];
    double beta_cboma = p[25];

    // vw +- a
    double alpha_vwpa = p[26];
    double beta_vwpa = p[27];
    double alpha_vwma = p[28];
    double beta_vwma = p[29];

    // vo +- a
    double alpha_vopa = p[30];
    double beta_vopa = p[31];
    double alpha_voma = p[32];
    double beta_voma = p[33];

    // vw-cbo +- a
    double alpha_vwmcbopa = p[34];
    double beta_vwmcbopa = p[35];
    double alpha_vwmcboma = p[36];
    double beta_vwmcboma = p[37];

    // vo +- cbo
    double alpha_vopcbo = p[38];
    double beta_vopcbo = p[39];
    double alpha_vomcbo = p[40];
    double beta_vomcbo = p[41];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - TMath::Exp(-2*time/tau_cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbo*TMath::Cos(omega_vwmcbo*time)+beta_vwmcbo*TMath::Sin(omega_vwmcbo*time))
                         - TMath::Exp(-time/tau_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         - TMath::Exp(-time/tau_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbopa*TMath::Cos((omega_vwmcbo+omega)*time)+beta_vwmcbopa*TMath::Sin((omega_vwmcbo+omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_50parasGPR(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_cbo = p[5];
    double alpha_cbo = p[46]*(*splines)["alpha_cbo"]->Eval(30.1384+p[47]*(time-30.1384));
    double beta_cbo = p[48]*(*splines)["beta_cbo"]->Eval(30.1384+p[49]*(time-30.1384));
    double omega_cbo = p[6];

    // kloss
    double kloss = p[7];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vw
    double tau_vw = p[8];
    double alpha_vw = p[9];
    double beta_vw = p[10];
    double omega_vw = p[11]*fvw;

    // vo
    double tau_vo = p[12];
    double alpha_vo = p[13];
    double beta_vo = p[14];
    double omega_vo = p[15]*fvo;

    // dcbo
    double alpha_2cbo = p[16];
    double beta_2cbo = p[17];

    // vw-cbo
    double tau_vwmcbo = p[18];
    double alpha_vwmcbo = p[19];
    double beta_vwmcbo = p[20];
    double omega_vwmcbo = p[21];

    // vw+cbo
    double tau_vwpcbo = p[22];
    double alpha_vwpcbo = p[23];
    double beta_vwpcbo = p[24];
    double omega_vwpcbo = p[25];

    // cbo +- a
    double alpha_cbopa = p[26];
    double beta_cbopa = p[27];
    double alpha_cboma = p[28];
    double beta_cboma = p[29];

    // vw +- a
    double alpha_vwpa = p[30];
    double beta_vwpa = p[31];
    double alpha_vwma = p[32];
    double beta_vwma = p[33];

    // vo +- a
    double alpha_vopa = p[34];
    double beta_vopa = p[35];
    double alpha_voma = p[36];
    double beta_voma = p[37];

    // vw-cbo +- a
    double alpha_vwmcbopa = p[38];
    double beta_vwmcbopa = p[39];
    double alpha_vwmcboma = p[40];
    double beta_vwmcboma = p[41];

    // vo +- cbo
    double alpha_vopcbo = p[42];
    double beta_vopcbo = p[43];
    double alpha_vomcbo = p[44];
    double beta_vomcbo = p[45];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - TMath::Exp(-2*time/tau_cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbo*TMath::Cos(omega_vwmcbo*time)+beta_vwmcbo*TMath::Sin(omega_vwmcbo*time))
                         - TMath::Exp(-time/tau_vwpcbo) * (alpha_vwpcbo*TMath::Cos(omega_vwpcbo*time)+beta_vwpcbo*TMath::Sin(omega_vwpcbo*time))
                         - TMath::Exp(-time/tau_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         - TMath::Exp(-time/tau_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbopa*TMath::Cos((omega_vwmcbo+omega)*time)+beta_vwmcbopa*TMath::Sin((omega_vwmcbo+omega)*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         - TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}