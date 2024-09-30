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

double func_9paras(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];
    double wiggle = norm * TMath::Exp(-time/life) * (1 + asym*TMath::Cos(omega*time-phi));

    // cbo
    double tau_cbo = p[5];
    double alpha_cbo = p[6];
    double beta_cbo = p[7];
    double omega_cbo = p[8];
    double cbo = 1 - TMath::Exp(-time/tau_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time));

    return wiggle * cbo;
}

double func_10paras(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];
    double wiggle = norm * TMath::Exp(-time/life) * (1 + asym*TMath::Cos(omega*time-phi));

    // cbo
    double tau_cbo = p[5];
    double alpha_cbo = p[6];
    double beta_cbo = p[7];
    double omega_cbo = p[8];
    double cbo = 1 - TMath::Exp(-time/tau_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time));

    // kloss
    double kloss = p[9];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    return wiggle * (1-kloss*aloss) * cbo;
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

    // vo
    double tau_vo = p[10];
    double alpha_vo = p[11];
    double beta_vo = p[12];
    double omega_vo = p[13]*fvo*1e-3;
    double vo = 1 - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time));

    // vw
    double tau_vw = p[14];
    double alpha_vw = p[15];
    double beta_vw = p[16];
    double omega_vw = p[17]*fvw*1e-3;

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
                         + TMath::Exp(-time/tau_cbo-time/tau_vw) * (alpha_vwpcbo*TMath::Cos((omega_vw+omega_cbo)*time)+beta_vwpcbo*TMath::Sin((omega_vw+omega_cbo)*time)
                         + alpha_vwmcbo*TMath::Cos((omega_vw-omega_cbo)*time)+beta_vwmcbo*TMath::Sin((omega_vw-omega_cbo)*time));

    return wiggle * (1-kloss*aloss) * expansion * vo * dcbo;
}

double func_28paras_phi(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_cbo = p[5];
    double A_cbo = p[6];
    double phi_cbo = p[7];
    double omega_cbo = p[8];

    // kloss
    double kloss = p[9];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vo
    double tau_vo = p[10];
    double A_vo = p[11];
    double phi_vo = p[12];
    double omega_vo = p[13]*fvo*1e-3;
    double vo = 1 - TMath::Exp(-time/tau_vo) * A_vo*TMath::Cos(omega_vo*time+phi_vo);

    // vw
    double tau_vw = p[14];
    double A_vw = p[15];
    double phi_vw = p[16];
    double omega_vw = p[17]*fvw*1e-3;

    // dcbo
    double A_dcbo = p[18];
    double phi_dcbo = p[19];
    double dcbo = 1 - TMath::Exp(-2*time/tau_cbo) * A_dcbo*TMath::Cos(2*omega_cbo*time+phi_dcbo);

    // At and phit
    double A_At = p[20];
    double phi_At = p[21];
    double At = 1 - TMath::Exp(-time/tau_cbo) * A_At*TMath::Cos(omega_cbo*time+phi_At);
    double A_phit = p[22];
    double phi_phit = p[23];
    double phit = 1 - TMath::Exp(-time/tau_cbo) * A_phit*TMath::Cos(omega_cbo*time+phi_phit);

    double wiggle = norm * TMath::Exp(-time/life) * (1 - asym*At*TMath::Cos(omega*time+phi*phit));

    // expansion
    double A_vwpcbo = p[24];
    double phi_vwpcbo = p[25];
    double A_vwmcbo = p[26];
    double phi_vwmcbo = p[27];

    double expansion = 1 - TMath::Exp(-time/tau_cbo) * A_cbo*TMath::Cos(omega_cbo*time+phi_cbo)
                         - TMath::Exp(-time/tau_vw) * A_vw*TMath::Cos(omega_vw*time+phi_vw)
                         + TMath::Exp(-time/tau_cbo-time/tau_vw) * (A_vwpcbo*TMath::Cos((omega_vw+omega_cbo)*time+phi_vwpcbo)
                         + A_vwmcbo*TMath::Cos((omega_vw-omega_cbo)*time+phi_vwmcbo));

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

    // vo
    double tau_vo = p[10];
    double alpha_vo = p[11];
    double beta_vo = p[12];
    double omega_vo = p[13]*fvo*1e-3;
    double vo = 1 - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time));

    // vw
    double tau_vw = p[14];
    double alpha_vw = p[15];
    double beta_vw = p[16];
    double omega_vw = p[17]*fvw*1e-3;

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
    double tau_cbovw = p[28];

    double expansion = 1 - TMath::Exp(-time/tau_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         + TMath::Exp(-time/tau_cbovw) * (alpha_vwpcbo*TMath::Cos((omega_vw+omega_cbo)*time)+beta_vwpcbo*TMath::Sin((omega_vw+omega_cbo)*time)
                         + alpha_vwmcbo*TMath::Cos((omega_vw-omega_cbo)*time)+beta_vwmcbo*TMath::Sin((omega_vw-omega_cbo)*time));

    return wiggle * (1-kloss*aloss) * expansion * vo * dcbo;
}

double func_14paras(double *x, double *p){
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

    // cbo +- a
    double alpha_cbopa = p[10];
    double beta_cbopa = p[11];
    double alpha_cboma = p[12];
    double beta_cboma = p[13];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - TMath::Exp(-time/tau_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_26paras(double *x, double *p){
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
    double omega_vw = p[13]*fvw*1e-3;

    // cbo +- a
    double alpha_cbopa = p[14];
    double beta_cbopa = p[15];
    double alpha_cboma = p[16];
    double beta_cboma = p[17];

    // vw +- a
    double alpha_vwpa = p[18];
    double beta_vwpa = p[19];
    double alpha_vwma = p[20];
    double beta_vwma = p[21];

    // vw +- cbo
    double alpha_vwpcbo = p[22];
    double beta_vwpcbo = p[23];
    double alpha_vwmcbo = p[24];
    double beta_vwmcbo = p[25];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - TMath::Exp(-time/tau_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         + TMath::Exp(-time/tau_vw-time/tau_cbo) * (alpha_vwpcbo*TMath::Cos((omega_vw+omega_cbo)*time)+beta_vwpcbo*TMath::Sin((omega_vw+omega_cbo)*time))
                         + TMath::Exp(-time/tau_vw-time/tau_cbo) * (alpha_vwmcbo*TMath::Cos((omega_vw-omega_cbo)*time)+beta_vwmcbo*TMath::Sin((omega_vw-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_38paras(double *x, double *p){
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
    double omega_vw = p[13]*fvw*1e-3;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo*1e-3;

    // cbo +- a
    double alpha_cbopa = p[18];
    double beta_cbopa = p[19];
    double alpha_cboma = p[20];
    double beta_cboma = p[21];

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

    // vw +- cbo
    double alpha_vwpcbo = p[30];
    double beta_vwpcbo = p[31];
    double alpha_vwmcbo = p[32];
    double beta_vwmcbo = p[33];

    // vo +- cbo
    double alpha_vopcbo = p[34];
    double beta_vopcbo = p[35];
    double alpha_vomcbo = p[36];
    double beta_vomcbo = p[37];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - TMath::Exp(-time/tau_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         + TMath::Exp(-time/tau_vw-time/tau_cbo) * (alpha_vwpcbo*TMath::Cos((omega_vw+omega_cbo)*time)+beta_vwpcbo*TMath::Sin((omega_vw+omega_cbo)*time))
                         + TMath::Exp(-time/tau_vw-time/tau_cbo) * (alpha_vwmcbo*TMath::Cos((omega_vw-omega_cbo)*time)+beta_vwmcbo*TMath::Sin((omega_vw-omega_cbo)*time))
                         + TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         + TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
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
    double omega_vw = p[13]*fvw*1e-3;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo*1e-3;

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
                         + TMath::Exp(-time/tau_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         + TMath::Exp(-time/tau_vw-time/tau_cbo) * (alpha_vwpcbo*TMath::Cos((omega_vw+omega_cbo)*time)+beta_vwpcbo*TMath::Sin((omega_vw+omega_cbo)*time))
                         + TMath::Exp(-time/tau_vw-time/tau_cbo) * (alpha_vwmcbo*TMath::Cos((omega_vw-omega_cbo)*time)+beta_vwmcbo*TMath::Sin((omega_vw-omega_cbo)*time))
                         + TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         + TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

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
    double omega_vw = p[13]*fvw*1e-3;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo*1e-3;

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
                         + TMath::Exp(-time/tau_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         + TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         + TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         + TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_46paras(double *x, double *p){
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
    double omega_vw = p[13]*fvw*1e-3;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo*1e-3;

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

    // vw-cbo-a
    double alpha_vwmcboma = p[40];
    double beta_vwmcboma = p[41];

    // vo +- cbo
    double alpha_vopcbo = p[42];
    double beta_vopcbo = p[43];
    double alpha_vomcbo = p[44];
    double beta_vomcbo = p[45];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - TMath::Exp(-time/tau_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - TMath::Exp(-2*time/tau_cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbo*TMath::Cos(omega_vwmcbo*time)+beta_vwmcbo*TMath::Sin(omega_vwmcbo*time))
                         - TMath::Exp(-time/tau_vwpcbo) * (alpha_vwpcbo*TMath::Cos(omega_vwpcbo*time)+beta_vwpcbo*TMath::Sin(omega_vwpcbo*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         + TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         + TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         + TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

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
    double omega_vw = p[13]*fvw*1e-3;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo*1e-3;

    // dcbo
    double tau_2cbo = p[18];
    double alpha_2cbo = p[19];
    double beta_2cbo = p[20];

    // vw-cbo
    double tau_vwmcbo = p[21];
    double alpha_vwmcbo = p[22];
    double beta_vwmcbo = p[23];
    double omega_vwmcbo = p[24];

    // cbo +- a
    double alpha_cbopa = p[25];
    double beta_cbopa = p[26];
    double alpha_cboma = p[27];
    double beta_cboma = p[28];

    // vw +- a
    double alpha_vwpa = p[29];
    double beta_vwpa = p[30];
    double alpha_vwma = p[31];
    double beta_vwma = p[32];

    // vo +- a
    double alpha_vopa = p[33];
    double beta_vopa = p[34];
    double alpha_voma = p[35];
    double beta_voma = p[36];

    // vw-cbo-a
    double alpha_vwmcboma = p[37];
    double beta_vwmcboma = p[38];

    // vo +- cbo
    double tau_vopmcbo = p[39];
    double alpha_vopcbo = p[40];
    double beta_vopcbo = p[41];
    double alpha_vomcbo = p[42];
    double beta_vomcbo = p[43];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - TMath::Exp(-time/tau_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - TMath::Exp(-time/tau_2cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbo*TMath::Cos(omega_vwmcbo*time)+beta_vwmcbo*TMath::Sin(omega_vwmcbo*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         + TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         + TMath::Exp(-time/tau_vopmcbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         + TMath::Exp(-time/tau_vopmcbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

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
    double omega_vw = p[13]*fvw*1e-3;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo*1e-3;

    // dcbo
    double tau_2cbo = p[18];
    double alpha_2cbo = p[19];
    double beta_2cbo = p[20];

    // vw-cbo
    double tau_vwmcbo = p[21];
    double alpha_vwmcbo = p[22];
    double beta_vwmcbo = p[23];
    double omega_vwmcbo = p[24];

    // vw+cbo
    double tau_vwpcbo = p[25];
    double alpha_vwpcbo = p[26];
    double beta_vwpcbo = p[27];
    double omega_vwpcbo = p[28];

    // cbo +- a
    double alpha_cbopa = p[29];
    double beta_cbopa = p[30];
    double alpha_cboma = p[31];
    double beta_cboma = p[32];

    // vw +- a
    double alpha_vwpa = p[33];
    double beta_vwpa = p[34];
    double alpha_vwma = p[35];
    double beta_vwma = p[36];

    // vo +- a
    double alpha_vopa = p[37];
    double beta_vopa = p[38];
    double alpha_voma = p[39];
    double beta_voma = p[40];

    // vw-cbo-a
    double alpha_vwmcboma = p[41];
    double beta_vwmcboma = p[42];

    // vo +- cbo
    double tau_vopmcbo = p[43];
    double alpha_vopcbo = p[44];
    double beta_vopcbo = p[45];
    double alpha_vomcbo = p[46];
    double beta_vomcbo = p[47];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - TMath::Exp(-time/tau_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - TMath::Exp(-time/tau_2cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbo*TMath::Cos(omega_vwmcbo*time)+beta_vwmcbo*TMath::Sin(omega_vwmcbo*time))
                         - TMath::Exp(-time/tau_vwpcbo) * (alpha_vwpcbo*TMath::Cos(omega_vwpcbo*time)+beta_vwpcbo*TMath::Sin(omega_vwpcbo*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         + TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         + TMath::Exp(-time/tau_vopmcbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         + TMath::Exp(-time/tau_vopmcbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

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
    double omega_vw = p[13]*fvw*1e-3;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo*1e-3;

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
    double tau_vopmcbo = p[38];
    double alpha_vopcbo = p[39];
    double beta_vopcbo = p[40];
    double alpha_vomcbo = p[41];
    double beta_vomcbo = p[42];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - TMath::Exp(-time/tau_cbo) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbo*TMath::Cos(omega_vwmcbo*time)+beta_vwmcbo*TMath::Sin(omega_vwmcbo*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         + TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         + TMath::Exp(-time/tau_vopmcbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         + TMath::Exp(-time/tau_vopmcbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

// GPR functions
double func_47paras(double *x, double *p){
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
    double omega_vw = p[13]*fvw*1e-3;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo*1e-3;

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

    // vw-cbo-a
    double alpha_vwmcboma = p[40];
    double beta_vwmcboma = p[41];

    // vo +- cbo
    double alpha_vopcbo = p[42];
    double beta_vopcbo = p[43];
    double alpha_vomcbo = p[44];
    double beta_vomcbo = p[45];

    // cbo envelope constant
    double c = p[46];

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - TMath::Exp(-time/tau_cbo+c) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - TMath::Exp(-2*time/tau_cbo) * (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbo*TMath::Cos(omega_vwmcbo*time)+beta_vwmcbo*TMath::Sin(omega_vwmcbo*time))
                         - TMath::Exp(-time/tau_vwpcbo) * (alpha_vwpcbo*TMath::Cos(omega_vwpcbo*time)+beta_vwpcbo*TMath::Sin(omega_vwpcbo*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         + TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         + TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         + TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_45paras_cboGPR(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double alpha_cbo = p[5]*(*splines)["alpha_cbo"]->Eval(time);
    double beta_cbo = p[6]*(*splines)["beta_cbo"]->Eval(time);
    double omega_cbo = p[7];

    // kloss
    double kloss = p[8];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vw
    double tau_vw = p[9];
    double alpha_vw = p[10];
    double beta_vw = p[11];
    double omega_vw = p[12]*fvw*1e-3;

    // vo
    double tau_vo = p[13];
    double alpha_vo = p[14];
    double beta_vo = p[15];
    double omega_vo = p[16]*fvo*1e-3;

    // dcbo
    double alpha_2cbo = p[17]*(*splines)["alpha_2cbo"]->Eval(time);
    double beta_2cbo = p[18]*(*splines)["beta_2cbo"]->Eval(time);

    // vw-cbo
    double tau_vwmcbo = p[19];
    double alpha_vwmcbo = p[20];
    double beta_vwmcbo = p[21];
    double omega_vwmcbo = p[22];

    // vw+cbo
    double tau_vwpcbo = p[23];
    double alpha_vwpcbo = p[24];
    double beta_vwpcbo = p[25];
    double omega_vwpcbo = p[26];

    // cbo +- a
    double alpha_cbopa = p[27]*(*splines)["alpha_cbopa"]->Eval(time);
    double beta_cbopa = p[28]*(*splines)["beta_cbopa"]->Eval(time);
    double alpha_cboma = p[29]*(*splines)["alpha_cboma"]->Eval(time);
    double beta_cboma = p[30]*(*splines)["beta_cboma"]->Eval(time);

    // vw +- a
    double alpha_vwpa = p[31];
    double beta_vwpa = p[32];
    double alpha_vwma = p[33];
    double beta_vwma = p[34];

    // vo +- a
    double alpha_vopa = p[35];
    double beta_vopa = p[36];
    double alpha_voma = p[37];
    double beta_voma = p[38];

    // vw-cbo-a
    double alpha_vwmcboma = p[39];
    double beta_vwmcboma = p[40];

    // vo +- cbo
    double alpha_vopcbo = p[41]*(*splines)["alpha_vopcbo"]->Eval(time);
    double beta_vopcbo = p[42]*(*splines)["beta_vopcbo"]->Eval(time);
    double alpha_vomcbo = p[43]*(*splines)["alpha_vomcbo"]->Eval(time);
    double beta_vomcbo = p[44]*(*splines)["beta_vomcbo"]->Eval(time);

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcbo*TMath::Cos(omega_vwmcbo*time)+beta_vwmcbo*TMath::Sin(omega_vwmcbo*time))
                         - TMath::Exp(-time/tau_vwpcbo) * (alpha_vwpcbo*TMath::Cos(omega_vwpcbo*time)+beta_vwpcbo*TMath::Sin(omega_vwpcbo*time))
                         + (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         + (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         + TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         + (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         + (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_41paras_GPR(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double alpha_cbo = p[5]*(*splines)["alpha_cbo"]->Eval(time);
    double beta_cbo = p[6]*(*splines)["beta_cbo"]->Eval(time);
    double omega_cbo = p[7];

    // kloss
    double kloss = p[8];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vw
    double alpha_vw = p[9]*(*splines)["alpha_vw"]->Eval(time);
    double beta_vw = p[10]*(*splines)["beta_vw"]->Eval(time);
    double omega_vw = p[11]*fvw*1e-3;

    // vo
    double alpha_vo = p[12]*(*splines)["alpha_vo"]->Eval(time);
    double beta_vo = p[13]*(*splines)["beta_vo"]->Eval(time);
    double omega_vo = p[14]*fvo*1e-3;

    // dcbo
    double alpha_2cbo = p[15]*(*splines)["alpha_2cbo"]->Eval(time);
    double beta_2cbo = p[16]*(*splines)["beta_2cbo"]->Eval(time);

    // vw-cbo
    double alpha_vwmcbo = p[17]*(*splines)["alpha_vwmcbo"]->Eval(time);
    double beta_vwmcbo = p[18]*(*splines)["beta_vwmcbo"]->Eval(time);
    double omega_vwmcbo = p[19];

    // vw+cbo
    double alpha_vwpcbo = p[20]*(*splines)["alpha_vwpcbo"]->Eval(time);
    double beta_vwpcbo = p[21]*(*splines)["beta_vwpcbo"]->Eval(time);
    double omega_vwpcbo = p[22];

    // cbo +- a
    double alpha_cbopa = p[23]*(*splines)["alpha_cbopa"]->Eval(time);
    double beta_cbopa = p[24]*(*splines)["beta_cbopa"]->Eval(time);
    double alpha_cboma = p[25]*(*splines)["alpha_cboma"]->Eval(time);
    double beta_cboma = p[26]*(*splines)["beta_cboma"]->Eval(time);

    // vw +- a
    double alpha_vwpa = p[27]*(*splines)["alpha_vwpa"]->Eval(time);
    double beta_vwpa = p[28]*(*splines)["beta_vwpa"]->Eval(time);
    double alpha_vwma = p[29]*(*splines)["alpha_vwma"]->Eval(time);
    double beta_vwma = p[30]*(*splines)["beta_vwma"]->Eval(time);

    // vo +- a
    double alpha_vopa = p[31]*(*splines)["alpha_vopa"]->Eval(time);
    double beta_vopa = p[32]*(*splines)["beta_vopa"]->Eval(time);
    double alpha_voma = p[33]*(*splines)["alpha_voma"]->Eval(time);
    double beta_voma = p[34]*(*splines)["beta_voma"]->Eval(time);

    // vw-cbo-a
    double alpha_vwmcboma = p[35]*(*splines)["alpha_vwmcboma"]->Eval(time);
    double beta_vwmcboma = p[36]*(*splines)["beta_vwmcboma"]->Eval(time);

    // vo +- cbo
    double alpha_vopcbo = p[37]*(*splines)["alpha_vopcbo"]->Eval(time);
    double beta_vopcbo = p[38]*(*splines)["beta_vopcbo"]->Eval(time);
    double alpha_vomcbo = p[39]*(*splines)["alpha_vomcbo"]->Eval(time);
    double beta_vomcbo = p[40]*(*splines)["beta_vomcbo"]->Eval(time);

    double expansion = 1 + asym*TMath::Cos(omega*time-phi)
                         - (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         - (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time))
                         - (alpha_2cbo*TMath::Cos(2*omega_cbo*time)+beta_2cbo*TMath::Sin(2*omega_cbo*time))
                         - (alpha_vwmcbo*TMath::Cos(omega_vwmcbo*time)+beta_vwmcbo*TMath::Sin(omega_vwmcbo*time))
                         - (alpha_vwpcbo*TMath::Cos(omega_vwpcbo*time)+beta_vwpcbo*TMath::Sin(omega_vwpcbo*time))
                         + (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         + (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         + (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         + (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         + (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         + (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         + (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         + (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         + (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_46paras_GPR(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_cbo = p[5];
    double alpha_cbo = p[6]*(*splines)["alpha_cbo"]->Eval(time);
    double beta_cbo = p[7]*(*splines)["beta_cbo"]->Eval(time);
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
    double omega_vw = p[13]*fvw*1e-3;

    // vo
    double tau_vo = p[14];
    double alpha_vo = p[15];
    double beta_vo = p[16];
    double omega_vo = p[17]*fvo*1e-3;

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

    // vw-cbo-a
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
                         + TMath::Exp(-time/tau_cbo) * (alpha_cbopa*TMath::Cos((omega_cbo+omega)*time)+beta_cbopa*TMath::Sin((omega_cbo+omega)*time))
                         + TMath::Exp(-time/tau_cbo) * (alpha_cboma*TMath::Cos((omega_cbo-omega)*time)+beta_cboma*TMath::Sin((omega_cbo-omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwpa*TMath::Cos((omega_vw+omega)*time)+beta_vwpa*TMath::Sin((omega_vw+omega)*time))
                         + TMath::Exp(-time/tau_vw) * (alpha_vwma*TMath::Cos((omega_vw-omega)*time)+beta_vwma*TMath::Sin((omega_vw-omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_vopa*TMath::Cos((omega_vo+omega)*time)+beta_vopa*TMath::Sin((omega_vo+omega)*time))
                         + TMath::Exp(-time/tau_vo) * (alpha_voma*TMath::Cos((omega_vo-omega)*time)+beta_voma*TMath::Sin((omega_vo-omega)*time))
                         + TMath::Exp(-time/tau_vwmcbo) * (alpha_vwmcboma*TMath::Cos((omega_vwmcbo-omega)*time)+beta_vwmcboma*TMath::Sin((omega_vwmcbo-omega)*time))
                         + TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vopcbo*TMath::Cos((omega_vo+omega_cbo)*time)+beta_vopcbo*TMath::Sin((omega_vo+omega_cbo)*time))
                         + TMath::Exp(-time/tau_vo-time/tau_cbo) * (alpha_vomcbo*TMath::Cos((omega_vo-omega_cbo)*time)+beta_vomcbo*TMath::Sin((omega_vo-omega_cbo)*time));

    return norm * TMath::Exp(-time/life) * (1-kloss*aloss) * expansion;
}

double func_29paras_c(double *x, double *p){
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

    // vo
    double tau_vo = p[10];
    double alpha_vo = p[11];
    double beta_vo = p[12];
    double omega_vo = p[13]*fvo*1e-3;
    double vo = 1 - TMath::Exp(-time/tau_vo) * (alpha_vo*TMath::Cos(omega_vo*time)+beta_vo*TMath::Sin(omega_vo*time));

    // vw
    double tau_vw = p[14];
    double alpha_vw = p[15];
    double beta_vw = p[16];
    double omega_vw = p[17]*fvw*1e-3;

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
    double c = p[28];

    double expansion = 1 - TMath::Exp(-time/tau_cbo+c) * (alpha_cbo*TMath::Cos(omega_cbo*time)+beta_cbo*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (alpha_vw*TMath::Cos(omega_vw*time)+beta_vw*TMath::Sin(omega_vw*time))
                         + TMath::Exp(-time/tau_cbo-time/tau_vw) * (alpha_vwpcbo*TMath::Cos((omega_vw+omega_cbo)*time)+beta_vwpcbo*TMath::Sin((omega_vw+omega_cbo)*time)
                         + alpha_vwmcbo*TMath::Cos((omega_vw-omega_cbo)*time)+beta_vwmcbo*TMath::Sin((omega_vw-omega_cbo)*time));

    return wiggle * (1-kloss*aloss) * expansion * vo * dcbo;
}