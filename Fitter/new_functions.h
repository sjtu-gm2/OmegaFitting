#include "TH1.h"
#include "TMath.h"

using namespace std;
using namespace blinding;

const double Pi = 3.14159265358979323846264338327950288419716939937510582;
const double pipy = 3.141592653589793;

double time_scale = 1.0;
Blinders::fitType ftype = Blinders::kOmega_a;

// #if data_version_major == 4
//     Blinders *getBlinded = new Blinders(ftype, "Unexpected virtue of ignorance.");
// #elif data_version_major == 2
//     Blinders *getBlinded = new Blinders(ftype, "stay home, stay healthy!");
// #elif data_version_major == 3
//     Blinders *getBlinded = new Blinders(ftype, "Bla Bla Bla!");
// #endif

string blindedString = "Unexpected virtue of ignorance.";
Blinders *getBlinded;

TH1 *lost_muon;

// fit function

double func_5paras(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];
    double wiggle = norm * TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time+phi));

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
    double wiggle = norm * TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time+phi));

    // cbo
    double tau_cbo = p[5];
    double Acbo_c = p[6];
    double Acbo_s = p[7];
    double omega_cbo = p[8];
    double cbo = 1 - TMath::Exp(-time/tau_cbo) * (Acbo_c*TMath::Cos(omega_cbo*time)+Acbo_s*TMath::Sin(omega_cbo*time));

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
    double wiggle = norm * TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time+phi));

    // cbo
    double tau_cbo = p[5];
    double Acbo_c = p[6];
    double Acbo_s = p[7];
    double omega_cbo = p[8];
    double cbo = 1 - TMath::Exp(-time/tau_cbo) * (Acbo_c*TMath::Cos(omega_cbo*time)+Acbo_s*TMath::Sin(omega_cbo*time));

    // kloss
    double kloss = p[9];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    return wiggle * (1-kloss*aloss) * cbo;
}

double func_14paras_vo(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];
    double wiggle = norm * TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time+phi));

    // cbo
    double tau_cbo = p[5];
    double Acbo_c = p[6];
    double Acbo_s = p[7];
    double omega_cbo = p[8];
    double cbo = 1 - TMath::Exp(-time/tau_cbo) * (Acbo_c*TMath::Cos(omega_cbo*time)+Acbo_s*TMath::Sin(omega_cbo*time));

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
    double Avo_c = p[11];
    double Avo_s = p[12];
    double omega_vo = p[13]*fvo;
    double vo = 1 - TMath::Exp(-time/tau_vo) * (Avo_c*TMath::Cos(omega_vo*time)+Avo_s*TMath::Sin(omega_vo*time));

    return wiggle * (1-kloss*aloss) * cbo * vo;
}

double func_14paras_vw(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];
    double wiggle = norm * TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time+phi));

    // cbo
    double tau_cbo = p[5];
    double Acbo_c = p[6];
    double Acbo_s = p[7];
    double omega_cbo = p[8];
    double cbo = 1 - TMath::Exp(-time/tau_cbo) * (Acbo_c*TMath::Cos(omega_cbo*time)+Acbo_s*TMath::Sin(omega_cbo*time));

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
    double Avw_c = p[11];
    double Avw_s = p[12];
    double omega_vw = p[13]*fvw;
    double vw = 1 - TMath::Exp(-time/tau_vw) * (Avw_c*TMath::Cos(omega_vw*time)+Avw_s*TMath::Sin(omega_vw*time));

    return wiggle * (1-kloss*aloss) * cbo * vw;
}

double func_18paras(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];
    double wiggle = norm * TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time+phi));

    // cbo
    double tau_cbo = p[5];
    double Acbo_c = p[6];
    double Acbo_s = p[7];
    double omega_cbo = p[8];
    double cbo = 1 - TMath::Exp(-time/tau_cbo) * (Acbo_c*TMath::Cos(omega_cbo*time)+Acbo_s*TMath::Sin(omega_cbo*time));

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
    double Avo_c = p[11];
    double Avo_s = p[12];
    double omega_vo = p[13]*fvo;
    double vo = 1 - TMath::Exp(-time/tau_vo) * (Avo_c*TMath::Cos(omega_vo*time)+Avo_s*TMath::Sin(omega_vo*time));

    // vw
    double tau_vw = p[14];
    double Avw_c = p[15];
    double Avw_s = p[16];
    double omega_vw = p[17]*fvw;
    double vw = 1 - TMath::Exp(-time/tau_vw) * (Avw_c*TMath::Cos(omega_vw*time)+Avw_s*TMath::Sin(omega_vw*time));

    return wiggle * (1-kloss*aloss) * cbo * vo * vw;
}

double func_20paras(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];
    double wiggle = norm * TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time+phi));

    // cbo
    double tau_cbo = p[5];
    double Acbo_c = p[6];
    double Acbo_s = p[7];
    double omega_cbo = p[8];
    double cbo = 1 - TMath::Exp(-time/tau_cbo) * (Acbo_c*TMath::Cos(omega_cbo*time)+Acbo_s*TMath::Sin(omega_cbo*time));

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
    double Avo_c = p[11];
    double Avo_s = p[12];
    double omega_vo = p[13]*fvo;
    double vo = 1 - TMath::Exp(-time/tau_vo) * (Avo_c*TMath::Cos(omega_vo*time)+Avo_s*TMath::Sin(omega_vo*time));

    // vw
    double tau_vw = p[14];
    double Avw_c = p[15];
    double Avw_s = p[16];
    double omega_vw = p[17]*fvw;
    double vw = 1 - TMath::Exp(-time/tau_vw) * (Avw_c*TMath::Cos(omega_vw*time)+Avw_s*TMath::Sin(omega_vw*time));

    // dcbo
    double Adcbo_c = p[18];
    double Adcbo_s = p[19];
    double dcbo = 1 - TMath::Exp(-2*time/tau_cbo) * (Adcbo_c*TMath::Cos(2*omega_cbo*time)+Adcbo_s*TMath::Sin(2*omega_cbo*time));

    return wiggle * (1-kloss*aloss) * cbo * vo * vw * dcbo;
}

double func_24paras(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_cbo = p[5];
    double Acbo_c = p[6];
    double Acbo_s = p[7];
    double omega_cbo = p[8];
    double cbo = 1 - TMath::Exp(-time/tau_cbo) * (Acbo_c*TMath::Cos(omega_cbo*time)+Acbo_s*TMath::Sin(omega_cbo*time));

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
    double Avo_c = p[11];
    double Avo_s = p[12];
    double omega_vo = p[13]*fvo;
    double vo = 1 - TMath::Exp(-time/tau_vo) * (Avo_c*TMath::Cos(omega_vo*time)+Avo_s*TMath::Sin(omega_vo*time));

    // vw
    double tau_vw = p[14];
    double Avw_c = p[15];
    double Avw_s = p[16];
    double omega_vw = p[17]*fvw;
    double vw = 1 - TMath::Exp(-time/tau_vw) * (Avw_c*TMath::Cos(omega_vw*time)+Avw_s*TMath::Sin(omega_vw*time));

    // dcbo
    double Adcbo_c = p[18];
    double Adcbo_s = p[19];
    double dcbo = 1 - TMath::Exp(-2*time/tau_cbo) * (Adcbo_c*TMath::Cos(2*omega_cbo*time)+Adcbo_s*TMath::Sin(2*omega_cbo*time));

    // At and phit
    double AAt_c = p[20];
    double AAt_s = p[21];
    double At = 1 - TMath::Exp(-time/tau_cbo) * (AAt_c*TMath::Cos(omega_cbo*time)+AAt_s*TMath::Sin(omega_cbo*time));
    double Aphit_c = p[22];
    double Aphit_s = p[23];
    double phit = 1 - TMath::Exp(-time/tau_cbo) * (Aphit_c*TMath::Cos(omega_cbo*time)+Aphit_s*TMath::Sin(omega_cbo*time));

    double wiggle = norm * TMath::Exp(-time/life) * (1 - asym*At*TMath::Cos(omega*time+phi*phit));

    return wiggle * (1-kloss*aloss) * cbo * vo * vw * dcbo;
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
    double Acbo_c = p[6];
    double Acbo_s = p[7];
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
    double Avo_c = p[11];
    double Avo_s = p[12];
    double omega_vo = p[13]*fvo;
    double vo = 1 - TMath::Exp(-time/tau_vo) * (Avo_c*TMath::Cos(omega_vo*time)+Avo_s*TMath::Sin(omega_vo*time));

    // vw
    double tau_vw = p[14];
    double Avw_c = p[15];
    double Avw_s = p[16];
    double omega_vw = p[17]*fvw;

    // dcbo
    double Adcbo_c = p[18];
    double Adcbo_s = p[19];
    double dcbo = 1 - TMath::Exp(-2*time/tau_cbo) * (Adcbo_c*TMath::Cos(2*omega_cbo*time)+Adcbo_s*TMath::Sin(2*omega_cbo*time));

    // At and phit
    double AAt_c = p[20];
    double AAt_s = p[21];
    double At = 1 - TMath::Exp(-time/tau_cbo) * (AAt_c*TMath::Cos(omega_cbo*time)+AAt_s*TMath::Sin(omega_cbo*time));
    double Aphit_c = p[22];
    double Aphit_s = p[23];
    double phit = 1 - TMath::Exp(-time/tau_cbo) * (Aphit_c*TMath::Cos(omega_cbo*time)+Aphit_s*TMath::Sin(omega_cbo*time));

    double wiggle = norm * TMath::Exp(-time/life) * (1 - asym*At*TMath::Cos(omega*time+phi*phit));

    // expansion
    double Avwpcbo_c = p[24];
    double Avwpcbo_s = p[25];
    double Avwmcbo_c = p[26];
    double Avwmcbo_s = p[27];

    double expansion = 1 - TMath::Exp(-time/tau_cbo) * (Acbo_c*TMath::Cos(omega_cbo*time)+Acbo_s*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (Avw_c*TMath::Cos(omega_vw*time)+Avw_s*TMath::Sin(omega_vw*time))
                         + TMath::Exp(-time/tau_cbo-time/tau_vw) * (Avwpcbo_c*TMath::Cos((omega_vw+omega_cbo)*time)+Avwpcbo_s*TMath::Sin((omega_vw+omega_cbo)*time)
                         + Avwmcbo_c*TMath::Cos((omega_vw-omega_cbo)*time)+Avwmcbo_s*TMath::Sin((omega_vw-omega_cbo)*time));

    return wiggle * (1-kloss*aloss) * expansion * vo * dcbo;
}

double func_6paras_kloss(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];
    double wiggle = norm * TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time+phi));

    // kloss
    double kloss = p[5];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    return wiggle * (1-kloss*aloss);
}

double func_11paras_cbo(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];
    double wiggle = norm * TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time+phi));

    // cbo
    double tau_cbo = p[5];
    double c_cbo = p[6];
    double Acbo_c = p[7];
    double Acbo_s = p[8];
    double omega_cbo = p[9];
    double cbo = 1 - (TMath::Exp(-time/tau_cbo) + c_cbo) * (Acbo_c*TMath::Cos(omega_cbo*time)+Acbo_s*TMath::Sin(omega_cbo*time));

    // kloss
    double kloss = p[10];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    return wiggle * (1-kloss*aloss) * cbo;
}

double func_13paras_res(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];
    double wiggle = norm * TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time+phi));

    // cbo
    double tau_cbo = p[5];
    double c_cbo = p[6];
    double Acbo_c = p[7];
    double Acbo_s = p[8];
    double omega_cbo = p[9];
    double cbo = 1 - (TMath::Exp(-time/tau_cbo) + c_cbo) * (Acbo_c*TMath::Cos(omega_cbo*time)+Acbo_s*TMath::Sin(omega_cbo*time));

    // kloss
    double kloss = p[10];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // residual
    double res_norm = p[11];
    double res_tau = p[12];
    double res = (1 + res_norm * TMath::Exp(-time/res_tau));

    return wiggle * (1-kloss*aloss) * cbo * res;
}