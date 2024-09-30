#include "TH1.h"
#include "TMath.h"
#include "TSpline.h"
// #include <Math/Interpolator.h>

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
TSpline3 *cboAmp;
TSpline3 *dcboAmp;
TSpline3 *AtAmp;
TSpline3 *phitAmp;
TSpline3 *vwpcboAmp;
TSpline3 *vwmcboAmp;

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

double func_9paras_phi(double *x, double *p){
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
    double A_cbo = p[6];
    double phi_cbo = p[7];
    double omega_cbo = p[8];
    double cbo = 1 - TMath::Exp(-time/tau_cbo) * A_cbo*TMath::Cos(omega_cbo*time+phi_cbo);

    return wiggle * cbo;
}

double func_10paras_phi(double *x, double *p){
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
    double A_cbo = p[6];
    double phi_cbo = p[7];
    double omega_cbo = p[8];
    double cbo = 1 - TMath::Exp(-time/tau_cbo) * A_cbo*TMath::Cos(omega_cbo*time+phi_cbo);

    // kloss
    double kloss = p[9];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    return wiggle * (1-kloss*aloss) * cbo;
}

double func_14paras_vo_phi(double *x, double *p){
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
    double A_cbo = p[6];
    double phi_cbo = p[7];
    double omega_cbo = p[8];
    double cbo = 1 - TMath::Exp(-time/tau_cbo) * A_cbo*TMath::Cos(omega_cbo*time+phi_cbo);

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
    double omega_vo = p[13]*fvo;
    double vo = 1 - TMath::Exp(-time/tau_vo) * A_vo*TMath::Cos(omega_vo*time+phi_vo);

    return wiggle * (1-kloss*aloss) * cbo * vo;
}

double func_14paras_vw_phi(double *x, double *p){
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
    double A_cbo = p[6];
    double phi_cbo = p[7];
    double omega_cbo = p[8];
    double cbo = 1 - TMath::Exp(-time/tau_cbo) * A_cbo*TMath::Cos(omega_cbo*time+phi_cbo);

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
    double A_vw = p[11];
    double phi_vw = p[12];
    double omega_vw = p[13]*fvw;
    double vw = 1 - TMath::Exp(-time/tau_vw) * A_vw*TMath::Cos(omega_vw*time+phi_vw);

    return wiggle * (1-kloss*aloss) * cbo * vw;
}

double func_18paras_phi(double *x, double *p){
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
    double A_cbo = p[6];
    double phi_cbo = p[7];
    double omega_cbo = p[8];
    double cbo = 1 - TMath::Exp(-time/tau_cbo) * A_cbo*TMath::Cos(omega_cbo*time+phi_cbo);

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
    double omega_vo = p[13]*fvo;
    double vo = 1 - TMath::Exp(-time/tau_vo) * A_vo*TMath::Cos(omega_vo*time+phi_vo);

    // vw
    double tau_vw = p[14];
    double A_vw = p[15];
    double phi_vw = p[16];
    double omega_vw = p[17]*fvw;
    double vw = 1 - TMath::Exp(-time/tau_vw) * A_vw*TMath::Cos(omega_vw*time+phi_vw);

    return wiggle * (1-kloss*aloss) * cbo * vo * vw;
}

double func_20paras_phi(double *x, double *p){
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
    double A_cbo = p[6];
    double phi_cbo = p[7];
    double omega_cbo = p[8];
    double cbo = 1 - TMath::Exp(-time/tau_cbo) * A_cbo*TMath::Cos(omega_cbo*time+phi_cbo);

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
    double omega_vo = p[13]*fvo;
    double vo = 1 - TMath::Exp(-time/tau_vo) * A_vo*TMath::Cos(omega_vo*time+phi_vo);

    // vw
    double tau_vw = p[14];
    double A_vw = p[15];
    double phi_vw = p[16];
    double omega_vw = p[17]*fvw;
    double vw = 1 - TMath::Exp(-time/tau_vw) * A_vw*TMath::Cos(omega_vw*time+phi_vw);

    // dcbo
    double A_dcbo = p[18];
    double phi_dcbo = p[19];
    double dcbo = 1 - TMath::Exp(-2*time/tau_cbo) * A_dcbo*TMath::Cos(2*omega_cbo*time+phi_dcbo);

    return wiggle * (1-kloss*aloss) * cbo * vo * vw * dcbo;
}

double func_24paras_phi(double *x, double *p){
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
    double cbo = 1 - TMath::Exp(-time/tau_cbo) * A_cbo*TMath::Cos(omega_cbo*time+phi_cbo);

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
    double omega_vo = p[13]*fvo;
    double vo = 1 - TMath::Exp(-time/tau_vo) * A_vo*TMath::Cos(omega_vo*time+phi_vo);

    // vw
    double tau_vw = p[14];
    double A_vw = p[15];
    double phi_vw = p[16];
    double omega_vw = p[17]*fvw;
    double vw = 1 - TMath::Exp(-time/tau_vw) * A_vw*TMath::Cos(omega_vw*time+phi_vw);

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

    return wiggle * (1-kloss*aloss) * cbo * vo * vw * dcbo;
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
    double omega_vo = p[13]*fvo;
    double vo = 1 - TMath::Exp(-time/tau_vo) * A_vo*TMath::Cos(omega_vo*time+phi_vo);

    // vw
    double tau_vw = p[14];
    double A_vw = p[15];
    double phi_vw = p[16];
    double omega_vw = p[17]*fvw;

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

// ***************** for cbo envelope study *****************

double func_29paras_phi(double *x, double *p){
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
    double omega_vo = p[13]*fvo;
    double vo = 1 - TMath::Exp(-time/tau_vo) * A_vo*TMath::Cos(omega_vo*time+phi_vo);

    // vw
    double tau_vw = p[14];
    double A_vw = p[15];
    double phi_vw = p[16];
    double omega_vw = p[17]*fvw;

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

    // A_cbo envelope
    double C = p[28];

    double expansion = 1 - (TMath::Exp(-time/tau_cbo)+C) * A_cbo*TMath::Cos(omega_cbo*time+phi_cbo)
                         - TMath::Exp(-time/tau_vw) * A_vw*TMath::Cos(omega_vw*time+phi_vw)
                         + TMath::Exp(-time/tau_cbo-time/tau_vw) * (A_vwpcbo*TMath::Cos((omega_vw+omega_cbo)*time+phi_vwpcbo)
                         + A_vwmcbo*TMath::Cos((omega_vw-omega_cbo)*time+phi_vwmcbo));

    return wiggle * (1-kloss*aloss) * expansion * vo * dcbo;
}

double func_12paras_dcbo_phi(double *x, double *p){
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
    double A_cbo = p[6];
    double phi_cbo = p[7];
    double omega_cbo = p[8];
    double cbo = 1 - TMath::Exp(-time/tau_cbo) * A_cbo*TMath::Cos(omega_cbo*time+phi_cbo);

    // kloss
    double kloss = p[9];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // dcbo
    double A_dcbo = p[10];
    double phi_dcbo = p[11];
    double dcbo = 1 - TMath::Exp(-2*time/tau_cbo) * A_dcbo*TMath::Cos(2*omega_cbo*time+phi_dcbo);

    return wiggle * (1-kloss*aloss) * cbo * dcbo;
}

double func_15paras_onlyCBO(double *x, double *p){
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
    double cbo = 1 - TMath::Exp(-time/tau_cbo) * A_cbo*TMath::Cos(omega_cbo*time+phi_cbo);

    // dcbo
    double A_dcbo = p[9];
    double phi_dcbo = p[10];
    double dcbo = 1 - TMath::Exp(-2*time/tau_cbo) * A_dcbo*TMath::Cos(2*omega_cbo*time+phi_dcbo);

    // At and phit
    double A_At = p[11];
    double phi_At = p[12];
    double At = 1 - TMath::Exp(-time/tau_cbo) * A_At*TMath::Cos(omega_cbo*time+phi_At);
    double A_phit = p[13];
    double phi_phit = p[14];
    double phit = 1 - TMath::Exp(-time/tau_cbo) * A_phit*TMath::Cos(omega_cbo*time+phi_phit);

    double wiggle = norm * TMath::Exp(-time/life) * (1 - asym*At*TMath::Cos(omega*time+phi*phit));

    return wiggle * cbo * dcbo;
}

double func_28paras_phi_forRun5(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_cbo = p[26];
    double A_cbo = p[27];
    double phi_cbo = p[5];
    double omega_cbo = p[6];

    // kloss
    double kloss = p[7];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vo
    double tau_vo = p[8];
    double A_vo = p[9];
    double phi_vo = p[10];
    double omega_vo = p[11]*fvo;
    double vo = 1 - TMath::Exp(-time/tau_vo) * A_vo*TMath::Cos(omega_vo*time+phi_vo);

    // vw
    double tau_vw = p[12];
    double A_vw = p[13];
    double phi_vw = p[14];
    double omega_vw = p[15]*fvw;

    // dcbo
    double A_dcbo = p[16];
    double phi_dcbo = p[17];
    double dcbo = 1 - TMath::Exp(-2*time/tau_cbo) * A_dcbo*TMath::Cos(2*omega_cbo*time+phi_dcbo);

    // At and phit
    double A_At = p[18];
    double phi_At = p[19];
    double At = 1 - TMath::Exp(-time/tau_cbo) * A_At*TMath::Cos(omega_cbo*time+phi_At);
    double A_phit = p[20];
    double phi_phit = p[21];
    double phit = 1 - TMath::Exp(-time/tau_cbo) * A_phit*TMath::Cos(omega_cbo*time+phi_phit);

    double wiggle = norm * TMath::Exp(-time/life) * (1 - asym*At*TMath::Cos(omega*time+phi*phit));

    // expansion
    double A_vwpcbo = p[22];
    double phi_vwpcbo = p[23];
    double A_vwmcbo = p[24];
    double phi_vwmcbo = p[25];

    double expansion = 1 - TMath::Exp(-time/tau_cbo) * A_cbo*TMath::Cos(omega_cbo*time+phi_cbo)
                         - TMath::Exp(-time/tau_vw) * A_vw*TMath::Cos(omega_vw*time+phi_vw)
                         + TMath::Exp(-time/tau_cbo-time/tau_vw) * (A_vwpcbo*TMath::Cos((omega_vw+omega_cbo)*time+phi_vwpcbo)
                         + A_vwmcbo*TMath::Cos((omega_vw-omega_cbo)*time+phi_vwmcbo));

    return wiggle * (1-kloss*aloss) * expansion * vo * dcbo;
}

// double func_26paras_phi_changingA(double *x, double *p){
//     double time = x[0] / time_scale;

//     // wiggle
//     double norm = p[0];
//     double life = p[1];
//     double asym = p[2];
//     double omega = getBlinded->paramToFreq(p[3]);
//     double phi = p[4];

//     // cbo
//     double A_cbo = itp->Eval(time);
//     double phi_cbo = p[5];
//     double omega_cbo = p[6];

//     // kloss
//     double kloss = p[7];
//     double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

//     // frequencies reference
//     double fcbo = 2.34;
//     double fc = 2*Pi/0.1492;
//     double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
//     double fvw = fc-2*fvo;

//     // vo
//     double tau_vo = p[8];
//     double A_vo = p[9];
//     double phi_vo = p[10];
//     double omega_vo = p[11]*fvo;
//     double vo = 1 - TMath::Exp(-time/tau_vo) * A_vo*TMath::Cos(omega_vo*time+phi_vo);

//     // vw
//     double tau_vw = p[12];
//     double A_vw = p[13];
//     double phi_vw = p[14];
//     double omega_vw = p[15]*fvw;

//     // dcbo
//     double A_dcbo = p[16];
//     double phi_dcbo = p[17];
//     double dcbo = 1 - A_dcbo*TMath::Cos(2*omega_cbo*time+phi_dcbo);

//     // At and phit
//     double A_At = p[18];
//     double phi_At = p[19];
//     double At = 1 - A_At*TMath::Cos(omega_cbo*time+phi_At);
//     double A_phit = p[20];
//     double phi_phit = p[21];
//     double phit = 1 - A_phit*TMath::Cos(omega_cbo*time+phi_phit);

//     double wiggle = norm * TMath::Exp(-time/life) * (1 - asym*At*TMath::Cos(omega*time+phi*phit));

//     // expansion
//     double A_vwpcbo = p[22];
//     double phi_vwpcbo = p[23];
//     double A_vwmcbo = p[24];
//     double phi_vwmcbo = p[25];

//     double expansion = 1 - A_cbo*TMath::Cos(omega_cbo*time+phi_cbo)
//                          - TMath::Exp(-time/tau_vw) * A_vw*TMath::Cos(omega_vw*time+phi_vw)
//                          + TMath::Exp(-time/tau_vw) * (A_vwpcbo*TMath::Cos((omega_vw+omega_cbo)*time+phi_vwpcbo)
//                          + A_vwmcbo*TMath::Cos((omega_vw-omega_cbo)*time+phi_vwmcbo));

//     return wiggle * (1-kloss*aloss) * expansion * vo * dcbo;
// }

double func_28paras_phi_minus(double *x, double *p){
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
    double omega_vo = p[13]*fvo;
    double vo = 1 - TMath::Exp(-time/tau_vo) * A_vo*TMath::Cos(omega_vo*time-phi_vo);

    // vw
    double tau_vw = p[14];
    double A_vw = p[15];
    double phi_vw = p[16];
    double omega_vw = p[17]*fvw;

    // dcbo
    double A_dcbo = p[18];
    double phi_dcbo = p[19];
    double dcbo = 1 - TMath::Exp(-2*time/tau_cbo) * A_dcbo*TMath::Cos(2*omega_cbo*time-phi_dcbo);

    // At and phit
    double A_At = p[20];
    double phi_At = p[21];
    double At = 1 - TMath::Exp(-time/tau_cbo) * A_At*TMath::Cos(omega_cbo*time-phi_At);
    double A_phit = p[22];
    double phi_phit = p[23];
    double phit = 1 - TMath::Exp(-time/tau_cbo) * A_phit*TMath::Cos(omega_cbo*time-phi_phit);

    double wiggle = norm * TMath::Exp(-time/life) * (1 + asym*At*TMath::Cos(omega*time-phi*phit));

    // expansion
    double A_vwpcbo = p[24];
    double phi_vwpcbo = p[25];
    double A_vwmcbo = p[26];
    double phi_vwmcbo = p[27];

    double expansion = 1 - TMath::Exp(-time/tau_cbo) * A_cbo*TMath::Cos(omega_cbo*time-phi_cbo)
                         - TMath::Exp(-time/tau_vw) * A_vw*TMath::Cos(omega_vw*time-phi_vw)
                         + TMath::Exp(-time/tau_cbo-time/tau_vw) * (A_vwpcbo*TMath::Cos((omega_vw+omega_cbo)*time-phi_vwpcbo)
                         + A_vwmcbo*TMath::Cos((omega_vw-omega_cbo)*time-phi_vwmcbo));

    return wiggle * (1-kloss*aloss) * expansion * vo * dcbo;
}

double func_5paras_CBO_residual(double *x, double *p){
    double time = x[0] / time_scale;

    // cbo
    double tau_cbo = p[0];
    double A_cbo = p[1];
    double phi_cbo = p[2];
    double omega_cbo = p[3];
    double offset = p[4];
    double cbo = offset + TMath::Exp(-time/tau_cbo) * A_cbo*TMath::Cos(omega_cbo*time-phi_cbo);

    return cbo;
}

double func_7paras_dCBO_residual(double *x, double *p){
    double time = x[0] / time_scale;

    // cbo
    double tau_cbo = p[0];
    double A_cbo = p[1];
    double phi_cbo = p[2];
    double omega_cbo = p[3];
    double offset = p[4];
    double cbo = offset + TMath::Exp(-time/tau_cbo) * A_cbo*TMath::Cos(omega_cbo*time-phi_cbo);

    // dcbo
    double A_dcbo = p[5];
    double phi_dcbo = p[6];
    double dcbo = TMath::Exp(-2*time/tau_cbo) * A_dcbo*TMath::Cos(2*omega_cbo*time-phi_dcbo);

    return cbo + dcbo;
}


// newCBO functions

double func_28paras_cboAmp(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double tau_cbo = p[5];
    double A_cbo = p[6]*cboAmp->Eval(time);
    double phi_cbo = p[7];
    double omega_cbo = p[8];

    // std::cout << "A_cbo: " << cboAmp->Eval(time) << std::endl;

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
    double omega_vo = p[13]*fvo;
    double vo = 1 - TMath::Exp(-time/tau_vo) * A_vo*TMath::Cos(omega_vo*time+phi_vo);

    // vw
    double tau_vw = p[14];
    double A_vw = p[15];
    double phi_vw = p[16];
    double omega_vw = p[17]*fvw;

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

    double expansion = 1 - A_cbo*TMath::Cos(omega_cbo*time+phi_cbo)
                         - TMath::Exp(-time/tau_vw) * A_vw*TMath::Cos(omega_vw*time+phi_vw)
                         + TMath::Exp(-time/tau_cbo-time/tau_vw) * (A_vwpcbo*TMath::Cos((omega_vw+omega_cbo)*time+phi_vwpcbo)
                         + A_vwmcbo*TMath::Cos((omega_vw-omega_cbo)*time+phi_vwmcbo));

    return wiggle * (1-kloss*aloss) * expansion * vo * dcbo;

    
}

double func_27paras_cboRelatedAmp(double *x, double *p){
    double time = x[0] / time_scale;

    // wiggle
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // cbo
    double A_cbo = p[5]*cboAmp->Eval(time);
    double phi_cbo = p[6];
    double omega_cbo = p[7];

    // std::cout << "A_cbo: " << cboAmp->Eval(time) << std::endl;

    // kloss
    double kloss = p[8];
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // frequencies reference
    double fcbo = 2.34;
    double fc = 2*Pi/0.1492;
    double fvo = TMath::Sqrt(fcbo*(2*fc-fcbo));
    double fvw = fc-2*fvo;

    // vo
    double tau_vo = p[9];
    double A_vo = p[10];
    double phi_vo = p[11];
    double omega_vo = p[12]*fvo;
    double vo = 1 - TMath::Exp(-time/tau_vo) * A_vo*TMath::Cos(omega_vo*time+phi_vo);

    // vw
    double tau_vw = p[13];
    double A_vw = p[14];
    double phi_vw = p[15];
    double omega_vw = p[16]*fvw;

    // dcbo
    double A_dcbo = p[17]*dcboAmp->Eval(time);
    double phi_dcbo = p[18];
    double dcbo = 1 - A_dcbo*TMath::Cos(2*omega_cbo*time+phi_dcbo);

    // At and phit
    double A_At = p[19]*AtAmp->Eval(time);
    double phi_At = p[20];
    double At = 1 - A_At*TMath::Cos(omega_cbo*time+phi_At);
    double A_phit = p[21]*phitAmp->Eval(time);
    double phi_phit = p[22];
    double phit = 1 - A_phit*TMath::Cos(omega_cbo*time+phi_phit);

    double wiggle = norm * TMath::Exp(-time/life) * (1 - asym*At*TMath::Cos(omega*time+phi*phit));

    // expansion
    double A_vwpcbo = p[23]*vwpcboAmp->Eval(time);
    double phi_vwpcbo = p[24];
    double A_vwmcbo = p[25]*vwmcboAmp->Eval(time);
    double phi_vwmcbo = p[26];

    double expansion = 1 - A_cbo*TMath::Cos(omega_cbo*time+phi_cbo)
                         - TMath::Exp(-time/tau_vw) * A_vw*TMath::Cos(omega_vw*time+phi_vw)
                         + (A_vwpcbo*TMath::Cos((omega_vw+omega_cbo)*time+phi_vwpcbo)
                         + A_vwmcbo*TMath::Cos((omega_vw-omega_cbo)*time+phi_vwmcbo));

    return wiggle * (1-kloss*aloss) * expansion * vo * dcbo;

    
}


double func_32paras(double *x, double *p){
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
    double Avopcbo_c = p[28];
    double Avopcbo_s = p[29];
    double Avomcbo_c = p[30];
    double Avomcbo_s = p[31];

    double expansion = 1 - TMath::Exp(-time/tau_cbo) * (Acbo_c*TMath::Cos(omega_cbo*time)+Acbo_s*TMath::Sin(omega_cbo*time))
                         - TMath::Exp(-time/tau_vw) * (Avw_c*TMath::Cos(omega_vw*time)+Avw_s*TMath::Sin(omega_vw*time))
                         - TMath::Exp(-time/tau_vo) * (Avo_c*TMath::Cos(omega_vo*time)+Avo_s*TMath::Sin(omega_vo*time))
                         + TMath::Exp(-time/tau_cbo-time/tau_vw) * (Avwpcbo_c*TMath::Cos((omega_vw+omega_cbo)*time)+Avwpcbo_s*TMath::Sin((omega_vw+omega_cbo)*time)
                         + Avwmcbo_c*TMath::Cos((omega_vw-omega_cbo)*time)+Avwmcbo_s*TMath::Sin((omega_vw-omega_cbo)*time))
                         + TMath::Exp(-time/tau_cbo-time/tau_vo) * (Avopcbo_c*TMath::Cos((omega_vo+omega_cbo)*time)+Avopcbo_s*TMath::Sin((omega_vo+omega_cbo)*time)
                         + Avomcbo_c*TMath::Cos((omega_vo-omega_cbo)*time)+Avomcbo_s*TMath::Sin((omega_vo-omega_cbo)*time));

    return wiggle * (1-kloss*aloss) * expansion * dcbo;
}