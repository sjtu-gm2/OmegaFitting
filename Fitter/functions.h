using namespace std;
using namespace blinding;

const double Pi = 3.14159265358979323846264338327950288419716939937510582;
const double pipy = 3.141592653589793;

double time_scale = 1.0;
Blinders::fitType ftype = Blinders::kOmega_a;

#if data_version_major == 4
Blinders * getBlinded = new Blinders(ftype, "Unexpected virtue of ignorance.");
#elif data_version_major == 2
Blinders * getBlinded = new Blinders(ftype, "stay home, stay healthy!");
#elif data_version_major == 3
Blinders * getBlinded = new Blinders(ftype, "Bla Bla Bla!");
#endif

TH1 * lost_muon;

//fit functions

double func_5paras(double *x,double *p) {
  double time = x[0] / time_scale;
  double norm = p[0];
  double life = p[1];
  double asym = p[2];  
  double omega = getBlinded->paramToFreq(p[3]);
  double phi = p[4];

  return norm * TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time + phi));
}

double func_9paras_cbo(double *x, double *p) {
    double time = x[0] / time_scale;

    // 5 paras
    double norm = p[0];
    double life = p[1];
    double asym = p[2];  
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // 4 paras: cbo
    double tau_cbo = p[5]; //fix to infinity
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];
    
    double cbo = 1-TMath::Exp(-time/tau_cbo)*asym_cbo*TMath::Cos(omega_cbo*time + phi_cbo);
    return  norm * TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time + phi)) *  cbo;
}


double func_10paras_cbo_lost(double *x, double *p) {
    double time = x[0] / time_scale;
    
    // 5 paras
    double norm = p[0];
    double life = p[1];
    double asym = p[2];  
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // 4 paras: cbo
    double tau_cbo = p[5];
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];

    // k_loss
    double k = p[9]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    double cbo = 1-TMath::Exp(-time/tau_cbo)*asym_cbo*TMath::Cos(omega_cbo*time + phi_cbo);

    return  norm *(1 - k*aloss)* TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time + phi)) *  cbo;
}


double func_14paras_cbo_lost_vw(double *x, double *p) {
    double time = x[0] / time_scale;
    
    // 5 paras
    double norm = p[0];
    double life = p[1];
    double asym = p[2];  
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // 4 paras: cbo
    double tau_cbo = p[5];
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];

    // k_loss
    double k = p[9]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // 4 paras: vw
    double fcbo = 2.32657;
    double fc = 2*M_PI/0.1492;
    double fvo = sqrt(fcbo*(2*fc - fcbo));
    double fvw = fc - 2*fvo;

    double tau_vw = p[10];
    double asym_vw = p[11];
    double omega_vw = p[12]*fvw;
    double phi_vw = p[13];




    double cbo = 1-TMath::Exp(-time/tau_cbo)*asym_cbo*TMath::Cos(omega_cbo*time + phi_cbo);
    double vw  = 1-TMath::Exp(-time/tau_vw)*asym_vw*TMath::Cos(omega_vw*time + phi_vw);
    return  norm *(1 - k*aloss)* TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time + phi)) *  cbo * vw;
}


double func_14paras_cbo_lost_vo(double *x, double *p) {
    double time = x[0] / time_scale;
    
    // 5 paras
    double norm = p[0];
    double life = p[1];
    double asym = p[2];  
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // 4 paras: cbo
    double tau_cbo = p[5];
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];

    // k_loss
    double k = p[9]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // 4 paras: vw
    double fcbo = 2.34;
    double fc = 2*M_PI/0.1492;
    double fvo = sqrt(fcbo*(2*fc - fcbo));
    double fvw = fc - 2*fvo;

    double tau_vo = p[10];
    double asym_vo = p[11];
    double omega_vo = p[12]*fvo;
    double phi_vo = p[13];


    

    double cbo = 1-TMath::Exp(-time/tau_cbo)*asym_cbo*TMath::Cos(omega_cbo*time + phi_cbo);
    double vo  = vo = 1 - exp(-time/tau_vo)*asym_vo*cos(omega_vo*time + phi_vo);
    return  norm *(1 - k*aloss)* TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time + phi)) *  cbo * vo;
}




double func_18paras_cbo_lost_vo_vw(double *x, double *p) {
    double time = x[0] / time_scale;
    // 5-par
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double R = p[3];
    double phi = p[4];
    double omega = getBlinded->paramToFreq(R);
    
    // cbo-par
    double tau_cbo = p[5];
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];
    double cbo = 1-TMath::Exp(-time/tau_cbo)*asym_cbo*TMath::Cos(omega_cbo*time + phi_cbo);
    
    double fcbo = 2.32657;
    double fc = 2*M_PI/0.1492;
    double fvo = sqrt(fcbo*(2*fc - fcbo));
    double fvw = fc - 2*fvo;
    
    // k_loss
    double k = p[9]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    
    // vo-par
    double tau_vo = p[10];
    double asym_vo = p[11];
    double omega_vo = p[12]*fvo;
    double phi_vo = p[13];
    double vo = 1 - exp(-time/tau_vo)*asym_vo*cos(omega_vo*time + phi_vo);

    // vw-par
    double tau_vw = p[14];
    double asym_vw = p[15];
    double omega_vw = p[16]*fvw;
    double phi_vw = p[17];
    double vw = 1-exp(-time/tau_vw)*asym_vw*cos(omega_vw*time + phi_vw);

    return (1 - k*aloss) * norm * exp(-time/life) * (1 - asym*cos(omega*time + phi)) * cbo * vw * vo;
}

//1.9MHz func, 28 paras
double func_28paras_cbo_lost_vw_expansion(double *x, double *p) {
    double time = x[0] / time_scale;
    // 5-par
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double R = p[3];
    double phi = p[4];
    double omega = getBlinded->paramToFreq(R);
    
    // cbo-par
    double tau_cbo = p[5];
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];
    
    double fcbo = 2.32657;
    double fc = 2*M_PI/0.1492;
    double fvo = sqrt(fcbo*(2*fc - fcbo));
    double fvw = fc - 2*fvo;
    
    // k_loss
    double k = p[9]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // vw-par
    double tau_vw = p[10];
    double asym_vw = p[11];
    double omega_vw = p[12]*fvw;
    double phi_vw = p[13];

    // expansion-par
    double asym_vwcbo = p[24];
    double phi_vwcbo = p[25];
    double asym_vw_cbo = p[26];
    double phi_vw_cbo = p[27];
    double expan = (1 - exp(-time/tau_cbo)*asym_cbo*cos(omega_cbo*time + phi_cbo) - exp(-time/tau_vw)*asym_vw*cos(omega_vw*time + phi_vw) + exp(-time/tau_cbo - time/tau_vw)*(asym_vwcbo*cos((omega_vw + omega_cbo)*time + phi_vwcbo) + asym_vw_cbo*cos((omega_vw - omega_cbo)*time + phi_vw_cbo)));
    
    // dcbo-par
    double asym_dcbo = p[14];
    double phi_dcbo = p[15];
    double dcbo = 1 - exp(-2*time/tau_cbo)*asym_dcbo*cos(2*omega_cbo*time + phi_dcbo);
    
    // vo-par
    double tau_vo = p[20];
    double asym_vo = p[21];
    double omega_vo = p[22]*fvo;
    double phi_vo = p[23];
    double vo = 1 - exp(-time/tau_vo)*asym_vo*cos(omega_vo*time + phi_vo);
    
    // modification of A and phi
    double A1_cbo = p[16];
    double phi1_cbo = p[17];
    double A2_cbo = p[18];
    double phi2_cbo = p[19];
    double At = 1 - A1_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi1_cbo);
    double phit = 1 - A2_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi2_cbo);
    return (1 - k*aloss) * norm * exp(-time/life) * (1 - asym*At*cos(omega*time + phi*phit)) * expan * dcbo * vo;
}

//1.9MHz func, 22 paras
double func_22paras_cbo_lost_vw_expansion_lite(double *x, double *p) {
    
    double time = x[0] / time_scale;

    int nvar = 0;
    // 5-par
    double norm = p[nvar++];
    double life = p[nvar++];
    double asym = p[nvar++];
    double R = p[nvar++];
    double phi = p[nvar++];
    double omega = getBlinded->paramToFreq(R);
    
    // cbo-par
    double tau_cbo = p[nvar++];
    double asym_cbo = p[nvar++];
    double omega_cbo = p[nvar++];
    double phi_cbo = p[nvar++];
    
    double fcbo = 2.32657;
    double fc = 2*M_PI/0.1492;
    double fvo = sqrt(fcbo*(2*fc - fcbo));
    double fvw = fc - 2*fvo;
    
    // k_loss
    double k = p[nvar++]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // vw-par
    double tau_vw = p[nvar++];
    double asym_vw = p[nvar++];
    double omega_vw = p[nvar++]*fvw;
    double phi_vw = p[nvar++];

    // expansion-par
    // double asym_vwcbo = p[24];
    // double phi_vwcbo = p[25];
    double asym_vw_cbo = p[nvar++];
    double phi_vw_cbo = p[nvar++];
    double expan = (1 - exp(-time/tau_cbo)*asym_cbo*cos(omega_cbo*time + phi_cbo) - exp(-time/tau_vw)*asym_vw*cos(omega_vw*time + phi_vw) + exp(-time/tau_cbo - time/tau_vw)*(asym_vw_cbo*cos((omega_vw - omega_cbo)*time + phi_vw_cbo)));
    
    // // dcbo-par
    // // double asym_dcbo = p[14];
    // // double phi_dcbo = p[15];
    double dcbo = 1; // - exp(-2*time/tau_cbo)*asym_dcbo*cos(2*omega_cbo*time + phi_dcbo);
    
    // vo-par
    double tau_vo = p[nvar++];
    double asym_vo = p[nvar++];
    double omega_vo = p[nvar++]*fvo;
    double phi_vo = p[nvar++];
    double vo = 1 - exp(-time/tau_vo)*asym_vo*cos(omega_vo*time + phi_vo);
    
    // modification of A and phi
    double A1_cbo = p[nvar++];
    double phi1_cbo = p[nvar++];
    // double A2_cbo = p[18]; 
    // double phi2_cbo = p[19];
    double At = 1 - A1_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi1_cbo);
    double phit = 1; //- A2_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi2_cbo);

    return (1 - k*aloss) * norm * exp(-time/life) * (1 - asym*At*cos(omega*time + phi*phit)) * expan * dcbo * vo;
}


//frequency changing cbo per calo fit
// 5 + 4 cbo + 2 exp + 1 kloss
double func_12paras_changing_cbo(double *x, double *p) {
    double time = x[0] / time_scale;

    int nvar = 0;
    // 5-par
    double norm = p[nvar++];
    double life = p[nvar++];
    double asym = p[nvar++];
    double R = p[nvar++];
    double phi = p[nvar++];
    double omega = getBlinded->paramToFreq(R);
    
    // cbo-par
    double tau_cbo = p[nvar++];
    double asym_cbo = p[nvar++];
    double omega_cbo = p[nvar++];
    double phi_cbo = p[nvar++];
    

    // k_loss
    double k = p[nvar++]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // time chaging cbo terms
    double a1 = p[nvar++];
    double tau1 = p[nvar++];

    // double fcbo = 2.32657;
    // double fc = 2*M_PI/0.1492;
    // double fvo = sqrt(fcbo*(2*fc - fcbo));
    // double fvw = fc - 2*fvo;
    double cbo = 1-TMath::Exp(-time/tau_cbo)*asym_cbo*TMath::Cos(omega_cbo*time + a1*exp(-time/tau1)+ phi_cbo);
    return (1 - k*aloss) * norm * exp(-time/life) * (1 - asym*cos(omega*time + phi)) *cbo;
}

double func_11paras_changing_cbo(double *x, double *p) {
    double time = x[0] / time_scale;

    int nvar = 0;
    // 5-par
    double norm = p[nvar++];
    double life = p[nvar++];
    double asym = p[nvar++];
    double R = p[nvar++];
    double phi = p[nvar++];
    double omega = getBlinded->paramToFreq(R);
    
    // cbo-par
    double tau_cbo = p[nvar++];
    double asym_cbo = p[nvar++];
    double omega_cbo = p[nvar++];
    double phi_cbo = p[nvar++];
    

    // k_loss
    // double k = p[nvar++]*1e-9;
    // double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // time chaging cbo terms
    double a1 = p[nvar++];
    double tau1 = p[nvar++];

    // double fcbo = 2.32657;
    // double fc = 2*M_PI/0.1492;
    // double fvo = sqrt(fcbo*(2*fc - fcbo));
    // double fvw = fc - 2*fvo;
    double cbo = 1-TMath::Exp(-time/tau_cbo)*asym_cbo*TMath::Cos(omega_cbo*time + a1*exp(-time/tau1)+ phi_cbo);
    return norm * exp(-time/life) * (1 - asym*cos(omega*time + phi)) *cbo;
}

double func_13paras_cbo_vo(double *x, double *p) {
    double time = x[0] / time_scale;
    
    // 5 paras
    double norm = p[0];
    double life = p[1];
    double asym = p[2];  
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // 4 paras: cbo
    double tau_cbo = p[5];
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];

    // // k_loss
    // double k = p[9]*1e-9;
    // double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // 4 paras: vw
    double fcbo = 2.34;
    double fc = 2*M_PI/0.1492;
    double fvo = sqrt(fcbo*(2*fc - fcbo));
    double fvw = fc - 2*fvo;

    double tau_vo = p[9];
    double asym_vo = p[10];
    double omega_vo = p[11]*fvo;
    double phi_vo = p[12];

    double cbo = 1-TMath::Exp(-time/tau_cbo)*asym_cbo*TMath::Cos(omega_cbo*time + phi_cbo);
    double vo  = vo = 1 - exp(-time/tau_vo)*asym_vo*cos(omega_vo*time + phi_vo);
    return  norm * TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time + phi)) *  cbo * vo;
}

double func_15paras_changing_cbo_vo(double *x, double *p) {
    double time = x[0] / time_scale;
    
    // 5 paras
    double norm = p[0];
    double life = p[1];
    double asym = p[2];  
    double omega = getBlinded->paramToFreq(p[3]);
    double phi = p[4];

    // 4 paras: cbo
    double tau_cbo = p[5];
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];

    // // k_loss
    // double k = p[9]*1e-9;
    // double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // 4 paras: vw
    double fcbo = 2.34;
    double fc = 2*M_PI/0.1492;
    double fvo = sqrt(fcbo*(2*fc - fcbo));
    double fvw = fc - 2*fvo;

    double tau_vo = p[9];
    double asym_vo = p[10];
    double omega_vo = p[11]*fvo;
    double phi_vo = p[12];

    // double cbo = 1-TMath::Exp(-time/tau_cbo)*asym_cbo*TMath::Cos(omega_cbo*time + a1*exp(-time/tau1)+ phi_cbo);
    double cbo_a1 = p[13];
    double cbo_tau1 = p[14];

    double cbo = 1-TMath::Exp(-time/tau_cbo)*asym_cbo*TMath::Cos(omega_cbo*time + cbo_a1*exp(-time/cbo_tau1)+ phi_cbo);
    double vo  = vo = 1 - exp(-time/tau_vo)*asym_vo*cos(omega_vo*time + phi_vo);
    return  norm * TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time + phi)) *  cbo * vo;
}

//1.9MHz func, 29 paras
double func_29paras_cbo_envelope_C(double *x, double *p) {
    double time = x[0] / time_scale;
    // 5-par
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double R = p[3];
    double phi = p[4];
    double omega = getBlinded->paramToFreq(R);
    
    // cbo-par
    double tau_cbo = p[5];
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];
    
    double fcbo = 2.32657;
    double fc = 2*M_PI/0.1492;
    double fvo = sqrt(fcbo*(2*fc - fcbo));
    double fvw = fc - 2*fvo;
    
    // k_loss
    double k = p[9]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // vw-par
    double tau_vw = p[10];
    double asym_vw = p[11];
    double omega_vw = p[12]*fvw;
    double phi_vw = p[13];

    // expansion-par
    double asym_vwcbo = p[24];
    double phi_vwcbo = p[25];
    double asym_vw_cbo = p[26];
    double phi_vw_cbo = p[27];
    double cbo_constant = p[28];

    double expan = (1 - (exp(-time/tau_cbo) + cbo_constant)*asym_cbo*cos(omega_cbo*time + phi_cbo) - exp(-time/tau_vw)*asym_vw*cos(omega_vw*time + phi_vw) + exp(-time/tau_cbo - time/tau_vw)*(asym_vwcbo*cos((omega_vw + omega_cbo)*time + phi_vwcbo) + asym_vw_cbo*cos((omega_vw - omega_cbo)*time + phi_vw_cbo)));
    
    // dcbo-par
    double asym_dcbo = p[14];
    double phi_dcbo = p[15];
    double dcbo = 1 - exp(-2*time/tau_cbo)*asym_dcbo*cos(2*omega_cbo*time + phi_dcbo);
    
    // vo-par
    double tau_vo = p[20];
    double asym_vo = p[21];
    double omega_vo = p[22]*fvo;
    double phi_vo = p[23];
    double vo = 1 - exp(-time/tau_vo)*asym_vo*cos(omega_vo*time + phi_vo);
    
    // modification of A and phi
    double A1_cbo = p[16];
    double phi1_cbo = p[17];
    double A2_cbo = p[18];
    double phi2_cbo = p[19];
    double At = 1 - A1_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi1_cbo);
    double phit = 1 - A2_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi2_cbo);
    return (1 - k*aloss) * norm * exp(-time/life) * (1 - asym*At*cos(omega*time + phi*phit)) * expan * dcbo * vo;
}

//1.9MHz func, 30 paras
double func_30paras_cbo_freq(double *x, double *p) {
    double time = x[0] / time_scale;
    // 5-par
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double R = p[3];
    double phi = p[4];
    double omega = getBlinded->paramToFreq(R);
    
    // cbo-par
    double tau_cbo = p[5];
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];
    
    double fcbo = 2.32657;
    double fc = 2*M_PI/0.1492;
    double fvo = sqrt(fcbo*(2*fc - fcbo));
    double fvw = fc - 2*fvo;
    
    // k_loss
    double k = p[9]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // vw-par
    double tau_vw = p[10];
    double asym_vw = p[11];
    double omega_vw = p[12]*fvw;
    double phi_vw = p[13];

    // expansion-par
    double asym_vwcbo = p[24];
    double phi_vwcbo = p[25];
    double asym_vw_cbo = p[26];
    double phi_vw_cbo = p[27];
    double a_cbo_residual = p[28];
    double tau_cbo_residual = p[29];

    double expan = (1 - exp(-time/tau_cbo)*asym_cbo*cos(omega_cbo*time + phi_cbo + a_cbo_residual*exp(-time/tau_cbo_residual)) - exp(-time/tau_vw)*asym_vw*cos(omega_vw*time + phi_vw) + exp(-time/tau_cbo - time/tau_vw)*(asym_vwcbo*cos((omega_vw + omega_cbo)*time + phi_vwcbo) + asym_vw_cbo*cos((omega_vw - omega_cbo)*time + phi_vw_cbo)));
    
    // dcbo-par
    double asym_dcbo = p[14];
    double phi_dcbo = p[15];
    double dcbo = 1 - exp(-2*time/tau_cbo)*asym_dcbo*cos(2*omega_cbo*time + phi_dcbo);
    
    // vo-par
    double tau_vo = p[20];
    double asym_vo = p[21];
    double omega_vo = p[22]*fvo;
    double phi_vo = p[23];
    double vo = 1 - exp(-time/tau_vo)*asym_vo*cos(omega_vo*time + phi_vo);
    
    // modification of A and phi
    double A1_cbo = p[16];
    double phi1_cbo = p[17];
    double A2_cbo = p[18];
    double phi2_cbo = p[19];
    double At = 1 - A1_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi1_cbo);
    double phit = 1 - A2_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi2_cbo);
    return (1 - k*aloss) * norm * exp(-time/life) * (1 - asym*At*cos(omega*time + phi*phit)) * expan * dcbo * vo;
}

//1.9MHz func, 31 paras
double func_31paras_cbo_time(double *x, double *p) {
    double time = x[0] / time_scale;
    // 5-par
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double R = p[3];
    double phi = p[4];
    double omega = getBlinded->paramToFreq(R);
    
    // cbo-par
    double tau_cbo = p[5];
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];
    
    double fcbo = 2.32657;
    double fc = 2*M_PI/0.1492;
    double fvo = sqrt(fcbo*(2*fc - fcbo));
    double fvw = fc - 2*fvo;
    
    // k_loss
    double k = p[9]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // vw-par
    double tau_vw = p[10];
    double asym_vw = p[11];
    double omega_vw = p[12]*fvw;
    double phi_vw = p[13];

    // expansion-par
    double asym_vwcbo = p[24];
    double phi_vwcbo = p[25];
    double asym_vw_cbo = p[26];
    double phi_vw_cbo = p[27];

    //cbo time constants
    double tau_cbo_a = p[28];
    double tau_cbo_phi = p[29];
    double tau_cbo_2a = p[30];
    

    double expan = (1 - exp(-time/tau_cbo)*asym_cbo*cos(omega_cbo*time + phi_cbo) - exp(-time/tau_vw)*asym_vw*cos(omega_vw*time + phi_vw) + exp(-time/tau_cbo - time/tau_vw)*(asym_vwcbo*cos((omega_vw + omega_cbo)*time + phi_vwcbo) + asym_vw_cbo*cos((omega_vw - omega_cbo)*time + phi_vw_cbo)));
    
    // dcbo-par
    double asym_dcbo = p[14];
    double phi_dcbo = p[15];
    double dcbo = 1 - exp(-2*time/tau_cbo_2a)*asym_dcbo*cos(2*omega_cbo*time + phi_dcbo);
    
    // vo-par
    double tau_vo = p[20];
    double asym_vo = p[21];
    double omega_vo = p[22]*fvo;
    double phi_vo = p[23];
    double vo = 1 - exp(-time/tau_vo)*asym_vo*cos(omega_vo*time + phi_vo);
    
    // modification of A and phi
    double A1_cbo = p[16];
    double phi1_cbo = p[17];
    double A2_cbo = p[18];
    double phi2_cbo = p[19];
    
    
    double At = 1 - A1_cbo*exp(-time/tau_cbo_a)*cos(omega_cbo*time + phi1_cbo);
    double phit = 1 - A2_cbo*exp(-time/tau_cbo_phi)*cos(omega_cbo*time + phi2_cbo);
    return (1 - k*aloss) * norm * exp(-time/life) * (1 - asym*At*cos(omega*time + phi*phit)) * expan * dcbo * vo;
}


//simplified func for calorimeter fit
//1.9MHz func, 22 paras
double func_simplified_22paras_calos(double *x, double *p) {
    double time = x[0] / time_scale;
    // 5-par
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double R = p[3];
    double phi = p[4];
    double omega = getBlinded->paramToFreq(R);
    
    // cbo-par
    double tau_cbo = p[5];
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];
    
    double fcbo = 2.32657;
    double fc = 2*M_PI/0.1492;
    double fvo = sqrt(fcbo*(2*fc - fcbo));
    double fvw = fc - 2*fvo;

    // k_loss
    double k = p[9]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // vw-par
    double tau_vw = p[10];
    double asym_vw = p[11];
    double omega_vw = p[12]*fvw;
    double phi_vw = p[13];

    // dcbo-par
    double asym_dcbo = p[14];
    double phi_dcbo = p[15];

    // vo-par
    double tau_vo = p[16];
    double asym_vo = p[17];
    double omega_vo = p[18]*fvo;
    double phi_vo = p[19];

    // expansion-par        
    //cbo+vw, consistent with 0
    // double asym_vwcbo = p[24];
    // double phi_vwcbo = p[25];
    //cbo-vw,

    double asym_vw_cbo = p[20];
    double phi_vw_cbo = p[21];
    
    double expan = (1 - exp(-time/tau_cbo)*asym_cbo*cos(omega_cbo*time + phi_cbo) - exp(-time/tau_vw)*asym_vw*cos(omega_vw*time + phi_vw) + exp(-time/tau_cbo - time/tau_vw)*( asym_vw_cbo*cos((omega_vw - omega_cbo)*time + phi_vw_cbo)));
    
    
    double dcbo = 1 - exp(-2*time/tau_cbo)*asym_dcbo*cos(2*omega_cbo*time + phi_dcbo);
    
    
    double vo = 1 - exp(-time/tau_vo)*asym_vo*cos(omega_vo*time + phi_vo);
    
    // modification of A and phi
    // double A1_cbo = p[16];
    // double phi1_cbo = p[17];
    // double A2_cbo = p[18];
    // double phi2_cbo = p[19];

    // double At = 1 - A1_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi1_cbo);
    // double phit = 1 - A2_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi2_cbo);
    return (1 - k*aloss) * norm * exp(-time/life) * (1 - asym*cos(omega*time + phi)) * expan * dcbo * vo;
}

double func_calos_cbo(double *x, double *p){
    double time = x[0] / time_scale;
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double R = p[3];
    double phi = p[4];
    double omega = getBlinded->paramToFreq(R);

    double life_cbo = p[5];
    double asym_cbo = p[6];

    //cbo changing freq.
    double a_cbo_ch = p[13];
    double tau_cbo_ch = p[14];

    // double omega_cbo = p[7]*(1+a_cbo_ch*exp(-time/tau_cbo_ch));
    double omega_cbo = p[7];
    double phi_cbo = p[8];
    

    double k=p[9]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    double A1_cbo = p[10];
    double phi1_cbo = p[11];
    double cbo_constant = p[12];
    
    double At = 1 - A1_cbo*exp(-time/life_cbo)*cos(omega_cbo*time + phi1_cbo);


    // double cbo=1-(exp(-time/life_cbo)+cbo_constant)*asym_cbo*cos(omega_cbo*time + phi_cbo);
    double cbo=1-(exp(-time/life_cbo)+cbo_constant)*asym_cbo*cos(omega_cbo*time + phi_cbo + a_cbo_ch*exp(-time/tau_cbo_ch));

    return (1-k*aloss)*norm*exp(-time/life)*(1-asym*At*cos(omega*time+phi))*cbo;
}

//IRMA style
//all substitute
double func_31paras_calos_1(double *x, double *p) {
    double time = x[0] / time_scale;
    // 5-par
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double R = p[3];
    double phi = p[4];
    double omega = getBlinded->paramToFreq(R);
    
    // cbo-par
    double tau_cbo = p[5];
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];
    
    double fcbo = 2.32657;
    double fc = 2*M_PI/0.1492;
    double fvo = sqrt(fcbo*(2*fc - fcbo));
    double fvw = fc - 2*fvo;
    
    // k_loss
    double k = p[9]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // vw-par
    double tau_vw = p[10];
    double asym_vw = p[11];
    double omega_vw = p[12]*fvw;
    double phi_vw = p[13];

    // expansion-par
    double asym_vwcbo = p[24];
    double phi_vwcbo = p[25];
    double asym_vw_cbo = p[26];
    double phi_vw_cbo = p[27];
    double cbo_constant = p[28];
    double cbo_a_residual = p[29];
    double cbo_tau_residual = p[30];


    double expan = (1 - (exp(-time/tau_cbo)+cbo_constant)*asym_cbo*cos(omega_cbo*time + phi_cbo + cbo_a_residual*exp(-time/cbo_tau_residual)) - exp(-time/tau_vw)*asym_vw*cos(omega_vw*time + phi_vw) + ((exp(-time/tau_cbo) + cbo_constant) * exp(- time/tau_vw))*(asym_vwcbo*cos((omega_vw + omega_cbo)*time + phi_vwcbo + cbo_a_residual*exp(-time/cbo_tau_residual)) + asym_vw_cbo*cos((omega_vw - omega_cbo)*time + phi_vw_cbo + cbo_a_residual*exp(-time/cbo_tau_residual))));
    
    // dcbo-par
    double asym_dcbo = p[14];
    double phi_dcbo = p[15];
    double dcbo = 1 - (exp(-2*time/tau_cbo)+cbo_constant)*asym_dcbo*cos(2*omega_cbo*time + phi_dcbo + cbo_a_residual*exp(-time/cbo_tau_residual));

    // vo-par
    double tau_vo = p[20];
    double asym_vo = p[21];
    double omega_vo = p[22]*fvo;
    double phi_vo = p[23];
    double vo = 1 - exp(-time/tau_vo)*asym_vo*cos(omega_vo*time + phi_vo);
    
    // modification of A and phi
    double A1_cbo = p[16];
    double phi1_cbo = p[17];
    double A2_cbo = p[18];
    double phi2_cbo = p[19];
    double At = 1 - A1_cbo*(exp(-time/tau_cbo)+cbo_constant)*cos(omega_cbo*time + phi1_cbo + cbo_a_residual*exp(-time/cbo_tau_residual));
    double phit = 1 - A2_cbo*(exp(-time/tau_cbo)+cbo_constant)*cos(omega_cbo*time + phi2_cbo + cbo_a_residual*exp(-time/cbo_tau_residual));
    return (1 - k*aloss) * norm * exp(-time/life) * (1 - asym*At*cos(omega*time + phi*phit)) * expan * dcbo * vo;

}

//old style
//just substitute the A_cbo^{N}
double func_31paras_calos_2(double *x, double *p) {
    double time = x[0] / time_scale;
    // 5-par
    double norm = p[0];
    double life = p[1];
    double asym = p[2];
    double R = p[3];
    double phi = p[4];
    double omega = getBlinded->paramToFreq(R);
    
    // cbo-par
    double tau_cbo = p[5];
    double asym_cbo = p[6];
    double omega_cbo = p[7];
    double phi_cbo = p[8];
    
    double fcbo = 2.32657;
    double fc = 2*M_PI/0.1492;
    double fvo = sqrt(fcbo*(2*fc - fcbo));
    double fvw = fc - 2*fvo;
    
    // k_loss
    double k = p[9]*1e-9;
    double aloss = lost_muon->GetBinContent((int)(time/0.1492)+1);

    // vw-par
    double tau_vw = p[10];
    double asym_vw = p[11];
    double omega_vw = p[12]*fvw;
    double phi_vw = p[13];

    // expansion-par
    double asym_vwcbo = p[24];
    double phi_vwcbo = p[25];
    double asym_vw_cbo = p[26];
    double phi_vw_cbo = p[27];
    double cbo_constant = p[28];
    double cbo_a_residual = p[29];
    double cbo_tau_residual = p[30];

    double cbo_residual = cbo_a_residual*exp(-time/cbo_tau_residual);


    double expan = (1 - (exp(-time/tau_cbo)+cbo_constant)*asym_cbo*cos(omega_cbo*time + phi_cbo + cbo_residual) - exp(-time/tau_vw)*asym_vw*cos(omega_vw*time + phi_vw) + exp(-time/tau_cbo - time/tau_vw)*(asym_vwcbo*cos((omega_vw + omega_cbo)*time + phi_vwcbo) + asym_vw_cbo*cos((omega_vw - omega_cbo)*time + phi_vw_cbo)));
    
    // dcbo-par
    double asym_dcbo = p[14];
    double phi_dcbo = p[15];
    double dcbo = 1 - exp(-2*time/tau_cbo)*asym_dcbo*cos(2*omega_cbo*time + phi_dcbo);
    
    // vo-par
    double tau_vo = p[20];
    double asym_vo = p[21];
    double omega_vo = p[22]*fvo;
    double phi_vo = p[23];
    double vo = 1 - exp(-time/tau_vo)*asym_vo*cos(omega_vo*time + phi_vo);
    
    // modification of A and phi
    double A1_cbo = p[16];
    double phi1_cbo = p[17];
    double A2_cbo = p[18];
    double phi2_cbo = p[19];
    double At = 1 - A1_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi1_cbo);
    double phit = 1 - A2_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi2_cbo);
    return (1 - k*aloss) * norm * exp(-time/life) * (1 - asym*At*cos(omega*time + phi*phit)) * expan * dcbo * vo;

}

//fit pseudo data
//5pars
double func_pseudo_5pars(double *x,double *p) {
  double time = x[0] / time_scale;
  double norm = p[0];
  double life = p[1];
  double asym = p[2];
  double omega = 2*pipy*0.2291*(1+p[3]*1e-6);
  double phi = p[4];

  return norm * TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time + phi));
}

//5pars+cbo+c
double func_pseudo_10pars(double *x,double *p) {
  double time = x[0] / time_scale;
  double norm = p[0];
  double life = p[1];
  double asym = p[2];
  double omega = 2*pipy*0.2291*(1+p[3]*1e-6);
  double phi = p[4];

  double tau_cbo = p[5]; //fix to infinity
  double asym_cbo = p[6];
  double omega_cbo = p[7];
  double phi_cbo = p[8];
  double c_cbo = p[9];

  double cbo = 1-(TMath::Exp(-time/tau_cbo)+c_cbo)*asym_cbo*TMath::Cos(omega_cbo*time - phi_cbo);

  return norm * TMath::Exp(-time/life) * (1 - asym*TMath::Cos(omega*time + phi)) * cbo;
}