//********************************************************************************
//
// 27-par fit with variable CBO function, variable vertical waist and variable fy; cross terms for 1.9 MHz
// with Expo + Constant CBO decoherence term
//
//********************************************************************************
Double_t fitFunc27parConstant(Double_t *x, Double_t *par){

    Double_t t = x[0];
    Double_t N = par[0];
    Double_t B = par[1];
    Double_t A = par[2];
    Double_t R = par[3];
    while (par[4] < 0)              par[4] += TMath::TwoPi();
    while (par[4] > TMath::TwoPi()) par[4] -= TMath::TwoPi();
    Double_t p = par[4];
    
    //CBO term
    Double_t ACBO = par[5];
    Double_t w0 = par[6];
    while (par[7] < 0)              par[7] += TMath::TwoPi();
    while (par[7] > TMath::TwoPi()) par[7] -= TMath::TwoPi();
    Double_t pCBO = par[7];
    Double_t tA = par[8];


    //two constant as fixed parameters
    Double_t t1sys = par[28];
    Double_t t2sys = par[29];


    Double_t varyingFreq = w0*t + A1sys * exp(-t / t1sys) + A2sys * exp(-t / t2sys);
    
    //amplitude and phase only, frequency determined by AV * fcbo
    Double_t AVW = par[9];
    Double_t factorfy = par[10];
    while (par[11] < 0)              par[11] += TMath::TwoPi();
    while (par[11] > TMath::TwoPi()) par[11] -= TMath::TwoPi();
    Double_t pVW = par[11];
    Double_t kLM = par[12];
    
    Double_t ampA = par[13];
    while (par[14] < 0)              par[14] += TMath::TwoPi();
    while (par[14] > TMath::TwoPi()) par[14] -= TMath::TwoPi();
    Double_t phaseA = par[14];
    Double_t ampP = par[15];
    while (par[16] < 0)              par[16] += TMath::TwoPi();
    while (par[16] > TMath::TwoPi()) par[16] -= TMath::TwoPi();
    Double_t phaseP = par[16];
    Double_t twiceA = par[17];
    while (par[18] < 0)              par[18] += TMath::TwoPi();
    while (par[18] > TMath::TwoPi()) par[18] -= TMath::TwoPi();
    Double_t twiceP = par[18];
    
    Double_t tauVW = par[19];
    
    //Add the fy term
    Double_t Afy = par[20];
    Double_t taufy = tauVW * 2.0;
    while (par[21] < 0)              par[21] += TMath::TwoPi();
    while (par[21] > TMath::TwoPi()) par[21] -= TMath::TwoPi();
    Double_t phasefy = par[21];
    Double_t tCBO = tA;
    
    Double_t w = getBlinded->paramToFreq(R);    
    
    Double_t C = par[27];
    Double_t deco = exp(-t/tCBO) + C;
    Double_t decoHighOrder = deco;
    
    Double_t NCBOTerm = 1.0 + (deco * ACBO * cos(varyingFreq - pCBO));
    NCBOTerm         += ( pow(decoHighOrder, 2) * twiceA * cos(2*varyingFreq - twiceP));
    Double_t ACBOTerm = 1.0 + ampA * decoHighOrder * cos(varyingFreq - phaseA);
    Double_t PCBOTerm = ampP * decoHighOrder * cos(varyingFreq - phaseP);

    Double_t fcbo = varyingFreq / (t*TMath::TwoPi()); //MHz
    Double_t fc = 1.0 / tCyclotron; //MHz
    Double_t fy = factorfy * fcbo * sqrt( (2*fc/(factorfy * fcbo)) - 1.0 );
    Double_t fvw = fc - (2.0 * fy);
    
    Double_t varyingVW = fvw * TMath::TwoPi() * t;
    Double_t VWTerm = 1.0 + exp(-t / tauVW) * AVW * cos(varyingVW - pVW);
    
    Double_t varY = fy * TMath::TwoPi() * t;
    Double_t varYTerm = 1.0 + exp(-t / taufy) * Afy * cos(varY - phasefy);
    
    Double_t LMTerm = 1.0 - kLM*lost_muon->GetBinContent((int)(t/0.1492)+1); 
    
    Double_t tCBOVW = par[22];
    Double_t ACBOplusVW = par[23];
    while (par[24] < 0)              par[24] += TMath::TwoPi();
    while (par[24] > TMath::TwoPi()) par[24] -= TMath::TwoPi();
    Double_t PCBOplusVW = par[24];
    Double_t ACBOminusVW = par[25];
    while (par[26] < 0)              par[26] += TMath::TwoPi();
    while (par[26] > TMath::TwoPi()) par[26] -= TMath::TwoPi();
    Double_t PCBOminusVW = par[26];
    Double_t CBOVWTerm = 1.0 + exp(-t / tCBOVW) * ( ACBOplusVW * cos( varyingFreq + varyingVW - PCBOplusVW ) + ACBOminusVW * cos( varyingFreq - varyingVW - PCBOminusVW ) );
    
    return N * exp(-t / B) * (1 + A*(ACBOTerm) * cos( w * t - p + PCBOTerm )) * NCBOTerm * VWTerm * LMTerm * varYTerm * CBOVWTerm;
}