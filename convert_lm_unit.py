import ROOT as R

inF=R.TFile('/home/chencheng/data/Run23/Run2_3_LM_spectrum_and_integral.root')
outF=R.TFile('/home/chencheng/data/Run23/Run2_LM_intergral_ns.root','recreate')


h =inF.Get('run2all_LM_integral_from30')
N = h.GetNbinsX()


hnew = R.TH1D('run2all_integral','',N,0,149.2*N)


for n in range(N):
    v = h.GetBinContent(n+1)
    e = h.GetBinError(n+1)
    hnew.SetBinContent(n+1,v)
    hnew.SetBinError(n+1,e)

hnew.SetDirectory(outF)
outF.Write()
outF.Close()


