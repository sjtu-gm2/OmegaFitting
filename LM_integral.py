import ROOT as R
from math import exp
def make_lm_int(hist_cso,outf,name):
  integral_hist_from30 = R.TH1D('%s_from30'%(name),'',4500,0,671.4e3)  
  integral_hist = R.TH1D('%s'%(name),'',4500,0,671.4e3)  
  N = hist_cso.GetNbinsX()

  currentsum = 0.
  currentsum_from30 = 0.

  binw=hist_cso.GetBinWidth(1)  
  for i in range(N):
    t = hist_cso.GetBinCenter(i+1)
    if t>30.:
      currentsum_from30 += hist_cso.GetBinContent(i+1)*exp((i*binw+binw/2)/64.4)
      integral_hist_from30.SetBinContent(i+1,currentsum_from30)
    currentsum += hist_cso.GetBinContent(i+1)*exp((i*binw+binw/2)/64.4)
    integral_hist.SetBinContent(i+1,currentsum)

  integral_hist.SetDirectory(outf)
  integral_hist_from30.SetDirectory(outf)

  return [integral_hist,integral_hist_from30]
  




def run2_individual():
  tfile = R.TFile('/home/chencheng/data/LM_Run23_SJTU/SJTU_LM_run2.root')
  outf = R.TFile('/home/chencheng/data/LM_Run23_SJTU/LM_Run2_integral.root','recreate')

  hss = []
  for p in ['B','C','D','E','F','G','H','all']:
    hist_cso = tfile.Get('run2%s_LM_spec'%(p))
    if p == 'F':
      hists = make_lm_int(hist_cso,outf,'LM_integral_run2F_skip_calo18')
    else:
      hists = make_lm_int(hist_cso,outf,'LM_integral_run2%s'%(p))

    hss.append(hists)
  outf.Write()
  outf.Close()

def run3_individual():
  tfile = R.TFile('/home/chencheng/data/LM_Run23_SJTU/SJTU_LM_run3.root')
  outf = R.TFile('/home/chencheng/data/LM_Run23_SJTU/LM_Run3_integral.root','recreate')

  hss = []
  for p in ['B','C','D','E','F','G','I','J','K','L','M','N','O','all']:
    hist_cso = tfile.Get('lm3%s_spec'%(p))
    hists = make_lm_int(hist_cso,outf,'LM_integral_run3%s'%(p))
    hss.append(hists)
  outf.Write()
  outf.Close()
  print ('generated %s'%(outf))


def main():
  run2_individual()

if __name__ == '__main__':
  main()