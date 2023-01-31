import ROOT as R
import json
import sys
from array import array
from copy import deepcopy


def main(res,root_path,array_name):
    
    ar = array('d',res)
    f = R.TFile(root_path,'recreate')
    ad = R.TArrayD(len(res),ar)
    f.WriteObject(ad,array_name)
    f.Close()


def extract(in_file,in_tag,out_file,out_tag):

    f = R.TFile(in_file)
    func = f.Get(in_tag)
    N = func.GetNpar()
    vals = []
    for n in range(N):
        vals.append(func.GetParameter(n))
    ar = array('d',vals)
    fout = R.TFile(out_file,'recreate')
    ad = R.TArrayD(N,ar)
    fout.WriteObject(ad,out_tag)
    fout.Close()


if __name__ == '__main__':
    cbo_init = {
        'envelope' : [0.],
        'freq' : [6.86,6.2],
        'timing' : ['5','5','5'],
        
    }

    fout = R.TFile('./values/cbo_syst_init.root','recreate')
    res = {}
    for d in ['run2all_skip_calo18','run3a','run3b']:
        res[d] = {}
        for m in ['T1700','A']:
            res[d][m] = []
            fname = '/home/chencheng/Fitter_wa/run/fetch_note/run23_seed0_calos/condor_{0:}_{1:}.root'.format(d,m)
            f = R.TFile(fname)
            for n in range(25):
                func = f.Get('func_28paras_cbo_lost_vw_expansion_calos_default_{0:}_{1:}_{2:}'.format(m,d,n))
                vs = []
                for nvar in range(28):
                    v = func.GetParameter(nvar)
                    vs.append(v)
                for k,v in cbo_init.items():
                    vs_ = deepcopy(vs)
                    for p in v:
                        if isinstance(p, str):
                            vs_.append(vs[int(p)])
                        else:
                            vs_.append(p)
                    N = len(vs_)
                    ar = array('d',vs_)
                    ad = R.TArrayD(N,ar)
                    name = 'vals_{3:}_{0:}_{1:}_{2:}'.format(d,m,n,k)
                    fout.WriteObject(ad,name)
                    print(name,vs_)
    fout.Close()

                
                
    
    # for m in ['A','T1700']:
    #     intf = '/home/chencheng/Fitter_wa/run/values/run2_cbo/run2cbo_timing/condor_run2all_skip_calo18_{0:}/run2cbo_timing_{0:}_run2all_skip_calo18_10.root'.format(m)
    #     f = R.TFile(intf)
    #     ad = f.Get('run2cbo_timing_{0:}_run2all_skip_calo18_10'.format(m))
    #     pars = []
    #     for p in ad:
    #         pars.append(p)

    #     for k in cbo_init:
    #         appendix = cbo_init[k]
    #         res = deepcopy(pars)
    #         for v in appendix:
    #             if isinstance(v, str):
    #                 res.append(pars[int(v)]) 
    #             else:
    #                 res.append(v)
    #         print (k,res)

    #         ar = array('d',res)

    #         fout = R.TFile('./values/cbo_syst/vals_{0:}_{1:}_calo.root'.format(m,k),'recreate')
    #         N = len(res)
    #         ad = R.TArrayD(N,ar)
    #         fout.WriteObject(ad,'vals_{0:}_{1:}'.format(m,k))
    #         fout.Close()






    
