import ROOT as R
import json
import sys
from array import array


def main(json_path,root_path,array_name):
    # root_path = json_path.rstrip('json') + 'root'

    with open(json_path) as fp:
        res = json.load(fp)['0']

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
    extract('./output/T1650_calo1_shadow_PerFill/result_13paras_cbo_vo_T1650_calo1_shadow_PerFill.root','func_13paras_cbo_vo_T1650_calo1_shadow_PerFill', './values/start_time_T1650_calo1_shadow_PerFill/time_0.root','T1650_calo1_shadow_PerFill_time0_13paras_cbo_vo')
    # main(sys.argv[1], sys.argv[2])
    # for n in range(1000,3000,100):
    #     slice_name = '%s_%s'%(n,n+100)
    #     json_name = '/Users/cheng/workspace/TestCMake_fitter/run/output_slices_values/Vals_E%s_10paras_cbo_lost.json'%(slice_name)
    #     root_name = './values/empirical_normal_slice%s.root'%(slice_name)
    #     array_name = 'empirical_normal_slice%s'%(slice_name)
    #     main(json_name,root_name,array_name)
    # main('./jsons/run4_baseline_T.json', './values/empirical_normal_A_nominal_22paras.root' ,'empirical_normal_A_nominal_22paras')
