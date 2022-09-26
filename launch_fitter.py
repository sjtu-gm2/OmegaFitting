import subprocess
import os
import json
import sys

sys.path.append('/home/chencheng/Fitter_Wa/run')
import trans

default_config = {
    'lm_file' : '~/data/Run4U/nominal_emp/gm2pro_daq_offline_dqc_run4U_5307A_makeBULostMuons_gridRun.root',
    'lm_name' : 'topDir/Iter0/LostMuons/BaseCuts/Triples/Losses/triple_losses_spectra_integral',
}

def fit_random_type_calo(method,random_type,calo,time_offset):

    wiggle_file = './data_random_type/wiggles_T.root'
    wiggle_name = 'wiggle_T1650_calo{0:}_{1:}_{2:}'.format(calo,method,random_type)


    lm_file = default_config['lm_file']
    lm_name = default_config['lm_name']
    
    output_dir = './output/start_time_T1650_calo{0:}_{1:}_{2:}'.format(calo,method,random_type)
    output_val = './values/start_time_T1650_calo{0:}_{1:}_{2:}'.format(calo,method,random_type)
    os.system('mkdir -p %s'%(output_dir))
    os.system('mkdir -p %s'%(output_val))

    mode = 10

    func_name = '13paras_cbo_vo'
    initial_file = '{0:}/time_{1:}.root'.format(output_val,time_offset-5)
    initial_name = 'T1650_calo{0:}_{1:}_{2:}_time{3:}_{4:}'.format(calo,method,random_type,time_offset-5,func_name)

    output_file = '{0:}/time_{1:}.root'.format(output_val,time_offset)
    output_name = 'T1650_calo{0:}_{1:}_{2:}_time{3:}_{4:}'.format(calo,method,random_type,time_offset,func_name)

    tag_out = 'T1650_calo{0:}_{1:}_{2:}_time{3:}'.format(calo,method,random_type,time_offset)

    max_try = 5

    os.system('../build/MAIN {0:} {1:} {2:} {3:} {4:} {5:} {6:} {7:} {8:} {9:} {10:}'.format(
        wiggle_file,wiggle_name,lm_file,lm_name,initial_file,initial_name,output_dir,tag_out,mode,max_try,time_offset))
    
    res_name = '{0:}/result_{1:}_{2:}.root'.format(output_dir,func_name,tag_out)
    res_func_name = 'func_{0:}_{1:}'.format(func_name,tag_out)


    trans.extract(res_name,res_func_name,output_file,output_name)
    
    # os.system('hadd -f {0:}/result_{1:}.root {0:}/result_*{1:}.root'.format(output_dir,tag_out))


if __name__ == '__main__':
    n = int(sys.argv[1])+1
    method = 'shadow'
    random_type = 'PerFill'
    calo = 1    
    fit_random_type_calo(method, random_type,calo,n*5)
