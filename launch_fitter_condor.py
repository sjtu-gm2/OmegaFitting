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
    
    output_dir = './output/condor_{3:}_T1650_calo{0:}_{1:}_{2:}'.format(calo,method,random_type,job)
    
    output_val_condor = './values/condor_{3:}_T1650_calo{0:}_{1:}_{2:}'.format(calo,method,random_type,job)
    os.system('mkdir -p %s'%(output_dir))    
    os.system('mkdir -p %s'%(output_val_condor))

    mode = 10

    func_name = '13paras_cbo_vo'

    output_file = '{0:}/calo_{1:}.root'.format(output_val_condor,calo)
    output_name = 'T1650_calo{0:}_{1:}_{2:}_time{3:}_{4:}'.format(calo,method,random_type,time_offset,func_name)

    # initial_file = './values/condor_start_time_T1650_calo1_shadow_PerFill/time_0.root'
    # initial_name = 'T1650_calo1_shadow_PerFill_time0_13paras_cbo_vo'
    initial_file = './values/condor_{3:}_T1650_calo{0:}_{1:}_{2:}/calo_{0:}.root'.format(21,method,random_type,job)
    initial_name = 'T1650_calo{0:}_{1:}_{2:}_time{3:}_{4:}'.format(21,method,random_type,time_offset,func_name)

    

    tag_out = 'T1650_calo{0:}_{1:}_{2:}_time{3:}'.format(calo,method,random_type,time_offset)

    max_try = 5

    os.system('../build/MAIN {0:} {1:} {2:} {3:} {4:} {5:} {6:} {7:} {8:} {9:} {10:}'.format(
        wiggle_file,wiggle_name,lm_file,lm_name,initial_file,initial_name,output_dir,tag_out,mode,max_try,time_offset))

    res_name = '{0:}/result_{1:}_{2:}.root'.format(output_dir,func_name,tag_out)
    res_func_name = 'func_{0:}_{1:}'.format(func_name,tag_out)

    trans.extract(res_name,res_func_name,output_file,output_name)


def fit_random_type_seed_scan(method,random_type,seed):

    wiggle_file = './data_random_type/wiggles_seed_scan_fix.root'
    # wiggle_file = './data_random_type/wiggles_T.root'
    wiggle_name = 'wiggle_T1650_{1:}_{2:}_seed{0:}'.format(seed,method,random_type)
    # wiggle_name = 'wiggle_T1650_{0:}_{1:}'.format(method,random_type)


    lm_file = default_config['lm_file']
    lm_name = default_config['lm_name']
    
    output_dir = './output/condor_{0:}_T1650_{1:}_{2:}'.format(job,method,random_type)
    
    output_val_condor = './values/condor_{0:}_T1650_{1:}_{2:}'.format(job,method,random_type)

    os.system('mkdir -p %s'%(output_dir))    
    os.system('mkdir -p %s'%(output_val_condor))

    mode = 7

    func_name = '22paras_cbo_lost_vw_expansion_lite'

    output_file = '{0:}/seed_{1:}.root'.format(output_val_condor,seed)
    output_name = 'T1650_{1:}_{2:}_{3:}_{0:}'.format(seed,method,random_type,func_name)

    # initial_file = './values/condor_start_time_T1650_calo1_shadow_PerFill/time_0.root'
    # initial_name = 'T1650_calo1_shadow_PerFill_time0_13paras_cbo_vo'
    initial_file = '/home/chencheng/Fitter_Wa/run/jsons/nominal_fit_T/Vals_shadow_normal_22paras_cbo_lost_vw_expansion_lite.root'
    initial_name = 'Vals_shadow_normal_22paras_cbo_lost_vw_expansion_lite'

    

    tag_out = 'T1650_{1:}_{2:}_seed{0:}'.format(seed,method,random_type)

    max_try = 5
    time_offset = 0
    os.system('../build/MAIN {0:} {1:} {2:} {3:} {4:} {5:} {6:} {7:} {8:} {9:} {10:}'.format(
        wiggle_file,wiggle_name,lm_file,lm_name,initial_file,initial_name,output_dir,tag_out,mode,max_try,time_offset))

    res_name = '{0:}/result_{1:}_{2:}.root'.format(output_dir,func_name,tag_out)
    res_func_name = 'func_{0:}_{1:}'.format(func_name,tag_out)

    trans.extract(res_name,res_func_name,output_file,output_name)

def fit_random_type_pileup_scan(method,random_type,scan_point):

    wiggle_file = './data_pileup_scan/wiggles_pileup_scan_fix.root'
    # wiggle_file = './data_random_type/wiggles_T.root'
    wiggle_name = 'wiggle_T1650_{1:}_{2:}_pu{0:.1f}'.format(scan_point*0.1,method,random_type)

    lm_file = default_config['lm_file']
    lm_name = default_config['lm_name']

    output_dir = './output/condor_{0:}_T1650_{1:}_{2:}'.format(job,method,random_type)

    output_val_condor = './values/condor_{0:}_T1650_{1:}_{2:}'.format(job,method,random_type)

    os.system('mkdir -p %s'%(output_dir))    
    os.system('mkdir -p %s'%(output_val_condor))

    mode = 7

    func_name = '22paras_cbo_lost_vw_expansion_lite'

    output_file = '{0:}/pu_{1:.1f}.root'.format(output_val_condor,scan_point*0.1)
    output_name = 'T1650_{1:}_{2:}_{3:}_{0:.1f}'.format(scan_point*0.1,method,random_type,func_name)

    # initial_file = './values/condor_start_time_T1650_calo1_shadow_PerFill/time_0.root'
    # initial_name = 'T1650_calo1_shadow_PerFill_time0_13paras_cbo_vo'
    initial_file = '/home/chencheng/Fitter_Wa/run/jsons/nominal_fit_T/Vals_shadow_normal_22paras_cbo_lost_vw_expansion_lite.root'
    initial_name = 'Vals_shadow_normal_22paras_cbo_lost_vw_expansion_lite'

    

    tag_out = 'T1650_{1:}_{2:}_pu{0:.1f}'.format(scan_point*0.1,method,random_type)

    max_try = 5
    time_offset = 0
    os.system('../build/MAIN {0:} {1:} {2:} {3:} {4:} {5:} {6:} {7:} {8:} {9:} {10:}'.format(
        wiggle_file,wiggle_name,lm_file,lm_name,initial_file,initial_name,output_dir,tag_out,mode,max_try,time_offset))

    res_name = '{0:}/result_{1:}_{2:}.root'.format(output_dir,func_name,tag_out)
    res_func_name = 'func_{0:}_{1:}'.format(func_name,tag_out)

    trans.extract(res_name,res_func_name,output_file,output_name)


if __name__ == '__main__':    
    method = sys.argv[1]
    random_type = sys.argv[2]
    job = sys.argv[4]
    if job == 'seed_scan':    
        seed = int(sys.argv[3])%200    
        fit_random_type_seed_scan(method, random_type,seed)
    elif job == 'pileup_scan':
        scan_point = int(sys.argv[3])%20
        fit_random_type_pileup_scan(method,random_type,scan_point)