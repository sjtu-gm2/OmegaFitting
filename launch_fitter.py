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

def fit_random_type_calo(method,random_type,calo):

    tag = 'T1650_calo{2:}_{0:}_{1:}'.format(method,random_type,calo)

    wiggle_file = './data_random_type/wiggles_T.root'
    wiggle_name = 'wiggle_{0:}'.format(tag)

    lm_file = default_config['lm_file']
    lm_name = default_config['lm_name']

    func = '22paras_cbo_lost_vw_expansion_lite'

    output_dir = './output/{0:}'.format(tag)
    output_val = './values/{0:}'.format(tag)

    mode = 9

    initial_json = '/home/chencheng/Fitter_Wa/run/values_calofit/tmethod_calo_fit.json'
    initial_file = '/home/chencheng/Fitter_Wa/run/values_calofit/tmethod_calo_fit.root'
    initial_name = 'initial_values'

    trans.main(initial_json,initial_file,initial_name)

    max_try = 5

    os.system('mkdir -p %s'%(output_dir))
    os.system('mkdir -p %s'%(output_val))

    os.system('../build/MAIN {0:} {1:} {2:} {3:} {4:} {5:} {6:} {7:} {8:} {9:}'.format(
        wiggle_file,wiggle_name,lm_file,lm_name,initial_file,initial_name,output_dir,tag,mode,max_try))

method = 'shadow'
random_type = 'PerFill'
calo = 1

fit_random_type_calo(method, random_type,calo)