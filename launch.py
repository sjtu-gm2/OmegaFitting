import subprocess,os,json,sys,argparse,re,inspect,numpy,random,string

sys.path.append('%s/../OmegaFitting'%(os.environ['PWD']))
# import trans

####################################
# Template for job submission file #
####################################

condor = '''\
Universe   = vanilla
Executable = ./submit_caches/{1:}/run.sh
Log        = ./submit_caches/{1:}/submit.log
{0:}
'''

queue = '''\
Arguments  = {0:} {1:} {2:} {3:} $(Process)
Output     = {5:}/condor_{2:}.{3:}.$(Process).out
Error      = {5:}/condor_{2:}.{3:}.$(Process).error
Queue {4:}
'''

executable = '''\
#! /bin/bash
source  /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh
python ./submit_caches/%s/launch.py --run ./submit_caches/%s/${1} ${2} ${3} ${4} ${5}
echo -e "END - ./submit_caches/%s/${1} ${2} ${3} ${4} ${5}"
'''

#####################
# Utility functions #
#####################

def digest_scan_list(config,process=0):
    'Get scan list for the specified process no.'

    scan = config['scan']
    scan_per_queue = config['scan_per_queue']

    if isinstance(scan[0],str):
        scan_list = [k for k in scan]
    elif len(scan)==3:
        start = scan[0]
        end   = scan[1]
        step  = scan[2]
        scan_list = []
        index = 0
        while start+index*step<end:
            scan_list.append(start+index*step)
            index += 1
    else:
        raise ValueError()

    if (scan_per_queue<0):
        scan_per_queue = len(scan_list)

    tot_scan = len(scan_list)
    tot_process = int(numpy.ceil(tot_scan / scan_per_queue))
    this_process = int(process) % tot_process
    this_scan_list = []
    for subprocess in range(scan_per_queue):    
        this_index = this_process * scan_per_queue + subprocess
        if this_index < tot_scan:
            this_scan_list.append(scan_list[this_index])

    dummy = type("",(),{})

    dummy.scan_list = this_scan_list
    dummy.nq = tot_process
    dummy.q = this_process    
    dummy.nsubq = len(this_scan_list)
    dummy.subq_per_q = scan_per_queue
    return dummy

def regist_tag(cfg,entry):
    'Create cache directory for this job.'

    fname = os.path.basename(cfg).split('.')[0]

    while True:
        stamp = ''.join([random.choice(string.ascii_uppercase + string.digits) \
            for _ in range(3)])
        tag = '{0:}.{1:}.{2:}'.format(fname,entry,stamp)
        status = subprocess.getstatusoutput('touch ./submit_caches/jobs.log')
        status = subprocess.getstatusoutput('grep {0:} ./submit_caches/jobs.log'.format(tag))
        if status[0] == 1:
            os.system('echo "{0:}" >> ./submit_caches/jobs.log'.format(tag))
            os.system('mkdir -p ./submit_caches/{0:}'.format(tag))
            break
        elif status[0] == 2:
            print (status)
            raise ValueError
    return tag    

def parse_config(config,entry,dataset,job,scan,scan_id):

    kw = {
        'job' : job,
        'dataset' : dataset,
        'scan' : scan,
        'process' : scan_id,
        'key' :entry
    }

    kw_extra = {}

    if 'keys' in config:
        for k,str_key in config['keys'].items():
            if str_key[:4] == 'def ':
                space = {}
                fname = re.split('def\ |\(',str_key)[1]
                exec(str_key,space)
                var_dict = dict.fromkeys(inspect.getfullargspec(space[fname])[0])
                for key in var_dict:
                    var_dict[key] = kw[key]
                kw_extra[k] = space[fname](**var_dict)                
            else:
                kw_extra[k] = str(str_key).format(**kw)

    kw.update(kw_extra)
    def f(key):
        return str(config[key]).format(**kw).replace(';','\\;')
    
    parsed_config = {}
    for key in ['wiggle_file', 'wiggle_name', 'lm_file', 'lm_name', 'initial_file', 'initial_name', 'start_bin', 'end_bin', 'output_dir', 'tag']:
        parsed_config[key] = f(key)
    
    for key in ['mode', 'max_try']:
        parsed_config[key] = config[key]

    return parsed_config

##################
# Main functions #
##################

def submit(args):
    '''
    Create run scripts and job submission files. \\
    Called by --submit.
    '''

    config = {
        'logDir' : './logs',
        'scan_per_queue' : 1,
    }
    with open(args[0],'r') as fp:
        config.update(json.load(fp)[args[1]])
    logDir = config['logDir'].format(key=args[1])

    os.system('mkdir -p %s'%(logDir))
    os.system('mkdir -p submit_caches')
    os.system('touch ./submit_caches/jobs.log')

    tag = regist_tag(args[0],args[1])

    N = digest_scan_list(config).nq

    cfg_base_name = os.path.basename(args[0])

    Queues = []
    for job in config['job']:
        for dataset in config['dataset']:
            q = queue.format(cfg_base_name,args[1],job,dataset,N,logDir)
            Queues.append(q)

    condor_str = condor.format('\n'.join(Queues),tag)
    condor_file = './submit_caches/{0:}/submit.condor'.format(tag)
    with open(condor_file,'w') as fp:
        fp.write(condor_str)
    run_file = './submit_caches/{0:}/run.sh'.format(tag)
    steering_file = './submit_caches/{0:}/launch.py'.format(tag)
    os.system('cp ./launch.py {0:}'.format(steering_file))
    os.system('cp {0:} ./submit_caches/{1:}/'.format(args[0],tag))

    with open(run_file,'w') as fp:
        fp.write(executable%(tag,tag,tag))
    os.system('chmod 777 {0:}'.format(run_file))
    os.system('condor_submit {0:}'.format(condor_file))

def run_queue(args):
    '''
    Invoke "run" function to run scan jobs in a queue (process).
    '''

    config = {
        'max_try' : 3,
        'start_bin' : 202,
        'end_bin' : 4357,
        'use_list' : 0,        
        'scan_per_queue' : 1,
        'auto_combine' : False,
        'refit_by_scan' : True
    }
    with open(args[0],'r') as fp:
        config.update(json.load(fp)[args[1]])

    entry = args[1]    
    job = args[2]
    dataset = args[3]
    process = args[4]

    dummy = digest_scan_list(config,process)
    out_files_match = [] # output files to be hadded
    out_files = {}       # output files to become initial value
    out_names = {}

    for subq, scan in enumerate(dummy.scan_list):
        print('\n','-'*30)

        tag = '{0:}_{1:}'.format(job,dataset)
        if config['use_list']:
            if (not tag in config) or (not scan in config[tag]):
                print ('skip {0:} {1:}'.format(tag,scan))
                continue
        print('Processing {0:}.{1:}.{2:}  q: {3:}/{4:} sub-q: {5:}/{6:} scan point: {7:}'.format(entry,dataset,job,dummy.q+1,dummy.nq,subq+1,dummy.nsubq,scan))

        scan_id = subq+dummy.subq_per_q*dummy.q
        new_config = parse_config(config,entry,dataset,job,scan,scan_id)

        # Refit by pluging in the previous output to the initial value
        now_fit = 0
        max_fit = 1
        if config['refit_by_scan']:
            max_fit = config['max_try']
            config['max_try'] = 1

        while now_fit < max_fit and subq >= now_fit:
            if (now_fit > 0):
                new_config['initial_file'] = out_files[subq - now_fit]
                new_config['initial_name'] = out_names[subq - now_fit]

            out_dir,out_name,full_tag = run(new_config)
            out_file_match = '{0:}/result_*{1:}.root'.format(out_dir,out_name)
            out_files[subq] = ('{0:}/result_{1:}.root'.format(out_dir,full_tag))
            out_names[subq] = ("func_{0:}".format(full_tag))

            now_fit += 1
        
        out_files_match.append(out_file_match)

    if config['auto_combine']:
        tag = '{0:}_{1:}_{2:}'.format(job,dataset,process)
        in_files = ' '.join(out_files_match)
        os.system('hadd -f {0:}/combined_{1:}.root {2:}'.format(out_dir,tag,in_files))

    # print(out_files)
    # print(out_names)


def run(config):
    '''
    Run a fitting for the specified scan id.
    '''

    wiggle_file  = config['wiggle_file']
    wiggle_name  = config['wiggle_name']
    lm_file      = config['lm_file']
    lm_name      = config['lm_name']
    initial_file = config['initial_file']
    initial_name = config['initial_name']

    start_bin    = config['start_bin']
    end_bin      = config['end_bin']

    output_dir   = config['output_dir']
    tag_out      = config['tag']

    mode         = config['mode']
    max_try      = config['max_try']


    os.system('mkdir -p {0:}'.format(output_dir))    
    

    fix = "None"
    if 'fix' in config:
        fix = ''
        for key,value in config['fix'].items():
            fix += '{0:} {1:} '.format(key,value)

    range_par = "None"
    if 'range' in config:
        range_par = ''
        for key,ranges in config['range'].items():
            range_par += '{0:} {1:} {2:} '.format(key,ranges[0],ranges[1])

    cmd = '../build/MAIN {0:} {1:} {2:} {3:} {4:} {5:} {6:} {7:} {8:} {9:} {10:} {11:} --fix {12:} --range {13:}'.format(
        wiggle_file,wiggle_name,lm_file,lm_name,initial_file,initial_name,output_dir,tag_out,mode,max_try,start_bin,end_bin,fix,range_par)

    print (cmd)
    # os.system(cmd)
    output = os.popen(cmd).read()
    print(output)

    pattern = r"function\s+:\s+func_(\w+)"
    matches = re.findall(pattern, output)
    if matches:
        full_tag = matches[-1]
    else:
        print("Error! Cannot find function name in output!")

    return output_dir,tag_out,full_tag
    
    # fetch values 
    # value_dir    = f('value_dir')
    # value_of_func = config['value_of_func']

    # os.system('mkdir -p {0:}'.format(value_dir))
    # res_name = '{0:}/result_{1:}_{2:}.root'.format(output_dir,value_of_func,tag_out)
    # res_func_name = 'func_{0:}_{1:}'.format(value_of_func,tag_out)
    # output_file = '{0:}/{1:}.root'.format(value_dir,tag_out)

    # trans.extract(res_name,res_func_name,output_file,tag_out)

def fetch_clean(args,mode='fetch'):
    config = {
        'auto_combine' : 0
    }

    with open(args[0],'r') as fp:
        config.update(json.load(fp)[args[1]])

    fetch_dir = './fetch/'
    if len(args)>2:
        fetch_dir = args[2]
    if config['auto_combine']==1:
        filter_str = 'combined_*'
    else:
        filter_str = '*'

    if len(args)>3:
        filter_str = args[3]

    if mode == 'fetch':
        os.system('mkdir -p {0:}'.format(fetch_dir))

    for job in config['job']:
        for dataset in config['dataset']:
            input_dir = config['output_dir'].format(job=job,dataset=dataset,key=args[1])
            basename = os.path.basename(input_dir.rstrip('/'))
            if mode == 'fetch':
                cmd = 'hadd -f {0:}/{4:}_{1:}.root {2:}/{3:}.root'.format(fetch_dir,basename,input_dir,filter_str,args[1])
                os.system(cmd)
            elif mode == 'clean':
                cmd = 'rm -rf {0:}'.format(input_dir)
                go = input('Delting: {0:}\n----\nPress y or [return] to delete, others to abort!'.format(input_dir))
                if len(go)==0 or go[0] == 'y':
                    os.system(cmd)
                    print('\ndeleted {0:}\n'.format(input_dir),'-'*20)
                else:                    
                    raise ValueError('\ndelete aborted!')
            else:
                raise ValueError(mode)            

def fetch(args):
    fetch_clean(args,'fetch')

def clean(args):
    fetch_clean(args,'clean')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--submit',nargs='+')
    parser.add_argument('--run',nargs='+')
    parser.add_argument('--fetch',nargs='+')
    parser.add_argument('--clean',nargs='+')    

    args = parser.parse_args()

    keys = []    
    for key in args.__dict__.keys():
        if args.__dict__[key]!=None:
            keys.append(key)

    if len(keys)!=1:
        parser.error('Should use 1 and only 1 job!')

    job = keys[0]

    if job == 'submit':
        submit(args.submit)
    elif job == 'run':
        run_queue(args.run)
    elif job == 'fetch':
        fetch(args.fetch)
    elif job == 'clean':
        clean(args.clean)
    else:
        raise ValueError('use "python ./launch.py --submit/run/fetch/clean ... "')