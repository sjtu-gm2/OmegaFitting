import subprocess,os,json,sys,argparse,re,inspect,numpy,random,string,copy

default_config = {
    'max_try' : 3,
    'start_bin' : 202,
    'end_bin' : 4357,
    'use_list' : 0,        
    'scan_per_queue' : 1,
    'auto_combine' : False,
    'refit_by_scan' : True
}

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

def digest_scan_list(config, process=0):
    'Get scan list for the specified process no.'

    scan = config['scan']
    scan_per_queue = config['scan_per_queue']

    if isinstance(scan[0], str):
        scan_list = copy.deepcopy(scan)
    elif isinstance(scan[0], (int, float)) and len(scan) == 3:
        scan_list = numpy.arange(scan[0], scan[1], scan[2])
    else:
        raise ValueError("Parameter scan should be list of string or 3 numbers representing start stop step of scanning numbers.")

    if (scan_per_queue < 0):
        scan_per_queue = len(scan_list)

    tot_scan = len(scan_list)
    tot_process = int(numpy.ceil(tot_scan / scan_per_queue))
    this_process = int(process) % tot_process
    this_scan_list = []
    for subprocess in range(scan_per_queue):    
        this_index = this_process * scan_per_queue + subprocess
        if this_index < tot_scan:
            this_scan_list.append(scan_list[this_index])

    scan_info = {"scan_list": this_scan_list,
                 "nq": tot_process,
                 "q": this_process,
                 "nsubq": len(this_scan_list),
                 "subq_per_q": scan_per_queue}

    return scan_info


def regist_tag(cfg, entry):
    'Create cache directory for this job.'

    fname = os.path.basename(cfg).split('.')[0]

    while True:
        stamp = ''.join([random.choice(string.ascii_uppercase + string.digits) \
            for _ in range(3)])
        tag = '{0:}.{1:}.{2:}'.format(fname, entry, stamp)
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
    '''
    Turn the origin config into formatted one.
    '''

    kw = {
        'job' : job,
        'dataset' : dataset,
        'scan' : scan,
        'process' : scan_id,
        'key' : entry
    }

    # Process "keys" in config
    kw_extra = {}

    if 'keys' in config:
        for key, str_value in config['keys'].items():
            if str_value[:4] == 'def ':
                # Turn the string in "keys" into function.
                space = {}
                funcname = re.split('def\ |\(',str_value)[1]
                exec(str_value, space)
                myfunc = space[funcname]

                # Perform the function and get extra keywords.
                args_list = inspect.getfullargspec(myfunc)[0]
                arg_dict = {}
                for arg_name in args_list:
                    arg_dict[arg_name] = kw[arg_name]
                kw_extra[key] = myfunc(**arg_dict)
            else:
                kw_extra[key] = str(str_value).format(**kw)

    kw.update(kw_extra)
    def get_parsed(arg):
        return str(config[arg]).format(**kw).replace(';','\\;')
    
    parsed_config = copy.deepcopy(config)
    for arg in ['wiggle_file', 'wiggle_name', 'lm_file', 'lm_name', 'initial_file', 'initial_name', 'start_bin', 'end_bin', 'output_dir', 'tag']:
        parsed_config[arg] = get_parsed(arg)

    # Process "fix" and "range" in config
    fix_pars = "None"
    if 'fix' in config:
        fix_pars = ''
        for key, value in config['fix'].items():
            fix_pars += '{0:} {1:} '.format(key, value)

    range_pars = "None"
    if 'range' in config:
        range_pars = ''
        for key, ranges in config['range'].items():
            range_pars += '{0:} {1:} {2:} '.format(key, ranges[0], ranges[1])

    parsed_config.update({"fix_pars": fix_pars, "range_pars": range_pars})

    return parsed_config


def execute_and_parse_output(cmd):
    '''
    Run the c++ fitting program and parse output to get useful information (tag name and fit valid)
    '''
    pattern_tag = r"function\s+:\s+func_(\w+)"
    matches_tag = []
    pattern_valid = r"Valid\s+=\s+(\d+)"
    matches_valid = []

    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    for line in process.stdout:
        print(line.strip())
        match_tag = re.search(pattern_tag, line)
        match_valid = re.search(pattern_valid, line)
        if match_tag:
            matches_tag.append(match_tag.group(1))
        if match_valid:
            matches_valid.append(match_valid.group(1))

    process.wait()

    if matches_tag:
        full_tag = matches_tag[-1]
    else:
        raise RuntimeError("Cannot find function name in output.")

    if matches_valid:
        fit_valid = matches_valid[-1]
    else:
        raise RuntimeError("Cannot find fit valid in output.")

    return full_tag, fit_valid

##################
# Main functions #
##################

def run(config):
    '''
    Run a fitting for the specified scan id. \\
    The config input should be the parsed one.
    '''
    output_dir = config['output_dir']
    tag_out = config['tag']

    os.system('mkdir -p {0:}'.format(output_dir))    

    cmd = '../build/MAIN {wiggle_file} {wiggle_name:} {lm_file:} {lm_name:} {initial_file:} {initial_name:} {output_dir:} {tag:} {mode:} {max_try:} {start_bin:} {end_bin:} --fix {fix_pars:}--range {range_pars:}'.format(**config)

    print(cmd)
    full_tag, fit_valid = execute_and_parse_output(cmd)

    return output_dir, tag_out, full_tag, fit_valid


def run_queue(args):
    '''
    Invoke "run" function to run scan jobs in a queue (scan_id). \\
    Called by --run.
    '''
    file_name, entry, job, dataset, scan_id = args

    config = copy.deepcopy(default_config)
    with open(file_name, 'r') as fp:
        config.update(json.load(fp)[args[1]])

    scan_info = digest_scan_list(config, scan_id)
    out_files_match = [] # match pattern of output files to be hadded
    out_files = {}       # output files to become initial value
    out_names = {}       # output names to become initial value

    max_fit = 1
    if config['refit_by_scan']:
        max_fit = config['max_try']
        config['max_try'] = 0

    for subq, scan in enumerate(scan_info["scan_list"]):
        print('\n','-'*30)

        tag = '{0:}_{1:}'.format(job, dataset)
        if config['use_list']:
            if (not tag in config) or (not scan in config[tag]):
                print ('skip {0:} {1:}'.format(tag, scan))
                continue

        scan_id = subq + scan_info["subq_per_q"] * scan_info["q"]
        new_config = parse_config(config, entry, dataset, job, scan, scan_id)

        # Refit by pluging in the previous output to the initial value
        now_fit = 0

        while now_fit < max_fit and subq >= now_fit:
            if (now_fit > 0):
                new_config['initial_file'] = out_files[subq - now_fit]
                new_config['initial_name'] = out_names[subq - now_fit]

            print('Processing {0:}.{1:}.{2:}  queue:({3:}/{4:}) sub-queue:({5:}/{6:}) scan point:({7:}) time of try:({8:})'.format(entry, dataset, job, scan_info["q"]+1, scan_info["nq"], subq+1, scan_info["nsubq"], scan, now_fit+1))

            out_dir, out_name, full_tag, fit_valid = run(new_config)
            out_file_match = '{0:}/result_*{1:}.root'.format(out_dir, out_name)
            out_files[subq] = ('{0:}/result_{1:}.root'.format(out_dir, full_tag))
            out_names[subq] = ("func_{0:}".format(full_tag))

            if fit_valid == "1":
                print("Fitting succeed.")
                break
            now_fit += 1
        
        out_files_match.append(out_file_match)

    if config['auto_combine']:
        tag = '{0:}_{1:}_{2:}'.format(job, dataset, scan_id)
        in_files = ' '.join(out_files_match)
        os.system('hadd -f {0:}/combined_{1:}.root {2:}'.format(out_dir, tag, in_files))


def submit(args):
    '''
    Create run scripts and job submission files. \\
    Called by --submit.
    '''
    filename, entry = args

    config = {
        'logDir' : './logs',
        'scan_per_queue' : 1,
    }
    with open(filename,'r') as fp:
        config.update(json.load(fp)[entry])
    logDir = config['logDir'].format(key=entry)

    os.system('mkdir -p %s'%(logDir))
    os.system('mkdir -p submit_caches')
    os.system('touch ./submit_caches/jobs.log')

    tag = regist_tag(filename, entry)

    N = digest_scan_list(config)["nq"]

    cfg_base_name = os.path.basename(filename)

    Queues = []
    for job in config['job']:
        for dataset in config['dataset']:
            q = queue.format(cfg_base_name, entry, job, dataset, N, logDir)
            Queues.append(q)

    condor_str = condor.format('\n'.join(Queues), tag)
    condor_file = './submit_caches/{0:}/submit.condor'.format(tag)
    with open(condor_file, 'w') as fp:
        fp.write(condor_str)
    run_file = './submit_caches/{0:}/run.sh'.format(tag)
    steering_file = './submit_caches/{0:}/launch.py'.format(tag)
    os.system('cp ./launch.py {0:}'.format(steering_file))
    os.system('cp {0:} ./submit_caches/{1:}/'.format(args[0], tag))

    with open(run_file, 'w') as fp:
        fp.write(executable%(tag, tag, tag))
    os.system('chmod 777 {0:}'.format(run_file))
    os.system('condor_submit {0:}'.format(condor_file))


def fetch_clean(args, mode='fetch'):
    config = {
        'auto_combine' : 0
    }

    with open(args[0], 'r') as fp:
        config.update(json.load(fp)[args[1]])

    fetch_dir = './fetch/'
    if len(args) > 2:
        fetch_dir = args[2]
    if config['auto_combine'] == 1:
        filter_str = 'combined_*'
    else:
        filter_str = '*'

    if len(args) > 3:
        filter_str = args[3]

    if mode == 'fetch':
        os.system('mkdir -p {0:}'.format(fetch_dir))

    for job in config['job']:
        for dataset in config['dataset']:
            input_dir = config['output_dir'].format(job=job, dataset=dataset, key=args[1])
            basename = os.path.basename(input_dir.rstrip('/'))
            if mode == 'fetch':
                cmd = 'hadd -f {0:}/{4:}_{1:}.root {2:}/{3:}.root'.format(fetch_dir, basename, input_dir, filter_str, args[1])
                os.system(cmd)
            elif mode == 'clean':
                cmd = 'rm -rf {0:}'.format(input_dir)
                go = input('Delting: {0:}\n----\nPress y or [return] to delete, others to abort!'.format(input_dir))
                if len(go) == 0 or go[0] == 'y':
                    os.system(cmd)
                    print('\ndeleted {0:}\n'.format(input_dir), '-'*20)
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
    parser.add_argument('--submit', nargs='+')
    parser.add_argument('--run', nargs='+')
    parser.add_argument('--fetch', nargs='+')
    parser.add_argument('--clean', nargs='+')    

    args = parser.parse_args()

    keys = []    
    for key in args.__dict__.keys():
        if args.__dict__[key] != None:
            keys.append(key)

    if len(keys) != 1:
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