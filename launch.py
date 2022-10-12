import subprocess,os,json,sys,argparse
sys.path.append('%s/../OmegaFitting'%(os.environ['PWD']))
import trans

condor = '''\
Universe   = vanilla
Executable = run_temp.sh
Log        = submit.log
{0:}
'''

queue = '''\
Arguments  = {0:} {1:} {2:} {3:} $(Process)
Output     = logs/condor_{2:}.{3:}.$(Process).out
Error      = logs/condor_{2:}.{3:}.$(Process).error
Queue {4:}
'''

executable = '''\
#! /bin/bash
source  /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh
python ./launch.py --run ${1} ${2} ${3} ${4} ${5}
echo -e "END - ${1} ${2} ${3} ${4} ${5}"
'''

def submit(args):
    with open(args[0],'r') as fp:
        config = json.load(fp)[args[1]]
    scan = config['scan']
    N = int((scan[1]-scan[0])/scan[2])
    Queues = []
    for job in config['job']:
        for dataset in config['dataset']:
            q = queue.format(args[0],args[1],job,dataset,N)
            Queues.append(q)
    condor_str = condor.format('\n'.join(Queues))
    condor_file = 'submit_{0:}.condor'.format(args[1])
    with open(condor_file,'w') as fp:
        fp.write(condor_str)
    run_file = 'run_{0:}.sh'.format(args[1])
    with open(run_file,'w') as fp:
        fp.write(executable)
    os.system('chmod 777 {0:}'.format(run_file))
    os.system('condor_submit {0:}'.format(condor_file))

def run(args):
    config = {
        'max_try' : 3,
        'time_offset' : 0
    }

    with open(args[0],'r') as fp:
        config.update(json.load(fp)[args[1]])

    job = args[2]
    dataset = args[3]
    N = int((config['scan'][1]-config['scan'][0])/config['scan'][2])
    
    process = int(args[4]) % N
    step = config['scan'][2]
    scan = process * step

    kw = {
        'job' : job,
        'dataset' : dataset,
        'scan' : scan,
        'process' : process,        
    }
    
    f = lambda x:str(config[x]).format(**kw)

    wiggle_file = f('wiggle_file')
    wiggle_name = f('wiggle_name')
    lm_file = f('lm_file')
    lm_name = f('lm_name')
    initial_file = f('initial_file')
    initial_name = f('initial_name')
    output_dir = f('output_dir')
    value_dir = f('value_dir')
    tag_out = f('tag')
    mode = config['mode']
    max_try = config['max_try']
    time_offset = f('time_offset')
    value_of_func = config['value_of_func']

    os.system('mkdir -p {0:}'.format(output_dir))
    os.system('mkdir -p {0:}'.format(value_dir))

    fix = "None"
    if 'fix' in config:
        fix = ''
        for key,value in config['fix'].iteritems():
            fix += '{0:} {1:}'.format(key,value)
    os.system('../build/MAIN {0:} {1:} {2:} {3:} {4:} {5:} {6:} {7:} {8:} {9:} {10:} --fix {11:}'.format(
        wiggle_file,wiggle_name,lm_file,lm_name,initial_file,initial_name,output_dir,tag_out,mode,max_try,time_offset,fix))

    res_name = '{0:}/result_{1:}_{2:}.root'.format(output_dir,value_of_func,tag_out)
    res_func_name = 'func_{0:}_{1:}'.format(value_of_func,tag_out)
    output_file = '{0:}/{1:}.root'.format(value_dir,tag_out)

    trans.extract(res_name,res_func_name,output_file,tag_out)

def fetch(args):
    with open(args[0],'r') as fp:
        config = json.load(fp)[args[1]]
    
    fetch_dir = './fetch/'
    if len(args)>2:
        fetch_dir = args[2]
    os.system('mkdir -p {0:}'.format(fetch_dir))

    for job in config['job']:
        for dataset in config['dataset']:
            input_dir = config['output_dir'].format(job=job,dataset=dataset)
            basename = os.path.basename(input_dir.rstrip('/'))            
            cmd = 'hadd -f {0:}/{1}.root {2:}/*.root'.format(fetch_dir,basename,input_dir)            
            os.system(cmd)
            
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--submit',nargs='+')
    parser.add_argument('--run',nargs='+')
    parser.add_argument('--fetch',nargs='+')

    args = parser.parse_args()    
    if args.fetch!=None:
        fetch(args.fetch)
    elif args.submit!=None and args.run==None:
        submit(args.submit)
    elif args.submit==None and args.run!=None:
        run(args.run)
    else:
        raise ValueError('use "python ./launch.py --submit/run/fetch ... "')
