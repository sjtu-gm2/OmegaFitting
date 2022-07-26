import os
import re
import json
import sys
from collections import OrderedDict

import ROOT as R
R.gROOT.SetBatch(1)

class Output:
    def __repr__(self):
        res = '-----------\n'
        res += 'Status = %s\n'%(self.status)
        res += 'Chi2/NDF = %.2f (%.2f/%d)\n'%(self.chi2/self.ndf,self.chi2,self.ndf)
        res += 'nvars = %s\n'%(self.nvars)
        if self.nvars!=None:
            for n in range(self.nvars):
                res+= '{0:<15} {1:>10.2E} +- {2:<2.1E}\n'.format(*self.paras[n])
        return res

    def __init__(self):
        self.status = None
        self.chi2 = None
        self.ndf = None
        self.paras = None

    def getvals(self):
        res = []
        for n in range(self.nvars):
            val = self.paras[n][1]
            if abs(val)<1e-40:
                val = 0.
            res.append(val)
        return res

def ParseTxt(entry):
    status = re.findall('STATUS=([^ ]+) ', entry[0])[0]
    chi2, ndf = re.findall('chi2=([^ ]+) ndf=([0-9]+)', entry[-1])[0]
    out = Output()
    out.chi2 = float(chi2)
    out.ndf = int(ndf)
    out.status = status
    paras = []
    for line in entry[4:-1]:
        res = re.split('[ ]+',line)
        num = int(res[1]) - 1
        name = res[2]
        val = float(res[3])
        err = float(res[4])
        paras.append([name,val,err])
    out.paras = paras
    out.nvars = len(paras)
    return out

def ExtractTxt(file_name):
    with open(file_name,'r') as f:
        lines = f.readlines()
    start = False
    end = False
    entries = []
    for line in lines:
        line = line[:-1]        
        if 'FCN=' in line:
            start = True
            entry = []
        if 'Chi2/NDF = ' in line:            
            end = True
        if start == True:
            entry.append(line)
        if end == True:
            entries.append(ParseTxt(entry))
            start = False
            end = False
    return entries

def ExtractROOT(file_name,func_name):
    f = R.TFile(file_name,'read')
    
    func = f.Get(func_name)    
    if func==None:
        return None
    chi2 = func.GetChisquare()
    ndf = func.GetNDF()
    nvars = func.GetNpar()
    paras = []
    for n in range(nvars):
        name = func.GetParName(n)
        val = func.GetParameter(n)
        err = func.GetParError(n)
        paras.append([name,val,err])
    entry = Output()
    entry.chi2 = chi2
    entry.ndf = ndf
    entry.nvars = nvars
    entry.paras = paras
    return entry


def GetInitValuesFromROOT(file_name,json_name):
    entries = OrderedDict()
    for n in range(20):
        e_low = 1000 + 100*n
        e_high = 1100 + 100*n
        func_name = 'func_10paras_cbo_lost_Slice_%s_%s_MeV'%(e_low,e_high)
        entry = ExtractROOT(file_name, func_name)        
        if entry == None:
            continue
        entries[n] = entry.getvals()        
    with open(json_name,'w') as json_f:
        json.dump(entries, json_f, sort_keys=True, indent=4)

def GetInitValuesFromROOT_TmethodFit(file_name,json_name):
    entries = OrderedDict()

    e_low = 1700
    e_high = 3100
    func_name = 'func_10paras_cbo_lost_Slice_%s_%s_MeV'%(e_low,e_high)
    entry = ExtractROOT(file_name, func_name)        
    
    entries[0] = entry.getvals()        
    with open(json_name,'w') as json_f:
        json.dump(entries, json_f, sort_keys=True, indent=4)


def GetInitValuesFromROOT_Run23Official(file_name,json_name):
    entries = OrderedDict()

    
    func_name = '%smethod_fitfunc_seed_%s'%(method,seed)
    entry = ExtractROOT(file_name, func_name)        
    
    entries[0] = entry.getvals()        
    with open(json_name,'w') as json_f:
        json.dump(entries, json_f, sort_keys=True, indent=4)

def GetInitValuesFromROOT_Run23SystScan(file_name,json_name,method,syst):
    entries = OrderedDict()

    
    func_name = 'func_28paras_run23_sjtu_%s_%smethod_%s'%(dataset,method,syst)
    entry = ExtractROOT(file_name, func_name)        
    
    entries[0] = entry.getvals()        
    with open(json_name,'w') as json_f:
        json.dump(entries, json_f, sort_keys=True, indent=4)

    print ('generated: %s'%(json_name))


# entries = ExtractTxt('./logs/newfitter_All.log')
# for entry in entries:
#     print ('-----------------')
#     print (entry)
    


# for n in range(20):
#     entry = ExtractROOT('./output/newfitter_All.root', 'func_9paras_cbo_Slice_2500_2600_MeV')
#     print (entry)


# version = sys.argv[1]
# cut = sys.argv[2]


# dataset = "Run2all"

# seed = "gain_A_10"


# GetInitValuesFromROOT_Run23Official('./data/%s_seedscan_result_Bingzhi_149ns.root'%(dataset), 'pre_fitted_values_%s_%smethod_seed_%s.json'%(dataset,method,seed))
# GetInitValuesFromROOT_Run23Official('./data/%s_seedscan_result_Bingzhi_149ns.root'%(dataset), 'pre_fitted_values_%s_%smethod_seed_%s.json'%(dataset,method,seed))
# GetInitValuesFromROOT_Run23SystScan('./output/Run2C_Amethod_gain_A_Scan/result_28paras_run23_sjtu_Run2C_Amethod_gain_A_10.root', 'fitter_28paras_run23_sjtu_Run2C_Amethod_gain_A_10.json')
# GetInitValuesFromROOT_Run23SystScan('./output/Run2C_Amethod_gain_A_Scan_random_adhoc.root', 'fitter_28paras_run23_sjtu_Run2C_Amethod_gain_A_random_adhoc_10.json','A')
# GetInitValuesFromROOT_Run23SystScan('./output/Run2C_Tmethod_gain_A_Scan_random_adhoc.root', 'fitter_28paras_run23_sjtu_Run2C_Tmethod_gain_A_random_adhoc_10.json','T')

def main():
    for method in ['T','A']:
        for gain in ['stdp','gain']:
            for sys in ['A','T']:
                for version in ['random','random_adhoc']:
                    scan = '%s_%s'%(gain,sys)
                    
                    outputDir = 'json/Run2C_%smethod_%s_%s_Scan_%s'%(method,gain,sys,version)
                    for m in range(20):
                        inputROOT = 'output/Run2C_%smethod_%s_Scan_%s/result_28paras_run23_sjtu_Run2C_%smethod_%s_%s.root'%(method,scan,version,method,scan,m)
                        json_name = 'Scan_%s.json'%(m)
                        syst = '%s_%s_%s'%(gain,sys,m)
                        os.system('mkdir -p %s'%(outputDir))
                        GetInitValuesFromROOT_Run23SystScan(inputROOT, '%s/%s'%(outputDir,json_name),method,syst)


if __name__ == '__main__':
    method = sys.argv[1]
    scan = sys.argv[2]
    m = sys.argv[3]
    version = sys.argv[4]
    dataset = sys.argv[5]

    # method = 'T'
    # scan = 'gain_A'
    # m = '10'
    # version = 'sys.argv[4]'
    # if version.startswith('Run4U'): 
    #     dataset = 'Run4U'
    # else:
    #     dataset = 

    # inputROOT = 'output/Run4U_Tmethod_gain_A_Scan_Run4U_preliminary.root'
    # outputDir = './'


    inputROOT = 'output/%s_%smethod_%s_Scan_%s/result_28paras_run23_sjtu_%s_%smethod_%s_%s.root'%(dataset,method,scan,version,dataset,method,scan,m)
    outputDir = 'json/%s_%smethod_%s_Scan_%s'%(dataset,method,scan,version)


    os.system('mkdir -p %s'%(outputDir))
    syst = '%s_%s'%(scan,m)
    json_name = 'Scan_%s.json'%(m)

    GetInitValuesFromROOT_Run23SystScan(inputROOT, '%s/%s'%(outputDir,json_name),method,syst)
    # main()
# GetInitValuesFromROOT('./output/%s_%s.root'%(version,cut), 'fitted_values_%s_%s.json'%(version,cut))
# GetInitValuesFromROOT_TmethodFit('./output/%s_%s.root'%(version,cut), 'fitted_values_%s_%s.json'%(version,cut))
