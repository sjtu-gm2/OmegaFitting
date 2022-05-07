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
        entries[n] = entry.getvals()        
    with open(json_name,'w') as json_f:
        json.dump(entries, json_f, sort_keys=True, indent=4)














# entries = ExtractTxt('./logs/newfitter_All.log')
# for entry in entries:
#     print ('-----------------')
#     print (entry)
    


# for n in range(20):
#     entry = ExtractROOT('./output/newfitter_All.root', 'func_9paras_cbo_Slice_2500_2600_MeV')
#     print (entry)


version = sys.argv[1]
cut = sys.argv[2]

GetInitValuesFromROOT('./output/%s_%s.root'%(version,cut), 'fitted_values_%s_%s.json'%(version,cut))

