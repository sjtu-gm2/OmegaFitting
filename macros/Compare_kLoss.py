import ROOT as R
import subprocess
import sys
R.gROOT.SetBatch(1)

R.gStyle.SetStatBorderSize(0)
R.gStyle.SetStatX(.94)
R.gStyle.SetStatY(.94)
R.gStyle.SetPalette(57)
R.gStyle.SetOptFit(1111)
R.gStyle.SetOptStat(0)

outpath = './plots_cmp_kLoss/'
subprocess.getstatusoutput('mkdir -p %s'%(outpath))

colors = [R.kRed+1, R.kAzure+2, R.kGreen+2]

cmps = {}

cmps['newfitter'] = [
    ['newfitter','All'],
    ['newfitter','Cut30'],
    ['newfitter','Cut40'],
]

cmps['oldfitter'] = [
    ['oldfitter','All'],
    ['oldfitter','Cut30'],
    ['oldfitter','Cut40'],
]

cmps['All'] = [
    ['newfitter','All'],
    ['oldfitter','All'],
]

cmps['Cut30'] = [
    ['newfitter','Cut30'],
    ['oldfitter','Cut30'],
]

cmps['Cut40'] = [
    ['newfitter','Cut40'],
    ['oldfitter','Cut40'],
]

paras = {
    5 : '5paras',
    9 : '9paras_cbo',
    10 : '10paras_cbo_lost'
}

names = {
    'newfitter' : 'kLoss of newfitter',
    'oldfitter' : 'kLoss of oldfitter',
    'All' : 'kLoss without hit energy cut',
    'Cut30' : 'Hit energy Cut @ 30 MeV',
    'Cut40' : 'Hit energy Cut @ 40 MeV',
}

def DrawCmp(name,versions):
    c = R.TCanvas()
    mg = R.TMultiGraph()    
    n=0
    for version in versions:
        graph = Get_kLoss(version)
        deco(graph, colors[n]);n+=1
        mg.Add(graph)
    title = '%s;Energy bin center [MeV];kLoss'%(names[name])
    mg.SetTitle(title)
    mg.Draw('ALP')
    c.BuildLegend()
    c.Draw()
    c.SaveAs(outpath + 'CMP_%s.png'%(name))

def Get_kLoss(version):
    f = R.TFile('./output/%s_%s.root'%(version[0],version[1]))
    central_graph_kloss = R.TGraphErrors()
    title = '%s %s MeV;Energy bin center [MeV];kLoss'%(version[0],version[1])
    central_graph_kloss.SetName('kLoss_graph_%s_%s_MeV'%(version[0],version[1]))
    central_graph_kloss.SetTitle(title)

    for n in range(20):
        func_name = form_name('func',10,n)
        func = f.Get(func_name)
        kloss = func.GetParameter(9)
        kloss_error = func.GetParError(9)
        e_start = (2*n+21 -1)*50
        e_end = (2*n+22)*50 
        bin_center = (e_start + e_end)*0.5
        
        
        central_graph_kloss.SetPoint(n,bin_center,kloss)
        central_graph_kloss.SetPointError(n,0.,kloss_error)    
    return central_graph_kloss



def form_name(tp,nvar,nslice):

    
    e_start = (2*nslice+21 -1)*50;
    e_end = (2*nslice+22) * 50;     
    
    return '%s_%s_Slice_%s_%s_MeV'%(tp,paras[nvar],e_start,e_end)

def deco(graph,color,title=None):
    graph.SetMarkerSize(0.9)
    graph.SetLineStyle(1)
    graph.SetLineWidth(1)
    graph.SetMarkerStyle(8)
    graph.SetLineColor(color)
    graph.SetMarkerColor(color)
    if title!=None:
        graph.SetTitle(title)

for name,versions in cmps.items():
    DrawCmp(name, versions)



