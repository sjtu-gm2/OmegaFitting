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

version = sys.argv[1] #'newfitter'
cut = sys.argv[2] #'Cut30'

f = R.TFile('./output/%s_%s.root'%(version,cut))
outpath = './plots_%s_%s/'%(version,cut)

subprocess.getstatusoutput('mkdir -p %s'%(outpath))

paras = {
    5 : '5paras',
    9 : '9paras_cbo',
    10 : '10paras_cbo_lost'
}

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

def GetNDF(nvar):
    graph = R.TGraph()
    graph.SetName('Chi2_%s'%(paras[nvar]))
    graph.SetTitle('Chi2 / NDF')
    for n in range(20):
        func_name = form_name('func',nvar,n)
        func = f.Get(func_name)
        # print (func)
        chi2 = func.GetChisquare()
        ndf = func.GetNDF()
        e_start = (2*n+21 -1)*50
        e_end = (2*n+22)*50 
        bin_center = (e_start + e_end)*0.5
        graph.SetPoint(n,bin_center,chi2/ndf)
    return graph


for nvar in [5,9,10]:
    # wiggle plot
    c = R.TCanvas('wiggle %sparas'%(nvar),'wiggle %sparas'%(nvar),2500,2000)
    c.Divide(5,4)
    for n in range(20):
        c.cd(n+1)
        name = form_name('wiggle',nvar,n)
        hist = f.Get(name)
        # print (name)
        hist.GetXaxis().SetRangeUser(30e3,100e3)
        hist.Draw()
    c.cd(0)
    c.Draw()
    c.SaveAs(outpath+'WigglePlot_%sParasFit.png'%(nvar))

    # fft     
    c = R.TCanvas('fft %sparas'%(nvar),'fft %sparas'%(nvar),2500,2000)
    c.Divide(5,4)
    for n in range(20):
        c.cd(n+1)
        name = form_name('fft',nvar,n)
        hist = f.Get(name)
        # print (name)
        hist.GetXaxis().SetRangeUser(0,1./0.1402/2.)
        hist.Draw()
    c.cd(0)
    c.Draw()
    c.SaveAs(outpath+'FFT_%sParasFit.png'%(nvar))


#kLoss
central_graph_kloss = R.TGraphErrors()
title = ';Energy bin center [MeV];kLoss'
kloss_array = []
kloss_error_array = []
centers = []
central_graph_kloss.SetName('kLoss graph')
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
    kloss_array.append(kloss)
    kloss_error_array.append(kloss_error)
    centers.append(bin_center)

    # print ('{0:<3} {3:<5.2f} {1:>10.2E} +- {2:<3.2E} '.format(n,kloss,kloss_error,bin_center))


central_graph_kloss.SetMarkerSize(0.9)
central_graph_kloss.SetLineStyle(1)
central_graph_kloss.SetLineWidth(1)
central_graph_kloss.SetMarkerStyle(8)

c = R.TCanvas()
central_graph_kloss.Draw('AP')
c.Draw()
c.SaveAs(outpath+'kLoss.png')


# Chi2 / NDF
ndf_5par = GetNDF(5)
ndf_9par = GetNDF(9)
ndf_10par = GetNDF(10)
c = R.TCanvas()
title = '%s parameters fit;Energy [MeV];Chi2/NDF'
deco(ndf_5par,R.kRed+1,title%(5))
deco(ndf_9par,R.kAzure+2,title%(9))
deco(ndf_10par,R.kGreen+2,title%(10))
mg = R.TMultiGraph()
mg.Add(ndf_5par)
mg.Add(ndf_9par)
mg.Add(ndf_10par)
mg.SetTitle('Chisquare in Energy slices fitting;Energy [MeV];Chi2 / NDF')
mg.Draw('ALP')

c.BuildLegend()
c.Draw()
c.SaveAs(outpath+'Chi2NDF.png')


