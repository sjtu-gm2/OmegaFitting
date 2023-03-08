import os,re
import subprocess
import ROOT as R
import time
import sys

weights_run2 = [
    0.0050,
    0.0278,
    0.0532,
    0.0812,
    0.1118,
    0.1445,
    0.1797,
    0.2172,
    0.2571,
    0.2993,
    0.3437,
    0.3905,
    0.4398,
    0.4914,
    0.5452,
    0.5999,
    0.6534,
    0.7035,
    0.7506,
    0.8018,
]

weights_run3a = [
    0.0047,
    0.0274,
    0.0531,
    0.0815,
    0.1124,
    0.1456,
    0.1813,
    0.2195,
    0.2602,
    0.3033,
    0.3485,
    0.3962,
    0.4465,
    0.4991,
    0.5538,
    0.6087,
    0.6615,
    0.7095,
    0.7508,
    0.7833,
]
weights_run3b = [
    0.0033,
    0.0254,
    0.0508,
    0.0787,
    0.1092,
    0.1419,
    0.1771,
    0.2146,
    0.2546,
    0.2968,
    0.3412,
    0.3880,
    0.4374,
    0.4889,
    0.5426,
    0.5966,
    0.6495,
    0.6984,
    0.7402,
    0.7741,
]
weights = {
    'run2' : weights_run2,
    'run3a' : weights_run3a,
    'run3b' : weights_run3b,

}

def MakeWiggles_T_method(file_path,out_path,tag):
    tf = R.TFile(file_path)
    outf = R.TFile('%s/wiggles_T_%s.root'%(out_path,tag),'recreate')
    wiggles = []
    for seed in range(200):
        dN = (seed // 20)*10
        dA = (seed % 20)*5
        raw = tf.Get('shadowperClusterdN{0:}dA{1:}/ET_raw'.format(dN,dA)).Clone();raw.SetDirectory(0)
        pud = tf.Get('shadowperClusterdN{0:}dA{1:}/ET_double'.format(dN,dA)).Clone();pud.SetDirectory(0)
        put = tf.Get('shadowperClusterdN{0:}dA{1:}/ET_triple'.format(dN,dA)).Clone();put.SetDirectory(0)

        corr = raw.Clone();corr.SetDirectory(0)
        corr.Add(pud,-1)
        corr.Add(put,-1)

        e_bin_start = corr.GetYaxis().FindBin(1700)
        e_bin_end   = corr.GetYaxis().FindBin(3100) - 1

        name = 'wiggle_T1700_dN%sdA%s'%(dN,dA)
        wiggle = corr.ProjectionX(name,e_bin_start,e_bin_end)


        wiggle.SetDirectory(outf)
        wiggles.append(wiggle)
        print( name,wiggle.Integral() )
        
    outf.Write()
    outf.Close()


def MakeWiggles_A_method(file_path,out_path,tag,period):
    

    tf = R.TFile(file_path)
    outf = R.TFile('%s/wiggles_A_%s.root'%(out_path,tag),'recreate')
    wiggles = []
    for seed in range(200):
        dN = (seed // 20)*10
        dA = (seed % 20)*5
        raw = tf.Get('shadowperClusterdN{0:}dA{1:}/ET_raw'.format(dN,dA)).Clone();raw.SetDirectory(0)
        pud = tf.Get('shadowperClusterdN{0:}dA{1:}/ET_double'.format(dN,dA)).Clone();pud.SetDirectory(0)
        put = tf.Get('shadowperClusterdN{0:}dA{1:}/ET_triple'.format(dN,dA)).Clone();put.SetDirectory(0)

        corr = raw.Clone();corr.SetDirectory(0)
        corr.Add(pud,-1)
        corr.Add(put,-1)
        # fill wiggle_A
        for n_slice in range(20):
            e_start = 1000+n_slice*100
            e_end = 1100+n_slice*100
            e_bin_start = corr.GetYaxis().FindBin(e_start)
            e_bin_end   = corr.GetYaxis().FindBin(e_end) - 1

            name = 'wiggle_E%sto%s_seed%s'%(e_start,e_end,seed)
            wiggle_slice = corr.ProjectionX(name,e_bin_start,e_bin_end);wiggle_slice.SetDirectory(0)
            if n_slice == 0:
                wiggle = wiggle_slice.Clone()
                wiggle.Reset()
                name_wiggle = 'wiggle_A_dN%sdA%s'%(dN,dA)
                wiggle.SetName(name_wiggle)
                wiggle.SetDirectory(outf)
                wiggles.append(wiggle)

            wiggle.Add(wiggle_slice,weights[period][n_slice])
        print (name_wiggle,wiggle.Integral())
    outf.Write()
    outf.Close()


def MakeWiggles_slices(file_path,out_path,tag):
    tf = R.TFile(file_path)
    outf = R.TFile('%s/wiggles_Eslices_%s.root'%(out_path,tag),'recreate')
    wiggles = []
    for seed in range(200):
        dN = (seed // 20)*10
        dA = (seed % 20)*5
        raw = tf.Get('shadowperClusterdN{0:}dA{1:}/ET_raw'.format(dN,dA)).Clone();raw.SetDirectory(0)
        pud = tf.Get('shadowperClusterdN{0:}dA{1:}/ET_double'.format(dN,dA)).Clone();pud.SetDirectory(0)
        put = tf.Get('shadowperClusterdN{0:}dA{1:}/ET_triple'.format(dN,dA)).Clone();put.SetDirectory(0)

        corr = raw.Clone();corr.SetDirectory(0)
        corr.Add(pud,-1)
        corr.Add(put,-1)

        for n_slice in range(20):            
            e_start = 1000+n_slice*100
            e_end = 1100+n_slice*100
            e_bin_start = corr.GetYaxis().FindBin(e_start)
            e_bin_end   = corr.GetYaxis().FindBin(e_end) - 1

            name = 'wiggle_E%s_dN%sdA%s'%(e_start,dN,dA)
            wiggle = corr.ProjectionX(name,e_bin_start,e_bin_end)


            wiggle.SetDirectory(outf)
            wiggles.append(wiggle)
            
            print (name,wiggle.Integral())

        
    outf.Write()
    outf.Close()

if __name__ == '__main__':

    # for p in ['2C','3D','3b']:
    for p in ['3b']:
        file_path = '/home/chencheng/data/Adhoc_study/Run%s_adhoc.root'%(p)
        out_path = '/home/chencheng/data/Adhoc_study/'
        tag = 'run%s'%(p)
        period = 'run%s'%(p)

        MakeWiggles_T_method(file_path,out_path,tag)
        MakeWiggles_A_method(file_path,out_path,tag,period)
        MakeWiggles_slices(file_path, out_path, tag)