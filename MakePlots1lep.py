import ROOT
import os
import sys
import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from ROOT import gROOT

from infofile import infos

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadLeftMargin(0.13)
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetGridStyle(2)
ROOT.gStyle.SetPadLeftMargin(0.13)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

channel = sys.argv[1]  # ee or uu


# Variables

variables = [ 'eta1', 'phi1', 'etcone', 'etconeRel', 'mT_QCDregion', 'mT', 'mT_QCDregionlog', 'mTLog', 'pT_QCDregion','met_QCDregion', 'met', 'pTLog_QCDregion','pTLog', 'pt1']

xtitles = {'pt1': 'Leading lepton p_{T} (GeV)', 'eta1': 'Leading lepton #eta', 'phi1': 'Leading lepton #phi', 'met': 'E_{T}^{miss} (GeV)', 'etcone':'isolation variable','etconeRel':'relative isolation', 'mT_QCDregion':'m (GeV)', 'mT': 'm_{T} (GeV)', 'mT_QCDregionlog':'m (GeV)', 'mTLog': 'm_{T} (GeV)', 'pT_QCDregion': 'Leading lepton p_{T} (GeV)', 'met_QCDregion': 'E_{T}^{miss} (GeV)','pTLog_QCDregion':'Leading lepton p_{T} (GeV)','pTLog':'Leading lepton p_{T} (GeV)'}


# Backgrounds

backgrounds = ['Wjets','Zjets', 'Top','Diboson']#put the higher one first because this is the order they are stacked in the histogram, if not, it will eclipse the others


Zjets = [364100, 364101, 364102, 364103, 364104, 364105, 364106, 364107, 364108, 364109, 364110, 364111, 364112, 364113, 364114, 364115, 364116, 364117, 364118, 364119, 364120, 364121, 364122, 364123, 364124,
         364125, 364126, 364127, 364128, 364129, 364130, 364131, 364132, 364133, 364134, 364135, 364136, 364137, 364138, 364139, 364140, 364141]

Wjets = [364156, 364157, 364158, 364159, 364160, 364161, 364162, 364163, 364164, 364165, 364166, 364167, 364168, 364169, 364170, 364171, 364172, 364173, 364174, 364175, 364176, 364177, 364178, 364179, 364180,
         364181, 364182, 364183, 364184, 364185, 364186, 364187, 364188, 364189, 364190, 364191, 364192, 364193, 364194, 364195, 364196, 364197]

Diboson = [363356, 363358, 363359, 363360, 363489, 363490, 363491, 363492, 363493] 

Top = [410000, 410011, 410012, 4100013, 410014, 410025, 410026]

signals = ['Wprime2000', 'Wprime3000', 'Wprime4000', 'Wprime5000']

Zprime2000 = [301329];

Zprime3000 = [301333];

Zprime4000 = [301337];

Zprime5000 = [301341];

Wprime2000 = [311215, 311220];

Wprime3000 = [311216, 311221];

Wprime4000 = [311217, 311222];

Wprime5000 = [311218, 311223];

fileIDs = {'Diboson':Diboson, 'Zjets':Zjets, 'Wjets':Wjets, 'Top':Top, 'Zprime2000':Zprime2000, 'Zprime3000':Zprime3000, 'Zprime4000':Zprime4000, 'Zprime5000':Zprime5000, 'Wprime2000':Wprime2000, 'Wprime3000':Wprime3000, 'Wprime4000':Wprime4000, 'Wprime5000':Wprime5000}

hist_bkg = {}#array of histograms of the background
hist_sig = {}#array for signal(s)
for var in variables:
    hist_bkg[var] = {}#for each variable to plot(invariant mass or whatever)
    hist_sig[var] = {}
    for bkg in backgrounds:
        hist_bkg[var][bkg] = ROOT.TH1F()#each of the backgrounds contributing
    for sig in signals:
        hist_sig[var][sig] = ROOT.TH1F()

colours = dict(Diboson=ROOT.kAzure+1, Top=ROOT.kRed+1, Zjets=ROOT.kOrange-2, Wjets=ROOT.kGray, Zprime2000=ROOT.kCyan, Zprime3000=ROOT.kMagenta, Zprime4000=ROOT.kYellow, Zprime5000=ROOT.kGreen, Wprime2000=ROOT.kBlack, Wprime3000=ROOT.kViolet, Wprime4000=ROOT.kOrange, Wprime5000=ROOT.kPink)#colours to plot variables and signal


# Extract info about cross section and sum of weights from infofile

info = {}
for key in infos.keys():
    ID = infos[key]['DSID']
    info[ID] = {}
    info[ID]['xsec'] = infos[key]['xsec']
    info[ID]['sumw'] = infos[key]['sumw']
    info[ID]['events'] = infos[key]['events']

L = 10.6  # integrated luminosity = 10.6 fb^-1


# Function for making histograms


def fill_hist(h, h_name, key, ID):

    h_midl = infile.Get(h_name).Clone("h_midl")

    xsec = 1000*info[ID]['xsec']
    nev = info[ID]['sumw']

    N_mc = xsec*L

    sf = N_mc/nev

    if not h.GetName():
        h = infile.Get(h_name)
        h.Scale(sf)
        n = h.GetNbinsX()
        for i in range(n):
            bc = h.GetBinContent(i)
            if bc < 0:
                h.SetBinContent(i, 0)
        h.SetFillColor(colours[key])
        h.SetLineColor(colours[key])
    else:
        h_midl.Scale(sf)
        n = h_midl.GetNbinsX()
        for i in range(n):
            bc = h_midl.GetBinContent(i)
            if bc < 0:
                h_midl.SetBinContent(i, 0)
        h.Add(h_midl)

    return h


# Loop over files in MC directory

for filename in os.listdir('Escritorio/CodeExample/Histograms1lep good one/MC/'):
    if '.root' in filename:
        filepath = 'Escritorio/CodeExample/Histograms1lep good one/MC/'+filename
        infile = ROOT.TFile(filepath)
        file_id = int(filename.split('.')[2])
        # print filename
        for var in variables:
            for bkg in backgrounds:
                if file_id in fileIDs[bkg]:
                    #print(var)
                    hist_bkg[var][bkg] = fill_hist(
                        hist_bkg[var][bkg], 'h_'+channel+'_'+var, bkg, file_id)
            for sig in signals:
                if file_id in fileIDs[sig]: 
                    hist_sig[var][sig] = fill_hist(hist_sig[var][sig], 'h_'+channel+'_'+var, sig, file_id)

# Get data

data = ROOT.TFile('Escritorio/CodeExample/Histograms1lep good one/Data/hist.Data.2016.root')
hist_d = {}
#variables2 adds mtvsIsolation and mTvsIsolationrel, to construct 2D histograms with the data only
variables2=['mTvsIsolation', 'mTvsIsolationrel']
for i in range(0, len(variables)):
    variables2.append(variables[i])
xtitles2= {'mTvsIsolation':'isolation variable','mTvsIsolationrel':'relative isolation'}
xtitles2.update(xtitles)
for var in variables2:
    hist_d[var] = data.Get('h_'+channel+'_'+var)
    hist_d[var].SetMarkerStyle(20)
    hist_d[var].SetMarkerSize(0.7)
    hist_d[var].SetLineColor(ROOT.kBlack)
    hist_d[var].GetYaxis().SetTitle("Events")
    hist_d[var].GetXaxis().SetTitle(xtitles2[var])
    hist_d[var].GetXaxis().SetTitleFont(43)
    hist_d[var].GetXaxis().SetTitleSize(16)
    hist_d[var].GetYaxis().SetTitleFont(43)
    hist_d[var].GetYaxis().SetTitleSize(16)
    hist_d[var].GetXaxis().SetLabelFont(43)
    hist_d[var].GetXaxis().SetLabelSize(16)
    hist_d[var].GetYaxis().SetLabelFont(43)
    hist_d[var].GetYaxis().SetLabelSize(16)
    hist_d[var].GetXaxis().SetTitleOffset(4)
    hist_d[var].GetYaxis().SetTitleOffset(1.5)


# Style histograms, make stack and histograms with full background

stack = {}
hist_r = {}#this does the data/bkg ratio
hist_mc = {}

#pure QCD histogram


hist_pureQCDlog=hist_d['mT_QCDregionlog'].Clone()
###histograms to add qcd as background
variables3=['mT_QCDregion', 'pT_QCDregion', 'met_QCDregion', 'pTLog_QCDregion']#important that these variables come before their respective normal ones in the "variables" array so that the loop works correctly
hist_pureQCD={}
addmedef={}
histogramtoadd={}
addme={}
for var in variables3:
    hist_pureQCD[var]=hist_d[var].Clone()#histpureqcd to get data of qdc region minus allbk
    addmedef[var]=hist_d[var].Clone()
    addmedef[var].Add(hist_d[var],-1.)
addmedef2=hist_d['mT_QCDregionlog'].Clone()#clone it to get same binning but with bins to 0
addmedef2.Add(hist_d['mT_QCDregionlog'],-1.)

escala=2.9

for var in variables:
    stack[var] = ROOT.THStack(var, "")
    hist_mc[var] = ROOT.TH1F()
    hist_r[var] = ROOT.TH1F()
    for bkg in reversed(backgrounds):
        hist_bkg[var][bkg].GetYaxis().SetTitle("Events")
        hist_bkg[var][bkg].GetXaxis().SetTitle(xtitles[var])
        hist_bkg[var][bkg].GetXaxis().SetTitleFont(43)
        hist_bkg[var][bkg].GetXaxis().SetTitleSize(16)
        hist_bkg[var][bkg].GetYaxis().SetTitleFont(43)
        hist_bkg[var][bkg].GetYaxis().SetTitleSize(16)
        hist_bkg[var][bkg].GetXaxis().SetLabelFont(43)
        hist_bkg[var][bkg].GetXaxis().SetLabelSize(16)
        hist_bkg[var][bkg].GetYaxis().SetLabelFont(43)
        hist_bkg[var][bkg].GetYaxis().SetLabelSize(16)
        hist_bkg[var][bkg].GetXaxis().SetTitleOffset(4)
        hist_bkg[var][bkg].GetYaxis().SetTitleOffset(1.5)
        stack[var].Add(hist_bkg[var][bkg])
        if var in variables3:
            hist_pureQCD[var].Add(hist_bkg[var][bkg],-1.)
        if var=='mT_QCDregionlog':
            hist_pureQCDlog.Add(hist_bkg[var][bkg],-1.)
        if not hist_mc[var].GetName():
            hist_mc[var] = hist_bkg[var][bkg].Clone()
        else:
            print("Here")
            hist_mc[var].Add(hist_bkg[var][bkg])
            #qcd aditions
            if var=='pt1' and bkg==backgrounds[0]:
                hist_mc[var].Add(addmedef['pT_QCDregion'])
            if var in ['mT','met','pTLog'] and bkg==backgrounds[0]:#we dont want to add the same thing multiple times
                hist_mc[var].Add(addmedef[var+'_QCDregion'])#addme gets data at the end of the loop
                print(var)
            if var=='mTLog' and bkg==backgrounds[0]:
                hist_mc[var].Add(addmedef2)

            print("Bg = "+bkg)
            print("Done")

        hist_r[var] = hist_d[var].Clone()
        hist_r[var].Divide(hist_mc[var])
        hist_r[var].SetTitle("")
        hist_r[var].GetXaxis().SetTitle(xtitles[var])
        hist_r[var].GetYaxis().SetTitle("Data/#SigmaMC")
        hist_r[var].GetYaxis().SetNdivisions(506)
        hist_r[var].SetMarkerStyle(20)
        hist_r[var].SetMarkerSize(0.7)

    if var in variables3:
        #filter negative values inside histQCD
        for i in range(1, hist_pureQCD[var].GetNbinsX()):
            if hist_pureQCD[var].GetBinContent(i)<0:
                hist_pureQCD[var].SetBinContent(i,0)

        addme[var]=hist_pureQCD[var].Clone()
        addme[var].Scale(escala)
        addmedef[var].Add(addme[var])
        
    if var=='mT_QCDregionlog':
        #filter negative values inside histQCD in case qcd "data"<backgrounds
        for i in range(1, hist_pureQCDlog.GetNbinsX()):
            if hist_pureQCDlog.GetBinContent(i)<0:
                hist_pureQCDlog.SetBinContent(i,0)
        addmelog=hist_pureQCDlog.Clone()
        addmelog.Scale(escala)
        addmedef2.Add(addmelog)
        addmedef2.SetFillColor(ROOT.kGreen)
        addmedef2.SetLineColor(ROOT.kGreen)
#Adding QCD pure to mT background stack
for var in variables3:
    histogramtoadd[var]=hist_pureQCD[var].Clone()
    histogramtoadd[var].Scale(escala)
    histogramtoadd[var].SetFillColor(ROOT.kGreen)
    histogramtoadd[var].SetLineColor(ROOT.kGreen)
    if var=='pT_QCDregion':
        stack['pt1'].Add(histogramtoadd[var])
    else:
        stack[var.split('_')[0]].Add(histogramtoadd[var])
histogramtoaddlog=hist_pureQCDlog.Clone()
histogramtoaddlog.Scale(escala)
stack['mTLog'].Add(addmedef2)




# Make plot legend

leg = ROOT.TLegend(0.70, 0.50, 0.88, 0.88)
leg.SetFillStyle(4000)
leg.SetFillColor(0)
leg.SetTextFont(42)
leg.SetBorderSize(0)

bkg_labels = {'Zjets': 'Z+jets', 'Top': 'Top',
              'Diboson': 'Diboson', 'Wjets': 'W+jets'}

sig_labels = {'Zprime2000':"Z' (2 TeV)", 'Zprime3000':"Z' (3 TeV)", 'Zprime4000':"Z' (4 TeV)", 'Zprime5000':"Z' (5 TeV)", 'Wprime2000':"W' (2 TeV)", 'Wprime3000':"W' (3 TeV)", 'Wprime4000':"W' (4 TeV)", 'Wprime5000':"W' (5 TeV)"}

for bkg in backgrounds:
    leg.AddEntry(hist_bkg['pt1'][bkg], bkg_labels[bkg], "f")

legQCD=leg.Clone()#not to include qcd as background in the qcd plot

leg.AddEntry(addmedef2, 'QCD jets', "f")#qcd analysis background

for sig in signals: 
    hist_sig['pt1'][sig].SetFillColor(0);

    leg.AddEntry(hist_sig['pt1'][sig], sig_labels[sig], "f")
    legQCD.AddEntry(hist_sig['pt1'][sig], sig_labels[sig], "f")

leg.AddEntry(hist_d['pt1'], "Data", "ple")
legQCD.AddEntry(hist_d['pt1'], "Data", "ple")

selection = ""
if channel == "ee":
    selection = "e #nu"
if channel == "uu":
    selection = "#mu#nu"

# Make plots
logvariables=['mTLog','mT_QCDregionlog','pTLog','pTLog_QCDregion']
qcdvariables=['mT_QCDregionlog','pTLog_QCDregion','pT_QCDregion']#variables we don't want showing the QCD jets legend
for var in variables:
    
    cnv = ROOT.TCanvas("cnv_"+var, "", 500, 500)
    cnv.SetTicks(1, 1)
    cnv.SetLeftMargin(0.13)

    p1 = ROOT.TPad("p1", "", 0, 0.35, 1, 1)
    p2 = ROOT.TPad("p2", "", 0, 0.0, 1, 0.35)

    p1.SetLogy()
    
    if var in logvariables:#logarithmic x axis plot
        p1.SetLogx()
    p1.SetBottomMargin(0.0)
    p1.Draw()
    p1.cd()

    stack[var].Draw("hist")
    stack[var].SetMinimum(10E-2)
    stack[var].GetYaxis().SetTitle("Events")
    stack[var].GetYaxis().SetTitleFont(43)
    stack[var].GetYaxis().SetTitleSize(16)
    stack[var].GetYaxis().SetLabelFont(43)
    stack[var].GetYaxis().SetLabelSize(16)
    stack[var].GetYaxis().SetTitleOffset(1.5)
    if var in ['eta1', 'eta2', 'phi1', 'phi2','etcone','etconeRel','pTLog'] or var in qcdvariables:
        maximum = stack[var].GetMaximum()
        stack[var].SetMaximum(maximum*10E4)

    if var=='mT_QCDregionlog':#this fixes the qcd plot y axis limit of the representation
        stack[var].SetMaximum(10000)
   
    hist_d[var].Draw("same e0")


    if var in qcdvariables:
        legQCD.Draw("same")#legQCD doesnt have the qcd jets legend
    else:
        leg.Draw("same")

    for sig in signals:
        hist_sig[var][sig].SetFillColor(0);
        hist_sig[var][sig].Draw("same hist");

    

    s = ROOT.TLatex()
    s.SetNDC(1)
    s.SetTextAlign(13)
    s.SetTextColor(ROOT.kBlack)
    s.SetTextSize(0.044)
    s.DrawLatex(0.46, 0.86, "#font[72]{ATLAS} Open Data")
    s.DrawLatex(0.46, 0.81, "#bf{#sqrt{s} = 13 TeV,^{}%.1f^{}fb^{-1}}" % (L))
    s.DrawLatex(0.46, 0.76, "#bf{"+selection+" selection}")

    p1.Update()
    p1.RedrawAxis()

    cnv.cd()

    p2.Draw()
    p2.cd()

    p2.SetGridy()
    if var in logvariables:
        p2.SetLogx()
    hist_r[var].SetMaximum(2.99)
    if var in qcdvariables or var=='etcone' or var=='etconeRel':
        hist_r[var].SetMaximum(600)
    hist_r['mT_QCDregionlog'].SetMaximum(20)
    hist_r['pTLog'].SetMaximum(8)
    hist_r['met'].SetMaximum(5)
    hist_r[var].SetMinimum(0.01)
    hist_r[var].Draw("0PZ")

    p2.SetTopMargin(0)
    p2.SetBottomMargin(0.35)
    p2.Update()
    p2.RedrawAxis()

    cnv.cd()
    cnv.Update()
    cnv.Print('Plots1lep/'+channel+'_'+var+'.png')
    cnv.Close()

####Pure QCD plot
cnv = ROOT.TCanvas("cnv_pureQCD", "", 500, 500)
cnv.SetTicks(1, 1)
cnv.SetLeftMargin(0.13)
hist_pureQCD['mT_QCDregion'].Draw()
cnv.cd()
cnv.Update()
cnv.Print('Plots1lep/pureQCD.png')
cnv.Close()
#pure qcd log plot
cnv = ROOT.TCanvas("cnv_pureQCDlog", "", 500, 500)
cnv.SetTicks(1, 1)
cnv.SetLeftMargin(0.13)
p1 = ROOT.TPad("p1", "", 0, 0.35, 1, 1)
p1.SetLogx()
p1.Draw()
p1.cd()
hist_pureQCDlog.Draw()
cnv.cd()
cnv.Update()
cnv.Print('Plots1lep/pureQCDlog.png')
cnv.Close()

############
#Not tight study for data: mT vs Isolation and mT vs relative isolation
cnv = ROOT.TCanvas("cnv_mTIsolation", "", 500, 500)
cnv.SetTicks(1, 1)
cnv.SetLeftMargin(0.13)
hist_d['mTvsIsolation'].Draw("COLZ")
cnv.cd()
cnv.Update()
cnv.Print('Plots1lep/mTvsIsolation.png')
cnv.Close()
cnv = ROOT.TCanvas("cnv_mTIsolationrel", "", 500, 500)
cnv.SetTicks(1, 1)
cnv.SetLeftMargin(0.13)
hist_d['mTvsIsolationrel'].Draw("COLZ")
cnv.cd()
cnv.Update()
cnv.Print('Plots1lep/mTvsIsolationrel.png')
cnv.Close()

###############
#histograms for magnar
#histogramdata=hist_d['mTLog']
#histogrampureQCD=addmelog.Clone()
#histogrampureQCD.Scale(1/escala)
#histogrammontecarlo=hist_mc['mTLog'].Clone()
#histogrammontecarlo.Add(addmelog,-1.)#pure montecarlo without QCD
#data.Close()
#datamagnar= ROOT.TFile("hist.data.root","RECREATE");
#histogramdata.Write();
#datamagnar.Close();
#pureQCD= ROOT.TFile("hist.pureQCD.root","RECREATE");
#histogrampureQCD.Write();
#pureQCD.Close()
#montecarlomagnar= ROOT.TFile("hist.montecarlo.root","RECREATE");
#histogrammontecarlo.Write();
#montecarlomagnar.Close()
########################################
#optimization and table
####Integration of the histogram. Need to integrate the background, the data and the signal for the invariant mass plot (mT).

totalbkgmT=stack['mT'].GetStack().Last()
datamT=hist_d['mT']
signalmT=hist_sig['mT']['Wprime2000']

#find optimal bin cut for which Z, significance, is maximum. The Integral method takes initial bin, last bin, so receives integer as a number. cant minimize with a normal minimizing method that uses double. for length of the histogram use GetSize(). If 100 bins, returns 102(100 bins+underflow+overflow). Size of data is 52
#print(datamT.Integral(0, datamT.GetSize()-2))

def minimizeme(cut):
    n_obs=datamT.Integral(cut, datamT.GetSize()-2)#also use integralanderror
    b=totalbkgmT.Integral(cut, totalbkgmT.GetSize()-2)
    if n_obs==0:#limit of xlog(x) when x-->0 is 0
        Zobs=math.sqrt(2*b)
    else:
        Zobs=math.sqrt(2*(n_obs*math.log(n_obs/b)-n_obs+b))
    return Zobs

def ZwithoutWpeak(cut):
    n_obs=datamT.Integral(cut, datamT.GetSize()-2)#also use integralanderror
    b=totalbkgmT.Integral(cut, totalbkgmT.GetSize()-49)
    if n_obs==0:#limit of xlog(x) when x-->0 is 0
        Zobswithoutpeak=math.sqrt(2*b)#but this is nonsense
    else:
        Zobswithoutpeak=math.sqrt(2*(n_obs*math.log(n_obs/b)-n_obs+b))
        return Zobswithoutpeak
#print(minimizeme(10))
bincuts=[]
possibleoptimalZ=[]
#search for local maximums and save them
for i in range(0,datamT.GetSize()-110):#dont put the cut too far, no data last bins, change limits in myselector, need to rerun (aprox 1hour)
    if minimizeme(i)>minimizeme(i+1):
        possibleoptimalZ.append(minimizeme(i))
        bincuts.append(i)
#get the absolute maximum, optimal significance Z
optimalZ=max(possibleoptimalZ)
#print(possibleoptimalZ)
#print(datamT.GetBin(bincuts[5]))
#print("the maximum observed significance is "+ str(optimalZ))

#expected significance
masscut=[1500, 2250, 3000, 3750]#75%mass cut. Divergence for 1500????

a=-1
arrayfortable=[]
scuts=[]
serrors=[]
for sig in signals:
    a=a+1
    signalmT=hist_sig['mT'][sig]
    def expectedsignificance(cut):
        b=totalbkgmT.Integral(cut, totalbkgmT.GetSize()-2)
        s=signalmT.Integral(cut, signalmT.GetSize()-2)
        if b==0:
            return 0
        else:
            Zexp=math.sqrt(2*((s+b)*math.log(1+s/b)-s))
            return Zexp

    def confidencelimit(cut):
        b=totalbkgmT.Integral(cut, totalbkgmT.GetSize()-2)
        s=signalmT.Integral(cut, signalmT.GetSize()-2)
        #return math.exp(-(b+s))
        return math.exp(-s)
    expectedZs=[]
    bincuts=[]
    bincutsexp=[]
    possibleoptimalZexp=[]
    #search for local maximums and save them
    for i in range(0,totalbkgmT.GetSize()-48):#tbins=150 so size=152 and last diboson at i=104, no more background after. Leads to log(1+s/0)
        expectedZs.append(expectedsignificance(i))
        bincuts.append(totalbkgmT.GetBinCenter(i))
        if expectedsignificance(i)>expectedsignificance(i+1):
            possibleoptimalZexp.append(expectedsignificance(i))
            bincutsexp.append(i)
#print(possibleoptimalZexp)
#get the absolute maximum, optimal significance Z

    optimalZexp=max(possibleoptimalZexp)
    indice=expectedZs.index(optimalZexp)
    optimalcut=bincuts[indice]
    #for any given cut, not the optimal one
    optimalcut=masscut[a]
    indice=0
    #search from which bin we have to start integrating knowing the mass cut
    for j in range(1, len(bincuts)-1):
        if bincuts[j]<=optimalcut and bincuts[j+1]>optimalcut:
            indice=j
            print(indice)

    berror=ROOT.Double()
    serror=ROOT.Double()
    bcut=totalbkgmT.IntegralAndError(indice, totalbkgmT.GetSize()-2, berror)
    #bcut=totalbkgmT.Integral(indice, totalbkgmT.GetSize()-2)
    scut=signalmT.IntegralAndError(indice, signalmT.GetSize()-2, serror)
    zobscut=minimizeme(indice)
    n_obs=datamT.Integral(indice, datamT.GetSize()-2)
    confidencelimitcut=1-confidencelimit(indice)


    print(bcut)
    arrayfortable.append([sig.split('e')[1],round(optimalcut,2),round(bcut,3),round(berror,5),round(scut,5), n_obs,round(optimalZexp,2), round(zobscut,2), round(confidencelimitcut, 2)])
    
    scuts.append(scut)
    serrors.append(serror)
    #print("the maximum expected significance is "+ str(optimalZexp))
    plt.plot(bincuts,expectedZs)
    plt.xlabel('bincut/GeV')
    plt.ylabel('Optimal Z_exp for '+channel+' channel')
    plt.legend(("mW'=2 TeV","mW'=3 TeV","mW'=4 TeV", "mW'=5 TeV"))
    plt.savefig('Plots1lep/optimalZexp_'+channel+sig+'.png')
#plt.show()
#write it in a table

fig, ax = plt.subplots()

fig.patch.set_visible(False)
ax.axis('off')
ax.axis('tight')

df=pd.DataFrame(arrayfortable, columns=["Wmass/GeV", "Cut/GeV", "b/counts", "berror/counts", "s/counts", "n_obs/counts", "Zexp", "Zobs", "Exclusion CL"])

ax.table(cellText=df.values, colLabels=df.columns, loc='center', colWidths=[0.115, 0.115,0.115, 0.115, 0.115, 0.115, 0.115, 0.115, 0.12])

fig.tight_layout()

fig.savefig('Plots1lep/table_'+channel+".png")


##calculating signal efficiency to combine channels
#cross section: xsec
#luminosity: L
#signal values: scuts[]
if channel=='ee':
    keys=['Wprime2000enu','Wprime3000enu','Wprime4000enu','Wprime5000enu']
else:
    keys=['Wprime2000munu','Wprime3000munu','Wprime4000munu','Wprime5000munu']
efficiency=[]
efficiencyerror=[]
for i in range(0,len(scuts)):
    efficiency.append(scuts[i]/(L*1000*infos[keys[i]]['xsec']))
    efficiencyerror.append(serrors[i]/(L*1000*infos[keys[i]]['xsec']))
print("signal error: ", serrors)
print("efficiencies for "+channel+" channel:",efficiency)
print("efficiencies error for "+channel+" channel:", efficiencyerror)


