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

variables = ['pt1', 'pt2', 'eta1', 'eta2', 'phi1', 'phi2', 'mll', 'met', 'mllLog']

xtitles = {'pt1': 'Leading lepton p_{T} (GeV)', 'pt2': 'Subleading lepton p_{T} (GeV)', 'eta1': 'Leading lepton #eta', 'eta2': 'Subleading lepton #eta',
           'phi1': 'Leading lepton #phi', 'phi2': 'Subleading lepton #phi', 'mll': 'm_{ll} (GeV)', 'met': 'E_{T}^{miss} (GeV)', 'mllLog': 'm_{ll} (GeV)'}


# Backgrounds

backgrounds = ['Zjets', 'Top','Diboson', 'Wjets']


Zjets = [364100, 364101, 364102, 364103, 364104, 364105, 364106, 364107, 364108, 364109, 364110, 364111, 364112, 364113, 364114, 364115, 364116, 364117, 364118, 364119, 364120, 364121, 364122, 364123, 364124,
         364125, 364126, 364127, 364128, 364129, 364130, 364131, 364132, 364133, 364134, 364135, 364136, 364137, 364138, 364139, 364140, 364141]

Wjets = [364156, 364157, 364158, 364159, 364160, 364161, 364162, 364163, 364164, 364165, 364166, 364167, 364168, 364169, 364170, 364171, 364172, 364173, 364174, 364175, 364176, 364177, 364178, 364179, 364180,
         364181, 364182, 364183, 364184, 364185, 364186, 364187, 364188, 364189, 364190, 364191, 364192, 364193, 364194, 364195, 364196, 364197]

Diboson = [363356, 363358, 363359, 363360, 363489, 363490, 363491, 363492, 363493] 

Top = [410000, 410011, 410012, 4100013, 410014, 410025, 410026]

signals = ['Zprime2000', 'Zprime3000', 'Zprime4000','Zprime5000']

Zprime2000 = [301215, 301220];

Zprime3000 = [301216, 301221];

Zprime4000 = [301217, 301222];

Zprime5000 = [301218, 301223];

fileIDs = {'Diboson':Diboson, 'Zjets':Zjets, 'Wjets':Wjets, 'Top':Top, 'Zprime2000':Zprime2000, 'Zprime3000':Zprime3000, 'Zprime4000':Zprime4000, 'Zprime5000':Zprime5000}

hist_bkg = {}#dictionary of histograms of the background
hist_sig = {}#dictionary for signal(s)
for var in variables:
    hist_bkg[var] = {}#for each variable to plot(invariant mass or whatever)
    hist_sig[var] = {}
    for bkg in backgrounds:
        hist_bkg[var][bkg] = ROOT.TH1F()#each of the backgrounds contributing, declare the histogram (empty)
    for sig in signals:
        hist_sig[var][sig] = ROOT.TH1F()

colours = dict(Diboson=ROOT.kAzure+1, Top=ROOT.kRed+1, Zjets=ROOT.kOrange-2, Wjets=ROOT.kGray, Zprime2000=ROOT.kCyan, Zprime3000=ROOT.kMagenta, Zprime4000=ROOT.kYellow, Zprime5000=ROOT.kGreen)#colours to plot variables and signal


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

for filename in os.listdir('Escritorio/CodeExample/Histograms/MC/'):
    if '.root' in filename:
        filepath = 'Escritorio/CodeExample/Histograms/MC/'+filename
        infile = ROOT.TFile(filepath)
        file_id = int(filename.split('.')[2])
        # print filename
        for var in variables:
            for bkg in backgrounds:
                if file_id in fileIDs[bkg]:
                    hist_bkg[var][bkg] = fill_hist(
                        hist_bkg[var][bkg], 'h_'+channel+'_'+var, bkg, file_id)
            for sig in signals:
                if file_id in fileIDs[sig]: 
                    hist_sig[var][sig] = fill_hist(hist_sig[var][sig], 'h_'+channel+'_'+var, sig, file_id)

##################################
#Smoothing background 03/05/2020
#Replace Z+jets with a fit
gROOT.LoadMacro("Escritorio/CodeExample/fitStuff.cpp");
var = 'mllLog';
bkg = 'Zjets';
new={}
contador=0
downlimit=[300, 410]
functions=[0,1,2]
highend=hist_bkg[var][bkg].GetBinLowEdge(hist_bkg[var][bkg].GetNbinsX()+1)
for i in downlimit:
    for j in functions:
        clonehist=hist_bkg[var][bkg].Clone()
        print(j,i)
        new[contador]=ROOT.fitAndReplace(clonehist, i, 3500.0, 411.0, j, highend)#ojo esto borra el histograma y hay que clonarlo. Ajusta la funcion j de i a 3500 y sustituye el nuevo valor desde 801 hasta highend.
        print("hi")
        contador=contador+1

cnv = ROOT.TCanvas("cnv_fits", "", 500, 500)
cnv.SetTicks(1, 1)
cnv.SetLeftMargin(0.13)
pfit = ROOT.TPad("pfit", "", 0, 0.35, 1, 1)
pfit.SetLogy()
pfit.SetLogx()
pfit.Draw()
pfit.cd()
color=[ROOT.kBlue,ROOT.kBlack,ROOT.kGreen,ROOT.kRed,ROOT.kOrange,ROOT.kPink,ROOT.kMagenta,ROOT.kYellow,ROOT.kCyan]
#Create the averaged one
hist_averagefit=new[0].Clone()#clone for same bins and structure
hist_averagefit.Add(new[0],-1.0)#add negative one so everything's at 0.
for i in range(0, len(new)):
    new[i].SetFillColor(0)
    new[i].SetLineColor(color[i])
    new[i].GetYaxis().SetTitle("Events")
    new[i].GetXaxis().SetTitle("m_{ll} [GeV]")
    new[i].GetXaxis().SetTitleOffset(1.5)
    new[i].Draw("same hist");
    pfit.Update()
    hist_averagefit.Add(new[i])
hist_bkg['mllLog']['Zjets'].Draw("samePE")
pfit.Update()
hist_averagefit.Scale(1.0/len(new))#divide the histogram by the number of hist used to average
for i in range(0, hist_averagefit.GetSize()-2):
    array=[]
    for j in range(0, len(new)):
        array.append(new[j].GetBinContent(i))
    maximum=max(array)
    uperror=maximum-hist_averagefit.GetBinContent(i)
    #print(uperror)
    hist_averagefit.SetBinError(i,uperror)
cnv.cd()
# Make plot legend

legfit = ROOT.TLegend(0.70, 0.7, 0.88, 0.91)
legfit.SetFillStyle(4000)
legfit.SetFillColor(0)
legfit.SetTextFont(42)
legfit.SetBorderSize(0)
fitlabels=['f1 limit=300','f2 limit=300','f3 limit=300','f1 limit=410','f2 limit=410','f3 limit=410']
for i in range(0,len(fitlabels)):
    legfit.AddEntry(new[i],fitlabels[i],"f")
legfit.AddEntry(hist_bkg['mllLog']['Zjets'], "Original Z+jets","ple")

legfit.Draw()
cnv.cd
cnv.Update()
cnv.Print('Plots/fitting'+channel+'.png')
cnv.Close()

#Replace the background histogram that goes then into the stacked plot loop
hist_bkg[var][bkg]=hist_averagefit.Clone()



#hist_bkg[var][bkg] = ROOT.fitAndReplace(hist_bkg[var][bkg], 400.0, 3500.0, 600.0, 0, hist_bkg[var][bkg].GetBinLowEdge(hist_bkg[var][bkg].GetNbinsX()+1));

var='mll'
bkg = 'Zjets';
new2={}
contador=0
downlimit=[]
if channel=="uu":
    downlimit=[310, 410]
else:
    downlimit=[300,340]#fit crashes in some points and they are different for ee or uu
functions=[0,1,2]
highend=hist_bkg[var][bkg].GetBinLowEdge(hist_bkg[var][bkg].GetNbinsX()+1)
for i in downlimit:
    for j in functions:
        clonehist=hist_bkg[var][bkg].Clone()
        print(j,i)
        new2[contador]=ROOT.fitAndReplace(clonehist, i, 3500.0, 801.0, j, highend)#ojo esto borra el histograma y hay que clonarlo. Ajusta la funcion j de i a 3500 y sustituye el nuevo valor desde 801 hasta highend.
        print("hi")
        contador=contador+1

#Create the averaged one
hist_averagefit=new2[0].Clone()#clone for same bins and structure
hist_averagefit.Add(new2[0],-1.0)#add negative one so everything's at 0.
for i in range(0, len(new2)):
    hist_averagefit.Add(new2[i])
print(len(new2))
hist_averagefit.Scale(1/len(new2))#divide the histogram by the number of hist used to average
#Replace the background histogram that goes then into the stacked plot loop
for i in range(0, hist_averagefit.GetSize()-2):
    array=[]
    for j in range(0, len(new2)):
        array.append(new2[j].GetBinContent(i))
    #print(array)
    maximum=max(array)
    uperror=maximum-hist_averagefit.GetBinContent(i)
    #print(uperror, maximum)
    hist_averagefit.SetBinError(i,uperror)
hist_bkg[var][bkg]=hist_averagefit.Clone()



#hist_bkg[var][bkg] = ROOT.fitAndReplace(hist_bkg[var][bkg], 400.0, 3500.0, 600.0, 0, hist_bkg[var][bkg].GetBinLowEdge(hist_bkg[var][bkg].GetNbinsX()+1));
#print(hist_bkg[var][bkg].GetBinContent(3))



##################################

# Get data

data = ROOT.TFile('Escritorio/CodeExample/Histograms/Data/hist.Data.2016.root')
hist_d = {}

for var in variables:
    hist_d[var] = data.Get('h_'+channel+'_'+var)
    hist_d[var].SetMarkerStyle(20)
    hist_d[var].SetMarkerSize(0.7)
    hist_d[var].SetLineColor(ROOT.kBlack)
    hist_d[var].GetYaxis().SetTitle("Events")
    hist_d[var].GetXaxis().SetTitle(xtitles[var])
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
hist_r = {}
hist_mc = {}

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
        if not hist_mc[var].GetName():
            hist_mc[var] = hist_bkg[var][bkg].Clone()
        else:
            print("Here")
            hist_mc[var].Add(hist_bkg[var][bkg])
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




# Make plot legend

leg = ROOT.TLegend(0.70, 0.50, 0.88, 0.88)
leg.SetFillStyle(4000)
leg.SetFillColor(0)
leg.SetTextFont(42)
leg.SetBorderSize(0)

bkg_labels = {'Zjets': 'Z+jets', 'Top': 'Top',
              'Diboson': 'Diboson', 'Wjets': 'W+jets'}

sig_labels = {'Zprime2000':"Z' (2 TeV)", 'Zprime3000':"Z' (3 TeV)", 'Zprime4000':"Z' (4 TeV)", 'Zprime5000':"Z' (5 TeV)" }

for bkg in backgrounds:
    leg.AddEntry(hist_bkg['pt1'][bkg], bkg_labels[bkg], "f")

for sig in signals: 
    leg.AddEntry(hist_sig['pt1'][sig], sig_labels[sig], "f")

leg.AddEntry(hist_d['pt1'], "Data", "ple")

selection = ""
if channel == "ee":
    selection = "ee"
if channel == "uu":
    selection = "#mu#mu"

# Make plots

for var in variables:

    cnv = ROOT.TCanvas("cnv_"+var, "", 500, 500)
    cnv.SetTicks(1, 1)
    cnv.SetLeftMargin(0.13)

    p1 = ROOT.TPad("p1", "", 0, 0.35, 1, 1)
    p2 = ROOT.TPad("p2", "", 0, 0.0, 1, 0.35)

    p1.SetLogy()
    if var=='mllLog':#logarithmic x axis plot
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
    if var in ['eta1', 'eta2', 'phi1', 'phi2']:
        maximum = stack[var].GetMaximum()
        stack[var].SetMaximum(maximum*10E4)

    hist_d[var].Draw("same e0")
    leg.Draw("same")

    for sig in signals:
        hist_sig[var][sig].SetFillColor(0);
        hist_sig[var][sig].Draw("same hist");

    s = ROOT.TLatex()
    s.SetNDC(1)
    s.SetTextAlign(13)
    s.SetTextColor(ROOT.kBlack)
    s.SetTextSize(0.044)
    s.DrawLatex(0.4, 0.86, "#font[72]{ATLAS} Open Data")
    s.DrawLatex(0.4, 0.81, "#bf{#sqrt{s} = 13 TeV,^{}%.1f^{}fb^{-1}}" % (L))
    s.DrawLatex(0.4, 0.76, "#bf{"+selection+" selection}")

    p1.Update()
    p1.RedrawAxis()

    cnv.cd()

    p2.Draw()
    p2.cd()

    p2.SetGridy()
    if var=='mllLog':
        p2.SetLogx()
    hist_r[var].SetMaximum(1.99)
    hist_r[var].SetMinimum(0.01)
    hist_r[var].Draw("0PZ")

    p2.SetTopMargin(0)
    p2.SetBottomMargin(0.35)
    p2.Update()
    p2.RedrawAxis()

    cnv.cd()
    cnv.Update()
    cnv.Print('Plots/'+channel+'_'+var+'.png')
    cnv.Close()

####Integration of the histogram. Need to integrate the background, the data and the signal for the invariant mass plot (mll).

totalbkgmll=stack['mll'].GetStack().Last()
datamll=hist_d['mll']
signalmll=hist_sig['mll']['Zprime2000']

#find optimal bin cut for which Z, significance, is maximum. The Integral method takes initial bin, last bin, so receives integer as a number. cant minimize with a normal minimizing method that uses double. for length of the histogram use GetSize(). If 100 bins, returns 102(100 bins+underflow+overflow). Size of data is 52
#print(datamll.Integral(0, datamll.GetSize()-2))

def minimizeme(cut):
    n_obs=datamll.Integral(cut, datamll.GetSize()-2)#also use integralanderror
    b=totalbkgmll.Integral(cut, totalbkgmll.GetSize()-2)
    if n_obs==0:#limit of xlog(x) when x-->0 is 0
        Zobs=math.sqrt(2*b)
    else:
        Zobs=math.sqrt(2*(n_obs*math.log(n_obs/b)-n_obs+b))
    return Zobs

def ZwithoutWpeak(cut):
    n_obs=datamll.Integral(cut, datamll.GetSize()-2)#also use integralanderror
    b=totalbkgmll.Integral(cut, totalbkgmll.GetSize()-49)
    if n_obs==0:#limit of xlog(x) when x-->0 is 0
        Zobswithoutpeak=math.sqrt(2*b)#but this is nonsense, it is not defined when b=0, this is just for it not to launch an error
    else:
        Zobswithoutpeak=math.sqrt(2*(n_obs*math.log(n_obs/b)-n_obs+b))
        return Zobswithoutpeak
#print(minimizeme(10))
bincuts=[]
possibleoptimalZ=[]
#search for local maximums and save them
for i in range(0,datamll.GetSize()-110):#dont put the cut too far, no data last bins, change limits in myselector, need to rerun (aprox 1hour)
    if minimizeme(i)>minimizeme(i+1):
        possibleoptimalZ.append(minimizeme(i))
        bincuts.append(i)
#get the absolute maximum, optimal significance Z
optimalZ=max(possibleoptimalZ)
#print(possibleoptimalZ)
#print(datamll.GetBin(bincuts[5]))
#print("the maximum observed significance is "+ str(optimalZ))

#expected significance
masscut=[1600, 2400, 3200, 4000]

a=-1
arrayfortable=[]

scuts=[]#to be used in the signal efficiency part
serrors=[]

for sig in signals:
    a=a+1
    signalmll=hist_sig['mll'][sig]
    def expectedsignificance(cut):
        b=totalbkgmll.Integral(cut, totalbkgmll.GetSize()-2)
        s=signalmll.Integral(cut, signalmll.GetSize()-2)
        if b==0:
            return 0
        else:
            Zexp=math.sqrt(2*((s+b)*math.log(1+s/b)-s))
            return Zexp

    def confidencelimit(cut):
        b=totalbkgmll.Integral(cut, totalbkgmll.GetSize()-2)
        s=signalmll.Integral(cut, signalmll.GetSize()-2)
        #return math.exp(-(b+s))
        return math.exp(-s)
    expectedZs=[]
    bincuts=[]
    bincutsexp=[]
    possibleoptimalZexp=[]
    #search for local maximums and save them
    for i in range(0,totalbkgmll.GetSize()-48):#tbins=150 so size=152 and last diboson at i=104, no more background after. Leads to log(1+s/0)
        expectedZs.append(expectedsignificance(i))
        bincuts.append(totalbkgmll.GetBinCenter(i))
        if expectedsignificance(i)>expectedsignificance(i+1):
            possibleoptimalZexp.append(expectedsignificance(i))
            bincutsexp.append(i)
#print(possibleoptimalZexp)
#get the absolute maximum, optimal significance Z

    optimalZexp=max(possibleoptimalZexp)
    indice=expectedZs.index(optimalZexp)
    optimalcut=bincuts[indice]#table shows optimalZexp, no the one corresponding to the arbitrary 80% mass cut
    #for any given cut, not the optimal one
    optimalcut=masscut[a]
    indice=0
    for j in range(1, len(bincuts)-1):
        if bincuts[j]<=optimalcut and bincuts[j+1]>optimalcut:
            indice=j
            print(indice)
    berror=ROOT.Double()
    serror=ROOT.Double()#to be used in the signal efficiency part
    bcut=totalbkgmll.IntegralAndError(indice, totalbkgmll.GetSize()-2, berror)
    #bcut=totalbkgmll.Integral(indice, totalbkgmll.GetSize()-2)
    scut=signalmll.IntegralAndError(indice, signalmll.GetSize()-2, serror)
    zobscut=minimizeme(indice)
    n_obs=datamll.Integral(indice, datamll.GetSize()-2)
    confidencelimitcut=1-confidencelimit(indice)


    print(bcut)
    arrayfortable.append([sig.split('e')[1],round(optimalcut,2),round(bcut,9),round(berror,9),round(scut,5), n_obs,round(optimalZexp,2), round(zobscut,2), round(confidencelimitcut, 2)])
    
    scuts.append(scut)
    serrors.append(serror)
    #print("the maximum expected significance is "+ str(optimalZexp))
    plt.plot(bincuts,expectedZs)
    plt.xlabel('bincut/GeV')
    plt.ylabel('Optimal Z_exp for '+channel+' channel')
    plt.legend(("mZ'=2 TeV","mZ'=3 TeV","mZ'=4 TeV", "mZ'=5 TeV"))
    plt.savefig('Plots/optimalZexp_'+channel+sig+'.png')
#plt.show()
#write it in a table

fig, ax = plt.subplots()

fig.patch.set_visible(False)
ax.axis('off')
ax.axis('tight')

df=pd.DataFrame(arrayfortable, columns=["Zmass/GeV", "Cut/GeV", "b/counts", "berror/counts", "s/counts", "n_obs/counts", "Zexp", "Zobs", "Exclusion CL"])

ax.table(cellText=df.values, colLabels=df.columns, loc='center', colWidths=[0.115, 0.115,0.115, 0.115, 0.115, 0.115, 0.115, 0.115, 0.12])

fig.tight_layout()

fig.savefig('Plots/table_'+channel+".png")


##calculating signal efficiency to combine channels
#cross section: xsec
#luminosity: L
#signal values: scuts[]
#sigma= 1000*info[ID]['xsec']
if channel=='ee':
    keys=['Zprime2000ee','Zprime3000ee','Zprime4000ee','Zprime5000ee']
else:
    keys=['Zprime2000mumu','Zprime3000mumu','Zprime4000mumu','Zprime5000mumu']
efficiency=[]
efficiencyerror=[]
for i in range(0,len(scuts)):
    efficiency.append(scuts[i]/(L*1000*infos[keys[i]]['xsec']))
    efficiencyerror.append(serrors[i]/(L*1000*infos[keys[i]]['xsec']))
print("signal error: ", serrors)
print("efficiencies for "+channel+" channel:",efficiency)
print("efficiencies error for "+channel+" channel:", efficiencyerror)




