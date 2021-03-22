import ROOT as R
import sys, os, time
#para usar esto escribir en terminal python ~/Escritorio/CodeExample/runSelector.py Data o lo mismo pero terminado en MC (montecarlo en vez de datos)
#codigo ya adaptado a python3, print con parentesis, ROOT como R en vez de from root import.
t0 = time.time() 

arg1 = sys.argv[1]  

if arg1 == 'Data': 
        input_dir = 'Escritorio/1lep/Data/'
elif arg1 == 'MC': 
        input_dir = 'Escritorio/1lep/MC/'


myChain = R.TChain('mini') 


for filename in os.listdir(input_dir):
        if not '.root' in filename: continue 
        print(filename)
        myChain.Add(input_dir+filename) 
if arg1 == 'MC': 
        for filename in os.listdir(input_dir):
                if not '.root' in filename: continue 
                print(filename)  
                myChain.Add(input_dir+filename) 


entries = myChain.GetEntries() 

print("-------------------------------------------")
if arg1 == 'Data': 
        print("Running on real data!")
else: 
        print("Running on Monte Carlo!")
print("Number of events to process: %d" %entries)
print("-------------------------------------------")

if arg1 == 'Data': 
        myChain.Process("Escritorio/CodeExample/MySelector1lep.C+", "Data")
else: 
        myChain.Process("Escritorio/CodeExample/MySelector1lep.C+", "MC") 

t = int( time.time()-t0 )/60  

print("Time spent: %d min" %t)
