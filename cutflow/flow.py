import os
import numpy as np 
from matplotlib import pyplot as plt



paths =  ["./data/decay/",
          "./data/production/",
          "./data/interference/"]

names = [
"all",
"Photon",
"Lepton",
"b tag",
"Kin. cuts",
"W reco.",
"Matching"]

for p in paths:
  files = os.listdir(p)
  effis = np.zeros(8)
  for data in files:
    print(data)
    f = open(p+data,"r")
    line = f.readline()
    lsp = line.split()
    effis[0] = 1
    #print(lsp)
    for i in range(1,8):
      line = f.readline()
      lsp = line.split()
      #print(lsp)
      effis[i] += float((lsp[-1])[:-1])/100/len(files)

  pspl = p.split("/")[-2]
  
  xloc = [0,1,2,3,4,5,6,7]
  #plt.hlines(effis[:-1],xloc[:-1],xloc[1:])
  effis = np.insert(effis,0,0)
  print(effis)
  plt.errorbar(xloc, effis[:-1],drawstyle='steps')
  plt.xlim(0,7)
  plt.ylim(0.01,1.1)
  plt.yscale("log")
  plt.ylabel("Efficiency")
  plt.xticks((xloc[:-1]+np.diff(xloc)-0.5),names)
  plt.axes().set_xticks((xloc[:-1]+np.diff(xloc)),minor=True)
  #plt.xticks((xloc[:-1]+np.diff(xloc)),minor=True)
  plt.grid(alpha=0.5,which='both',axis='y')
  plt.grid(alpha=0.5,which='minor',axis='x')
  plt.savefig(pspl+"_eff.pdf")
  plt.clf()