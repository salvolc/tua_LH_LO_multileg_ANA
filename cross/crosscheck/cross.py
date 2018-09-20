import numpy as np 
import os
from uncertainties import ufloat
import uncertainties
from uncertainties import unumpy as unp 

names = ["decay","interference","production"]

files = os.listdir("./data/")

cross = [[],[],[]]

for i in range(0,3):
  print(i)
  for f in files:
    if names[i] in f:
      data = open("./data/"+f,"r")
      line = data.readline()
      lsp = line.split()
      print(f)
      print(lsp)
      cross[i].append(ufloat(float(lsp[0]),float(lsp[1])))

  print(cross[i])


for i in range(0,3):
  print("Crossection für "+str(names[i])+" mode ist:")
  #print(str(np.mean(cross[i]))+" +/- "+str(np.std(cross[i])))
  print(str(sum(cross[i])/len(cross[i])))