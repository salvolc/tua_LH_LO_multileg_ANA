import numpy as np
#from scipy.integrate import quad
from matplotlib import pyplot as plt
import os,inspect
import matplotlib.ticker as ticker

def parameters_kin(var,par):
	topmassup=800
	topmassdown=0
	bmassup=60
	bmassdown=0
	umassup=1
	umassdown=0
	wmassup=300
	wmassdown=0
	jmassup=60
	jmassdown=0

	default_range_up=800
	default_range_down=0

	default_nbins = 64
	mass_nbins = 64
	nbin = default_nbins
	lower_range=default_range_down
	upper_range=default_range_up
	up_l = 0.1
	#plt.style.use("atlas.mplstyle")
	if(par == "Photon" and var=="PT"):
		upper_range=default_range_up
		lower_range=default_range_down
		nbin=default_nbins
	if(par == "TopQuark" and var=="PT"):
		upper_range=default_range_up
		lower_range=default_range_down
		nbin=64
	if (var == "Eta" or var == "Phi"):
		lower_range=-5
		upper_range=5
		nbin = 32
		up_l = 0
	if (var == "Eta" and (par=="Photon" or par=="bJet" or par=="Jet")):
		lower_range=-2.5
		upper_range=2.5
		nbin = 32
		up_l = 0
	if(var=="M"):
		up_l=0
		nbin=mass_nbins
		if("W" in par):
			lower_range=wmassdown
			upper_range=wmassup
		if("Top" in par):
			lower_range=topmassdown
			upper_range=topmassup
		if("Jet" in par):
				lower_range=jmassdown
				upper_range=jmassup
		if("bJet" in par):
			lower_range=bmassdown
			upper_range=bmassup
		if("UQuark" in par):
			lower_range=umassdown
			upper_range=umassup
	if(var=="PT"):
		if("bJet" in par):
			upper_range=400
	return [lower_range,upper_range,nbin,up_l]


sam = ["dec","int","pro"]
var = ["PT","Eta","Phi","M"]
par = ["Photon","Jet","bJet","WBoson","TopQuark"]

PREFIX = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/"

datas = os.listdir(PREFIX+"data")


if("corr" not in os.listdir(PREFIX+"plots/")):
	os.chdir(PREFIX+"plots/")
	os.mkdir("corr")
	os.chdir(PREFIX)

ratiox= 15
ratioy= 15

bfmpd = np.array([])
bfmi = np.array([])
names = []

i=0;j=0;
for pa1 in par:
	for va1 in var:
		if(pa1 == "Photon" and va1 == "M"):
			continue

		evp = np.genfromtxt("data/pro_"+pa1+"_"+va1+".txt");evd = np.genfromtxt("data/dec_"+pa1+"_"+va1+".txt");evi = np.genfromtxt("data/int_"+pa1+"_"+va1+".txt")
		evp = evp[(np.abs(evp)>0) & (evp!=999.9)];evd = evd[(np.abs(evd)>0) & (evd!=999.9)];evi = evi[(np.abs(evi)>0) & (evi!=999.9)]
		evpd = np.concatenate((evp, evd), axis=0)

		if bfmpd.size == 0:
			bfmpd = np.concatenate((bfmpd,evpd))
		else:
			bfmpd = np.vstack((bfmpd,evpd))

		if bfmi.size == 0:
			bfmi = np.concatenate((bfmi,evi))
		else:
			bfmi = np.vstack((bfmi,evi))

		n = pa1+va1
		names.append(n)	
		i+=1			
			
		

fig = plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')

corrpd = np.corrcoef(bfmpd)

plt.imshow(corrpd)
plt.colorbar()
plt.xticks(np.arange(i),names, rotation=45)
plt.yticks(np.arange(i),names, rotation=45)

import itertools
thresh = corrpd.max() / 2.

for i, j in itertools.product(range(corrpd.shape[0]), range(corrpd.shape[1])):
	plt.text(j, i,format(corrpd[i, j], '.2f'),fontsize=13,
	horizontalalignment="center",
	color="white" if corrpd[i, j] > thresh else "black")
plt.title("Correlation Matrix Production + Decay")
plt.savefig("plots/CorrMatrixpd.pdf",bbox_inches='tight')

plt.close()
fig.clear()

plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')

corri = np.corrcoef(bfmi)

plt.imshow(corri)
plt.colorbar()
plt.xticks(np.arange(i),names, rotation=45)
plt.yticks(np.arange(i),names, rotation=45)

thresh = corri.max() / 2.

for i, j in itertools.product(range(corri.shape[0]), range(corri.shape[1])):
	plt.text(j, i,format(corri[i, j], '.2f'),fontsize=13,
	horizontalalignment="center",
	color="white" if corri[i, j] > thresh else "black")
plt.title("Correlation Matrix Interfernce")
plt.savefig("plots/CorrMatrixi.pdf",bbox_inches='tight')






