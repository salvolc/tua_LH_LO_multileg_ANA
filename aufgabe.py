import numpy as np
#from scipy.integrate import quad
from matplotlib import pyplot as plt
import os,inspect

def plot_error_region(yvalues,yerrors,bins,color='orange'):

	for i in range(len(bins)-1):
		plt.fill_between(x=[bins[i],bins[i+1]],
			y1=yvalues[i]-yerrors[i],
			y2=yvalues[i]+yerrors[i],
			alpha=0.7,
			color=color)

def plot_error_region2(yvalues,yerrors,bins,color='orange'):
	bins_a = bins[:-1]
	bins_e = bins[1:]
	bins_m = (bins_a+bins_e)/2
	plt.errorbar(bins_m,yvalues,yerr=yerrors,drawstyle = 'steps-mid',color=color)


###################################################################################################
########################################TRUTH PLOTS################################################
###################################################################################################

samples = ["dec","int","pro"]
vari = ["PT","Eta","Phi","M"]
part = ["TopQuark","Photon","BQuark","WBoson","UQuark"]
#vgle = ["LH vs RH" , "$tuγ$ vs. $tcγ$"]
#utyp = ["tua","tca"]
#hand = ["LH","RH"]
PREFIX = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/"

datas = os.listdir(PREFIX+"data")

ntruth = 0
for i in range(len(datas)):
	if("truth" in datas[i]):
		ntruth += 1

for par in part:
	if(par+"_truth" not in os.listdir(PREFIX+"plots/")):
		os.chdir(PREFIX+"plots/")
		os.mkdir(par+"_truth/")
		os.chdir(PREFIX)

topmassup=180
topmassdown=165
bmassup=5.2
bmassdown=3.2
umassup=1
umassdown=0
wmassup=95
wmassdown=65

default_range_up=500
default_range_down=0

default_nbins = 128
mass_nbins = 64

ratiox=10
ratioy=4

for sam in samples:
	for par in part:
		for var in vari:
			nbin = default_nbins
			lower_range=default_range_down
			upper_range=default_range_up

			if(par == "Photon" and var=="M"):
				continue
			if(par == "Photon" and var=="PT"):
				upper_range=default_range_up
				lower_range=default_range_down
				nbin=default_nbins
			if (var == "Eta" or var == "Phi"):
				lower_range=-5
				upper_range=5
				nbin = 32
			if(var=="M"):
				nbin=mass_nbins
				if("W" in par):
					lower_range=wmassdown
					upper_range=wmassup
				if("Top" in par):
					lower_range=topmassdown
					upper_range=topmassup
				if("BQuark" in par):
					lower_range=bmassdown
					upper_range=bmassup
				if("UQuark" in par):
					lower_range=umassdown
					upper_range=umassup

			ev = np.genfromtxt("data/"+sam+"_"+par+"_"+var+"_truth.txt")
			ev = ev[(ev!=0) & (ev!=999.9)]
			print(len(ev[(ev!=0) & (ev!=999.9)]))
			print("data/"+sam+"_"+par+"_"+var+"_truth.txt")
			plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')
			n,bins,a = plt.hist(ev,label=sam,bins=nbin,lw=0.5,color="blue",fill=False,normed=False,range=(lower_range,upper_range),histtype='step')
			plot_error_region2(n, np.sqrt(n), bins,"blue")
			plt.xlabel(r"$"+var+"("+par+")"+r"$")
			plt.ylabel("N")
			plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			plt.legend(loc="best")
			if ("dec" in sam):
				plt.title(var+"("+par+") decaymode")
				plt.savefig("plots/"+par+"_truth/"+"decaymode"+"_"+par+"_"+var+".pdf",bbox_inches='tight')
				plt.close()
			if ("int" in sam):
				plt.title(var+"("+par+") interference")
				plt.savefig("plots/"+par+"_truth/"+"interference"+"_"+par+"_"+var+".pdf",bbox_inches='tight')
				plt.close()
			if ("pro" in sam):
				plt.title(var+"("+par+") production")
				plt.savefig("plots/"+par+"_truth/"+"production"+"_"+par+"_"+var+".pdf",bbox_inches='tight')
				plt.close()
			# plt.title(var+"("+par+") "+sam[0:2]+r"$γ$"+sam[3:])
			# plt.savefig("plots/"+par+"_truth/"+sam+"_"+par+"_"+var+".pdf",bbox_inches='tight')
			# plt.close()
		



topmassup=180
topmassdown=165
bmassup=5.2
bmassdown=3.2
umassup=1
umassdown=0
wmassup=95
wmassdown=65

default_range_up=500
default_range_down=0

default_nbins = 32
mass_nbins = 64

ratiox=10
ratioy=4
for par in part:
	for var in vari:
		nbin = default_nbins
		lower_range=default_range_down
		upper_range=default_range_up

		if(par == "Photon" and var=="M"):
			continue
		if(par == "Photon" and var=="PT"):
			upper_range=default_range_up
			lower_range=default_range_down
			nbin=default_nbins
		if (var == "Eta" or var == "Phi"):
			lower_range=-5
			upper_range=5
			nbin = 32
		if(var=="M"):
			nbin=mass_nbins
			if("W" in par):
				lower_range=wmassdown
				upper_range=wmassup
			if("Top" in par):
				lower_range=topmassdown
				upper_range=topmassup
			if("BQuark" in par):
				lower_range=bmassdown
				upper_range=bmassup
			if("UQuark" in par):
				lower_range=umassdown
				upper_range=umassup

		evp = np.genfromtxt("data/pro_"+par+"_"+var+"_truth.txt")
		evd = np.genfromtxt("data/dec_"+par+"_"+var+"_truth.txt")
		evi = np.genfromtxt("data/int_"+par+"_"+var+"_truth.txt")
		evp = evp[evp!=0]
		evd = evd[evd!=0]
		evi = evi[evi!=0]
		ev = np.concatenate((evp, evd), axis=0)
		print(par+var)
		vn,vbins,va = plt.hist(ev,label=r"production+decay",bins=nbin,lw=0.5,color="blue",fill=False,normed=False,range=(lower_range,upper_range),histtype='step')
		vnI,vbinsI,vaI = plt.hist(evi,label=r"interference",bins=nbin,lw=0.5,color="red",fill=False,normed=False,range=(lower_range,upper_range),histtype='step')

		plt.close()
		plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')
		n,bins,a = plt.hist(ev,label=r"production+decay",bins=nbin,lw=0.5,color="blue",fill=False,normed=True,range=(lower_range,upper_range),histtype='step')
		plot_error_region2(n,1/np.sqrt(vn)*n, bins,"blue")
		nI,binsI,aI = plt.hist(evi,label=r"interference",bins=nbin,lw=0.5,color="red",fill=False,normed=True,range=(lower_range,upper_range),histtype='step')
		plot_error_region2(nI,1/np.sqrt(vnI)*nI, binsI,"red")
		plt.xlabel(r"$"+var+"("+par+")"+r"$")
		plt.ylabel("N")
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.legend(loc="best")
		plt.title(var+"("+par+") all")
		plt.savefig("plots/"+par+"_truth/"+"all"+"_"+par+"_"+var+".pdf",bbox_inches='tight')
		plt.close()