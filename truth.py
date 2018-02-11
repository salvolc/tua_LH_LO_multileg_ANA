import numpy as np
#from scipy.integrate import quad
from matplotlib import pyplot as plt
import os,inspect
import matplotlib.ticker as ticker

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
	plt.errorbar(bins_m,yvalues,yerr=yerrors,drawstyle = 'steps-mid',color=color,lw=0.5)


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

default_range_up=800
default_range_down=0

default_nbins = 128
mass_nbins = 64

ratiox=10
ratioy=4

ratiox= 8
ratioy= 6


for sam in samples:
	for par in part:
		for var in vari:
			nbin = default_nbins
			lower_range=default_range_down
			upper_range=default_range_up
			up_l = 10
			#plt.style.use("atlas.mplstyle")
			if(par == "Photon" and var=="M"):
				continue
			if(par == "Photon" and var=="PT"):
				upper_range=default_range_up
				lower_range=default_range_down
				nbin=default_nbins
			if((sam=="pro") | (sam=="int")):
				upper_range=800
			if (var == "Eta" or var == "Phi"):
				lower_range=-5
				upper_range=5
				nbin = 32
				up_l = 0
			if(var=="M"):
				nbin=mass_nbins
				up_l = 0
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
			if(var=="PT"):
				if("BQuark" in par):
					upper_range=500

			ev = np.genfromtxt("data/"+sam+"_"+par+"_"+var+"_truth.txt")
			ev = ev[(np.abs(ev)>up_l) & (ev!=999.9)]
			print(len(ev[(np.abs(ev)>up_l) & (ev!=999.9)]))
			print("data/"+sam+"_"+par+"_"+var+"_truth.txt")
			fig = plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')

			binning = np.arange(lower_range,upper_range+0.001,(upper_range-lower_range)/nbin)
			binning = np.array([-1000,*binning,1000])

			n,bins,a = plt.hist(ev,label=sam,bins=binning,lw=0.5,color="blue",fill=False,normed=False,range=(lower_range,upper_range),histtype='step')
			plt.xlim(lower_range,upper_range)
			plot_error_region2(n, np.sqrt(n), bins,"blue")

			ax = plt.gca()
			ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
			plt.grid(alpha=0.5)
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


topmassup=180
topmassdown=165
bmassup=5.2
bmassdown=3.2
umassup=1
umassdown=0
wmassup=95
wmassdown=65

default_range_up=800
default_range_down=0

default_nbins = 32
mass_nbins = 64

ratiox=8
ratioy=6

up_l=0.01

for par in part:
	for var in vari:
		nbin = default_nbins
		lower_range=default_range_down
		upper_range=default_range_up
		up_l = 10

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
			if("BQuark" in par):
				lower_range=bmassdown
				upper_range=bmassup
			if("UQuark" in par):
				lower_range=umassdown
				upper_range=umassup
		if(var=="PT"):
			if("BQuark" in par):
				upper_range=500


		evp = np.genfromtxt("data/pro_"+par+"_"+var+"_truth.txt")
		evd = np.genfromtxt("data/dec_"+par+"_"+var+"_truth.txt")
		evi = np.genfromtxt("data/int_"+par+"_"+var+"_truth.txt")
		evp = evp[(np.abs(evp)>up_l) & (evp!=999.9)]
		evd = evd[(np.abs(evd)>up_l) & (evd!=999.9)]
		evi = evi[(np.abs(evi)>up_l) & (evi!=999.9)]
		ev = np.concatenate((evp, evd), axis=0)
		#print(par+var)

		binning = np.arange(lower_range,upper_range+0.001,(upper_range-lower_range)/nbin)
		binning = np.array([-1000,*binning,1000])
		vn,vbins,va = plt.hist(ev,label=r"production+decay",bins=binning,lw=0.5,color="blue",fill=False,normed=False,range=(lower_range,upper_range),histtype='step')
		vnI,vbinsI,vaI = plt.hist(evi,label=r"interference",bins=binning,lw=0.5,color="red",fill=False,normed=False,range=(lower_range,upper_range),histtype='step')
		plt.close()

		plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')
		
		n,bins,a = plt.hist(ev,label=r"production+decay",bins=binning,lw=0.5,color="blue",fill=False,normed=True,range=(lower_range,upper_range),histtype='step')
		plot_error_region2(n,1/np.sqrt(vn)*n, bins,"blue")
		nI,binsI,aI = plt.hist(evi,label=r"interference",bins=binning,lw=0.5,color="red",fill=False,normed=True,range=(lower_range,upper_range),histtype='step')
		plot_error_region2(nI,1/np.sqrt(vnI)*nI, binsI,"red")
		
		ax = plt.gca()
		ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
		plt.grid(alpha=0.5)
		plt.xlim(lower_range,upper_range+0.0001*upper_range)
		plt.xlabel(r"$"+var+"("+par+")"+r"$")
		plt.ylabel("N")
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.legend(loc="best")
		plt.title(var+"("+par+") all")
		plt.savefig("plots/"+par+"_truth/"+"all"+"_"+par+"_"+var+".pdf",bbox_inches='tight')
		plt.close()




if("R_truth" not in os.listdir(PREFIX+"plots/")):
	os.chdir(PREFIX+"plots/")
	os.mkdir("R_truth/")
	os.chdir(PREFIX)

if("M_truth" not in os.listdir(PREFIX+"plots/")):
	os.chdir(PREFIX+"plots/")
	os.mkdir("M_truth/")
	os.chdir(PREFIX)

		
var=["R","M"]
RMPart=["Photon","TopQuark"]
RMPartP=["Photon","TopQuark","BQuark","WBoson"]

up_l=0.01

for sam in samples:
	for va in var:
		for p1 in RMPart:
			for p2 in RMPartP:
				if(p1==p2):
					continue

				if va == "R":
					lower_range = 0
					upper_range = 7
					nbin = 32
				if va == "M":
					lower_range = 0
					if p1 =="Photon":
						lower_range = 0
						upper_range = 1400
						if p2 == "TopQuark":
							lower_range = 150
							upper_range = 800
					if p1 =="TopQuark":
						lower_range = 150
						upper_range = 500
						if p2 == "WBoson":
							lower_range = 200
						if p2 == "Photon":
							lower_range = 150
							upper_range = 800
					nbin = 64

				ev = np.genfromtxt("data/"+sam+"_"+p1+"_"+p2+"_"+va+"_truth.txt")
				ev = ev[(np.abs(ev)>up_l) & (ev!=999.9)]

				print(len(ev[(np.abs(ev)>up_l) & (ev!=999.9)]))
				print("data/"+sam+"_"+p1+"_"+p2+"_"+va+"_truth.txt")

				fig = plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')

				binning = np.arange(lower_range,upper_range+0.001,(upper_range-lower_range)/nbin)
				binning = np.array([-10000,*binning,10000])

				n,bins,a = plt.hist(ev,label=sam,lw=0.5,color="blue",bins=binning,fill=False,normed=False,histtype='step',range=(lower_range,upper_range))
				plot_error_region2(n, np.sqrt(n), bins,"blue")

				ax = plt.gca()
				ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
				plt.xlim(lower_range,upper_range+0.01*upper_range)
				plt.grid(alpha=0.5)
				plt.xlabel(r"$"+va+"("+p1+" "+p2+")"+r"$")
				plt.ylabel("N")
				plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
				plt.legend(loc="best")
				if ("dec" in sam):
					plt.title(va+"("+p1+"_"+p2+") decaymode")
					plt.savefig("plots/"+va+"_truth/"+"decaymode"+"_"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')
					plt.close()
				if ("int" in sam):
					plt.title(va+"("+p1+"_"+p2+") interference")
					plt.savefig("plots/"+va+"_truth/"+"interference"+"_"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')
					plt.close()
				if ("pro" in sam):
					plt.title(va+"("+p1+"_"+p2+") production")
					plt.savefig("plots/"+va+"_truth/"+"production"+"_"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')
					plt.close()





for va in var:
	for p1 in RMPart:
		for p2 in RMPartP:
			if(p1==p2):
				continue

			if va == "R":
				lower_range = 0
				upper_range = 7
				nbin = 32
			if va == "M":
				lower_range = 0
				if p1 =="Photon":
					lower_range = 0
					upper_range = 1400
					if p2 == "TopQuark":
						lower_range = 150
						upper_range = 800
				if p1 =="TopQuark":
					lower_range = 150
					upper_range = 500
					if p2 == "WBoson":
						lower_range = 200
					if p2 == "Photon":
						lower_range = 150
						upper_range = 800
				nbin = 64

			evp = np.genfromtxt("data/pro_"+p1+"_"+p2+"_"+va+"_truth.txt")
			evd = np.genfromtxt("data/dec_"+p1+"_"+p2+"_"+va+"_truth.txt")
			evi = np.genfromtxt("data/int_"+p1+"_"+p2+"_"+va+"_truth.txt")
			evp = evp[(np.abs(evp)>up_l) & (evp!=999.9)]
			evd = evd[(np.abs(evd)>up_l) & (evd!=999.9)]
			evi = evi[(np.abs(evi)>up_l) & (evi!=999.9)]
			ev = np.concatenate((evp, evd), axis=0)

			fig = plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')

			binning = np.arange(lower_range,upper_range+0.1,(upper_range-lower_range)/nbin)
			binning = np.array([-10000,*binning,10000])

			vn2,vbins2,va2 = plt.hist(ev,label=r"production+decay",bins=binning,lw=0.5,color="blue",fill=False,normed=False,histtype='step')
			#vn,vbins,va = plt.hist(ev,label=r"production+decay",bins=binning,lw=0.5,color="blue",fill=False,normed=False,histtype='step')#,range=(lower_range,upper_range))
			vnI2,vbinsI2,vaI2 = plt.hist(evi,label=r"interference",bins=binning,lw=0.5,color="red",fill=False,normed=False,histtype='step')#,range=(lower_range,upper_range))
			plt.close()

			plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')
			
			n,bins,a = plt.hist(ev,label=r"production+decay",bins=binning,lw=0.5,color="blue",fill=False,normed=True,histtype='step')#,range=(lower_range,upper_range))
			plot_error_region2(n,1/np.sqrt(vn2)*n, bins,"blue")
			nI,binsI,aI = plt.hist(evi,label=r"interference",bins=binning,lw=0.5,color="red",fill=False,normed=True,histtype='step')#,range=(lower_range,upper_range))
			plot_error_region2(nI,1/np.sqrt(vnI2 )*nI, binsI,"red")

			ax = plt.gca()
			ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
			plt.xlim(lower_range,upper_range+0.01*upper_range)
			plt.legend(loc='best')
			plt.grid(alpha=0.5)
			plt.xlabel(r"$"+va+"("+p1+" "+p2+")"+r"$")
			plt.ylabel("N")
			plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			plt.title(va+"("+p1+"_"+p2+") interference")
			plt.savefig("plots/"+va+"_truth/"+"all_"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')
			plt.close()