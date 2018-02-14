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

def parameters_kin(sam,var,par):
	topmassup=250
	topmassdown=110
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
	if((sam=="pro") | (sam=="int")):
		upper_range=800
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



def parameters_RM(p1,p2,va):
	up_l=0.001
	if va == "R":
		lower_range = 0
		upper_range = 5
		nbin = 32
		if p1 == "TopQuark" and p2 == "WBoson":
			lower_range = 0
			upper_range = 4
	
	if va == "M":
		lower_range = 0
		if p1 =="Photon":
			lower_range = 0
			upper_range = 1400
			if p2 == "TopQuark":
				lower_range = 0
				upper_range = 1400
			if p2 == "LeadingJet":
				lower_range = 0
				upper_range = 400
			if p2 == "bJet":
				lower_range = 0
				upper_range = 400
		if p1 =="TopQuark":
			lower_range = 0
			upper_range = 800
			if p2 == "WBoson":
				lower_range = 50
			if p2 == "Photon":
				lower_range = 0
				upper_range = 1400
			if p2 == "LeadingJet":
				upper_range = 1200
		nbin = 64
	return [lower_range,upper_range,nbin,up_l]


###################################################################################################
########################################DELPH PLOTS################################################
###################################################################################################

samples = ["dec","int","pro"]
vari = ["PT","Eta","Phi","M"]
part = ["TopQuark","Photon","bJet","Jet","WBoson"]
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
	if(par not in os.listdir(PREFIX+"plots/")):
		os.chdir(PREFIX+"plots/")
		os.mkdir(par)
		os.chdir(PREFIX)

ratiox= 8
ratioy= 6


for sam in samples:
	for par in part:
		for var in vari:
			if(par == "Photon" and var=="M"):
				continue
			para = parameters_kin(sam, var, par)
			lower_range=para[0]
			upper_range=para[1]
			nbin=para[2]
			up_l=para[3]

			ev = np.genfromtxt("data/"+sam+"_"+par+"_"+var+".txt")
			ev = ev[(np.abs(ev)>up_l) & (ev!=999.9)]
			ev = np.clip(ev, lower_range, upper_range)
			
			print(len(ev[(np.abs(ev)>up_l) & (ev!=999.9)]))
			print("data/"+sam+"_"+par+"_"+var+".txt")
			
			fig = plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')

			n,bins,a = plt.hist(ev,label=sam,bins=nbin,lw=0.5,color="blue",fill=False,normed=False,range=(lower_range,upper_range),histtype='step')
			plot_error_region2(n, np.sqrt(n), bins,"blue")

			plt.xlim(lower_range,upper_range)
			ax = plt.gca()
			ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
			plt.grid(alpha=0.5)
			plt.xlabel(r"$"+var+"("+par+")"+r"$")
			plt.ylabel("N")
			plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			plt.legend(loc="best")
			if ("dec" in sam):
				plt.title(var+"("+par+") decaymode")
				plt.savefig("plots/"+par+"/"+"decaymode"+"_"+par+"_"+var+".pdf",bbox_inches='tight')
				plt.close()
			if ("int" in sam):
				plt.title(var+"("+par+") interference")
				plt.savefig("plots/"+par+"/"+"interference"+"_"+par+"_"+var+".pdf",bbox_inches='tight')
				plt.close()
			if ("pro" in sam):
				plt.title(var+"("+par+") production")
				plt.savefig("plots/"+par+"/"+"production"+"_"+par+"_"+var+".pdf",bbox_inches='tight')
				plt.close()


for par in part:
	for var in vari:
		if(par == "Photon" and var=="M"):
			continue
		para = parameters_kin(sam, var, par)
		lower_range=para[0]
		upper_range=para[1]
		nbin=para[2]
		up_l=para[3]
		evp = np.genfromtxt("data/pro_"+par+"_"+var+".txt")
		evd = np.genfromtxt("data/dec_"+par+"_"+var+".txt")
		evi = np.genfromtxt("data/int_"+par+"_"+var+".txt")
		evp = evp[(np.abs(evp)>up_l) & (evp!=999.9)]
		evd = evd[(np.abs(evd)>up_l) & (evd!=999.9)]
		evi = evi[(np.abs(evi)>up_l) & (evi!=999.9)]
		ev = np.concatenate((evp, evd), axis=0)
		#print(par+var)
		ev = np.clip(ev, lower_range, upper_range)
		evi = np.clip(evi, lower_range, upper_range)


		vn,vbins,va = plt.hist(ev,label=r"production+decay",bins=nbin,lw=0.5,color="blue",fill=False,normed=False,range=(lower_range,upper_range),histtype='step')
		vnI,vbinsI,vaI = plt.hist(evi,label=r"interference",bins=vbins,lw=0.5,color="red",fill=False,normed=False,range=(lower_range,upper_range),histtype='step')
		plt.close()

		plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')
		n,bins,a = plt.hist(ev,label=r"production+decay",bins=nbin,lw=0.5,color="blue",fill=False,weights=np.ones_like(ev)/float(len(ev)),range=(lower_range,upper_range),histtype='step')
		plot_error_region2(n,1/np.sqrt(vn)*n, bins,"blue")
		nI,binsI,aI = plt.hist(evi,label=r"interference",bins=nbin,lw=0.5,color="red",fill=False,weights=np.ones_like(evi)/float(len(evi)),range=(lower_range,upper_range),histtype='step')
		plot_error_region2(nI,1/np.sqrt(vnI)*nI, binsI,"red")
		
		ax = plt.gca()
		ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
		plt.grid(alpha=0.5)
		plt.xlim(lower_range,upper_range)
		plt.xlabel(r"$"+var+"("+par+")"+r"$")
		plt.ylabel("N")
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.legend(loc="best")
		plt.title(var+"("+par+") all")
		plt.savefig("plots/"+par+"/"+"all"+"_"+par+"_"+var+".pdf",bbox_inches='tight')
		plt.close()


if("R" not in os.listdir(PREFIX+"plots/")):
	os.chdir(PREFIX+"plots/")
	os.mkdir("R/")
	os.chdir(PREFIX)

if("M" not in os.listdir(PREFIX+"plots/")):
	os.chdir(PREFIX+"plots/")
	os.mkdir("M")
	os.chdir(PREFIX)

		
var=["R","M"]
RMPart=["Photon","TopQuark"]
RMPartP=["Photon","TopQuark","bJet","WBoson","LeadingJet"]



for sam in samples:
	for va in var:
		for p1 in RMPart:
			for p2 in RMPartP:
				if(p1==p2):
					continue
				para = parameters_RM(p1, p2, va)
				lower_range=para[0]
				upper_range=para[1]
				nbin=para[2]
				up_l=para[3]

				ev = np.genfromtxt("data/"+sam+"_"+p1+"_"+p2+"_"+va+".txt")
				ev = ev[(np.abs(ev)>up_l) & (ev!=999.9)]

				print(len(ev[(np.abs(ev)>up_l) & (ev!=999.9)]))
				print("data/"+sam+"_"+p1+"_"+p2+"_"+va+".txt")

				fig = plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')

				ev = np.clip(ev, lower_range, upper_range)
				n,bins,a = plt.hist(ev,label=sam,lw=0.5,color="blue",bins=nbin,fill=False,normed=False,histtype='step',range=(lower_range,upper_range))
				plot_error_region2(n, np.sqrt(n), bins,"blue")

				ax = plt.gca()
				ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
				plt.xlim(lower_range,upper_range)
				plt.grid(alpha=0.5)
				plt.xlabel(r"$"+va+"("+p1+" "+p2+")"+r"$")
				plt.ylabel("N")
				plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
				plt.legend(loc="best")
				if ("dec" in sam):
					plt.title(va+"("+p1+"_"+p2+") decaymode")
					plt.savefig("plots/"+va+"/"+"decaymode"+"_"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')
					plt.close()
				if ("int" in sam):
					plt.title(va+"("+p1+"_"+p2+") interference")
					plt.savefig("plots/"+va+"/"+"interference"+"_"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')
					plt.close()
				if ("pro" in sam):
					plt.title(va+"("+p1+"_"+p2+") production")
					plt.savefig("plots/"+va+"/"+"production"+"_"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')
					plt.close()


for va in var:
	for p1 in RMPart:
		for p2 in RMPartP:
			if(p1==p2):
				continue

			para = parameters_RM(p1, p2, va)
			lower_range=para[0]
			upper_range=para[1]
			nbin=para[2]
			up_l=para[3]
			evp = np.genfromtxt("data/pro_"+p1+"_"+p2+"_"+va+".txt")
			evd = np.genfromtxt("data/dec_"+p1+"_"+p2+"_"+va+".txt")
			evi = np.genfromtxt("data/int_"+p1+"_"+p2+"_"+va+".txt")
			evp = evp[(np.abs(evp)>up_l) & (evp!=999.9)]
			evd = evd[(np.abs(evd)>up_l) & (evd!=999.9)]
			evi = evi[(np.abs(evi)>up_l) & (evi!=999.9)]
			ev = np.concatenate((evp, evd), axis=0)

			fig = plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')

			ev = np.clip(ev, lower_range, upper_range)
			evi = np.clip(evi, lower_range, upper_range)


			vn2,vbins2,va2 = plt.hist(ev,label=r"production+decay",bins=nbin,lw=0.5,color="blue",fill=False,normed=False,histtype='step',range=(lower_range,upper_range))
			vnI2,vbinsI2,vaI2 = plt.hist(evi,label=r"interference",bins=vbins2,lw=0.5,color="red",fill=False,normed=False,histtype='step',range=(lower_range,upper_range))
			plt.close()

			plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')
			
			n,bins,a = plt.hist(ev,label=r"production+decay",bins=nbin,lw=0.5,color="blue",fill=False,weights=np.ones_like(ev)/float(len(ev)),histtype='step')#,range=(lower_range,upper_range))
			plot_error_region2(n,1/np.sqrt(vn2)*n, bins,"blue")
			nI,binsI,aI = plt.hist(evi,label=r"interference",bins=bins,lw=0.5,color="red",fill=False,weights=np.ones_like(evi)/float(len(evi)),histtype='step')#,range=(lower_range,upper_range))
			plot_error_region2(nI,1/np.sqrt(vnI2)*nI, binsI,"red")

			ax = plt.gca()
			ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
			plt.xlim(lower_range,upper_range)

			plt.legend(loc='best')
			plt.grid(alpha=0.5)
			plt.xlabel(r"$"+va+"("+p1+" "+p2+")"+r"$")
			plt.ylabel("N")
			plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			plt.title(va+"("+p1+"_"+p2+") interference")
			plt.savefig("plots/"+va+"/"+"all_"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')
			plt.close()