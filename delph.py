import numpy as np
#from scipy.integrate import quad
from matplotlib import pyplot as plt
import os,inspect
import matplotlib.ticker as ticker
from matplotlib import gridspec

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
	plt.errorbar(bins_m,yvalues,yerr=yerrors,drawstyle = 'steps-mid',color=color,lw=0.8)

def plot_error_region2ax(ax,yvalues,yerrors,bins,color='orange'):
	bins_a = bins[:-1]
	bins_e = bins[1:]
	bins_m = (bins_a+bins_e)/2
	ax.errorbar(bins_m,yvalues,yerr=yerrors,drawstyle = 'steps-mid',color=color,lw=0.8)

def parameters_kin(sam,var,par):
	topmassup=400
	topmassdown=0
	bmassup=60
	bmassdown=0
	umassup=1
	umassdown=0
	wmassup=100
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
	if((par == "Jet") and (var=="PT")):
		upper_range=500
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
		if("W" in par):
			nbin = 32
		if("Top" in par):
			nbin = 32
		if("Jet" in par):
			nbin = 32
		if("bJet" in par):
			nbin = 32
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
		nbin = 64
		lower_range = 0
		if p1 =="Photon":
			lower_range = 0
			upper_range = 1400
			if p2 == "TopQuark":
				lower_range = 0
				upper_range = 1400
				nbin = 32
			if p2 == "LeadingJet":
				lower_range = 0
				upper_range = 700
			if p2 == "bJet":
				lower_range = 0
				upper_range = 500
				nbin = 32
			if p2 == "WBoson":
				nbin = 32
		if p1 =="TopQuark":
			lower_range = 0
			upper_range = 800
			if p2 == "WBoson":
				lower_range = 50
				nbin = 32
			if p2 == "Photon":
				lower_range = 0
				upper_range = 1400
				nbin = 32
			if p2 == "LeadingJet":
				upper_range = 1200
				nbin = 32
	return [lower_range,upper_range,nbin,up_l]


def labelkin(va,par):
	label = r"$"
	if va == "PT":
		label += r"p_T"
	if va == "Eta":
		label += r"\eta"
	if va == "Phi":
		label += r"\Phi"
	if va == "M":
		label += r"m"
	if par == "Photon":
		label += r"(\gamma)$"
	if par == "TopQuark":
		label += r"(t)$"
	if par == "bJet":
		label += r"(bJet)$"
	if par == "WBoson":
		label += r"(W)$"
	if par == "Jet":
		label += r"(Jet)$"
	if va == "PT" or va == "M":
		label += "/GeV"
	plt.xlabel(label)
	plt.ylabel("number of entries (normalized)",rotation=90)
	ax = plt.gca()
	ax.yaxis.set_label_coords(-0.08,0.5)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	plt.grid(alpha=0.5)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

	if va == "M":
		plt.legend(loc="upper right")
	else:
		plt.legend(loc="best")
	if va == "Eta":
		plt.legend(loc="best",fontsize="small")	

def labelkinax(ax,ax2,va,par):
	label = r"$"
	if va == "PT":
		label += r"p_T"
	if va == "Eta":
		label += r"\eta"
	if va == "Phi":
		label += r"\Phi"
	if va == "M":
		label += r"m"
	if par == "Photon":
		label += r"(\gamma)$"
	if par == "TopQuark":
		label += r"(t)$"
	if par == "bJet":
		label += r"(bJet)$"
	if par == "WBoson":
		label += r"(W)$"
	if par == "Jet":
		label += r"(Jet)$"
	if va == "PT" or va == "M":
		label += "/GeV"
	#ax.xlabel(label)
	ax.set_ylabel("number of entries (normalized)",rotation=90)
	ax.yaxis.set_label_coords(-0.08,0.5)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	ax.grid(alpha=0.5)
	ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	ax2.set_xlabel(label)
	ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	ax2.grid(alpha=0.5)
	ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

	if va == "M":
		ax.legend(loc="upper right")
	else:
		ax.legend(loc="best")
	if va == "Eta":
		ax.legend(loc="best",fontsize="small")		


def labelRM(va,p1,p2):
	vari = r"$"
	if va == "R":
		vari += r"\Delta R"
	if va == "M":
		vari += r"m"

	parts = ""

	if p1 == "Photon":
		parts += r"\gamma,"
	if p1 == "TopQuark":
		parts += r"t,"

	if p2 == "Photon":
		parts += r"\gamma"
	if p2 == "TopQuark":
		parts += r"t"
	if p2 == "bJet":
		parts += r"bJet"
	if p2 == "WBoson":
		parts += r"W"
	if p2 == "LeadingJet":
		parts += r"Jet"

	einh = ""
	if va == "M":
		einh += "/GeV"

	plt.xlabel(vari+"("+parts+")"+r"$"+einh)
	plt.ylabel("number of entries (normalized)",rotation=90)
	ax = plt.gca()
	ax.yaxis.set_label_coords(-0.08,0.5)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	plt.grid(alpha=0.5)
	#if p1 == "Photon" and p2 == "LeadingJet" and va=="M":
	#	pass
	#else:
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")


def labelRMax(ax,ax2,va,p1,p2):
	vari = r"$"
	if va == "R":
		vari += r"\Delta R"
	if va == "M":
		vari += r"m"

	parts = ""

	if p1 == "Photon":
		parts += r"\gamma,"
	if p1 == "TopQuark":
		parts += r"t,"

	if p2 == "Photon":
		parts += r"\gamma"
	if p2 == "TopQuark":
		parts += r"t"
	if p2 == "bJet":
		parts += r"bJet"
	if p2 == "WBoson":
		parts += r"W"
	if p2 == "LeadingJet":
		parts += r"Jet"

	einh = ""
	if va == "M":
		einh += "/GeV"

	#ax.xlabel(vari+"("+parts+")"+r"$"+einh)
	ax.set_ylabel("number of entries (normalized)",rotation=90)
	ax.yaxis.set_label_coords(-0.08,0.5)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	ax.grid(alpha=0.5)
	#if p1 == "Photon" and p2 == "LeadingJet" and va=="M":
	#	pass
	#else:
	ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	ax.legend(loc="best")
	ax2.set_xlabel(vari+"("+parts+")"+r"$"+einh)
	ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	ax2.grid(alpha=0.5)
	ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))


###################################################################################################
########################################DELPH PLOTS################################################
###################################################################################################

samples = ["dec","int","pro"]
vari = ["PT","Eta","Phi","M"]
part = ["TopQuark","Photon","bJet","Jet","WBoson"]

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

decay_fraction = 2266./7084.
prod_fraction = 4818./7084.


ratiox=6
ratioy=4


for par in part:
	for var in vari:
		if(par == "Photon" and var=="M"):
			continue
		para = parameters_kin("dec", var, par)
		lower_range=para[0]
		upper_range=para[1]
		nbin=para[2]
		up_l=para[3]

		evi = np.genfromtxt("data/"+"int"+"_"+par+"_"+var+".txt")
		evd = np.genfromtxt("data/"+"dec"+"_"+par+"_"+var+".txt")
		evp = np.genfromtxt("data/"+"pro"+"_"+par+"_"+var+".txt")
		evi = evi[(np.abs(evi)>up_l) & (evi!=999.9)];evd = evd[(np.abs(evd)>up_l) & (evd!=999.9)];evp = evp[(np.abs(evp)>up_l) & (evp!=999.9)]
		evi = np.clip(evi, lower_range, upper_range);evd = np.clip(evd, lower_range, upper_range);evp = np.clip(evp, lower_range, upper_range)

		npr, bpr, ar = plt.hist(evp,label="production mode",range=(lower_range,upper_range),bins=nbin,histtype="step",lw=0.5,color="blue")
		ndr, bdr, ar = plt.hist(evd,label="decay mode",range=(lower_range,upper_range),bins=bpr,histtype="step",lw=0.5,color="red")
		plt.close()
		npN, bpN, a1 = plt.hist(evp,label="production mode",weights=np.ones_like(evp)/len(evp),range=(lower_range,upper_range),bins=bpr,histtype="step",lw=0.8,color="blue")
		ndN, bdN, a2 = plt.hist(evd,label="decay mode",weights=np.ones_like(evd)/len(evd),range=(lower_range,upper_range),bins=bpr,histtype="step",lw=0.8,color="red")

		plot_error_region2(npN,1/np.sqrt(npr)*npN, bpr,"blue")
		plot_error_region2(ndN,1/np.sqrt(ndr)*ndN, bpr,"red")

		plt.xlim(lower_range,upper_range)
		labelkin(var,par)
		
		plt.savefig("plots/"+par+"/"+"decpro"+"_"+par+"_"+var+".pdf",bbox_inches='tight')

		npNV, bpN, a1 = plt.hist(evp,label="production mode",weights=np.ones_like(evp)*prod_fraction,range=(lower_range,upper_range),bins=bpr,histtype="step",lw=0.8,color="blue")
		ndNV, bdN, a2 = plt.hist(evd,label="decay mode",weights=np.ones_like(evd)*decay_fraction,range=(lower_range,upper_range),bins=bpr,histtype="step",lw=0.8,color="red")
		
		niV, biV, aV = plt.hist(evi,range=(lower_range,upper_range),bins=bpr,histtype="step",lw=0.8,color="green")
		nV,bV,aV = plt.hist(bpr[1:]-np.diff(bpr)/2,bins=bpr,weights=npNV+ndNV,histtype="step",lw=0.8,color="black",range=(lower_range,upper_range))

		plt.close()
		gs = gridspec.GridSpec(2, 1, width_ratios=[1],height_ratios=[3,1]) 
		f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
		ax1 = plt.subplot(gs[0])
		ax2 = plt.subplot(gs[1])

		npN, bpN, a1 = ax1.hist(evp,label="production mode",weights=np.ones_like(evp)/len(evp)*prod_fraction,range=(lower_range,upper_range),bins=bpr,histtype="step",lw=0.8,color="blue")
		ndN, bdN, a2 = ax1.hist(evd,label="decay mode",weights=np.ones_like(evd)/len(evd)*decay_fraction,range=(lower_range,upper_range),bins=bpr,histtype="step",lw=0.8,color="red")
		nir, bir, a = ax1.hist(evi,weights=np.ones_like(evi)/len(evi),label="interference sample",bins=bpr,histtype="step",lw=0.8,color="black",range=(lower_range,upper_range))
		n,b,a = ax1.hist(bpr[1:]-np.diff(bpr)/2,bins=bpr,label="decay+production mode",weights=npN+ndN,histtype="step",lw=0.8,color="green")

		plot_error_region2ax(ax1,npN,1/np.sqrt(npNV)*npN, bpr,"blue")
		plot_error_region2ax(ax1,ndN,1/np.sqrt(ndNV)*ndN, bpr,"red")
		plot_error_region2ax(ax1,nir,1/np.sqrt(niV)*nir, bpr,"black")
		plot_error_region2ax(ax1,n,1/np.sqrt(nV)*n, bpr,"green")

		ax2.hist(bpr[1:]-np.diff(bpr)/2,bins=bpr,weights=np.ones_like(bpr[1:]),histtype="step",lw=0.8,color="black")
		w=n/nir;w[w==np.inf] = np.nan;w[np.isnan(w)] = 0

		nn,bb,aa = ax2.hist(bpr[1:]-np.diff(bpr)/2,alpha=0,bins=bpr,weights=w,histtype="step",lw=0.8,color="green")

		error = np.sqrt((np.sqrt(nV)/nV)**2 + (np.sqrt(niV)/niV)**2)
		ax2.errorbar(bpr[1:]-np.diff(bpr)/2,nn,xerr=np.diff(bpr)/2,yerr=error,fmt="none",color="green",label="(Pro+Dec)/Int",lw=0.8)
		ax2.set_ylabel("Dec+Pro/Int")
		#ax2.legend(loc="best")
		ax2.grid(alpha=0.6)

		ax1.set_xlim(lower_range,upper_range)
		ax2.set_xlim(lower_range,upper_range)


		labelkinax(ax1,ax2,var,par)
		
		plt.savefig("plots/"+par+"/"+"decproint"+"_"+par+"_"+var+".pdf",bbox_inches='tight')


for par in part:
	for var in vari:
		if(par == "Photon" and var=="M"):
			continue

		para = parameters_kin("dec", var, par)
		lower_range=para[0]
		upper_range=para[1]
		nbin=para[2]
		up_l=para[3]

		evp = np.genfromtxt("data/pro_"+par+"_"+var+".txt")
		evd = np.genfromtxt("data/dec_"+par+"_"+var+".txt")
		evi = np.genfromtxt("data/int_"+par+"_"+var+".txt")
		evp = evp[(np.abs(evp)>up_l) & (evp!=999.9)];evd = evd[(np.abs(evd)>up_l) & (evd!=999.9)];evi = evi[(np.abs(evi)>up_l) & (evi!=999.9)]
		ev = np.concatenate((evp, evd), axis=0)
		#print(par+var)
		ev = np.clip(ev, lower_range, upper_range)
		evi = np.clip(evi, lower_range, upper_range)
		evp = np.clip(evp, lower_range, upper_range)
		evd = np.clip(evd, lower_range, upper_range)


		nbpV,binspV,apV = plt.hist(evp,bins=nbin,lw=0.5,alpha=0.0,color="blue",fill=False,weights=np.ones_like(evp)*prod_fraction,range=(lower_range,upper_range),histtype='step')
		nbdV,binsdV,adV = plt.hist(evd,bins=binspV,lw=0.5,alpha=0.0,color="blue",fill=False,weights=np.ones_like(evd)*decay_fraction,range=(lower_range,upper_range),histtype='step')
		nV,binsV,aV 	= plt.hist(binspV[1:]-np.diff(binspV)/2,label=r"production+decay",bins=binspV,lw=0.5,color="blue",fill=False,weights=nbpV+nbdV,range=(lower_range,upper_range),histtype='step')

		vnI,vbinsI,vaI = plt.hist(evi,label=r"interference",bins=nbin,lw=0.5,color="red",fill=False,normed=False,range=(lower_range,upper_range),histtype='step')
		plt.close()

		plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')

		nbp,binsp,ap = plt.hist(evp,bins=nbin,lw=0.5,alpha=0.0,color="blue",fill=False,weights=np.ones_like(evp)*prod_fraction/len(evp),range=(lower_range,upper_range),histtype='step')
		nbd,binsd,ad = plt.hist(evd,bins=binsp,lw=0.5,alpha=0.0,color="blue",fill=False,weights=np.ones_like(evd)*decay_fraction/len(evd),range=(lower_range,upper_range),histtype='step')
		n,bins,a 	 = plt.hist(binsp[1:]-np.diff(binsp)/2,label=r"production + decay mode",bins=binsp,lw=0.8,color="blue",fill=False,weights=nbp+nbd,range=(lower_range,upper_range),histtype='step')

		plot_error_region2(n,1/np.sqrt(nV)*n, bins,"blue")
		nI,binsI,aI = plt.hist(evi,label=r"interference sample",bins=bins,lw=0.8,color="red",fill=False,weights=np.ones_like(evi)/float(len(evi)),range=(lower_range,upper_range),histtype='step')
		plot_error_region2(nI,1/np.sqrt(vnI)*nI, binsI,"red")
		

		plt.xlim(lower_range,upper_range)
		labelkin(var,par)

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




for va in var:
	for p1 in RMPart:
		for p2 in RMPartP:
			if(p1==p2):
				continue

			para = parameters_RM(p1, p2, va)
			lower_range=para[0];upper_range=para[1];nbin=para[2];up_l=para[3]

			evi = np.genfromtxt("data/int_"+p1+"_"+p2+"_"+va+".txt")
			evd = np.genfromtxt("data/dec_"+p1+"_"+p2+"_"+va+".txt")
			evp = np.genfromtxt("data/pro_"+p1+"_"+p2+"_"+va+".txt")
			evi = evi[(np.abs(evi)>up_l) & (evi!=999.9)];evd = evd[(np.abs(evd)>up_l) & (evd!=999.9)];evp = evp[(np.abs(evp)>up_l) & (evp!=999.9)]
			evi = np.clip(evi, lower_range, upper_range);evd = np.clip(evd, lower_range, upper_range);evp = np.clip(evp, lower_range, upper_range)

			npr, bpr, ar = plt.hist(evp,label="production mode",range=(lower_range,upper_range),bins=nbin,histtype="step",lw=0.5,color="blue")
			ndr, bdr, ar = plt.hist(evd,label="decay mode",range=(lower_range,upper_range),bins=bpr,histtype="step",lw=0.5,color="red")
			plt.close()
			npN, bpN, a1 = plt.hist(evp,label="production mode",weights=np.ones_like(evp)/len(evp),range=(lower_range,upper_range),bins=bpr,histtype="step",lw=0.8,color="blue")
			ndN, bdN, a2 = plt.hist(evd,label="decay mode",weights=np.ones_like(evd)/len(evd),range=(lower_range,upper_range),bins=bpr,histtype="step",lw=0.8,color="red")

			plot_error_region2(npN,1/np.sqrt(npr)*npN, bpr,"blue")
			plot_error_region2(ndN,1/np.sqrt(ndr)*ndN, bpr,"red")

			plt.xlim(lower_range,upper_range)
			labelRM(va,p1,p2)
			
			plt.savefig("plots/"+va+"/"+"decpro"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')

			npNV, bpN, a1 = plt.hist(evp,label="production mode",weights=np.ones_like(evp)*prod_fraction,range=(lower_range,upper_range),bins=bpr,histtype="step",lw=0.8,color="blue")
			ndNV, bdN, a2 = plt.hist(evd,label="decay mode",weights=np.ones_like(evd)*decay_fraction,range=(lower_range,upper_range),bins=bpr,histtype="step",lw=0.8,color="red")
			
			niV, biV, aV = plt.hist(evi,range=(lower_range,upper_range),bins=bpr,histtype="step",lw=0.8,color="green")
			nV,bV,aV = plt.hist(bpr[1:]-np.diff(bpr)/2,bins=bpr,weights=npNV+ndNV,histtype="step",lw=0.8,color="black",range=(lower_range,upper_range))

			plt.close()
			gs = gridspec.GridSpec(2, 1, width_ratios=[1],height_ratios=[3,1]) 
			f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
			ax1 = plt.subplot(gs[0])
			ax2 = plt.subplot(gs[1])
			gs.update(wspace=0.03)

			npN, bpN, a1 = ax1.hist(evp,label="production mode",weights=np.ones_like(evp)/len(evp)*prod_fraction,range=(lower_range,upper_range),bins=bpr,histtype="step",lw=0.8,color="blue")
			ndN, bdN, a2 = ax1.hist(evd,label="decay mode",weights=np.ones_like(evd)/len(evd)*decay_fraction,range=(lower_range,upper_range),bins=bpr,histtype="step",lw=0.8,color="red")
			nir, bir, a = ax1.hist(evi,weights=np.ones_like(evi)/len(evi),label="interference sample",bins=bpr,histtype="step",lw=0.8,color="black",range=(lower_range,upper_range))
			n,b,a = ax1.hist(bpr[1:]-np.diff(bpr)/2,bins=bpr,label="decay+production mode",weights=npN+ndN,histtype="step",lw=0.8,color="green")

			plot_error_region2ax(ax1,npN,1/np.sqrt(npNV)*npN, bpr,"blue")
			plot_error_region2ax(ax1,ndN,1/np.sqrt(ndNV)*ndN, bpr,"red")
			plot_error_region2ax(ax1,nir,1/np.sqrt(niV)*nir, bpr,"black")
			plot_error_region2ax(ax1,n,1/np.sqrt(nV)*n, bpr,"green")


			ax2.hist(bpr[1:]-np.diff(bpr)/2,bins=bpr,weights=np.ones_like(bpr[1:]),histtype="step",lw=0.8,color="black")
			w=n/nir;w[w==np.inf] = np.nan;w[np.isnan(w)] = 0

			nn,bb,aa = ax2.hist(bpr[1:]-np.diff(bpr)/2,alpha=0,bins=bpr,weights=w,histtype="step",lw=0.8,color="green")

			error = np.sqrt((np.sqrt(nV)/nV)**2 + (np.sqrt(niV)/niV)**2)
			ax2.errorbar(bpr[1:]-np.diff(bpr)/2,nn,xerr=np.diff(bpr)/2,yerr=error,fmt="none",color="green",label="(Pro+Dec)/Int",lw=0.8)
			ax2.set_ylabel("Dec+Pro/Int")
			#ax2.legend(loc="best")
			ax2.grid(alpha=0.6)

			ax1.set_xlim(lower_range,upper_range)
			ax2.set_xlim(lower_range,upper_range)

			labelRMax(ax1,ax2,va,p1,p2)
			
			plt.savefig("plots/"+va+"/"+"decproint_"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')




for va in var:
	for p1 in RMPart:
		for p2 in RMPartP:
			if(p1==p2):
				continue

			para = parameters_RM(p1, p2, va)
			lower_range=para[0];upper_range=para[1];nbin=para[2];up_l=para[3]

			evp = np.genfromtxt("data/pro_"+p1+"_"+p2+"_"+va+".txt")
			evd = np.genfromtxt("data/dec_"+p1+"_"+p2+"_"+va+".txt")
			evi = np.genfromtxt("data/int_"+p1+"_"+p2+"_"+va+".txt")
			evp = evp[(np.abs(evp)>up_l) & (evp!=999.9)];evd = evd[(np.abs(evd)>up_l) & (evd!=999.9)];evi = evi[(np.abs(evi)>up_l) & (evi!=999.9)]
			evi = np.clip(evi, lower_range, upper_range);evp = np.clip(evp, lower_range, upper_range);evd = np.clip(evd, lower_range, upper_range)
			#ev = np.concatenate((evp, evd), axis=0)
			#ev = np.clip(ev, lower_range, upper_range)

			fig = plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')


			nbpV2,binspV2,apV2 = plt.hist(evp,bins=nbin,lw=0.5,alpha=0.0,color="blue",fill=False,weights=np.ones_like(evp)*prod_fraction,range=(lower_range,upper_range),histtype='step')
			nbdV2,binsdV2,adV2 = plt.hist(evd,bins=binspV2,lw=0.5,alpha=0.0,color="blue",fill=False,weights=np.ones_like(evd)*decay_fraction,range=(lower_range,upper_range),histtype='step')
			nV2,binsV2,aV2 	= plt.hist(binspV2[1:]-np.diff(binspV2)/2,label=r"production+decay",bins=binspV2,lw=0.5,color="blue",fill=False,weights=nbpV2+nbdV2,range=(lower_range,upper_range),histtype='step')

			vnI2,vbinsI2,vaI2 = plt.hist(evi,label=r"interference",bins=nbin,lw=0.8,color="red",fill=False,normed=False,histtype='step')#,range=(lower_range,upper_range))
			
			plt.close()
			plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')
			
			nbp2,binsp2,ap2 = plt.hist(evp,bins=binspV2,lw=0.5,alpha=0.0,color="blue",fill=False,weights=np.ones_like(evp)*prod_fraction/len(evp),range=(lower_range,upper_range),histtype='step')
			nbd2,binsd2,ad2 = plt.hist(evd,bins=binspV2,lw=0.5,alpha=0.0,color="blue",fill=False,weights=np.ones_like(evd)*decay_fraction/len(evd),range=(lower_range,upper_range),histtype='step')
			n2,bins2,a2 	= plt.hist(binspV2[1:]-np.diff(binspV2)/2,label=r"production + decay mode",bins=binspV2,lw=0.8,color="blue",fill=False,weights=nbp2+nbd2,range=(lower_range,upper_range),histtype='step')
			
			nI,binsI,aI = plt.hist(evi,label=r"interference sample",bins=bins2,lw=0.8,color="red",fill=False,weights=np.ones_like(evi)/float(len(evi)),histtype='step')#,range=(lower_range,upper_range))
			
			plot_error_region2(n2,1/np.sqrt(nV2)*n2, bins2,"blue")
			plot_error_region2(nI,1/np.sqrt(vnI2)*nI, binsI,"red")

			plt.xlim(lower_range,upper_range)
			labelRM(va,p1,p2)
			plt.savefig("plots/"+va+"/"+"all_"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')
			plt.close()





		#plt.errorbar(binsp[:-1]+np.diff(binsp)/2,nbp+nbd,drawstyle = 'steps-mid',color="blue",lw=0.5)
		#plt.bar(binsp[:-1],(nbp+nbd),align="edge",label=r"production+decay",width=np.diff(binsp),lw=0.5,fill=False,edgecolor="blue")
