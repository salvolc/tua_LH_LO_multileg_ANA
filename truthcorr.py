import numpy as np
#from scipy.integrate import quad
from matplotlib import pyplot as plt
import os,inspect
import matplotlib.ticker as ticker

def parameters_kin(var,par):
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
	default_nbins = 25

	mass_nbins = 25

	nbin = default_nbins
	lower_range=default_range_down
	upper_range=default_range_up
	up_l = 0.1
	if(par == "Photon" and var=="PT"):
		up_l=20
		upper_range=500
		lower_range=20
		nbin=default_nbins
	if((sam=="pro") | (sam=="int")):
		upper_range=500
	if (var == "Eta"):
		lower_range=-5
		upper_range=5
		nbin = 16
		up_l = 0
	if (var == "Phi"):
		lower_range=-3.14
		upper_range=3.14
		nbin = 16
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
			upper_range=400
	return [lower_range,upper_range,nbin,up_l]


sam = ["dec","int","pro"]
var = ["PT","Eta","Phi","M"]
par = ["BQuark","Photon","WBoson","TopQuark"]

PREFIX = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/"

datas = os.listdir(PREFIX+"data")


if("corr" not in os.listdir(PREFIX+"plots/")):
	os.chdir(PREFIX+"plots/")
	os.mkdir("corr")
	os.chdir(PREFIX)

ratiox= 8
ratioy= 10

i=0;j=0;
for va1 in var:
	i +=1
	for pa1 in par:
		j +=1
		for va2 in var:
			for pa2 in par:
				if (pa1 == pa2):
					continue

				if((pa1 == "Photon" and va1 == "M") or (pa2 == "Photon" and va2 =="M")):
					continue

				ev1p = np.genfromtxt("data/pro_"+pa1+"_"+va1+"_truth.txt");ev1d = np.genfromtxt("data/dec_"+pa1+"_"+va1+"_truth.txt");ev1i = np.genfromtxt("data/int_"+pa1+"_"+va1+"_truth.txt")
				ev1p = ev1p[(np.abs(ev1p)>0) & (ev1p!=999.9)];ev1d = ev1d[(np.abs(ev1d)>0) & (ev1d!=999.9)];ev1i = ev1i[(np.abs(ev1i)>0) & (ev1i!=999.9)]
				ev1pd = np.concatenate((ev1p, ev1d), axis=0)

				ev2p = np.genfromtxt("data/pro_"+pa2+"_"+va2+"_truth.txt");ev2d = np.genfromtxt("data/dec_"+pa2+"_"+va2+"_truth.txt");ev2i = np.genfromtxt("data/int_"+pa2+"_"+va2+"_truth.txt")
				ev2p = ev2p[(np.abs(ev2p)>0) & (ev2p!=999.9)];ev2d = ev2d[(np.abs(ev2d)>0) & (ev2d!=999.9)];ev2i = ev2i[(np.abs(ev2i)>0) & (ev2i!=999.9)]
				ev2pd = np.concatenate((ev2p, ev2d), axis=0)

				plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')


				min1, max1, nbin1, upl1 = parameters_kin(va1, pa1)
				min2, max2, nbin2, upl2 = parameters_kin(va2, pa2)

				#if len(ev1i) > len(ev2i):
				#	help = np.ones(len(ev1i)-len(ev2i))*9999
				#	ev2i = np.concatenate((ev2i,help), axis=0)
				#if len(ev1i) < len(ev2i):
				#	help = np.ones(len(ev2i)-len(ev1i))*9999
				#	ev1i = np.concatenate((ev1i,help), axis=0)

				w = np.ones(len(ev1i))
				N = len(ev1i)

				f, (ax1, ax2) = plt.subplots(2, 1, num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')
				"""bins=[nbin1,nbin2],range=([min1,max1],[min2,max2]),"""
				h, xedg, yedg, img = ax1.hist2d(ev1i,ev2i,bins=[nbin1,nbin2],range=([min1,max1],[min2,max2]),weights=w/N)

				"""if len(ev1pd) > len(ev2pd):
					help = np.ones(len(ev1pd)-len(ev2pd))*9999
					ev2pd = np.concatenate((ev2pd,help), axis=0)
				if len(ev1pd) < len(ev2pd):
					help = np.ones(len(ev2pd)-len(ev1pd))*9999
					ev1pd = np.concatenate((ev1pd,help), axis=0)"""

				w = np.ones(len(ev1pd))
				N = len(ev1pd)
				"""bins=[nbin1,nbin2],range=([min1,max1],[min2,max2]),"""
				h, xedg, yedg, img = ax2.hist2d(ev1pd,ev2pd,bins=[nbin1,nbin2],range=([min1,max1],[min2,max2]),weights=w/N)
				ax1.set_title("Decay + production mode")
				ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator())
				ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
				ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
				#ax1.title(va1+"("+pa1+") "+va2+"("+pa2+")"+" interference")
				ax1.set(xlabel=r"$"+va1+"("+pa1+")"+r"$",ylabel=r"$"+va2+"("+pa2+")"+r"$")
				ax2.set_title("Interference mode")
				ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator())
				ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
				ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
				#ax2.title(va1+"("+pa1+") "+va2+"("+pa2+")"+" prod+decay")
				ax2.set(xlabel=r"$"+va1+"("+pa1+")"+r"$",ylabel=r"$"+va2+"("+pa2+")"+r"$")


				plt.savefig("plots/corr/"+pa1+"_"+va1+"_"+pa2+"_"+va2+".png",dpi=200,bbox_inches='tight')
				plt.close()












