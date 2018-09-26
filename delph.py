from dfunk import *

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

sdp = 1571.+3370.
decay_fraction = 1571./sdp
prod_fraction = 3370./sdp


ratiox=6
ratioy=4
eps=0.001

sample_numbers=[1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]


for par in part:
	for var in vari:
		if(par == "Photon" and var=="M"):
			continue
		para = parameters_kin("dec", var, par)
		lower_range=para[0]
		upper_range=para[1]
		nbin=para[2]
		up_l=para[3]

		evi = np.genfromtxt("data/"+"int"+"_"+par+"_"+var+"_1.txt")
		evd = np.genfromtxt("data/"+"dec"+"_"+par+"_"+var+"_1.txt")
		evp = np.genfromtxt("data/"+"pro"+"_"+par+"_"+var+"_1.txt")
		evi = evi[(np.abs(evi)>up_l) & (evi!=999.9)];evd = evd[(np.abs(evd)>up_l) & (evd!=999.9)];evp = evp[(np.abs(evp)>up_l) & (evp!=999.9)]
		evi = np.clip(evi, lower_range+eps, upper_range-eps);evd = np.clip(evd, lower_range+eps, upper_range-eps);evp = np.clip(evp, lower_range+eps, upper_range-eps)

		binning = np.linspace(lower_range, upper_range, nbin)
		hd = rplot.Hist(binning)
		hp = rplot.Hist(binning)

		for sample_number in sample_numbers:
			evi = np.genfromtxt("data/"+"int"+"_"+par+"_"+var+"_"+str(sample_number)+".txt")
			evd = np.genfromtxt("data/"+"dec"+"_"+par+"_"+var+"_"+str(sample_number)+".txt")
			evp = np.genfromtxt("data/"+"pro"+"_"+par+"_"+var+"_"+str(sample_number)+".txt")

			evi = evi[(np.abs(evi)>up_l) & (evi!=999.9)];evd = evd[(np.abs(evd)>up_l) & (evd!=999.9)];evp = evp[(np.abs(evp)>up_l) & (evp!=999.9)]
			evi = np.clip(evi, lower_range+eps, upper_range-eps);evd = np.clip(evd, lower_range+eps, upper_range-eps);evp = np.clip(evp, lower_range+eps, upper_range-eps)
			
			decay_fraction = np.genfromtxt("cross/tua_LH_decay_LO_multileg_tt_wbua_2_"+str(sample_number)+"_cross")
			prod_fraction = np.genfromtxt("cross/tua_LH_production_LO_multileg_ta_taj_2_"+str(sample_number)+"_cross")

			sdp = decay_fraction+prod_fraction
			decay_fraction = decay_fraction/sdp
			prod_fraction = prod_fraction/sdp

			map(hd.Fill, evd,np.ones_like(evd)*decay_fraction)
			map(hp.Fill, evp,np.ones_like(evp)*prod_fraction)

		hp.Scale(1/(hp.Integral(0,hp.GetNbinsX()+1)))
		hp.linecolor = "blue";hp.linewidth = 1
		rplt.hist(hp,fmt="none",lw=0.4,color="blue",label="production mode")
		rplt.errorbar(hp,fmt="none",lw=0.4,color="blue",label="_nolegend_")

		hd.Scale(1/(hd.Integral(0,hd.GetNbinsX()+1)))
		hd.linecolor = "red";hd.linewidth = 1
		rplt.hist(hd,fmt="none",lw=0.4,color="red",label="decay mode")
		rplt.errorbar(hd,fmt="none",lw=0.4,color="red",label="_nolegend_")

		plt.xlim(lower_range,upper_range)
		labelkin(var,par)
		
		left = lower_range
		right = upper_range
		ax = plt.gca()
		bottom, top = ax.get_ylim()
		ax.text(0.7*(left+right),0.4*(bottom+top),r"$\sqrt{s}=13\,$ TeV"+"\n"+"reco level")

		plt.savefig("plots/"+par+"/"+"decpro"+"_"+par+"_"+var+".pdf",bbox_inches='tight')
		plt.close()

		gs = gridspec.GridSpec(2, 1, width_ratios=[1],height_ratios=[3,1]) 
		f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
		ax1 = plt.subplot(gs[0])
		ax2 = plt.subplot(gs[1])

		#binning = np.linspace(lower_range, upper_range, nbin)
		binning = set_dyn_binning(evi, lower_range, upper_range, nbin)
		hi		= rplot.Hist(binning)
		hdp		= rplot.Hist(binning)
		hd = rplot.Hist(binning)
		hp = rplot.Hist(binning)

		hi.Sumw2(True)
		hdp.Sumw2(True)
		hd.Sumw2(True)
		hp.Sumw2(True)

		for sample_number in sample_numbers:
			evi = np.genfromtxt("data/"+"int"+"_"+par+"_"+var+"_"+str(sample_number)+".txt")
			evd = np.genfromtxt("data/"+"dec"+"_"+par+"_"+var+"_"+str(sample_number)+".txt")
			evp = np.genfromtxt("data/"+"pro"+"_"+par+"_"+var+"_"+str(sample_number)+".txt")

			evi = evi[(np.abs(evi)>up_l) & (evi!=999.9)];evd = evd[(np.abs(evd)>up_l) & (evd!=999.9)];evp = evp[(np.abs(evp)>up_l) & (evp!=999.9)]
			evi = np.clip(evi, lower_range+eps, upper_range-eps);evd = np.clip(evd, lower_range+eps, upper_range-eps);evp = np.clip(evp, lower_range+eps, upper_range-eps)
			
			decay_fraction = np.genfromtxt("cross/tua_LH_decay_LO_multileg_tt_wbua_2_"+str(sample_number)+"_cross")
			prod_fraction = np.genfromtxt("cross/tua_LH_production_LO_multileg_ta_taj_2_"+str(sample_number)+"_cross")

			sdp = decay_fraction+prod_fraction
			decay_fraction = decay_fraction/sdp
			prod_fraction = prod_fraction/sdp

			map(hd.Fill, evd,np.ones_like(evd)*decay_fraction)
			map(hp.Fill, evp,np.ones_like(evp)*prod_fraction)
			map(hi.Fill, evi)
			map(hdp.Fill, evd,np.ones_like(evd)*decay_fraction)
			map(hdp.Fill, evp,np.ones_like(evp)*prod_fraction)

		hp.Scale(1/(hdp.Integral(0,hdp.GetNbinsX()+1)),"width")
		hp.linecolor = "blue";hp.linewidth = 1
		rplt.hist(hp,fmt="none",axes=ax1,lw=0.4,color="blue",label="production mode")
		rplt.errorbar(hp,fmt="none",axes=ax1,lw=0.4,color="blue",label="_nolegend_")

		hd.Scale(1/(hdp.Integral(0,hdp.GetNbinsX()+1)),"width")
		hd.linecolor = "red";hd.linewidth = 1
		rplt.hist(hd,fmt="none",axes=ax1,lw=0.4,color="red",label="decay mode")
		rplt.errorbar(hd,fmt="none",axes=ax1,lw=0.4,color="red",label="_nolegend_")

		hi.Scale(1/(hi.Integral(0,hi.GetNbinsX()+1)),"width")
		hi.linecolor = "black";hi.linewidth = 1
		rplt.hist(hi,fmt="none",axes=ax1,lw=0.4,color="black",label="interference mode")
		rplt.errorbar(hi,fmt="none",axes=ax1,lw=0.4,color="black",label="_nolegend_")

		hdp.Scale(1/(hdp.Integral(0,hdp.GetNbinsX()+1)),"width")
		hdp.linecolor = "green";hdp.linewidth = 1
		rplt.hist(hdp,fmt="none",axes=ax1,lw=0.4,color="green",label="decay + production mode")
		rplt.errorbar(hdp,fmt="none",axes=ax1,lw=0.4,color="green",label="_nolegend_")

		print(par + " " + var + " KS Test")
		ksval = hi.KolmogorovTest(hdp)
		print(hi.KolmogorovTest(hdp))
		print()


		ax2.hist(binning[1:]-np.diff(binning)/2,bins=binning,weights=np.ones_like(binning[1:]),histtype="step",lw=0.8,color="black")
		hdp.Divide(hi)
		rplt.errorbar(hdp,fmt="none",axes=ax2,lw=0.6,color="green",label="_nolegend_")
		
		ax2.set_ylabel("(Dec.+Pro.)/Int.")
		ax2.grid(alpha=0.6)

		ax1.set_xlim(lower_range,upper_range)
		ax2.set_xlim(lower_range,upper_range)
		ax2.set_ylim(0.5,1.5)

		labelkinax(ax1,ax2,var,par)
		leg = ax1.legend()#title="KS Test: "+"{:.7f}".format(ksval),fontsize="small")
		#leg.get_title().set_fontsize('small')
		
		plt.savefig("plots/"+par+"/"+"decproint"+"_"+par+"_"+var+".pdf",bbox_inches='tight')
		plt.close()


"""for par in part:
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
		nI,binsI,aI = plt.hist(evi,label=r"interference mode",bins=bins,lw=0.8,color="red",fill=False,weights=np.ones_like(evi)/float(len(evi)),range=(lower_range,upper_range),histtype='step')
		plot_error_region2(nI,1/np.sqrt(vnI)*nI, binsI,"red")
		

		plt.xlim(lower_range,upper_range)
		labelkin(var,par)

		plt.savefig("plots/"+par+"/"+"all"+"_"+par+"_"+var+".pdf",bbox_inches='tight')
		plt.close()"""



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

			evi = np.genfromtxt("data/int_"+p1+"_"+p2+"_"+va+"_1.txt")
			evd = np.genfromtxt("data/dec_"+p1+"_"+p2+"_"+va+"_1.txt")
			evp = np.genfromtxt("data/pro_"+p1+"_"+p2+"_"+va+"_1.txt")
			evi = evi[(np.abs(evi)>up_l) & (evi!=999.9)];evd = evd[(np.abs(evd)>up_l) & (evd!=999.9)];evp = evp[(np.abs(evp)>up_l) & (evp!=999.9)]
			evi = np.clip(evi, lower_range+eps, upper_range-eps);evd = np.clip(evd, lower_range+eps, upper_range-eps);evp = np.clip(evp, lower_range+eps, upper_range-eps)

			binning = np.linspace(lower_range, upper_range, nbin)
			hd = rplot.Hist(binning)
			hp = rplot.Hist(binning)

			for sample_number in sample_numbers:
				evi = np.genfromtxt("data/"+"int"+"_"+p1+"_"+p2+"_"+va+"_"+str(sample_number)+".txt")
				evd = np.genfromtxt("data/"+"dec"+"_"+p1+"_"+p2+"_"+va+"_"+str(sample_number)+".txt")
				evp = np.genfromtxt("data/"+"pro"+"_"+p1+"_"+p2+"_"+va+"_"+str(sample_number)+".txt")

				evi = evi[(np.abs(evi)>up_l) & (evi!=999.9)];evd = evd[(np.abs(evd)>up_l) & (evd!=999.9)];evp = evp[(np.abs(evp)>up_l) & (evp!=999.9)]
				evi = np.clip(evi, lower_range+eps, upper_range-eps);evd = np.clip(evd, lower_range+eps, upper_range-eps);evp = np.clip(evp, lower_range+eps, upper_range-eps)
				
				decay_fraction = np.genfromtxt("cross/tua_LH_decay_LO_multileg_tt_wbua_2_"+str(sample_number)+"_cross")
				prod_fraction = np.genfromtxt("cross/tua_LH_production_LO_multileg_ta_taj_2_"+str(sample_number)+"_cross")

				sdp = decay_fraction+prod_fraction
				decay_fraction = decay_fraction/sdp
				prod_fraction = prod_fraction/sdp

				map(hd.Fill, evd,np.ones_like(evd)*decay_fraction)
				map(hp.Fill, evp,np.ones_like(evp)*prod_fraction)

			hp.Scale(1/(hp.Integral(0,hp.GetNbinsX()+1)))
			hp.linecolor = "blue";hp.linewidth = 1
			rplt.hist(hp,fmt="none",lw=0.4,color="blue",label="production mode")
			rplt.errorbar(hp,fmt="none",lw=0.4,color="blue",label="_nolegend_")

			hd.Scale(1/(hd.Integral(0,hd.GetNbinsX()+1)))
			hd.linecolor = "red";hd.linewidth = 1
			rplt.hist(hd,fmt="none",lw=0.4,color="red",label="decay mode")
			rplt.errorbar(hd,fmt="none",lw=0.4,color="red",label="_nolegend_")

			plt.xlim(lower_range,upper_range)
			labelRM(va,p1,p2)

			left = lower_range
			right = upper_range
			ax = plt.gca()
			bottom, top = ax.get_ylim()
			ax.text(0.7*(left+right),0.4*(bottom+top),r"$\sqrt{s}=13\,$ TeV"+"\n"+"truth level")
				
			plt.savefig("plots/"+va+"/"+"decpro"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')
			plt.close()

			gs = gridspec.GridSpec(2, 1, width_ratios=[1],height_ratios=[3,1]) 
			f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
			ax1 = plt.subplot(gs[0])
			ax2 = plt.subplot(gs[1])
			gs.update(wspace=0.03)

			#binning = np.linspace(lower_range, upper_range, nbin)
			binning = set_dyn_binning(evi, lower_range, upper_range, nbin)
			hi		= rplot.Hist(binning)
			hdp		= rplot.Hist(binning)
			hd = rplot.Hist(binning)
			hp = rplot.Hist(binning)

			for sample_number in sample_numbers:
				evi = np.genfromtxt("data/"+"int"+"_"+p1+"_"+p2+"_"+va+"_"+str(sample_number)+".txt")
				evd = np.genfromtxt("data/"+"dec"+"_"+p1+"_"+p2+"_"+va+"_"+str(sample_number)+".txt")
				evp = np.genfromtxt("data/"+"pro"+"_"+p1+"_"+p2+"_"+va+"_"+str(sample_number)+".txt")

				evi = evi[(np.abs(evi)>up_l) & (evi!=999.9)];evd = evd[(np.abs(evd)>up_l) & (evd!=999.9)];evp = evp[(np.abs(evp)>up_l) & (evp!=999.9)]
				evi = np.clip(evi, lower_range+eps, upper_range-eps);evd = np.clip(evd, lower_range+eps, upper_range-eps);evp = np.clip(evp, lower_range+eps, upper_range-eps)
				
				decay_fraction = np.genfromtxt("cross/tua_LH_decay_LO_multileg_tt_wbua_2_"+str(sample_number)+"_cross")
				prod_fraction = np.genfromtxt("cross/tua_LH_production_LO_multileg_ta_taj_2_"+str(sample_number)+"_cross")

				sdp = decay_fraction+prod_fraction
				decay_fraction = decay_fraction/sdp
				prod_fraction = prod_fraction/sdp

				map(hd.Fill, evd,np.ones_like(evd)*decay_fraction)
				map(hp.Fill, evp,np.ones_like(evp)*prod_fraction)
				map(hi.Fill, evi)
				map(hdp.Fill, evd,np.ones_like(evd)*decay_fraction)
				map(hdp.Fill, evp,np.ones_like(evp)*prod_fraction)

			hp.Scale(1/(hdp.Integral(0,hdp.GetNbinsX()+1)),"width")
			hp.linecolor = "blue";hp.linewidth = 1
			rplt.hist(hp,fmt="none",axes=ax1,lw=0.4,color="blue",label="production mode")
			rplt.errorbar(hp,fmt="none",axes=ax1,lw=0.4,color="blue",label="_nolegend_")

			hd.Scale(1/(hdp.Integral(0,hdp.GetNbinsX()+1)),"width")
			hd.linecolor = "red";hd.linewidth = 1
			rplt.hist(hd,fmt="none",axes=ax1,lw=0.4,color="red",label="decay mode")
			rplt.errorbar(hd,fmt="none",axes=ax1,lw=0.4,color="red",label="_nolegend_")

			hi.Scale(1/(hi.Integral(0,hi.GetNbinsX()+1)),"width")
			hi.linecolor = "black";hi.linewidth = 1
			rplt.hist(hi,fmt="none",axes=ax1,lw=0.4,color="black",label="interference mode")
			rplt.errorbar(hi,fmt="none",axes=ax1,lw=0.4,color="black",label="_nolegend_")

			hdp.Scale(1/(hdp.Integral(0,hdp.GetNbinsX()+1)),"width")
			hdp.linecolor = "green";hdp.linewidth = 1
			rplt.hist(hdp,fmt="none",axes=ax1,lw=0.4,color="green",label="decay + production mode")
			rplt.errorbar(hdp,fmt="none",axes=ax1,lw=0.4,color="green",label="_nolegend_")

			print(p1 + " " + p2 + " " + va + " KS Test")
			ksval = hi.KolmogorovTest(hdp)
			print(hi.KolmogorovTest(hdp))
			print()

			ax2.hist(binning[1:]-np.diff(binning)/2,bins=binning,weights=np.ones_like(binning[1:]),histtype="step",lw=0.8,color="black")
			hdp.Divide(hi)
			rplt.errorbar(hdp,fmt="none",axes=ax2,lw=0.6,color="green",label="_nolegend_")

			ax2.set_ylabel("(Dec.+Pro.)/Int.")

			#ax2.legend(loc="best")
			ax2.grid(alpha=0.6)

			ax1.set_xlim(lower_range,upper_range)
			ax2.set_xlim(lower_range,upper_range)
			ax2.set_ylim(0.5,1.5)

			labelRMax(ax1,ax2,va,p1,p2)
			leg = ax1.legend()#title="KS Test: "+"{:.7f}".format(ksval),fontsize="small")
			#leg.get_title().set_fontsize('small')
			
			plt.savefig("plots/"+va+"/"+"decproint_"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')
			plt.close()



"""for va in var:
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
			
			nI,binsI,aI = plt.hist(evi,label=r"interference mode",bins=bins2,lw=0.8,color="red",fill=False,weights=np.ones_like(evi)/float(len(evi)),histtype='step')#,range=(lower_range,upper_range))
			
			plot_error_region2(n2,1/np.sqrt(nV2)*n2, bins2,"blue")
			plot_error_region2(nI,1/np.sqrt(vnI2)*nI, binsI,"red")

			plt.xlim(lower_range,upper_range)
			labelRM(va,p1,p2)
			plt.savefig("plots/"+va+"/"+"all_"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')
			plt.close()"""





		#plt.errorbar(binsp[:-1]+np.diff(binsp)/2,nbp+nbd,drawstyle = 'steps-mid',color="blue",lw=0.5)
		#plt.bar(binsp[:-1],(nbp+nbd),align="edge",label=r"production+decay",width=np.diff(binsp),lw=0.5,fill=False,edgecolor="blue")
