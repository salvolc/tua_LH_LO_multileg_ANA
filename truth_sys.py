from tfunk import *

samples = ["dec","int","pro"]
vari = ["PT","Eta","Phi","M"]
part = ["TopQuark","Photon","BQuark","WBoson","UQuark"]

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


sdp = 1571.+3370.
decay_fraction = 1571./sdp
prod_fraction = 3370./sdp



ratiox=6
ratioy=4
eps=0.001

tfl = ROOT.TFile("tua_crosscheck_salvatore.root","recreate")
sample_numbers=[1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

for par in part:
	for var in vari:
		if(par == "Photon" and var=="M"):
			continue
		para = parameters_kin("dec", var, par)
		lower_range=para[0];upper_range=para[1];nbin=para[2];up_l=para[3]

		evi = np.genfromtxt("data/"+"int"+"_"+par+"_"+var+"_truth_1.txt")
		evd = np.genfromtxt("data/"+"dec"+"_"+par+"_"+var+"_truth_1.txt")
		evp = np.genfromtxt("data/"+"pro"+"_"+par+"_"+var+"_truth_1.txt")
		evi = evi[(np.abs(evi)>up_l) & (evi!=999.9)];evd = evd[(np.abs(evd)>up_l) & (evd!=999.9)];evp = evp[(np.abs(evp)>up_l) & (evp!=999.9)]
		evi = np.clip(evi, lower_range+eps, upper_range-eps);evd = np.clip(evd, lower_range+eps, upper_range-eps);evp = np.clip(evp, lower_range+eps, upper_range-eps)

		binning = np.linspace(lower_range, upper_range, nbin)

		gs = gridspec.GridSpec(2, 1, width_ratios=[1],height_ratios=[3,1]) 
		f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
		ax1 = plt.subplot(gs[0])
		ax2 = plt.subplot(gs[1])

		binning = set_dyn_binning(evi, lower_range, upper_range, nbin)

		hi		= rplot.Hist(binning)
		hisys	= rplot.Hist(binning)
		hdp		= rplot.Hist(binning)
		hd = rplot.Hist(binning)
		hp = rplot.Hist(binning)

		hi.Sumw2(True)
		hisys.Sumw2(True)
		hdp.Sumw2(True)
		hd.Sumw2(True)
		hp.Sumw2(True)

		for sample_number in sample_numbers:
			evi = np.genfromtxt("data/"+"int"+"_"+par+"_"+var+"_truth_"+str(sample_number)+".txt")
			evd = np.genfromtxt("data/"+"dec"+"_"+par+"_"+var+"_truth_"+str(sample_number)+".txt")
			evp = np.genfromtxt("data/"+"pro"+"_"+par+"_"+var+"_truth_"+str(sample_number)+".txt")

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
			map(hisys.Fill, evi)


		hp.Scale(1/(hp.Integral(0,hp.GetNbinsX()+1))*prod_fraction,"width")
		hp.linecolor = "blue";hp.linewidth = 1
		rplt.hist(hp,fmt="none",axes=ax1,lw=0.4,color="blue",label="production sample")
		rplt.errorbar(hp,fmt="none",axes=ax1,lw=0.4,color="blue",label="_nolegend_")

		hd.Scale(1/(hd.Integral(0,hd.GetNbinsX()+1))*decay_fraction,"width")
		hd.linecolor = "red";hd.linewidth = 1
		rplt.hist(hd,fmt="none",axes=ax1,lw=0.4,color="red",label="decay sample")
		rplt.errorbar(hd,fmt="none",axes=ax1,lw=0.4,color="red",label="_nolegend_")

		hi.Scale(1/(hi.Integral(0,hi.GetNbinsX()+1)),"width")
		hi.linecolor = "black";hi.linewidth = 1
		rplt.hist(hi,fmt="none",axes=ax1,lw=0.4,color="black",label="interference sample")
		rplt.errorbar(hi,fmt="none",axes=ax1,lw=0.4,color="black",label="_nolegend_")

		hdp.Scale(1/(hdp.Integral(0,hdp.GetNbinsX()+1)),"width")
		hdp.linecolor = "green";hdp.linewidth = 1
		rplt.hist(hdp,fmt="none",axes=ax1,lw=0.4,color="green",label="decay + production sample")
		rplt.errorbar(hdp,fmt="none",axes=ax1,lw=0.4,color="green",label="_nolegend_")

		f = open("../scalescan/scalevarexp/scale_"+var+"_"+par+"_truth.txt")
		ledge = []
		uperr = []
		downerr = []

		while True:
			line = f.readline()
			if line == '':
				break
			lsp = line.split()
			ledge.append(float(lsp[0]))
			line = f.readline()
			lsp = line.split()
			uperr.append(float(lsp[0]))
			downerr.append(float(lsp[1]))

		errbin = 0
		for i in range(0,hisys.GetNbinsX()+1):
			while hisys.GetBinCenter(i) > ledge[errbin+1]:
				errbin += 1
			cont = hisys.GetBinContent(i)
			errup = 	np.fabs((uperr[errbin]*cont)-cont)
			errdown = 	np.fabs((downerr[errbin]*cont)-cont)

			if errup > errdown:
				syserr = errup
			else:
				syserr = errdown

			hisys.SetBinError(i,hisys.GetBinError(i)+syserr)

			cont = hi.GetBinContent(i)
			errup = np.fabs(uperr[errbin]-1)*cont
			errdown = np.fabs(downerr[errbin]-1)*cont
			#errup = np.fabs(uperr[errbin]-1)+hi.GetBinError(i)
			#errdown = np.fabs(downerr[errbin]-1)+hi.GetBinError(i)

			if errup > errdown:
				syserr = errup
			else:
				syserr = errdown

			y1pl = [cont+hi.GetBinError(i)+syserr,cont+hi.GetBinError(i)+syserr]
			y2pl = [cont-hi.GetBinError(i)-syserr,cont-hi.GetBinError(i)-syserr]
			xpl  = [hisys.GetBinLowEdge(i),hisys.GetBinLowEdge(i)+(hisys.GetBinCenter(i)-hisys.GetBinLowEdge(i))*2]
			if i == 0:
				ax1.fill_between(xpl,y1pl,y2pl,color="magenta",alpha=0.5,label="scale variation uncertainty",lw=0)
				continue
			ax1.fill_between(xpl,y1pl,y2pl,color="magenta",alpha=0.5,lw=0)
			
		ksval = hisys.KolmogorovTest(hdp)

		hisys.Scale(1/(hisys.Integral(0,hisys.GetNbinsX()+1)),"width")
		hisys.Divide(hi)
		for i in range(0,hisys.GetNbinsX()+1):
			cont = hisys.GetBinContent(i)
			xpl  = [hisys.GetBinLowEdge(i),hisys.GetBinLowEdge(i)+(hisys.GetBinCenter(i)-hisys.GetBinLowEdge(i))*2]
			y1pl = [cont+hisys.GetBinError(i),cont+hisys.GetBinError(i)]
			y2pl = [cont-hisys.GetBinError(i),cont-hisys.GetBinError(i)]
			ax2.fill_between(xpl,y1pl,y2pl,color="magenta",alpha=0.5,lw=0)



		hisys.linecolor = "magenta";hisys.linewidth = 1
		#rplt.hist(hisys,fmt="none",axes=ax1,lw=0.4,color="black",label="interference sample")
		#rplt.errorbar(hisys,fmt="none",axes=ax1,lw=0.4,color="magenta",label="_nolegend_")




		#rplt.errorbar(hisys,fmt="none",axes=ax2,lw=0.6,color="black",label="_nolegend_")



		ax2.hist(binning[1:]-np.diff(binning)/2,bins=binning,weights=np.ones_like(binning[1:]),histtype="step",lw=0.8,color="black")

		hdp.Divide(hi)
		rplt.errorbar(hdp,fmt="none",axes=ax2,lw=0.6,color="green",label="_nolegend_")


		ax2.set_ylabel("(Dec.+Pro.)/Int.")
		ax2.grid(alpha=0.6)

		ax1.set_xlim(lower_range,upper_range)
		ax2.set_xlim(lower_range,upper_range)
		ax2.set_ylim(0.5,1.5)

		if par == "Photon" and var=="PT":
			ax1.set_ylim(ymin=0)
			ax2.set_ylim(0.75,1.25)

		labelkinax(ax1,ax2,var,par)
		leg = ax1.legend()#title="KS Test: "+"{:.7f}".format(ksval),fontsize="small")
		#leg.get_title().set_fontsize('small')



		plt.savefig("plots/"+par+"_truth/"+"decproint"+"_"+par+"_"+var+".pdf",bbox_inches='tight')
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
RMPartP=["Photon","TopQuark","BQuark","WBoson","UQuark"]

up_l=0.01

for va in var:
	for p1 in RMPart:
		for p2 in RMPartP:
			if(p1==p2):
				continue

			para = parameters_RM(p1, p2, va)
			lower_range=para[0];upper_range=para[1];nbin=para[2];up_l=para[3]

			evi = np.genfromtxt("data/int_"+p1+"_"+p2+"_"+va+"_truth_1.txt")
			evd = np.genfromtxt("data/dec_"+p1+"_"+p2+"_"+va+"_truth_1.txt")
			evp = np.genfromtxt("data/pro_"+p1+"_"+p2+"_"+va+"_truth_1.txt")
			evi = evi[(np.abs(evi)>up_l) & (evi!=999.9)];evd = evd[(np.abs(evd)>up_l) & (evd!=999.9)];evp = evp[(np.abs(evp)>up_l) & (evp!=999.9)]
			evi = np.clip(evi, lower_range+eps, upper_range-eps);evd = np.clip(evd, lower_range+eps, upper_range-eps);evp = np.clip(evp, lower_range+eps, upper_range-eps)

			binning = np.linspace(lower_range, upper_range, nbin)
			hd = rplot.Hist(binning)
			hp = rplot.Hist(binning)

			for sample_number in sample_numbers:
				evi = np.genfromtxt("data/"+"int"+"_"+p1+"_"+p2+"_"+va+"_truth_"+str(sample_number)+".txt")
				evd = np.genfromtxt("data/"+"dec"+"_"+p1+"_"+p2+"_"+va+"_truth_"+str(sample_number)+".txt")
				evp = np.genfromtxt("data/"+"pro"+"_"+p1+"_"+p2+"_"+va+"_truth_"+str(sample_number)+".txt")

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
			rplt.hist(hp,fmt="none",lw=0.4,color="blue",label="production sample")
			rplt.errorbar(hp,fmt="none",lw=0.4,color="blue",label="_nolegend_")

			hd.Scale(1/(hd.Integral(0,hd.GetNbinsX()+1)))
			hd.linecolor = "red";hd.linewidth = 1
			rplt.hist(hd,fmt="none",lw=0.4,color="red",label="decay sample")
			rplt.errorbar(hd,fmt="none",lw=0.4,color="red",label="_nolegend_")

			plt.xlim(lower_range,upper_range)
			labelRM(va,p1,p2)
			
			plt.savefig("plots/"+va+"_truth/"+"decpro"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')
			plt.close()

			gs = gridspec.GridSpec(2, 1, width_ratios=[1],height_ratios=[3,1]) 
			f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
			ax1 = plt.subplot(gs[0])
			ax2 = plt.subplot(gs[1])
			gs.update(wspace=0.03)
			#nbin=128
			#binning = np.linspace(lower_range, upper_range, nbin)
			binning = set_dyn_binning(evi, lower_range, upper_range, nbin)
			hi		= rplot.Hist(binning)
			hdp		= rplot.Hist(binning)
			hisys	= rplot.Hist(binning)
			hd = rplot.Hist(binning)
			hp = rplot.Hist(binning)


			for sample_number in sample_numbers:
				evi = np.genfromtxt("data/"+"int"+"_"+p1+"_"+p2+"_"+va+"_truth_"+str(sample_number)+".txt")
				evd = np.genfromtxt("data/"+"dec"+"_"+p1+"_"+p2+"_"+va+"_truth_"+str(sample_number)+".txt")
				evp = np.genfromtxt("data/"+"pro"+"_"+p1+"_"+p2+"_"+va+"_truth_"+str(sample_number)+".txt")

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
				map(hisys.Fill, evi)

			hp.Scale(1/(hp.Integral(0,hp.GetNbinsX()+1))*prod_fraction,"width")
			hp.linecolor = "blue";hp.linewidth = 1
			rplt.hist(hp,fmt="none",axes=ax1,lw=0.4,color="blue",label="production sample")
			rplt.errorbar(hp,fmt="none",axes=ax1,lw=0.4,color="blue",label="_nolegend_")

			hd.Scale(1/(hd.Integral(0,hd.GetNbinsX()+1))*decay_fraction,"width")
			hd.linecolor = "red";hd.linewidth = 1
			rplt.hist(hd,fmt="none",axes=ax1,lw=0.4,color="red",label="decay sample")
			rplt.errorbar(hd,fmt="none",axes=ax1,lw=0.4,color="red",label="_nolegend_")

			hi.Scale(1/(hi.Integral(0,hi.GetNbinsX()+1)),"width")
			hi.linecolor = "black";hi.linewidth = 1
			rplt.hist(hi,fmt="none",axes=ax1,lw=0.4,color="black",label="interference sample")
			rplt.errorbar(hi,fmt="none",axes=ax1,lw=0.4,color="black",label="_nolegend_")

			hdp.Scale(1/(hdp.Integral(0,hdp.GetNbinsX()+1)),"width")
			hdp.linecolor = "green";hdp.linewidth = 1
			rplt.hist(hdp,fmt="none",axes=ax1,lw=0.4,color="green",label="decay + production sample")
			rplt.errorbar(hdp,fmt="none",axes=ax1,lw=0.4,color="green",label="_nolegend_")



			f = open("../scalescan/scalevarexp/scale_"+p1+"_"+p2+"_"+va+"_truth.txt")
			ledge = []
			uperr = []
			downerr = []

			while True:
				line = f.readline()
				if line == '':
					break
				lsp = line.split()
				ledge.append(float(lsp[0]))
				line = f.readline()
				lsp = line.split()
				uperr.append(float(lsp[0]))
				downerr.append(float(lsp[1]))

			errbin = 0
			for i in range(0,hisys.GetNbinsX()+1):
				while hisys.GetBinCenter(i) > ledge[errbin+1]:
					errbin += 1
				cont = hisys.GetBinContent(i)
				errup = 	np.fabs((uperr[errbin]*cont)-cont)
				errdown = 	np.fabs((downerr[errbin]*cont)-cont)

				if errup > errdown:
					syserr = errup
				else:
					syserr = errdown

				hisys.SetBinError(i,hisys.GetBinError(i)+syserr)

				cont = hi.GetBinContent(i)
				errup = np.fabs(uperr[errbin]-1)*cont
				errdown = np.fabs(downerr[errbin]-1)*cont
				#errup = np.fabs(uperr[errbin]-1)+hi.GetBinError(i)
				#errdown = np.fabs(downerr[errbin]-1)+hi.GetBinError(i)

				if errup > errdown:
					syserr = errup
				else:
					syserr = errdown

				y1pl = [cont+hi.GetBinError(i)+syserr,cont+hi.GetBinError(i)+syserr]
				y2pl = [cont-hi.GetBinError(i)-syserr,cont-hi.GetBinError(i)-syserr]
				xpl  = [hisys.GetBinLowEdge(i),hisys.GetBinLowEdge(i)+(hisys.GetBinCenter(i)-hisys.GetBinLowEdge(i))*2]
				if i == 0:
					ax1.fill_between(xpl,y1pl,y2pl,color="magenta",alpha=0.5,label="scale variation uncertainty",lw=0)
					continue
				ax1.fill_between(xpl,y1pl,y2pl,color="magenta",alpha=0.5,lw=0)
				
			ksval = hisys.KolmogorovTest(hdp)

			hisys.Scale(1/(hisys.Integral(0,hisys.GetNbinsX()+1)),"width")
			hisys.Divide(hi)
			for i in range(0,hisys.GetNbinsX()+1):
				cont = hisys.GetBinContent(i)
				xpl  = [hisys.GetBinLowEdge(i),hisys.GetBinLowEdge(i)+(hisys.GetBinCenter(i)-hisys.GetBinLowEdge(i))*2]
				y1pl = [cont+hisys.GetBinError(i),cont+hisys.GetBinError(i)]
				y2pl = [cont-hisys.GetBinError(i),cont-hisys.GetBinError(i)]
				ax2.fill_between(xpl,y1pl,y2pl,color="magenta",alpha=0.5,lw=0)



			hisys.linecolor = "magenta";hisys.linewidth = 1


			print(p1 + " " + p2 + " " + va + " KS Test")
			ksval = hi.KolmogorovTest(hdp)
			print(hi.KolmogorovTest(hdp))
			print()

			ax2.hist(binning[1:]-np.diff(binning)/2,bins=binning,weights=np.ones_like(binning[1:]),histtype="step",lw=0.8,color="black")
			hdp.Divide(hi)
			rplt.errorbar(hdp,fmt="none",axes=ax2,lw=0.6,color="green",label="_nolegend_")
			#hi.Divide(hi)
			#rplt.hist(hi,fmt="none",axes=ax2,lw=0.6,color="black",label="_nolegend_")


			ax2.set_ylabel("(Dec.+Pro.)/Int.")
			#ax2.legend(loc="best")
			ax2.grid(alpha=0.6)

			ax1.set_xlim(lower_range,upper_range)
			ax2.set_xlim(lower_range,upper_range)
			ax2.set_ylim(0.5,1.5)

			if p1 == "Photon":
				if p2 == "UQuark":
					if va == "M":
						ax2.set_ylim(0,3)


			labelRMax(ax1,ax2,va,p1,p2)
			leg = ax1.legend()#title="KS Test: "+"{:.7f}".format(ksval),fontsize="small")
			#leg.get_title().set_fontsize('small')
			
			plt.savefig("plots/"+va+"_truth/"+"decproint_"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')
			plt.close()
			hdp.write("tua_ratio_"+va+"_"+p1+"_"+p2)




