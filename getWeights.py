import numpy as np
from matplotlib import pyplot as plt
import os, inspect

PATH = "/run/user/1000/gvfs/sftp:host=arrakis,user=slacagni/media/slacagni/Elements/tua_onlyextra/crosssections/"

PREFIX = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/"

os.chdir(PATH)

XSECS = os.listdir(".")


for XSEC in XSECS:
	f = open(XSEC,"r")
	w = open(PREFIX+"cross/"+XSEC[:-5]+"_cross","w")
	two = False
	for line in f:
		if "Cross-section" in line:
			if two == True:
				lsp = line.split()
				w.write(lsp[2])
			two = True