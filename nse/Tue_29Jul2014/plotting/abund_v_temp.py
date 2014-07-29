import sys
sys.path.append('/opt/scipy/lib/python2.7/site-packages')
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import itertools as iter
import numpy
from matplotlib import rc
rc('text', usetex=True)

f1 = open('bbn.dat', 'r')
f1_lines = f1.readlines()

f1.close()

obh2 = f1_lines[0][20:29]

fltidx = 8 + int(f1_lines[5][5:11])

f1_lines = f1_lines[9:fltidx]
f1_floats = []
T1_floats = []
Yp_floats = []

for line in f1_lines:
	#print line
	T1_floats+=[float(line[1:10])]
	list_floats=[float(line[13:22])]
	Yp_floats = Yp_floats + [float(line[25:34])]
	for i in iter.islice(iter.count(0),7):
		#print line[37+12*i:47+12*i]
		list_floats=list_floats+[float(line[37+12*i:47+12*i])]
	f1_floats = f1_floats+[list_floats]

size1 = len(f1_floats)
f1_elems = [[],[],[],[],[],[],[],[],[]]

for i in iter.islice(iter.count(0),8):
	for j in iter.islice(iter.count(0),size1):
		f1_elems[i]=f1_elems[i]+[f1_floats[j][i]]

for i in iter.islice(iter.count(0),size1):
	f1_elems[4][i] = f1_elems[4][i]/4.0/Yp_floats[i]

f2 = open('equil.dat', 'r')
f2_lines = f2.readlines()

f2.close()

f2_lines = f2_lines[4:-1]
f2_floats = []
T2_floats = []

for line in f2_lines:
	T2_floats+=[float(line[3:12])]
	list_floats = []
	for i in iter.islice(iter.count(0),4):
		list_floats=list_floats+[float(line[15+12*i:24+12*i])]
	f2_floats = f2_floats+[list_floats]


size2 = len(f2_floats)
f2_elems = [[],[],[],[],[]]

for i in iter.islice(iter.count(0),4):
	for j in iter.islice(iter.count(0),size2):
		f2_elems[i]=f2_elems[i]+[f2_floats[j][i]]

maxT = T1_floats[0]
minT = T1_floats[-1]

fig = plt.figure()
#ax = fig.add_plot(rect)
#plt.xlim(maxT,minT)
plt.xlim(30.0,0.008)
plt.ylim(1.0e-24,1.0)
#plt.ylim(1.0e-35,1.0)
plt.loglog(T1_floats,f1_elems[0],'b',T1_floats,f1_elems[1],'g',
	T1_floats,f1_elems[2],'r',T1_floats,f1_elems[3],'c',
	T1_floats,f1_elems[4],'m',T1_floats,f1_elems[5],'y',
	T1_floats,f1_elems[6],'k',T1_floats,f1_elems[7],'0.5')
plt.loglog(T2_floats,f2_elems[0],'g--',T2_floats,f2_elems[3],'m--',
	T2_floats,f2_elems[1],'r--',T2_floats,f2_elems[2],'c--')
#plt.loglog(T2_floats,f2_elems[4],'0.5',linestyle='--')
plt.ylabel(r'$Y_i/Y_H$')
plt.xlabel('$T_\gamma\,(\mathrm{MeV})$')
titlestr1 = r'$\mathrm{Relative\,abundances\,wrt\,}Y_H\mathrm{\,vs.\,Plasma\,Temp.\,}$' + \
	r'$(\Omega_bh^2 = $' + obh2 + r'$)$'
fig.suptitle(titlestr1,ha='center',fontsize=14)
plt.legend((r'$\mathrm{N}$',r'$\mathrm{D}$',r'$\mathrm{T}$',r'$^3\mathrm{He}$', \
	r'$^4\mathrm{He}$',r'$^6\mathrm{Li}$',r'$^7\mathrm{Li}$',r'$^7\mathrm{Be}$'),
	'upper left')
plt.savefig('plot_abund_v_temp.pdf')


