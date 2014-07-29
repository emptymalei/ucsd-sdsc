import sys
sys.path.append('/opt/scipy/lib/python2.7/site-packages')
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import itertools as iter
from matplotlib import rc

rc('text', usetex=True)

f1 = open('ye.dat', 'r')
f1_lines = f1.readlines()
f1.close

obh2 = f1_lines[1][20:29]
f1_lines = f1_lines[4:-1]
f1_floats = []

for line in f1_lines:
	list_floats = []
	for i in iter.islice(iter.count(0),3):
		list_floats = list_floats + [float(line[3+12*i:12*(i+1)])]
	f1_floats = f1_floats + [list_floats]

size = len(f1_floats)
f1_elems = [[],[],[]]
for i in iter.islice(iter.count(0),3):
	for j in iter.islice(iter.count(0),size):
		f1_elems[i] = f1_elems[i] + [f1_floats[j][i]]

large = f1_elems[0][0]
small = f1_elems[0][-1]

fig = plt.figure()
plt.xlim(large,small)
#plt.ylim(0.218,0.255)
line = plt.semilogx(f1_elems[0],f1_elems[1],'b-'
	,f1_elems[0],f1_elems[2],'g-')
plt.ylabel(r'$Y_e$')
plt.xlabel(r'$T_\gamma\,(\mathrm{MeV})$')
plt.legend([r'$\mathrm{Equil}$',r'$\mathrm{Actual}$'],loc='upper left')
titlestr = r'$\mathrm{Electron\,Fraction\,vs.\,Plasma\,Temperature\,}(\Omega_bh^2=$' + obh2 + r'$)$'
fig.suptitle(titlestr)

plt.savefig('plot_ye_v_temp.pdf')
plt.clf()
