
PlotXe = 1
LineWidth = 2
FontSize = 18

File_1 = 'tmp1.dat'
File_2 = 'tmp2.dat'
X_limit = [10, 3000]

import matplotlib.pyplot as plt
import numpy as np

d1 = np.genfromtxt(File_1)
d2 = np.genfromtxt(File_2)
zp1 = 1 + d1[:,0]
zp2 = 1 + d2[:,0]
if PlotXe == 1:
	y1 = d1[:,1]
	y2 = d2[:,1]
else:
	y1 = d1[:,2]
	y2 = d2[:,2]

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
plt.plot(zp1,y1,'k',linewidth=LineWidth,label = 'd1')
plt.plot(zp2,y2,'--b',linewidth=LineWidth,label = 'd2')
plt.xlabel('$1+z$',fontsize=FontSize,fontname='Times New Roman')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.legend(['d1','d2'],fontsize=FontSize,loc = 'lower right')
plt.xlim(X_limit)

if PlotXe == 1:
	plt.ylabel('$x_{\\rm{e}}$',fontsize=FontSize,fontname='Times New Roman')
	plt.xscale('log')
	plt.yscale('log')
else:
	plt.ylabel('$T$',fontsize=FontSize,fontname='Times New Roman')
	plt.xscale('log')
	plt.yscale('log')

plt.savefig('/Users/cangtao/Desktop/tmp.png',bbox_inches='tight')
plt.show()
