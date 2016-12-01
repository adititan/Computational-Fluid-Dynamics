import matplotlib.pyplot as plt
import numpy as np
import copy
from copy import copy, deepcopy

def error(T,T_old,alpha,dTc,L):
	error = L*L/(alpha*dTc)*max(abs(T_old[j][i]-T[j][i]) for i in range(1,imax-1) for j in range(1,jmax-1))
	return error

def twod_cond(T_init,qgen):
	T = deepcopy(T_init)
	T_old = [[0.0 for x in range(imax)] for y in range(jmax)]
	while(error(T,T_old,alpha,dTc,L1)>epsilon):
		T_old = deepcopy(T)
		qx =  [[0.0 for x in range(imax)] for y in range(jmax)]
		qy =  [[0.0 for x in range(imax)] for y in range(jmax)]	
		Qcond =[[0.0 for x in range(imax)] for y in range(jmax)]
		for j in range(1,jmax-1):
			for i in range(0,imax-1):
				if(i == 0 or i == imax-2):
					qx[j][i] = -1* (T_old[j][i+1] - T_old[j][i])/(1/2.0)
				else:
					qx[j][i] = -1* (T_old[j][i+1] - T_old[j][i])
		for j in range(0,jmax-1):
			for i in range(1,imax-1):
				if(j==0 or j==jmax-2):
					qy[j][i] = -1*(T_old[j+1][i] -  T_old[j][i])/(1/2.0)
				else:
					qy[j][i] = -1*(T_old[j+1][i] -  T_old[j][i])/1
		for i in range(1,imax-1):
			for j in range(1,imax-1):
				Qcond[j][i] = (qx[j][i-1] -qx[j][i]) + (qy[j-1][i] - qy[j][i]) 
				Qgen = qgen*dx*dy
				T[j][i] = T_old[j][i] + dt* (Qcond[j][i] + Qgen)
	return T

def plot_fig(T):
	plt.figure()
	plt.clf()
	plt.axes([0.025, 0.025, 0.95, 0.95])
	plt.contourf(T , 40, alpha=.75, cmap=plt.cm.hot)
	C = plt.contour(T, 40, colors='black', linewidth=.5)
	plt.gcf().subplots_adjust(top =0.10,bottom=0.05)
	plt.clabel(C, inline=1, fontsize=10)
	plt.xlabel('x(m)')
	plt.ylabel('y(m)')
	plt.xticks([0.0,0.5, 1.0],[r'$0$', r'$0.5$',r'$1.0$'])
	plt.yticks([0.5,1.0],[r'$0.5$', r'$1.0$'])

if __name__ == '__main__':
	T0=30.0
	imax, jmax =12,12
	L1, L2 = 1.0,1.0
	k = 16.2
	epsilon = 0.1
	rho, cp = 7750.0, 500.0
	alpha= k/(rho*cp)
	T_wb ,T_eb = 100.0, 300.0	#BC's
	T_nb , T_sb = 200.0,400.0
	dx = L1/(imax-2)
	dy = L2/(jmax-2)
	dt = 0.001#0.99*0.5* (dx*dx/alpha)
	dTc = T_nb - T_wb
	xi = np.linspace(0,L1,imax-1)
	yi = np.linspace(0,L2,imax-1)
	xc = [ 0.0 for x in range(imax)]
	for i in range(1,imax-1):
		xc[i] = (xi[i]+ xi[i-1])/2.0
	xc[0] = xi[0]
	xc[imax-1] = xi[imax-2]
	yc = [ 0.0 for y in range(jmax)]
	for i in range(1,imax-1):
		yc[i] = (yi[i]+ yi[i-1])/2.0
	yc[0] = yi[0]
	yc[imax-1] = yi[imax-2]
	T_init = [[T0 for x in range(imax)] for y in range(jmax)]
	for j in range(0,jmax):
		T_init[j][0] = T_wb
		T_init[j][imax-1] = T_eb
	for i in range(0,imax):
		T_init[jmax-1][i] = T_sb
		T_init[0][i] = T_nb
	T1 = twod_cond(T_init,0.0)
	plot_fig(T1)
	plt.title('2-D Heat Conduction with Qgen = 0 W/m^3')
	plt.savefig('seven.png')
	T2 = twod_cond(T_init,50000.0)
	plot_fig(T2)
	plt.title('2-D Heat Conduction with Qgen = 50 KW/m^3')
	plt.savefig('eight.png')
