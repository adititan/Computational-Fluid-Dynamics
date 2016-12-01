import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy

def error(T,T_old,alpha,dTc,L):
	error = 0.01*0.01/(alpha*dTc)*max(abs(T_old[j][i]-T[j][i]) for i in range(1,imax-1) for j in range(1,jmax-1))
	print(error)	
	return error

def two_d(T_init,Q_vol_gen):
	Q_gen = [[0.0 for x in range(imax-1)] for y in range(jmax-1)]
	aE = [[0.0 for x in range(imax-1)] for y in range(jmax-1)]
	aN = [[0.0 for x in range(imax-1)] for y in range(jmax-1)]
	aP0 = [[0.0 for x in range(imax-1)] for y in range(jmax-1)]
	aP = [[0.0 for x in range(imax-1)] for y in range(jmax-1)]
	for i in range(0,imax-1):
		for j in range(0,imax-1):
			Q_gen[j][i]=Q_vol_gen*dx[i]*dy[j]	
			aE[j][i] = k*dy[j]/dx[i]
			aN[j][i] = k*dx[i]/dy[j]
			aP0[j][i]=rho*cp*dx[i]*dy[j]/Dt
			aP[j][i]=aP0[j][i]+aE[j][i]+aE[j][i-1] +aN[j][i] +aN[j-1][i]
 	DTc=T_wb-T_inf
	T = deepcopy(T_init)
	T_old =  [[0.0 for x in range(imax)] for y in range(jmax)]
	while(error(T,T_old,alpha,DTc,L1)>epsilon_st):
		for j in range(0,jmax):
			T[j][0] = T_wb
			T[j][imax-1] = qw*dx[imax-2]/k + T[j][imax-2] 
		for i in range(0,imax):	
			T[jmax-1][i] = ((k*T[jmax-2][i])+(h*dy[jmax-2]*T_inf/2)) /(k+h*dy[jmax-2]/2)
			T[0][i] = T[1][i]
		T_old = deepcopy(T)
		b =  [[0.0 for x in range(imax-1)] for y in range(jmax-1)]
		for i in range(0,imax-1):
			for j in range(0,imax-1):
				b[j][i]=aP0[j][i]*T_old[j][i]+Q_gen[j][i]
		Error=1
		while (Error>=epsilon):
			T_old_iter=deepcopy(T)
			for i in range(1,imax-1):
				for j in range(0,imax-1):
					T[j][i] = (aE[j][i]*T[j][i+1] + aE[j][i-1]*T[j][i-1] + aN[j-1][i]*T[j-1][i] + aN[j][i]*T[j+1][i] + b[j][i])/aP[j][i]
			Error=max(abs(T[j][i]-T_old_iter[j][i]) for i in range(1,imax-1) for j in range(1,jmax-1))
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
	epsilon = 0.0001
	rho, cp = 7750.0, 500.0
	alpha= k/(rho*cp)
	T_wb=100.0
	T_inf=30.0
	h=100.0
	qw = 10000.0
	epsilon_st=0.0001 
	Dt=1
	psi = np.linspace(0,1,imax-1)
	beta=1.2
	xi = np.linspace(0,L1,imax-1)
	for i in range(0,imax-1):
		a = ((beta+1)/(beta-1))**(2*psi[i]-1)
		xi[i]=L1*((beta + 1)*(a)-(beta-1))/(2*(1+a))
	xc = [ 0.0 for x in range(imax)]
	for i in range(1,imax-1):
		xc[i] = (xi[i]+ xi[i-1])/2.0
	xc[0] = xi[0]
	xc[imax-1] = xi[imax-2]
	dx = [ 0.0 for x in range(imax-1)]
	for i in range(0,imax-1):
		dx[i] = xc[i+1]-xc[i]
	dy = [ 0.0 for x in range(imax-1)]
	dy = deepcopy(dx)
	T_init = [[T0 for x in range(imax)] for y in range(jmax)]
	for j in range(0,jmax):
		T_init[j][0] = T_wb
		T_init[j][imax-1] = qw*dx[imax-2]/k + T_init[j][imax-2] 
	for i in range(0,imax):	
		T_init[jmax-1][i] = ((k*T_init[jmax-2][i])+(h*dy[jmax-2]*T_inf/2)) /(k+h*dy[jmax-2]/2)
		T_init[0][i] = T_init[1][i]
	T1 = two_d(T_init,0.0)
	plot_fig(T1)
	plt.title('2-D Heat Conduction with Qgen = 0 W/m^3')
	plt.savefig('five.png')
	T2 = two_d(T_init,50000.0)	
	plot_fig(T2)
	plt.title('2-D Heat Conduction with Qgen = 50 KW/m^3')
	plt.savefig('six.png')
