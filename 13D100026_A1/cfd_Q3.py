import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy

def exact_solution1(Q_vol_gen):
	if (Q_vol_gen == 0):
		T = (h*T_inf/(L*h+k))*xi + T_wb
	elif(Q_vol_gen>0.0):
		T = -Q_vol_gen/(2*k)*(xi*xi) + ((T_inf - T_wb)*h/(h*L + k) + Q_vol_gen/(2*k)*L*(h*L + 2*k)/(h*L+k))*xi + T_wb
	plt.plot(xi,T, color="blue", linewidth=1.0, linestyle="-",label = "exact")

def error(T,T_old,alpha,dTc,L):
	error = L*L/(alpha*dTc)*max(abs(T_old[i]-T[i]) for i in range(1,imax-1))
	return error

def one_d(T_init,Q_vol_gen):
	Q_gen = [ 0.0 for x in range(imax-1)]	
	aE = [ 0.0 for x in range(imax-1)]
	aP0 = [ 0.0 for x in range(imax-1)]
	aP = [ 0.0 for x in range(imax-1)]
	for i in range(0,imax-1):
	 	Q_gen[i]=Q_vol_gen*dx[i]
	for i in range(0,imax-1):	
	     aE[i]=k/dx[i]
	for i in range(0,imax-1):
		aP0[i]=rho*cp*dx[i]/Dt
		aP[i]=aP0[i]+aE[i]+aE[i-1]
	DTc=T_inf-T_wb
	T = deepcopy(T_init)
	T_old =  [0.0 for x in range(imax)]
	while(error(T,T_old,alpha,DTc,L)>epsilon_st):
		T[imax-1] = ((k*T[imax-2])+(h*dx[imax-2]*T_inf/2)) /(k+h*dx[imax-2]/2)
		T_old = deepcopy(T)
		b =  [0.0 for x in range(imax-1)]
		for i in range(1,imax-1):
			b[i]=aP0[i]*T_old[i]+Q_gen[i]
		epsilon=0.0001
		Error=1
		while (Error>=epsilon):
			T_old_iter=deepcopy(T)
			for i in range(1,imax-1):
				T[i]=(aE[i]*T[i+1]+aE[i-1]*T[i-1]+b[i])/aP[i]
			Error=max(abs(T[i]-T_old_iter[i]) for i in range(1,imax-1))
	return T

def plot_fig(T,Q_vol_gen):
	plt.figure()
	plt.clf()
	plt.axes([0.07, 0.07, 0.95, 0.88])
	plt.plot(xc,T,color="red", linewidth=1.0, linestyle="-", label ="Approximate")
	exact_solution1(Q_vol_gen)
	plt.gcf().subplots_adjust(top =0.10,bottom=0.05)
	plt.xlabel('x')
	plt.ylabel('Temperature(Celcius)')
	plt.legend(loc='lower right')
	plt.xticks([0, 0.01],[r'$0$', r'$0.01$'])
	plt.yticks([0.0, 100.0],[r'$0.0$', r'$100.0$'])
		


if __name__ == '__main__':
	T0=30.0
	imax =12
	L=0.01
	psi = np.linspace(0,1,imax-1)
	beta=1.2
	xi = np.linspace(0,L,imax-1)
	for i in range(0,imax-1):
		a = ((beta+1)/(beta-1))**(2*psi[i]-1)
		xi[i]=L*((beta + 1)*(a)-(beta-1))/(2*(1+a))
	xc = [ 0.0 for x in range(imax)]
	for i in range(1,imax-1):
		xc[i] = (xi[i]+ xi[i-1])/2.0
	xc[0] = xi[0]
	xc[imax-1] = xi[imax-2]
	dx = [ 0.0 for x in range(imax)]
	for i in range(0,imax-1):
		dx[i] = xc[i+1]-xc[i]
	k = 16.2
	epsilon = 0.0001
	rho, cp = 7750.0, 500.0
	alpha= k/(rho*cp)
	T0=30.0
	T_wb=0.0
	T_inf=100.0
	h=1000.0
	epsilon_st=0.0001 
	Dt=1
	T_init = [T0 for x in range(imax)]
	T_init[0] = T_wb
	T_init[imax-1] = ((k*T_init[imax-2])+(h*dx[imax-2]*T_inf/2)) /(k+h*dx[imax-2]/2)
	T1 = one_d(T_init,0.0)
	plot_fig(T1,0.0)
	plt.title('1-D Heat Conduction with Qgen = 0 W/m^3')
	plt.savefig('three.png')
	T2 = one_d(T_init,100000000.0)
	plot_fig(T2,100000000.0)
	plt.title('1-D Heat Conduction with Qgen = 100 MW/m^3')
	plt.savefig('four.png')
