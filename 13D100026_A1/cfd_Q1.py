import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy

def exact_solution1(Qgen):
	if(Qgen == 0.0):
		T = (T_eb-T_wb)*x/0.01+T_wb
	elif(Qgen > 0.0):
		T = -Qgen/(2*k)*(x**2) + (T_eb - T_wb + (Qgen/(2*k)*L*L))*x/L + T_wb
	plt.plot(x,T, color="blue", linewidth=1.0, linestyle="-",label = "exact")	 

def error(T,T_old,alpha,dTc,L):
	error = L*L/(alpha*dTc)*max(abs(T_old[i]-T[i]) for i in range(1,imax-1))
	return error

def oned_cond(Qgen,T_init):
	T = deepcopy(T_init)
	T_old =  [0.0 for x in range(imax)]
	while(error(T,T_old,alpha,dTc,L)>epsilon):
		T_old = deepcopy(T)
		qx =  [0.0 for x in range(imax-1)]
		Qcond = [0.0 for x in range(imax-1)]
		for i in range(0,imax-1):
			if (i == 0 or i == imax-2):
				qx[i] = k*(T_old[i]- T_old[i+1])/(dx/2)
			else:
				qx[i] = k*(T_old[i]- T_old[i+1])/dx
		for i in range(1,imax-1):
			Qcond[i] =(qx[i-1] - qx[i])
			T[i] = T_old[i] + dt*( Qcond[i]+ Qgen*dx)/(rho*cp*dx)
	return T

def plot_function(T,Qgen):
	plt.figure()
	plt.clf()
	plt.axes([0.07, 0.07, 0.95, 0.88])
	plt.plot(x,T,color="red", linewidth=1.0, linestyle="-", label ="Approximate")
	exact_solution1(Qgen)
	plt.gcf().subplots_adjust(top =0.10,bottom=0.05)
	plt.xlabel('x')
	plt.ylabel('Temperature(Celcius)')
	plt.legend(loc='upper right')
	plt.xticks([0, 0.01],[r'$0$', r'$0.01$'])
	plt.yticks([0.0, 100.0],[r'$0.0$', r'$100.0$'])

if __name__ == '__main__':
	T0=30.0
	imax =12
	L=0.01
	dx = L/(imax-2)
	T_init = [T0 for x in range(imax)]
	T_wb ,T_eb = 0.0 , 100.0
	T_init[0] = T_wb
	T_init[imax-1] = T_eb
	dTc = T_eb - T_wb
	k = 16.2
	epsilon = 0.0001
	rho, cp = 7750.0, 500.0
	alpha= k/(rho*cp)
	dt =0.99*0.5*dx*dx/alpha
	L = 0.01
	x=np.linspace(0,L,imax)
	T1 =oned_cond(0.0,T_init)
	plot_function(T1,0.0)
	plt.title('1-D Heat Conduction with Qgen = 0 W/m^3')
	plt.savefig('1d_1.png')
	T2 = oned_cond(100000000.0,T_init)
	plot_function(T2,100000000.0)	
	plt.title('1-D Heat Conduction with Qgen = 100 MW/m^3')
	plt.savefig('1d_2.png')
