import matplotlib.pyplot as plt
import numpy as np
import copy
from copy import copy, deepcopy

def error(T,T_old,alpha,dTc,L):
	sum0 = 0.0
	for j in range(0,jmax):
		for i in range(0,imax):
			sum0 += (T_old[j][i] - T[j][i])**2
	error = np.sqrt(sum0/((imax-2)*(jmax-2)))
	print(error)
	return error

def temp_at_x(Temp):
	temp1 = [0.0]*jmax
	temp2 = [0.0]*jmax
	temp3 = [0.0]*jmax
	temp4 = [0.0]*jmax
	temp5 = [0.0]*jmax
	for j in range(jmax):	
			temp1[j] = Temp[j][12]
			temp2[j] = Temp[j][24]
			temp3[j] = Temp[j][36]
			temp4[j] = Temp[j][48]
			temp5[j] = Temp[j][61]
	return temp1,temp2,temp3,temp4,temp5

def twod_conv(T_init,func):
	T = deepcopy(T_init)
	T_old = [[0.0 for x in range(imax)] for y in range(jmax)]
	if func == 'fou':
		w1,w2,w3 = 0.0,1.0,0.0
	elif func == 'quick':
		w1,w2,w3 = 3.0/8.0, 6.0/8.0 ,-1.0/8.0
	while(error(T,T_old,alpha,dTc,L1)>epsilon):	
		T_old = deepcopy(T)
		Te = deepcopy(T)
		Tn = deepcopy(T)
		qx =  [[0.0 for x in range(imax)] for y in range(jmax)]
		qy =  [[0.0 for x in range(imax)] for y in range(jmax)]	
		Qcond =[[0.0 for x in range(imax)] for y in range(jmax)]
		hx =  [[0.0 for x in range(imax)] for y in range(jmax)]
		hy =  [[0.0 for x in range(imax)] for y in range(jmax)]	
		Qadv =[[0.0 for x in range(imax)] for y in range(jmax)]
		for j in range(1,jmax-1):
			for i in range(1,imax-2):
				Te[j][i] =  w1*T[j][i+1] +w2*T[j][i] + w3*T[j][i-1] 
		for j in range(1,jmax-2):
			for i in range(1,imax-1):
				Tn[j][i] = w1*T[j-1][i] + w2*T[j][i] + w3*T[j+1][i]  
		for j in range(1,jmax-1):
			for i in range(0,imax-1):
				if(i == 0 or i == imax-2):
					qx[j][i] = -k* (T_old[j][i+1] - T_old[j][i])/(dx/2.0)
				else:
					qx[j][i] = -k* (T_old[j][i+1] - T_old[j][i])/dx
		for j in range(0,jmax-1):
			for i in range(1,imax-1):
				if(j==0 or j==jmax-2):
					qy[j][i] = -k*(T_old[j+1][i] -  T_old[j][i])/(dy/2.0)
				else:
					qy[j][i] = -k*(T_old[j+1][i] -  T_old[j][i])/dy
		for j in range(1,jmax-1):
			for i in range(0,imax-1):
				hx[j][i] = rho*u*cp*Te[j][i]
		for i in range(1,imax-1):
			for j in range(0,jmax-1):
				hy[j][i] = rho*v*cp*Tn[j][i]
		for i in range(1,imax-1):
			for j in range(1,jmax-1):
				Qcond[j][i] = (qx[j][i-1] -qx[j][i])*dy + (qy[j-1][i] - qy[j][i])*dx 
				Qadv[j][i] = (hx[j][i] - hx[j][i-1]) *dy + (hy[j-1][i] - hy[j][i])*dx
				T[j][i] = T_old[j][i] + dt/(rho*cp*dx*dy)*( Qcond[j][i] - Qadv[j][i] )
				for j in range(0,jmax):		
					T[j][imax-1] = T[j][imax-2]
	phi = [[0.0 for x in range(imax)] for y in range(jmax)]
	for i in range(0,imax):
			for j in range(0,jmax):
				phi[j][i] = (T[j][i] -T_nb) / (T_wb - T_nb)
	return phi

def plot_fig(T):
	plt.figure()
	plt.clf()
	plt.imshow(T)
	plt.colorbar(shrink=.92)

if __name__ == '__main__':
	T0=0.0
	imax, jmax =62,22
	L1, L2 = 6.0,1.0
	k = 10
	epsilon = 0.0001
	rho, cp = 7.0, 100.0
	alpha= k/(rho*cp)
	T_wb  = 100.0	#BC's
	T_nb,T_sb  = 0.0 , 0.0
	dx = L1/(imax-2)
	dy = L2/(jmax-2)
	u = 1.0
	v = 0.0
	dt1 = 0.5/alpha*(1/(dx**2) + 1/(dy**2))
	dt2 = 0.5*dx/u
	dt = 0.01 *min(dt1,dt2)
	dTc = T_wb - T_nb
	T_init = [[T0 for x in range(imax)] for y in range(jmax)]
	for i in range(0,imax):
		T_init[0][i] = T_nb
		T_init[jmax-1][i] = T_sb
	for j in range(0,jmax):
		T_init[j][0] = T_wb
		T_init[j][imax-1] = T_init[j][imax-2]
	T2 = twod_conv(T_init,'fou')
	plot_fig(T2)
	plt.title('2-D Heat Convection using FOU')
	plt.savefig('Q3_fou.png')
	T1 = twod_conv(T_init,'quick')
	plot_fig(T1)
	plt.title('2-D Heat Convection using QUICK')
	plt.savefig('Q3_quick.png')
	y = np.linspace(L1,0 ,jmax)
	Ta,Tb,Tc,Td,Te = temp_at_x(T1)
	plt.figure()
	plt.title('Variation of Nondimensional Temperature along Vertical lines(FOU)')
	plt.plot(Ta,y,label = "x/H = 0.2")
	plt.plot(Tb,y,label = "x/H = 0.4")
	plt.plot(Tc,y,label = "x/H = 0.6")
	plt.plot(Td,y,label = "x/H = 0.8")
	plt.plot(Te,y,label = "x/H = 1.0")
	plt.legend(loc='upper right')
	plt.savefig('q4_1.png')		
	Ta1,Tb1,Tc1,Td1,Te1 = temp_at_x(T2)
	plt.figure()
	plt.title('Variation of Nondimensional Temperature along Vertical lines(QUICK)')
	plt.plot(Ta1,y,label = "x/H = 0.2")
	plt.plot(Tb1,y,label = "x/H = 0.4")
	plt.plot(Tc1,y,label = "x/H = 0.6")
	plt.plot(Td1,y,label = "x/H = 0.8")
	plt.plot(Te1,y,label = "x/H = 1.0")
	plt.legend(loc='upper right')
	plt.savefig('q4_2.png')			
