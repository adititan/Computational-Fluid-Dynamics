import matplotlib.pyplot as plt
import numpy as np
import copy
from copy import copy, deepcopy

def error(T,T_old):
	sum0 = 0.0
	for j in range(0,jmax):
		for i in range(0,imax):
			sum0 += (T_old[j][i] - T[j][i])**2
	error = np.sqrt(sum0/((imax-2)*(jmax-2)))
	print("a2",error)
	return error

def temp_at_centreline(Temp):
	temp = [0.0]*jmax
	for j in range(jmax):	
			temp[j] = (Temp[j][16] + Temp[j][15]) /2.0
	return temp

def flux_based_adv(T_init,func):
	T = deepcopy(T_init)
	T_old = [[0.0 for x in range(imax)] for y in range(jmax)]
	if func == 'fou':
		w1,w2,w3 = 0.0,1.0,0.0
	elif func == 'sou':
		w1,w2,w3 = 0.0,1.5,-0.5
	elif func == 'quick':
		w1,w2,w3 = 3.0/8.0, 6.0/8.0 ,-1.0/8.0
	hx =  [[0.0 for x in range(imax)] for y in range(jmax)]
	hy =  [[0.0 for x in range(imax)] for y in range(jmax)]	
	Qadv =[[0.0 for x in range(imax)] for y in range(jmax)]
	while(error(T,T_old)>epsilon):
		T_old = deepcopy(T)
		Te = deepcopy(T)
		Tn = deepcopy(T)
		for j in range(1,jmax-1):
			for i in range(1,imax-2):
				Te[j][i] =  w1*T[j][i+1] +w2*T[j][i] + w3*T[j][i-1]
		for i in range(1,imax-1):
			Tn[jmax-2][i] = T[jmax-1][i] 
		for j in range(0,jmax-2):
			for i in range(1,imax-1):
				Tn[j][i] = w1*T[j][i] + w2*T[j+1][i] + w3*T[j+2][i] 
		for j in range(1,jmax-1):
			for i in range(0,imax-1):
				hx[j][i] = u*Te[j][i]
		for i in range(1,imax-1):
			for j in range(0,jmax-1):
				hy[j][i] = v*Tn[j][i]
		for j in range(1,jmax-1):
			for i in range(1,imax-1):
				Qadv[j][i] = (hx[j][i] - hx[j][i-1])*dy + (hy[j-1][i] - hy[j][i])*dx
				T[j][i] = T_old[j][i] - (dt/(dx*dy))*Qadv[j][i]
		for i in range(imax):
			T[0][i] = T[1][i]
		for j in range(jmax):	
			T[j][imax-1] = T[j][imax-2]
	phi = [[0.0 for x in range(imax)] for y in range(jmax)]
	for i in range(0,imax):
			for j in range(0,jmax):
				phi[j][i] = (T[j][i] -Tb) / (Tl - Tb)
	return phi

def plot_fig(T):
	plt.figure()
	plt.clf()
	plt.axes([0.060, 0.060, 0.95, 0.90])
	plt.imshow(T)
	plt.colorbar(shrink=.92)

if __name__ == '__main__':
	L= 1.0
	H =1.0
	u =1.0
	v =1.0
	imax =32
	jmax =32
	T0 = 0.0
	Tb = 0.0
	Tl = 100.0
	dTc =100.0
	epsilon =0.0001
	dx = L/(imax-2)
	dy = H/(jmax-2)
	dt = 0.01*0.5*dx/u 
	T_init = [[T0 for x in range(imax)] for y in range(jmax)]
	for i in range(0,imax):
		T_init[jmax-1][i] = Tb
	for j in range(0,jmax):
		T_init[j][0] = Tl
	T1 = flux_based_adv(T_init,'fou')
	plot_fig(T1)
	plt.title('2-D Heat Advection using FOU')
	plt.savefig('q1_fou.png')
	T2 = flux_based_adv(T_init,'sou')
	plot_fig(T2)
	plt.title('2-D Heat Advection using SOU')
	plt.savefig('q1_sou.png')
	T3 = flux_based_adv(T_init,'quick')
	plot_fig(T3)
	plt.title('2-D Heat Advection using QUICK')	
	plt.savefig('q1_quick.png')
	fou = temp_at_centreline(T1)
	sou = temp_at_centreline(T2)
	quick = temp_at_centreline(T3)
	y = np.linspace(H,0,jmax)
	plt.figure()
	plt.plot(fou,y)
	plt.title('Variation of Nondimensional Temperature along Vertical Centerline(FOU)')
	plt.savefig('q2_1.png')
	plt.figure()
	plt.plot(sou,y)
	plt.title('Variation of Nondimensional Temperature along Vertical Centerline(SOU)')
	plt.savefig('q2_2.png')
	plt.figure()
	plt.plot(quick,y)
	plt.title('Variation of Nondimensional Temperature along Vertical Centerline(QUICK)')
	plt.savefig('q2_3.png')
