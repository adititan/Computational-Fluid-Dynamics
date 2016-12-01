import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy

def error(T,T_old):
	sum0 = 0.0
	for j in range(0,jmax):
		for i in range(0,imax):
			sum0 += (T_old[j][i] - T[j][i])**2
	error = np.sqrt(sum0/((imax-2)*(jmax-2)))
	print("a3",error)
	return error

def temp_at_centreline(Temp):
	temp = [0.0]*jmax
	for j in range(jmax):	
			temp[j] = (Temp[j][16] + Temp[j][15]) /2.0
	return temp

def plot_fig(T):
	plt.figure()
	plt.clf()
	plt.imshow(T)
	plt.colorbar(shrink=.92)

def lae_based_adv(T_init,func):
	T = deepcopy(T_init)
	T_old = [[0.0 for x in range(imax)] for y in range(jmax)]	
	aP0 = [[0.0 for x in range(imax-1)] for y in range(jmax-1)]
	rho,cp =1.0,1.0
	for i in range(0,imax-1):
		for j in range(0,jmax-1):
			aP0[j][i] = rho*cp*dx[i]*dy[j]/dt
	while(error(T,T_old)>epsilon):
		T_old = deepcopy(T)
		Te = deepcopy(T)
		Tn = deepcopy(T)
		for j in range(1,jmax-1):
			for i in range(2,imax-2):
				if func == 'fou':
					we1,we2,we3 = 0.0,1.0,0.0
				elif func == 'sou':
					we1,we2,we3 = 0.0, (2*dx[i-1] + dx[i-2])/(dx[i-1] + dx[i-2]),-dx[i-1]/(dx[i-1] + dx[i-2])
				elif func == 'quick':
					we1 = dx[i-1]*(2*dx[i-1] + dx[i-2])/(dx[i-1] + dx[i])/(2*dx[i-1] + dx[i-2] + dx[i])
					we2 = dx[i]*(2*dx[i-1] + dx[i-2])/(dx[i-1] + dx[i])/(dx[i-1] + dx[i-2]) 
					we3 = -1.0* dx[i]*dx[i-1]/(dx[i-1] + dx[i-2])/(2*dx[i-1] + dx[i-2] + dx[i])
				Te[j][i] =  we1*T[j][i+1] +we2*T[j][i] + we3*T[j][i-1] 
		for i in range(1,imax-1):
			Tn[jmax-2][i] = T[jmax-1][i] 
		for j in range(0,jmax-2):
			for i in range(1,imax-1):
				if func == 'fou':
					wn1,wn2,wn3 = 0.0,1.0,0.0
				elif func == 'sou':
					wn1,wn2,wn3 = 0.0, (2*dy[j] + dy[j+1])/(dy[j] + dy[j+1]),-dy[j]/(dy[j] + dy[j+1])
				elif func == 'quick':
					wn1 = dy[j]*(2*dy[j] + dy[j+1])/(dy[j] + dy[j-1])/(2*dy[j] + dy[j+1] + dy[j-1])
					wn2 = dy[j-1]*(2*dy[j] + dy[j+1])/(dy[j] + dy[j-1])/(dy[j] + dy[j+1])
					wn3 = -1.0* dy[j-1]*dy[j]/(dy[j] + dy[j+1])/(2*dy[j] + dy[j+1] + dy[j-1])
				Tn[j][i] = wn1*T[j][i] + wn2*T[j+1][i] + wn3*T[j+2][i]
		b =  [[0.0 for x in range(imax-1)] for y in range(jmax-1)]
		for i in range(1,imax-1):
			for j in range(1,jmax-1):	
				b[j][i]=aP0[j][i]*T_old[j][i] - (u*Te[j][i] - u*Te[j][i-1]) *dy[j] - (v*Tn[j-1][i] - v*Tn[j][i])*dx[i]
		Error=1.0
		while (Error>=epsilon):
			T_old_iter=deepcopy(T)
			for i in range(1,imax-1):
				for j in range(1,jmax-1):
					T[j][i] = b[j][i]/aP0[j][i]
			for i in range(imax):
				T[0][i] = T[1][i]
			for j in range(jmax):	
				T[j][imax-1] = T[j][imax-2]
			Error=max(abs(T[j][i]-T_old_iter[j][i]) for i in range(1,imax-1) for j in range(1,jmax-1))
	phi = [[0.0 for x in range(imax)] for y in range(jmax)]
	for i in range(0,imax):
			for j in range(0,jmax):
				phi[j][i] = (T[j][i] -Tb) / (Tl - Tb)		
	return phi

if __name__ == '__main__':
	L= 1.0
	H =1.0
	u =1.0
	v =1.0
	imax =32 ######
	jmax =32 #####
	T0 = 0.0
	Tb = 0.0
	Tl = 100.0
	dTc =100.0
	epsilon =0.01
	T_init = [[T0 for x in range(imax)] for y in range(jmax)]
	for i in range(0,imax):
		T_init[jmax-1][i] = Tb
	for j in range(0,jmax):
		T_init[j][0] = Tl
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
	dx = [ 0.0 for x in range(imax-1)]
	for i in range(0,imax-1):
		dx[i] = xc[i+1]-xc[i]
	psi2 = np.linspace(0,1,jmax-1)
	yj = np.linspace(0,H,jmax-1)
	for j in range(0,jmax-1):
		a = ((beta+1)/(beta-1))**(2*psi2[j]-1)
		yj[j]=H*((beta + 1)*(a)-(beta-1))/(2*(1+a))
	yc = [ 0.0 for y in range(jmax)]
	for j in range(1,jmax-1):
		yc[j] = (yj[j]+ yj[j-1])/2.0
	yc[0] = yj[0]
	yc[jmax-1] = yj[jmax-2]
	dy = [ 0.0 for y in range(jmax-1)]
	for j in range(0,jmax-1):
		dy[j] = yc[j+1]-yc[j]
	dt = 0.1*0.5*min((dx[i]/u + dy[i]/v) for i in range(0,imax-1))
	phi1 = lae_based_adv(T_init,'fou')
	plot_fig(phi1)	
	plt.savefig('imp1.png')
	phi3 = lae_based_adv(T_init,'sou')
	plot_fig(phi3)
	plt.savefig('imp3.png')
	phi2 = lae_based_adv(T_init,'quick')
	plot_fig(phi2)
	plt.savefig('imp2.png')
	fou = temp_at_centreline(phi1)
	sou = temp_at_centreline(phi3)
	quick = temp_at_centreline(phi2)
	y = np.linspace(H,0,jmax)
	plt.figure()
	plt.plot(fou,y)
	plt.title('Variation of Nondimensional Temperature along Vertical Centerline(FOU)')
	plt.savefig('q3_1.png')
	plt.figure()
	plt.plot(sou,y)
	plt.title('Variation of Nondimensional Temperature along Vertical Centerline(SOU)')
	plt.savefig('q3_2.png')
	plt.figure()
	plt.plot(quick,y)
	plt.title('Variation of Nondimensional Temperature along Vertical Centerline(QUICK)')
	plt.savefig('q3_3.png')
