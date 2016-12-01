import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy

def error(T,T_old):
	sum0 = 0.0
	for j in range(0,jmax):
		for i in range(0,imax):
			sum0 += (T_old[j][i] - T[j][i])**2
	error = np.sqrt(sum0/((imax-2)*(jmax-2)))
	print("a1",error)
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

def lae_based_conv(T_init,func):
	T = deepcopy(T_init)
	T_old = [[0.0 for x in range(imax)] for y in range(jmax)]	
	aE = [[0.0 for x in range(imax-1)] for y in range(jmax-1)]
	aN = [[0.0 for x in range(imax-1)] for y in range(jmax-1)]
	aP0 = [[0.0 for x in range(imax-1)] for y in range(jmax-1)]
	aP = [[0.0 for x in range(imax-1)] for y in range(jmax-1)]
	for i in range(0,imax-1):
		for j in range(0,jmax-1):
			aE[j][i] = k*dy[j]/dx[i]
			aN[j][i] = k*dx[i]/dy[j]
			aP0[j][i] = rho*cp*dx[i]*dy[j]/dt
			aP[j][i] = aP0[j][i] + aE[j][i] + aE[j][i-1] +aN[j][i] +aN[j-1][i]
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
					we2 = dx[i]*(2*dx[i-1] + dx[i-2])/(dx[i-1] + dx[i])/(dx[i-1] + dx[i-2]) -1.0
					we3 = -1.0* dx[i]*dx[i-1]/(dx[i-1] + dx[i-2])/(2*dx[i-1] + dx[i-2] + dx[i])
				Te[j][i] =  we1*T[j][i+1] +we2*T[j][i] + we3*T[j][i-1]
		for i in range(1,imax-1):
			Tn[jmax-2][i] = T[jmax-1][i]  
		for j in range(1,jmax-2):
			for i in range(1,imax-1):
				if func == 'fou':
					wn1,wn2,wn3 = 0.0,1.0,0.0
				elif func == 'sou':
					wn1,wn2,wn3 = 0.0, (2*dy[j-1] + dy[j-2])/(dy[j-1] + dy[j-2]),-dy[j-1]/(dy[j-1] + dy[j-2])
				elif func == 'quick':
					wn1 = dy[j]*(2*dy[j] + dy[j+1])/(dy[j] + dy[j-1])/(2*dy[j] + dy[j+1] + dy[j-1])
					wn2 = dy[j-1]*(2*dy[j] + dy[j+1])/(dy[j] + dy[j-1])/(dy[j] + dy[j+1]) -1.0
					wn3 = -1.0* dy[j-1]*dy[j]/(dy[j] + dy[j+1])/(2*dy[j] + dy[j+1] + dy[j-1])
				Tn[j][i] = wn1*T[j][i] + wn2*T[j+1][i] + wn3*T[j+2][i]
		b =  [[0.0 for x in range(imax-1)] for y in range(jmax-1)]
		for i in range(0,imax-1):
			for j in range(0,jmax-1):
				b[j][i]=aP0[j][i]*T_old[j][i] - (u*Te[j][i] - u*Te[j][i-1]) *dy[j] - (v*Tn[j][i] - v*Tn[j-1][i])*dx[i]
		Error=1.0
		while (Error>=epsilon):
			T_old_iter=deepcopy(T)
			for i in range(1,imax-1):
				for j in range(1,jmax-1):
					T[j][i] = (aE[j][i]*T[j][i+1] + (aE[j][i-1]+ rho*u*dy[j])*T[j][i-1] + (aN[j-1][i]+ rho*v*dx[i])*T[j-1][i] + aN[j][i]*T[j+1][i] + b[j][i])/aP[j][i]
			Error=max(abs(T[j][i]-T_old_iter[j][i]) for i in range(1,imax-1) for j in range(1,jmax-1))
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
	L= 6.0
	H =1.0
	u =1.0
	v =0.0
	imax =62
	jmax =22
	T_wb  = 100.0	#BC's
	T_nb,T_sb  = 0.0 , 0.0
	k = 10.0
	rho, cp = 7.0, 100.0
	alpha= k/(rho*cp)
	epsilon =0.0001
	T_init = [[0.0 for x in range(imax)] for y in range(jmax)]
	for i in range(0,imax):
		T_init[0][i] = T_nb
		T_init[jmax-1][i] = T_sb
	for j in range(0,jmax):
		T_init[j][0] = T_wb
		T_init[j][imax-1] = T_init[j][imax-2]
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
	dt2 = 0.5*min(dx[i]/u for i in range(0,imax-1))
	dt1 = 0.5/alpha*(min(1/(dx[i]**2) for i in range(0,imax-1)) + + min(1/(dy[j]**2) for j in range(0,jmax-1))) 
	dt =  min(dt1,dt2) 
	phi1 = lae_based_conv(T_init,'fou')
	plot_fig(phi1)
	plt.title('2-D Heat Convection using Implicit Method in FOU scheme')	
	plt.savefig('impl_con_f0ou.png')
	phi2 = lae_based_conv(T_init,'quick')
	plot_fig(phi2)
	plt.title('2-D Heat Convection using Implicit Method in QUICK scheme')
	plt.savefig('impl_con_qui0ck.png')
	y = np.linspace(L,0 ,jmax)
	Ta,Tb,Tc,Td,Te = temp_at_x(phi1)
	plt.figure()
	plt.title('Variation of Nondimensional Temperature along Vertical lines(FOU)')
	plt.plot(Ta,y,label = "x/H = 0.2")
	plt.plot(Tb,y,label = "x/H = 0.4")
	plt.plot(Tc,y,label = "x/H = 0.6")
	plt.plot(Td,y,label = "x/H = 0.8")
	plt.plot(Te,y,label = "x/H = 1.0")
	plt.legend(loc='upper right')
	plt.savefig('q3_1.png')		
	Ta1,Tb1,Tc1,Td1,Te1 = temp_at_x(phi2)
	plt.figure()
	plt.title('Variation of Nondimensional Temperature along Vertical lines(QUICK)')
	plt.plot(Ta1,y,label = "x/H = 0.2")
	plt.plot(Tb1,y,label = "x/H = 0.4")
	plt.plot(Tc1,y,label = "x/H = 0.6")
	plt.plot(Td1,y,label = "x/H = 0.8")
	plt.plot(Te1,y,label = "x/H = 1.0")
	plt.legend(loc='upper right')
	plt.savefig('q3_2.png')		
