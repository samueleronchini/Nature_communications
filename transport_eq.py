from numpy import *
import numpy as np
import matplotlib.pyplot as plt


p=2.7

gamma_min=5.0
gamma_max=7.0


k=1e-18
t0=300
alfa=1.0
lam=1.0
y=0.0
t_inj=300.0
jmax=10000

h=1e-1
dt=1.e0


def F_0(gamma):

	if gamma>gamma_max or gamma<gamma_min:

		return 0.0

	else:
		return 10**50*10**(-p*gamma)

gamma_arr=np.arange(gamma_min-2.5,gamma_max+0.5,h)


N=[]
Np=[]

def gamma_dot(gamma,t):
	return -k*(1+t/t0)**(-2*lam)*10**(2*gamma)-alfa*2/(3*t0)*(10**gamma)/(1+t/t0)


N0=[]

for n in range(0,len(gamma_arr)):
	Np.append(0.0)
	N.append(0.0)



var=1.0

def N_inj(gamma,t):
	if t<t_inj:
		return F_0(gamma)*10**gamma*log(10)*h
	else:
		return var*(F_0(gamma)*10**gamma*log(10)*h)*(t/t_inj)**y

plt.yscale('log')
plt.plot(gamma_arr,N)


N_app=Np


count_0=t_inj
count=count_0
time_step=1000
time_step_peak=1000
count_max=int(jmax/time_step)

gamma_peak=[]
time_peak=[]


N_step=N
N_evo=[]
for n in range(0,count_max):
	N_evo.append(open('N_evo_%s.txt' %n, 'w'))
m=0
j=0
der=[]
B0=1
while j<jmax:
	t=j*dt

	if j==count:
		print(j)
		count=(m+1)*time_step+count_0
		for n in range(0,len(gamma_arr)):
			N_step[n]=Np[n]/(10**gamma_arr[n]*log(10)*h)
			N_evo[m].write('%s\t%s\t%s\n' %(10**gamma_arr[n],N_step[n], B0*(1+t/t0)**(-lam)))
		m+=1
		plt.plot(gamma_arr,N_step)

	for n in range(0,len(N)):
		delta=-gamma_dot(gamma_arr[n],t)*dt/(log(10)*10**gamma_arr[n])

		if n<len(N)-1:
			delta1=-gamma_dot(gamma_arr[n+1],t)*dt/(log(10)*10**gamma_arr[n+1])
			N_app[n]=(1-delta/h)*(Np[n]+N_inj(gamma_arr[n],t)*dt)+delta1/h*(Np[n+1]+N_inj(gamma_arr[n+1],t)*dt)
		else:
			N_app[n]=(1-delta/h)*(Np[n]+N_inj(gamma_arr[n],t)*dt)			


	Np=N_app

	if j==jmax-1:
		for n in range(0,len(gamma_arr)):
			N[n]=Np[n]/(10**gamma_arr[n]*log(10)*h)

	if j==time_step_peak :
		if j*dt>0:
			index=N.index(max(N))
			gamma_peak.append(gamma_arr[index])
			time_peak.append(j*dt)
		time_step_peak+=200



	j+=1

plt.plot(gamma_arr,N)
plt.title('t=%.1e s, tc=%.1e s, t0=%.1e s'%(j*dt,1/(k*10**(gamma_min)),t0))
plt.yscale('log')
plt.ylim(1e10,1e40)

plt.show()


